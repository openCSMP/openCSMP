// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DFM_PermeabilityAndFluxRatioAnalyzer.cpp
//  CSMP_DFM_Upscaling
//
//
#include "DFM_PermeabilityAndFluxRatioAnalyzer.h"
#include "ModelTime.h"
#include "CSMP_definitions.h"
#include "compareFloats.h"

// the model
#include "Region.h"
#include "Boundary.h"
#include "Model.h"
#include "BoundaryInterface.h"

#include "SteadyStateDiffusor.h"

// Interrelations
#include "ConstantFactor.h"

#include "Standard_IO_Handler.h"
#include "CSMP_highLevelUtilities.h"
#include "ModelTime.h"

// Analysis
#include "VTK_Interface.h"
#include "VTU_Interface.h"
// to analyze boundary fluxes in this linear model
#include "NodeCenteredFiniteVolumeTransport.h"
#include "EquivalentPermeabilityTensor.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
DFM_PermeabilityAndFluxRatioAnalyzer<dim>::DFM_PermeabilityAndFluxRatioAnalyzer( const char* model_name, bool verbose )
 : matrix_region_("matrix"), fracture_region_("fractures"),
   xsect_area_(0.),
   model_(string{model_name}),
   monitor_( model_, string("volume flux"), string("volume flux"), false ),
   stdio_( model_name ),
   flux_histogram_( model_ ),
   model_dimensions_(dim,0.),
   k_equivalent_(3,0.),
   flux_ratio_(3,0.),
   matrix_permeability_(3,0.),
   verbose_(verbose)
{
  cout <<"\n\nCSMP++ FRACTURE - MATRIX EQUIVALENT PERMEABILITY & FRACTURE-MATRIX FLUX RATIO ANALYZER"<< endl;
  cout <<"copyright (c) Stephan K. Matthai (2013)."<< endl << endl;
//  PrintSpecifications();
  // fluid volume source
  model_.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
  model_.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
  const Region<dim>& model_domain(model_.Region("Model"));
  model_volume_      = model_domain.Volume();
  model_pore_volume_ = model_domain.Volume(true);

  // model dimensions
  cout <<"\n\nModel: '"<< model_name <<"', physical dimensions: ";
  printModelDimensions( model_ );
  cout <<"\n\tmodel_ volume: "<< model_volume_ <<" and pore volume: "<< model_pore_volume_ << endl;
  cout <<"\n\tmodel_ contains the unique regions: ";
  for ( auto rit=model_.UniqueRegionsBegin(); rit!=model_.UniqueRegionsEnd(); rit++ )
    cout << (*rit).first <<"  ";
  cout <<"\n\n";
  Point<dim> xyz_min, xyz_max;
  model_.MinMaxCoordinates( xyz_min, xyz_max );
  model_dimensions_[0] = xyz_max[0] - xyz_min[0];
  model_dimensions_[1] = xyz_max[1] - xyz_min[1];
  if constexpr ( dim == 3 ) model_dimensions_[2] = xyz_max[2] - xyz_min[2];

  // screening for surface bumps greater than 1% of the model_ dimensions
  if ( dim == 3 && DetectOutOfPlaneNodesAtBoxBoudaries( model_name, 0.01 ) )
    throw csmp::Exception( FATAL_ERROR, "DFM_PermeabilityAndFluxRatioAnalyzer::(constructor):",
                          "model boundaries contain irregularities; check model before continueing." );

  // this will later be used in search operations involving finding fractures in non-unique regions
  InitializeRegionIdentifiers("region identifier");

// ------------------------------------------------------------------------------------
// 1. Establishing initial & essential conditions for testing
// ------------------------------------------------------------------------------------
  list<string> integral_properties;  integral_properties.push_back("volume flux");
  list<string> range_properties;     range_properties.push_back("fluid pressure");


// ------------------------------------------------------------------------------------
// 2. Analysis of results (analyser is used in loop below)
// ------------------------------------------------------------------------------------
   flux_histogram_.DefineBins( "flux_histogram.bins", velo_bins_ );
  

// ------------------------------------------------------------------------------------
// 3. SAMG-Solver settings for steady-state fluid pressure algorithm [K]{p} = {Q}
// ------------------------------------------------------------------------------------
/*
  // targeting 2012' SAMG DLL 1 for this pressure solver
  fluid_pressure_.GetSolverSettings().SetSolverInstance(1);
  // iout
  fluid_pressure_.GetSolverSettings().Set_iout1( 0 );
  fluid_pressure_.GetSolverSettings().Set_iout2( 0 );
  if ( !verbose_ ) fluid_pressure_.GetSolverSettings().Set_idmp( -1 );
  // one-time solver set-up: nothing is remembered for next try
  fluid_pressure_.GetSolverSettings().Set_iswit(5);
  // SAMG solution criteria
  fluid_pressure_.GetSolverSettings().Set_eps(0.); // absolute criterion
  // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
  fluid_pressure_.GetSolverSettings().Set_a_cmplx(2);
  // Pre-adjust SAMG mesh complexity, based on solver output
  fluid_pressure_.GetSolverSettings().Set_g_cmplx(1.5);
  fluid_pressure_.GetSolverSettings().Set_w_avrge(2);

#ifdef SAMG_OUTPUT_TO_FILE
  // trigger output to file
  fluid_pressure.GetSolverSettings().Set_idmp( 8 );
  // // define SAMG file output format for reduced file size, idmp > 1 is required
  fluid_pressure.GetSolverSettings().Set_ioform( "f" );
  // set filename for SAMG file output other than default "level", idmp > 1 is required
  fluid_pressure.GetSolverSettings().Set_filnam_dump( "ReservoirSimulator_steady_state_p" );
#endif
*/

} // end constructor





/**
    initialize element variable integers that shall identify unique model regions by number
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::InitializeRegionIdentifiers( const char* region_identifying_variable )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const csmp::Index idx_key(model_.Database().StorageKey(region_identifying_variable));
    if ( idx_key.place != ELEMENT ) {
         csmp_error.Note( ERROR, "DFM_PermeabilityAndFluxRatioAnalyzer<dim>::InitializeRegionIdentifiers",
                                   "placement of region identifier key must be on the element.");
         return;
      }
   
    double region_identifier(1.);

    for ( auto it=model_.UniqueRegionsBegin(); it!=model_.UniqueRegionsEnd(); ++it ) {
         for ( auto eit=(*it).second.CellsBegin(); eit!=(*it).second.CellsEnd(); ++eit )
           (*eit)->Store( idx_key, makeScalar(PLAIN, region_identifier) );
         region_identifier += 1.;
      }
 
 } // end InitializeRegionIdentifiers
  




 
/**
    Creates and writes histograms of the "volume flux" to a text file.
    This output can be pasted into Maple.
    The file name is appended the direction in which the velocity is measured
    translating the coordinate axis into x, y, z
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::OutputVelocityHistogramToMaple( size_t i ) const
 {
    map<string,pair<vector<pair<double,double> >,size_t> >  results;
	  flux_histogram_.RegionPropertyHistogramsElement( "flow speed", velo_bins_, results );
   
	  if      ( i == 0 ) flux_histogram_.OutputRegionPropertyHistogramsMaple( "flow-speed-x", velo_bins_, results, true );
	  else if ( i == 1 ) flux_histogram_.OutputRegionPropertyHistogramsMaple( "flow-speed-y", velo_bins_, results, true );
	  else               flux_histogram_.OutputRegionPropertyHistogramsMaple( "flow-speed-z", velo_bins_, results, true );
 }
	   
//cout <<"\nflux histogram vector.\n";
//for ( map<string,pair<vector<pair<double,double> >,size_t> >::iterator it=results.begin(); it!=results.end(); it++ )
//  {
//     cout <<"\n"<< (*it).first <<": ";
//     for ( vector<pair<double,double> >::iterator j=(*it).second.first.begin(); j!=(*it).second.first.end(); j++ )
//       cout << (*j).first <<","<< (*j).second <<"  ";
//  }
   



template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::PrintSpecifications() const
{
  cout <<"Determination of equivalent permeability through flow-based upscaling and\
  calculation of fracture - matrix flux ratio (Nelson's 'excess permeability' of NFRs)\
  \n\nInput:  box-shaped CSMP model_ with target fracture apertures expressed in terms of the 'thickness'\
  of lower-dimensional fracture elements (surfaces in this 3D model_):\
  \
  \n\n - fractures / fracture sets must be identified by the string 'FRAC' or 'SET' in their names\
  \n\n - rock matrix must contain 'MATRIX' in their names\
  \n\n - principal directions of permeability should be aligned with the coordinate axes of the model_\
  \n\n - matrix and fracture permeabilities must be scalars\
  \n\n - thickness must be defined on all elements; volumetric ones must have value of 1\
  \n\n - variables needed: 'thickness', 'permeability', 'hydraulic conductivity', 'velocity', volume flux'\
  \
  \n\n Approach: LFEM, flux integration using FVM\
  \
  \n\n Setup of virtual flow experiment: pressure - pressure, steady-state\
  \
  \n\n Output: diagonal values of k_equiv and qf/qm ratio tensors\
  \
  (formerly CSP FRACTURE - MATRIX SINGLE PHASE FLOW ANALYZER)\
  \
  \n tested: O.K. SKM 22/7/2006\
  \n revised by SKM 30/10/2013";
  
} // end




  // display effective properties to screen
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::PrintResultsToScreen( const char* model_name ) const
{
  cout <<"\n\n________________________________________________________________________________\n\n";
  cout <<"\n                              ANALYSIS RESULTS";
  cout <<"\n\n________________________________________________________________________________\n\n";
  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: The permeability range of the model_ was:"<< endl;
  printRangeOfVariable( model_, "permeability" );

  cout <<"Model: '"<< model_name <<", physical dimensions: ";
  printModelDimensions( model_ );
  cout <<"\nP32: "<< fracture_surface_area_ / model_volume_;
  cout <<"\nfracture-matrix interface area: "<< fracture_surface_area_ * 2. <<" (m2).";
  cout <<"\nfraction of void space due to fractures: "<< fracture_matrix_void_ratio_;
  cout <<"\ntotal fracture volume of model: "<< FractureVolume() <<" (m3).\n";

  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: The equivalent (flow-based upscaling) permeability (m2) of the model is: ";
  cout <<"\nin the x-direction: " << k_equivalent_[0];
  cout <<"\nin the y-direction: " << k_equivalent_[1];
  cout <<"\nin the z-direction: " << k_equivalent_[2];
  cout << endl;
  
  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: Contribution of the fractures relative to the matrix: ";
  cout <<"\nin the x-direction: " << flux_ratio_[0];
  cout <<"\nin the y-direction: " << flux_ratio_[1];
  cout <<"\nin the z-direction: " << flux_ratio_[2];
  cout << endl;

  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: The matrix-only permeability (m2) is: ";
  cout <<"\nin the x-direction: " << matrix_permeability_[0];
  cout <<"\nin the y-direction: " << matrix_permeability_[1];
  cout <<"\nin the z-direction: " << matrix_permeability_[2];
  cout << endl;
  
} // end printResultsToScreen




template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::WriteResultsToFile( const char* model_name ) const
{
  string effective_properties_file_name(model_name);
  effective_properties_file_name += "-equivalent-properties.text";
  ofstream  ofs;
  ofs.open( effective_properties_file_name.c_str(), ios::out|ios::trunc );

  double dim_x, dim_y, dim_z;
  if constexpr ( dim == 3 )
    ofs <<"Model: '"<< model_name <<", physical dimensions: "<< boundingBox( model_, dim_x, dim_y, dim_z );
  else if constexpr ( dim == 2 ) {
     Point<2> xy_min, xy_max;
     model_.MinMaxCoordinates( xy_min, xy_max );
     ofs <<"Model: '"<< model_name <<", physical dimensions: ";
     ofs << xy_max[0]-xy_min[0] <<"m (x=width), "<< xy_max[1]-xy_min[1] <<"m (y=height)."<< endl;
  }
  ofs <<"\nP32: "<< fracture_surface_area_ / model_volume_;
  ofs <<"\nfracture-matrix interface area: "<< fracture_surface_area_ * 2.;
  ofs <<"\nfraction of void space due to fractures: "<< fracture_matrix_void_ratio_ << endl;

  ofs <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: The equivalent permeability (m2) of the model is: ";
  ofs <<"\nin the x-direction: " << k_equivalent_[0];
  ofs <<"\nin the y-direction: " << k_equivalent_[1];
  ofs <<"\nin the z-direction: " << k_equivalent_[2];
  ofs << endl;
  
  ofs <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: Contribution of the fractures relative to the matrix: ";
  ofs <<"\nin the x-direction: " << flux_ratio_[0];
  ofs <<"\nin the y-direction: " << flux_ratio_[1];
  ofs <<"\nin the z-direction: " << flux_ratio_[2];
  ofs << endl;

  ofs <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: The matrix-only permeability (m2) is: ";
  ofs <<"\nin the x-direction: " << matrix_permeability_[0];
  ofs <<"\nin the y-direction: " << matrix_permeability_[1];
  ofs <<"\nin the z-direction: " << matrix_permeability_[2];
  
  ofs.close();
  
} // end WriteResultsToFile




/**
     For the 3 spatial directions i=0..2, VTK files are written with pressure distributions etc.
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::WriteResultsToVTK( const char* model_name, size_t i )
 {
   static VTK_Interface<dim>  vtk_output;
   static bool first_call(true);
   if ( first_call ) vtk_output.OutputDataToVTK( model_, "permeability", "permeability", 0 );

   // dependent variable ranges in groups
   for ( typename map<string,Region<dim> >::const_iterator
         grit=model_.RegionsBegin(); grit!=model_.RegionsEnd(); grit++ ) {
        printRangeOfVariable( model_, stdio_,(*grit).first.c_str(),  "fluid pressure" );
        printRangeOfVariable( model_, stdio_, (*grit).first.c_str(), "velocity" );
        printRangeOfVariable( model_, stdio_, (*grit).first.c_str(), "volume flux" );
        // output results for visualization using VTK
        vtk_output.OutputDataToVTK( model_, (*grit).first.c_str(), "fluid-pressure", "fluid pressure", static_cast<size_t>(i), true );
        vtk_output.OutputDataToVTK( model_, (*grit).first.c_str(), "velocity",  "velocity", static_cast<size_t>(i), true );
        vtk_output.OutputDataToVTK( model_, (*grit).first.c_str(), "volume-flux", "volume flux", static_cast<size_t>(i), true );
     }

   // output results for visualization using VTK
   vtk_output.OutputDataToVTK( model_, string{model_name} + "_permeability",   "permeability",    i, true );
   vtk_output.OutputDataToVTK( model_, string{model_name} + "_velocity",       "velocity",        i, true );
   vtk_output.OutputDataToVTK( model_, string{model_name} + "_fluid-pressure", "fluid pressure",  i, true );
   
   first_call = false;
  
 } // end WriteResultsToVTK








/**
    Converts BOX_BOUNDARY enumeration values o variable values 
    and outputs these to VTK for visualisation.
    
    @attention use this to test boundary flagging.
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::OutputBoxBoundaryFlagsAsNumbers( const char* variable )
 {
    const csmp::Index p_key = model_.Database().StorageKey(variable);
    // are the boundary nodes flagged as expected for box-shaped model?
    model_.InputPropertyValue( variable, makeScalar(ANY,0.) );
    Region<dim>& model_domain(model_.Region("Model"));
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      (*nit)->Store( p_key, makeScalar(ANY,static_cast<double>((*nit)->AtBoundary())) );

    VTK_Interface<dim>  vtk_output;
    vtk_output.OutputDataToVTK( model_, "box_boundary_flags", "fluid pressure", 999 );
   
 } // end TestBoundaryIntegrity






/**
     Computes the absolute volumes of the supplied domain using the thickness attribute
     which ascertains that elements in the domain with different dimensionality are counted correctly.
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::RegionVolume( const char* region ) const
 {
    const Region<dim>&  ref(model_.Region(region));
    const csmp::Index   thi_key = model_.Database().StorageKey("thickness");
    double volume(0.);
    for ( typename vector<Element<dim>*>::const_iterator
          eit=ref.CellsBegin(); eit!=ref.CellsEnd(); ++eit ) {
         volume += (*eit)->Volume() * (*eit)->Read( thi_key );
      }
    return volume;

} // end regionVolume


/**
   If the side boundaries of the model_ have defects, the permeability calculation becomes inaccurate.
   This methods detect this and reports boundaries and locations of spurious nodes.
   
   Assumption: the model_ is box-shaped and only has such boundaries.
   
   @test O.K. SKM 1/11/2013
*/
template<uint32_t dim>
bool DFM_PermeabilityAndFluxRatioAnalyzer<dim>::DetectOutOfPlaneNodesAtBoxBoudaries( const char* model__name,
                                                                                     double relative_tolerance ) const
 {
    assert( dim == 3 );
    bool with_boudary_defects(false);
    VTU_Interface<dim>  vtu_output( model_, "box_boundary_defects" );
   
    // For each boundary of the box-shaped model_, find the bounding box which should be a plane.
    // If there the smallest dimension is larger than 0.1% of the largest one the boundary will be output
    // and an error will be reported.
    set<Point<dim> > bpoints;
    for ( typename Model<dim>::boundaryConstIterator bit=model_.BoundariesBegin(); bit!=model_.BoundariesEnd(); bit++ )
      {
         for ( typename vector<Node<dim>*>::const_iterator
               it=(*bit).second.NodesBegin(); it!=(*bit).second.NodesEnd(); ++it )
           bpoints.insert( (*it)->Coordinate() );
        
         // checking points in the set
         Point<dim> bdiag = (*bpoints.begin()) - (*bpoints.rbegin());
         // reporting faulty boundary to VTU file
         const double box_max = std::max( std::max( fabs(bdiag[0]), fabs(bdiag[1])), fabs(bdiag[2]) );
         const double max_undulation = box_max * relative_tolerance;
         if ( fabs(bdiag[0]) > max_undulation && fabs(bdiag[1]) > max_undulation && fabs(bdiag[2]) > max_undulation ) {
              with_boudary_defects = true;
              cout <<"\ndetectOutOfPlaneNodesAtBoxBoudaries: Problem: bounding box dimensions of boundary: ";
              cout << (*bit).first <<"\n";
              cout <<"\t"<< fabs(bdiag[0]) <<"  "<< fabs(bdiag[1]) <<"  "<< fabs(bdiag[2]) << endl;
              vtu_output.OutputDataToVTU( model__name, "fluid pressure", (*bit).second, 0 );
           }
         bpoints.clear();
      }
   
   return with_boudary_defects;

 } // end DetectOutOfPlaneNodesAtBoxBoudaries





/**
     Groups all regions containing MATRIX or FRAC or SET  in their names 
     into the corresponding regions
 
     'fractures' and 'matrix'.
     
     Assuming that the fractures have a lower-dimensional representation, 
     the method returns the surface area of all fractures combined.
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::IdentifyFractureAndMatrixRegions()
 {
  bool         with_matrix(false), with_fractures(false);
  set<string>  matrix_regions;
  double     total_fracture_area(0.);
  
  // FRACTURES:  looping over all the unique regions that contain fractures and combining these
  // but only if this has not already been done
  if ( !model_.ContainsRegion("fractures") ) {
      set<string>  fracture_regions;
      for ( typename Model<dim>::regionConstIterator rit=model_.UniqueRegionsBegin(); rit!=model_.UniqueRegionsEnd(); rit++ )
        {
           // if the region name contains 'FRACTURE' it is pooled into an non-unique region 'fractures'
           if ( (*rit).first.find("FRAC") != std::string::npos || (*rit).first.find("SET") != std::string::npos) {
                // recording the region name
                fracture_regions.insert((*rit).first);
                total_fracture_area += (*rit).second.Volume();
                if ( with_fractures == false ) {
                     model_.CopyRegion( (*rit).first.c_str(), "fractures" );
                     with_fractures = true;
                  }
                else model_.AssimilateRegion( (*rit).first.c_str(), "fractures" );
             }
        }

      if ( fracture_regions.empty() )
        throw csmp::Exception( FATAL_ERROR, "IdentifyFractureAndMatrixRegions", "there is no region tagged 'FRACTURE' or 'SET' in the model; quitting.");
    }
  else total_fracture_area = model_.Region("fractures").Volume(); // = surface area given that fractures are lower dimensional


  // MATRIX:  looping over all the unique regions that contain fractures and or matrix and combining these
  for ( typename Model<dim>::regionConstIterator rit=model_.UniqueRegionsBegin(); rit!=model_.UniqueRegionsEnd(); rit++ )
    {
       // same for matrix
       if ( (*rit).first.find("MATRIX") != std::string::npos ) {
            // recording the region name
            matrix_regions.insert((*rit).first);
            if ( with_matrix == false ) {
                 model_.CopyRegion( (*rit).first.c_str(), "matrix" );
                 with_matrix = true;
              }
            else model_.AssimilateRegion( (*rit).first.c_str(), "matrix" );
         }
    }

  if ( matrix_regions.empty() )
    throw csmp::Exception( FATAL_ERROR, "IdentifyFractureAndMatrixRegions", "there is no region tagged 'MATRIX' in the model; quitting.");
  
  if ( matrix_regions.size() > 1 ) {
       cout <<"\nIdentifyFractureAndMatrixRegions: model_ contains multiple matrix regions: ";
       for ( set<string>::const_iterator it=matrix_regions.begin(); it!=matrix_regions.end(); it++ )
         cout << (*it) <<" ";
       cout <<"\n\n";
       throw csmp::Exception( ERROR, "IdentifyFractureAndMatrixRegions", "unfortunately, this tool can only handle single-valued matrix permeability");
    }
  
  return total_fracture_area;

 } // end IdentifyFractureAndMatrixRegions







/**
    Calculate fluxes across box boundaries, taking lower-dimensional elements into account.
    
    @attention since the normals of boundary faces are always outward pointing, 
    outgoing fluxes are reported as positive and incoming ones as negative.
    
    @test OK for volumetric-only meshes.
 
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::BoxBoundaryFluxFEM( BOX_BOUNDARY target_boundary ) const
 {
    const Region<dim>&  mref(model_.Region("Model"));
    const csmp::Index   v_key = model_.Database().StorageKey("velocity");
    const csmp::Index   t_key = model_.Database().StorageKey("thickness");
    double              boundary_flux{0.}, face_area{0.};
    VectorVariable<dim> velo;
    ScalarVariable      thickness;
    vector<size_t>      face_nids;
    vector<double>      nrml;
    Box().UnitNormalTo( target_boundary, dim, nrml );
    const Point<dim>    unrml(nrml);

    // for all boundary faces of the region model
    if constexpr ( dim == 3 ) {
        for ( auto eidx=mref.InteriorCells(); eidx<mref.Cells(); ++eidx )
          {
             // retrieve the velocity that must already be scaled by thickness if lower dimensional elements are used
             mref.E(eidx)->Read( v_key, velo );
             mref.E(eidx)->Read( t_key, thickness );
             // loop over the boundary faces of the elements
             for ( uint32_t j=0U; j<mref.PerimeterFaces(eidx); ++j )
               {  // accessing the one or multiple faces of the element that lie on the model boundary
                  const uint32_t face = mref.PerimeterFace(eidx,j);
                  // getting the normal to the face and its area
                  auto fnids = mref.E(eidx)->FE()->CornerNodesOfFace( face );
                  switch( fnids.size() ) {
                     case 2: // line
                         face_area = thickness(); // this is a point or an aperture
                       break;
                     case 3: // triangular face
                         face_area = triangleArea( mref.E(eidx)->N(fnids[0])->Coordinate(),
                                                   mref.E(eidx)->N(fnids[1])->Coordinate(),
                                                   mref.E(eidx)->N(fnids[2])->Coordinate() );
                       break;
                     case 4: // quadrilateral face
                         face_area = facetArea4( mref.E(eidx)->N(fnids[0])->Coordinate(),
                                                 mref.E(eidx)->N(fnids[1])->Coordinate(),
                                                 mref.E(eidx)->N(fnids[2])->Coordinate(),
                                                 mref.E(eidx)->N(fnids[3])->Coordinate() );
                       break;
                     default: throw csmp::Exception( ERROR, "DFM_PermeabilityAndFluxRatioAnalyzer<dim>::BoxBoundaryFluxFEM",
                                                    "type of boundary face was not recognised");
                  }
              }

             // computing the face flux Aj v . n
             boundary_flux += face_area * velo.DotProduct( unrml );
          }
       } // 3D

    // the corresponding boundaries are found
    if constexpr ( dim == 2 ) {
         const Boundary<2>& side = model_.Boundary( parseBoundary(target_boundary) );
         for ( const auto& face : side.CellVector() ) {
              // retrieve the velocity that must already be scaled by thickness if lower dimensional elements are used
              face->InnerParent()->Read( v_key, velo );
              boundary_flux += face->Area() * face->InnerParent()->Read(t_key) * velo.DotProduct( unrml );
           }
      } // end 2D

    return boundary_flux;

 } // end BoxBoundaryFluxFEM




/*  LOOPS OVER ENTIRE BOUNDARY OF MODEL GETTING A FLUX BALANCE
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::BoxBoundaryFluxFEM( BOX_BOUNDARY target_boundary ) const
 {
    const Region<dim>&   mref(model_.Region("Model"));
    const csmp::Index    v_key = model_.Database().StorageKey("velocity");
    const csmp::Index    t_key = model_.Database().StorageKey("thickness");
    double             boundary_flux(0.);
    VectorVariable<dim>  velo;

    // for all boundary faces of the region model
    for ( size_t eidx=mref.InteriorCells(); eidx<mref.Cells(); ++eidx )
      {
         // retrieve the velocity that must already be scaled by thickness if lower dimensional elements are used
         mref.E(eidx)->Read( v_key, velo );
         // loop over the boundary faces of the elements
         for ( uint32_t j=0U; j<mref.PerimeterFaces(eidx); ++j )
           {  // accessing the one or multiple faces of the element that lie on the model boundary
              const uint32_t face = mref.PerimeterFace(eidx,j);
              // getting the unit normal and area of target face
              Point<dim>  unrml = mref.E(eidx)->UnitNormalToFace( face );
              double    Aface = mref.E(eidx)->FaceArea( face ) * mref.E(eidx)->Read( t_key );
              // computing the face flux Aj v . n
              boundary_flux += Aface * velo.DotProduct( unrml );
           }
      }

    
    return boundary_flux;

 } // end BoxBoundaryFluxFEM
*/






/**
    Computes the fraction of the total pore space that is occupied by fractures.
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::FractureMatrixVoidRatio( bool verbose ) const
 {
    const Region<dim>&  fgref(model_.Region(fracture_region_.c_str()));
    const Region<dim>&  mgref(model_.Region(matrix_region_.c_str()));
    
    const csmp::Index phi_key = model_.Database().StorageKey("porosity");
    double matrix_pore_volume(0.);
    for ( auto eit=mgref.CellsBegin(); eit!=mgref.CellsEnd(); ++eit )
      matrix_pore_volume +=  (*eit)->Volume() * (*eit)->Read( phi_key );
      
    // fracture pore volume (of the lower-dimensional fractures) that is substracted from matrix volume
    const csmp::Index thi_key = model_.Database().StorageKey("thickness");
    double fracture_pore_volume(0.);
    for ( auto eit=fgref.CellsBegin(); eit!=fgref.CellsEnd(); ++eit )
      fracture_pore_volume += (*eit)->Volume() * (*eit)->Read( thi_key) * (*eit)->Read( phi_key );

    // fracture volume must be added since fractures and matrix overlap and fractures have higher porosity than matrix
    const double fracture_fraction_of_void_space = fracture_pore_volume / (matrix_pore_volume + fracture_pore_volume);
    if ( verbose ) {
         cout <<"\nFractureMatrixVoidRatio: fraction of the total void space that is due to fractures: ";
         cout << fracture_fraction_of_void_space <<"\n\n";
      }
   
    return fracture_fraction_of_void_space;

} // end FractureMatrixVoidRatio



/// computes the fracture volume
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::FractureVolume( bool verbose ) const
 {
    const Region<dim>&  fgref(model_.Region(fracture_region_.c_str()));
   
    const csmp::Index phi_key = model_.Database().StorageKey("porosity");
    // fracture pore volume (of the lower-dimensional fractures) that is substracted from matrix volume
    const csmp::Index thi_key = model_.Database().StorageKey("thickness");
    double fracture_pore_volume(0.);
    for ( auto eit=fgref.CellsBegin(); eit!=fgref.CellsEnd(); ++eit )
      fracture_pore_volume += (*eit)->Volume() * (*eit)->Read( thi_key) * (*eit)->Read( phi_key );

    if ( verbose ) {
         cout <<"\nFractureVolume of model: ";
         cout << fracture_pore_volume <<" (m3).\n\n";
      }
   
    return fracture_pore_volume;

} // end FractureMatrixVoidRatio








/**
    Computes matrix ensemble permeability by flow-based upscaling: works also for multilayer aggregates
    
    @attention Method does not need to know which way the flow is going it just measures how much 
    comes into the domain and how much goes out. However, the cross-sectional area has to be right.
    
    @test OK, SKM 29/5/14
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::EquivalentPermeability( const char* region,
                                                                          double& flux_through_model )
 {
    NodeCenteredFiniteVolumeTransport<dim>  advector( region, model_, "porosity", "fluid pressure",
                                                     "velocity", "nodal fluid volume source",
                                                      false, false, "thickness" );

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

	  // equivalent permeability computed from fluid throughput through model as determined by FVM
    double  inflow, outflow;
    // whichever shape of form this region has
    const bool ASSUME_BOX_SHAPED_MODEL(true);
    advector.BoundaryFluxes( inflow, outflow, ASSUME_BOX_SHAPED_MODEL, false, "no advected variable", "fluid pressure", PLAIN );
    // computing the maximum of the in- and outflow
    flux_through_model = (fabs(inflow) + fabs(outflow)) / 2.;
    // making sure that the difference between in- and outflux from the model is less than 1%
    const double flux_discrepancy_tolerance(0.03); // 3 percent.
    if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > flux_discrepancy_tolerance ) {
         cerr <<"\ninflux versus outflux: "<< inflow <<" vs. "<< outflow <<" m3/s.";
         csmp_error.Note( ERROR, "EquivalentPermeability", "flux analysis is inaccurate; check integrity of boundary mesh." );
         // using the maximum flux estimate in this case
         flux_through_model = std::max( fabs(inflow), fabs(outflow) );
      }
   
	  cout <<"\nEquivalentPermeability: x-sec. area="<< xsect_area_ <<" m2, total- and specific flux through model: "<< flux_through_model;
	  cout <<", "<< flux_through_model / xsect_area_ << endl;
    cout <<"\n\t"<<"estimated error of calculation (%): <= "<< 100. * (fabs(fabs(inflow)-fabs(outflow)) / flux_through_model) << endl << endl;
   
    // equivalent permeability
	  return (flux_through_model * fluid_viscosity_) / (xsect_area_ * farfield_pf_gradient_);

 } // end EquivalentPermeability




/**
    Generic version for any kind of inflow and outflow boundaries.
    
    @attention in- and out-fluxes are based on finite element approximation and are therefore not as accurate as FV based estimates.
*/
template<uint32_t dim>
double DFM_PermeabilityAndFluxRatioAnalyzer<dim>::EquivalentPermeability( BOX_BOUNDARY inflow_boundary,
                                                                          BOX_BOUNDARY outflow_boundary,
                                                                          double& flux_through_model )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

	  // fluid throughput through model as determined by velocity projection on model boundaries
    const double  inflow  = fabs(BoxBoundaryFluxFEM( inflow_boundary ));
    const double  outflow = fabs(BoxBoundaryFluxFEM( outflow_boundary ));

    // computing the maximum of the in- and outflow
    flux_through_model = (fabs(inflow) + fabs(outflow)) / 2.;
    // making sure that the difference between in- and outflux from the model is less than 1%
    const double flux_discrepancy_tolerance(0.03); // 3 percent.
    if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > flux_discrepancy_tolerance ) {
         cerr <<"\ninflux versus outflux: "<< inflow <<" vs. "<< outflow <<" m3/s.";
         csmp_error.Note( WARNING, "EquivalentPermeability", "flux analysis is inaccurate; check integrity of boundary mesh." );
         // using the maximum flux estimate in this case
         flux_through_model = std::max( fabs(inflow), fabs(outflow) );
      }
   
	  cout <<"\nEquivalentPermeability: x-sec. area="<< xsect_area_ <<" m2, total- and specific flux through model: "<< flux_through_model;
	  cout <<", "<< flux_through_model / xsect_area_ << endl;
    cout <<"\n\testimated error of calculation (%): <= "<< 100. * (fabs(fabs(inflow)-fabs(outflow)) / flux_through_model) << endl << endl;
   
    // equivalent permeability
	  return (flux_through_model * fluid_viscosity_) / (xsect_area_ * farfield_pf_gradient_);

 } // end EquivalentPermeability (using Boundaries)








/**
    taking into account the thickness attribute correctly
    @test: OK
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::ComputeTransmissivity()
{
  const csmp::Index k_key   = model_.Database().StorageKey("permeability");
  const csmp::Index K_key   = model_.Database().StorageKey("transmissivity");
  const csmp::Index thi_key = model_.Database().StorageKey("thickness");
  printRangeOfVariable( model_, stdio_, "permeability" );
  printRangeOfVariable( model_, stdio_, "thickness" );
  cout <<"\nComputeTransmissivity: note: 'thickness' here refers to 'fracture aperture'\n\n";
  Region<dim>&  model_domain(model_.Region("Model"));
  for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); it++ ) {
        const double k = ((*it)->Read(k_key) * (*it)->Read(thi_key)) / fluid_viscosity_;
       (*it)->Store( K_key, makeScalar(PLAIN,k) );
    }
} // end ComputeTransmissivity





/** 
    Computes Darcy- interstitial (pore-) velocity and the (scalar) volume flux
    taking into account the thickness attribute correctly.
 
    The flow speed is the magnitude of the interstitial velocity and is computed for
    the generation of flow velocity histograms.
    
    @attention the interstitial velocity is equivalent to the actual transport velocity 
    with which colloid particles would be transported by the flow.
    It is the one that has to be analyzed in the velocity histograms.
*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::ComputeVelocityAndVolumeFlux( const char* region )
{
   Region<dim>&  ref(model_.Region(region));
   const csmp::Index pf_key  = model_.Database().StorageKey("fluid pressure");
   const csmp::Index K_key   = model_.Database().StorageKey("transmissivity"); // contains conductivity * thickness product
   const csmp::Index phi_key = model_.Database().StorageKey("porosity");
   const csmp::Index thi_key = model_.Database().StorageKey("thickness");
   const csmp::Index vD_key  = model_.Database().StorageKey("velocity");
   const csmp::Index vi_key  = model_.Database().StorageKey("pore velocity");
   const csmp::Index vf_key  = model_.Database().StorageKey("volume flux");
   const csmp::Index fs_key  = model_.Database().StorageKey("flow speed");
   VectorVariable<dim>   velo;
   DenseMatrix<DM_MIN>  DERIV;
  
   for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
         // computing the total velocity: vt = -K grad P
         // (taking into account the thickness attribute of 2D elements)
         const double K = (*it)->Read( K_key );
         velo = 0.;
         (*it)->dN_AtBaryCenter( DERIV, 1U );
         for ( uint32_t i=0U; i<(*it)->Nodes(); i++ ) {
              double pf = (*it)->N(i)->Read( pf_key );
              velo(0)  += pf  * -DERIV(0,i) * K;
              velo(1)  += pf  * -DERIV(1,i) * K;
              if constexpr ( dim == 3 ) velo(2)  += pf  * -DERIV(2,i) * K;
           }
        
         // storing the thickness-weighted velocity
         (*it)->Store( vD_key, velo );
        
         // storing the (scalar) volume flux = velocity * thickness product
         const double volume_flux = velo.Length();
         (*it)->Store( vf_key, makeScalar(PLAIN,volume_flux) );
        
         // pore velocity (actual velocity with which the fluid moves through the pore/fracture space)
         double thickness = (*it)->Read( thi_key );
         velo /= (*it)->Read( phi_key ) * thickness;
         (*it)->Store( vi_key, velo );
        
         // flow speed (=scalar value of pore velocity) for the generation of flow velocity histograms
         const double flow_speed = velo.Length();
         (*it)->Store( fs_key, makeScalar(PLAIN,flow_speed) );
      }

} // end computeVelocityAndVolumeFlux





/**
    fracture - matrix equivalent permeability analysis
    
    Computes equivalent permeability using pressure - pressure boundary conditions.
    
    @attention: regions must consist out of fracture subregions that
    do not include multiple sets. The fracture regions must not be discontiguous.
    
    @attention: for lower-dimensional elements fracture aperture must an attribute stored
                in the variable called 'thickness'

    @attention  If there are bumps in the boundaries, then the calculation will not be correct; this is reported.
    
    @test SKM 29/5/14: in progress
    
    TODO: reduce the number of output files to a manageable amount

*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::DiagonalTensorAnalysis( const char* model_name, bool with_vtk_output )
{ 
  ErrorHandler&    csmp_error( ErrorHandler::Instance() );
  const double&  model_time( ModelTime::Instance().modelTime );

// --------------------------------------------------------------------------------------
// 1. Combine configuration-file defined regions and get some diagnostics
// --------------------------------------------------------------------------------------
  fracture_surface_area_ = IdentifyFractureAndMatrixRegions();
  
  cout <<"\nDiagonalTensorAnalysis: model_ has a:\n";
  cout <<"\tfracture-matrix interface area, Af, of: "<< (fracture_surface_area_ * 2.) <<" (m2)\n";
  cout <<"\tP32 of: "<< (fracture_surface_area_*2.) / model_volume_ <<"\n";
  cout <<"\tfraction of void space due to fractures: ";
  fracture_matrix_void_ratio_ = FractureMatrixVoidRatio( false );
  cout <<  fracture_matrix_void_ratio_ <<"\n\n";
  // print the fracture aperture range
  printRangeOfVariable( model_, stdio_, "fractures", "thickness" );


// ------------------------------------------------------------------------------------
// 2. Calculating transmissivity (= hydraulic conductivity x thickness product)
//    (thickness attribute is used to scale lower dimensional fractures)
// ------------------------------------------------------------------------------------
  ComputeTransmissivity();
  printRangeOfVariable( model_, stdio_, "transmissivity" );


// ------------------------------------------------------------------------------------
// 3. Loop over the three coordinate directions (0=x, 1=y, 2=z), calculation pressure
// ------------------------------------------------------------------------------------
  Region<dim>&    model_domain(model_.Region("Model"));
  const double    pmin(100325.); // (Pa), atmospheric pressure
  string          monitor_file(model_name);
  monitor_file += "-monitored-regions";

  VTK_Interface<dim>  vtk_output;

  // x, y, and z direction
  for ( uint32_t i=0; i<dim; i++ )
    {
      // releasing any pre-assigned boundary conditions
      model_.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
      model_domain.ChangePropertyStatus( "fluid pressure", PLAIN, COMPLETE );
      // calculate fluid pressure in x direction
      if ( i == 0 ) {
           cout << "\n\n\nmain: Computing 'fluid pressure' and 'velocity' in x-direction" << endl;
           const double pmax = pmin + farfield_pf_gradient_ * model_dimensions_[0];
           xsect_area_ = ( dim == 3 ) ? fabs(model_dimensions_[1] * model_dimensions_[2]) : fabs(model_dimensions_[1]);
           model_.Boundary("LEFT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmax) );
           model_.Boundary("RIGHT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmin) );
        }
      // calculate fluid pressure in y direction
      else if ( i == 1 ) {
           cout << "\n\n\nmain: Computing 'fluid pressure' and 'velocity' in y-direction" << endl;
           const double pmax = pmin + farfield_pf_gradient_ * model_dimensions_[1];
           xsect_area_ = ( dim == 3 ) ? fabs(model_dimensions_[0] * model_dimensions_[2]) : fabs(model_dimensions_[0]);
           model_.Boundary("BOTTOM").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmax) );
           model_.Boundary("TOP").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmin) );
        }
      // calculate fluid pressure in z direction (will not be called if 2D)
      else if ( i == 2 ) {
           cout << "\n\n\nmain: Computing 'fluid pressure' and 'velocity' in z-direction" << endl;
           const double pmax = pmin + farfield_pf_gradient_ * model_dimensions_[2];
           xsect_area_ = fabs(model_dimensions_[0] * model_dimensions_[1]); // XY plane always exists
           model_.Boundary("FRONT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmax) );
           model_.Boundary("BACK").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, pmin) );
        }

      // fluid pressure and velocity calculation
      SteadyStateDiffusor<dim,Element> fluid_pressure( model_, "transmissivity", "fluid pressure", "fluid volume source" );
      fluid_pressure.IntegrateOver( model_, model_domain );
      ComputeVelocityAndVolumeFlux( "Model" );
      printRangeOfVariable( model_, stdio_, "fluid pressure" );
      printRangeOfVariable( model_, stdio_, "flow speed" );

      // VTK output of pressure, velocity and volume flux
      if ( with_vtk_output && stdio_.YesNo("Do you want to output results to VTK") )
        WriteResultsToVTK( model_name, i );
      
      // monitor velocity spectra and region properties
      OutputVelocityHistogramToMaple(i);
      monitor_.ScalarPropertyIntegrals( model_, model_time );
      monitor_.ScalarPropertyRanges( model_, model_time );
      if ( i == 0 )      monitor_.Out( (monitor_file + "-x").c_str() );
      else if ( i == 1 ) monitor_.Out( (monitor_file + "-y").c_str() );
      else if ( i == 2 ) monitor_.Out( (monitor_file + "-z").c_str() );
       
       
    // ------------------------------------------------------------------------------------
    // 4. calculate the relative proportions of the fracture and the matrix flux and
    //    equivalent permeability computed from fluid throughput through model as determined by FVM
    // ------------------------------------------------------------------------------------
      double total_flux{0.}, matrix_flux{0.};
      if      ( i == 0 ) k_equivalent_[i] = EquivalentPermeability( LEFT,   RIGHT, total_flux );
      else if ( i == 1 ) k_equivalent_[i] = EquivalentPermeability( BOTTOM, TOP,   total_flux );
      else if ( i == 2 ) k_equivalent_[i] = EquivalentPermeability( BACK,   FRONT, total_flux );

      // computing matrix-only permeability (to cover the case where the matrix is heterogeneous)
      // evaluating matrix permeability
      Region<dim>&  mref(model_.Region("matrix"));
      double  matrix_k1, matrix_k2;
      mref.MinMaxOf( "permeability", matrix_k1, matrix_k2 );
      // if matrix permeability is not single-valued, a computation is needed
      if ( fabs(matrix_k1 - matrix_k2) > numeric_limits<double>::epsilon() ) {
          VectorVariable<dim> zero_velo; zero_velo=0.;
          model_.InputPropertyValue( "velocity", zero_velo );
          SteadyStateDiffusor<dim,Element> matrix_pressure( model_, "transmissivity", "fluid pressure", "fluid volume source" );
          matrix_pressure.IntegrateOver( mref );
          ComputeVelocityAndVolumeFlux( "matrix" );
          printRangeOfVariable( model_, stdio_, "fluid pressure" );
          printRangeOfVariable( model_, stdio_, "velocity" );
          matrix_permeability_[i] = EquivalentPermeability( "matrix", matrix_flux );
        }
      else {
           matrix_flux = xsect_area_ * (matrix_k1 / fluid_viscosity_) * farfield_pf_gradient_;
           matrix_permeability_[i] = matrix_k1;
        }
      if ( matrix_flux > total_flux ) {
           cerr <<"\n\ttotal flux is smaller than calculated matrix flux: "<< total_flux <<" vs. "<< matrix_flux <<" m3/s.\n";
           k_equivalent_[i] = matrix_permeability_[i];
           total_flux       = matrix_flux;
           csmp_error.Note( ERROR, "DFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis:",
                                     "flux computation problem in 'Model' domain; check boundaries in mesh.");
        }

      // fracture-matrix flux ratio
      flux_ratio_[i] = (total_flux - matrix_flux) / matrix_flux;
    
  } // end of loop over three coordinate directions
  
  
  // display effective properties to screen
  PrintResultsToScreen( model_name );
  // and write them to a file
  string effective_properties_file_name(model_name);
  effective_properties_file_name += "-equivalent-properties.text";
  WriteResultsToFile( effective_properties_file_name.c_str() );
  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::DiagonalTensorAnalysis: That's it..."<< endl;
  
} // end DiagonalTensorAnalysis





/**
    fracture - matrix equivalent permeability analysis by flow-based upscaling.
    
    Computes full tensor permeability using approach of Durlofsky (200x).

*/
template<uint32_t dim>
void DFM_PermeabilityAndFluxRatioAnalyzer<dim>::FullTensorAnalysis( const char* model_name )
{ 
// --------------------------------------------------------------------------------------
// 1. Combine configuration-file defined regions and get some diagnostics
// --------------------------------------------------------------------------------------
  fracture_surface_area_ = IdentifyFractureAndMatrixRegions();
  
  cout <<"\nFullTensorAnalysis: model_ has a:\n";
  cout <<"\tfracture-matrix interface area, Af, of: "<< (fracture_surface_area_ * 2.) <<" (m2)\n";
  cout <<"\tP32 of: "<< (fracture_surface_area_*2.) / model_volume_ <<"\n";
  cout <<"\tfraction of void space due to fractures: ";
  fracture_matrix_void_ratio_ = FractureMatrixVoidRatio( false );
  cout <<  fracture_matrix_void_ratio_ <<"\n\n";


// ------------------------------------------------------------------------------------
// 2. Calculating transmissivity (= hydraulic conductivity x thickness product)
//    (thickness attribute is used to scale lower dimensional fractures)
// ------------------------------------------------------------------------------------
  ComputeTransmissivity();
  printRangeOfVariable( model_, stdio_, "transmissivity" );


// ------------------------------------------------------------------------------------
// 3. SIROOS MAIN PROGRAM STARTS HERE
// ------------------------------------------------------------------------------------
   EquivalentPermeabilityTensor<dim>  tensor_k_calculator(model_);
 
// ------------------------------------------------------------------------------------
//              equivalent permeability tensor: whole model
// ------------------------------------------------------------------------------------
    // assuming that the permeability of the model is single-valued
    const double epsilon = 1e-1 * printRangeOfVariable( model_, stdio_, "permeability" );
    tensor_k_calculator.PermeabilityTolerance(epsilon);
    tensor_k_calculator.ComputeEquivalentPermeabilityTensor( model_, "Model" );

    Region<dim>& mref(model_.Region("Model"));
  
    ofstream file1( (string(model_name)+"PermTensor.txt").c_str() );
    file1 << "====================================================================================================================== \n";
    file1 << setw(50) << "Pemeability Tensor" << setw(57) << "Principal Components" << endl;
    file1 << setw(8) << "No." << setw(16) << "kxx" << setw(21) << "kxy" << setw(20) << "kxz" << endl;
    file1 << setw(10) << "Elements" << setw(14) << "kyx" << setw(21) << "kyy" << setw(20) << "kyz" << setw(24) << "kmin" << setw(20) << "kmax" << endl;
    file1 << setw(24) << "kzx" << setw(21) << "kzy" << setw(20) << "kzz" << setw(20) << "Trend" << setw(10) << "Plunge" << setw(10) << "Trend" << setw(10) << "Plunge" << endl;
    file1 << "====================================================================================================================== \n";
    // writing the results in output file
    // 1st line
    file1 << setw(28) << setprecision(3) << tensor_k_calculator(0,0);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(0,1);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(0,2) << endl;
    // 2nd line
    file1 << setw(8)  << mref.Cells();
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(1,0);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(1,1);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(1,2);
    // 3rd line
    file1 << setw(28) << setprecision(3) << tensor_k_calculator(2,0);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(2,1);
    file1 << setw(20) << setprecision(3) << tensor_k_calculator(2,2);
    file1 << "---------------------------------------------------------------------------------------------------------------------- \n";
  
  cout <<"\nDFM_PermeabilityAndFluxRatioAnalyzer::FullTensorAnalysis: That's it..."<< endl;
  
} // end FullTensorAnalysis






/**
    Highlighting the perimeter nodes of the model in VTK output file.

    Testing
 
    - is the flux sampled at the correct boundary nodes ?
*/
template<uint32_t dim>
bool DFM_PermeabilityAndFluxRatioAnalyzer<dim>::TestBoundaryIntegrity()
 {
    const Region<dim>&  model_domain(model_.Region("Model"));
    model_domain.RenumberNodes();
    // go over all boundaries collect the nodes found there and then compare
    // these with the perimeter nodes recovered from the Model.
    // If these are less an error is reported and the nodes are output to VTK
    // (NOTE: this model must only have the external boundaries as boundaries)
    // test model aperture_20_frac has 1102 nodes on box boundary
    set<size_t>  actual_bnodes;
    for ( auto bit=model_.BoundariesBegin(); bit!=model_.BoundariesEnd(); ++bit )
      for ( auto it=(*bit).second.CellsBegin(); it!=(*bit).second.CellsEnd(); ++it )
        for ( uint32_t i=0U; i<(*it)->Nodes(); ++i )
          actual_bnodes.insert( (*it)->N(i)->Idx() );

    set<size_t>  model_bnodes;
    for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      model_bnodes.insert( (*nit)->Idx() );

    // from box boundary (not included as boundary flags will be deprecated)
    set<size_t>  box_bnodes;
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      if ( (*nit)->AtBoundary() != NOT )
        box_bnodes.insert( (*nit)->Idx() );
   
    if ( actual_bnodes != model_bnodes ) {
          // getting the node idx of the nodes that were not found
          vector<size_t> missing_nodes;
          set_symmetric_difference( actual_bnodes.begin(), actual_bnodes.end(),
                                    model_bnodes.begin(), model_bnodes.end(),
                                    back_inserter( missing_nodes ) );
      
          cerr <<"\nDFM_PermeabilityAndFluxRatioAnalyzer<dim>::TestBoundaryIntegrity:";
          cerr <<"\n\tnodes derived from boundaries: "<< actual_bnodes.size();
          cerr <<"\n\tnodes from perimeter-flagging: "<< model_bnodes.size();
          cerr <<"\n\tnodes from box-boundary flags: "<< box_bnodes.size();
          cerr <<"\n\nmismatching nodes:";
          out( missing_nodes );
    
          const csmp::Index p_key = model_.Database().StorageKey("fluid pressure");
          // are the boundary nodes where expected ?
          model_.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
          const Region<dim>& modeldomain(model_.Region("Model"));
          for ( auto nit=modeldomain.PerimeterNodesBegin(); nit!=modeldomain.NodesEnd(); ++nit )
            (*nit)->Store( p_key, makeScalar(INIT_COND,5.0e5) );

          VTK_Interface<dim>  vtk_output;
          vtk_output.OutputDataToVTK( model_, "pf", "fluid pressure", 0 );
    
          return false;
       }
   
    return true;
   
 } // end TestBoundaryIntegrity




// explicit instantiation of the template
template class DFM_PermeabilityAndFluxRatioAnalyzer<2U>;
template class DFM_PermeabilityAndFluxRatioAnalyzer<3U>;

} // end csmp
