// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CSMP_definitions.h"

// the model
#include "Region.h"
#include "Boundary.h"
#include "Model.h"

#include "SAMG_Settings.h"
#include "SAMG_Solver.h"

// PDE integration
#include "SteadyStateDiffusor.h"
#include "NumIntegral_NT_op_N_dV.h"

// Interrelations
#include "ConstantFactor.h"

#include "Standard_IO_Handler.h"
#include "CSMP_highLevelUtilities.h"
#include "ModelTime.h"

// Analysis
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "StatisticalAnalyzer.h"
#include "RegionMonitor.h"
// to analyze boundary fluxes in this linear model
#include "StencilProcessor.h"
#include "NodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {
  void computeApertureDistribution( Model<3U>&, const char* region_name );
	                          
	void  removeAuxiliaryRegions( Model<3U>&, set<string>& extra_groups );

  template<size_t dim>
	double fractureMatrixFluxRatio( const Model<dim>& sg,
	                                   const char* matrix_region, const char* fracture_region, 
	                                   double fluid_viscosity, 
	                                   double grad_pf, double xsect_area, double total_flux );

	double equivalentMatrixPermeability( const Model<3>& sg,
                                       const char* matrix_region, const char* fracture_region,
                                       double fluid_viscosity,
                                       double grad_pf, double xsect_area, double total_flux );

  bool detectOutOfPlaneNodesAtBoxBoudaries( const Model<3U>& model, const char* model_name, double relative_tolerance );

static void printSpecifications()
{
cout <<"Determination of equivalent permeability through flow-based upscaling and\
calculation of fracture - matrix flux ratio (Nelson's 'excess permeability' of NFRs)\
\n\nInput:  box-shaped CSMP model with target fracture apertures expressed in terms of the 'thickness'\
of lower-dimensional fracture elements (surfaces in this 3D model):\
\
\n\n - fractures / fracture sets must be identified by the string 'FRAC' or 'SET' in their names\
\n\n - rock matrix must contain 'MATRIX' in their names\
\n\n - principal directions of permeability should be aligned with the coordinate axes of the model\
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
}


/**
     Computes the absolute volumes of the supplied domain using the thickness attribute
     which ascertains that elements in the domain with different dimensionality are counted correctly.
*/
template<uint32_t dim>
double regionVolume( const Model<dim>& model, const char* region )
 {
    const Region<dim>&  ref(model.Region(region));
    const csmp::Index   thi_key = model.Database().StorageKey("thickness");
    double volume(0.);
    for ( auto eit=ref.CellsBegin(); eit!=ref.CellsEnd(); ++eit ) {
         volume += (*eit)->Volume() * (*eit)->Read( thi_key );
      }
   
    return volume;

} // end regionVolume

template double regionVolume( const Model<2U>&, const char* );
template double regionVolume( const Model<3U>&, const char* );






/**
   If the side boundaries of the model have defects, the permeability calculation becomes inaccurate.
   This methods detect this and reports boundaries and locations of spurious nodes.
   
   Assumption: the model is box-shaped and only has such boundaries.
   
   @test O.K. SKM 1/11/2013
*/
bool detectOutOfPlaneNodesAtBoxBoudaries( const Model<3U>& model, const char* model_name, double relative_tolerance )
 {
    bool with_boudary_defects(false);
    VTU_Interface<3U>  vtu_output( model, "box_boundary_defects" );
   
    // For each boundary of the box-shaped model, find the bounding box which should be a plane.
    // If there the smallest dimension is larger than 0.1% of the largest one the boundary will be output
    // and an error will be reported.
    set<Point<3U> > bpoints;
    for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); bit++ )
      {
         for ( auto it=(*bit).second.NodesBegin(); it!=(*bit).second.NodesEnd(); ++it )
           bpoints.insert( (*it)->Coordinate() );
        
         // checking points in the set
         Point<3U> bdiag = (*bpoints.begin()) - (*bpoints.rbegin());
         // reporting faulty boundary to VTU file
         const double box_max = std::max( std::max( fabs(bdiag[0]), fabs(bdiag[1])), fabs(bdiag[2]) );
         const double max_undulation = box_max * relative_tolerance;
         if ( fabs(bdiag[0]) > max_undulation && fabs(bdiag[1]) > max_undulation && fabs(bdiag[2]) > max_undulation ) {
              with_boudary_defects = true;
              cout <<"\ndetectOutOfPlaneNodesAtBoxBoudaries: Problem: bounding box dimensions of boundary: ";
              cout << (*bit).first <<"\n";
              cout <<"\t"<< fabs(bdiag[0]) <<"  "<< fabs(bdiag[1]) <<"  "<< fabs(bdiag[2]) << endl;
              vtu_output.OutputDataToVTU( model_name, "fluid pressure", (*bit).second, static_cast<size_t>(0) );
           }
         bpoints.clear();
      }
   
   return with_boudary_defects;

 } // end detectOutOfPlaneNodesAtBoxBoudaries




/**
    Computes the fraction of the total pore space that is occupied by fractures.
*/
template<uint32_t dim>
double fractureMatrixVoidRatio( const Model<dim>& sg, const char* matrix_region, const char* fracture_region, bool verbose=false )
 {
    const Region<dim>&  fgref(sg.Region(fracture_region)); 
    const Region<dim>&  mgref(sg.Region(matrix_region));
    
    const csmp::Index phi_key = sg.Database().StorageKey("porosity");
    double matrix_pore_volume(0.);
    for ( auto eit=mgref.CellsBegin(); eit!=mgref.CellsEnd(); ++eit ) {
         assert( (*eit)->IsVolume() );
         matrix_pore_volume +=  (*eit)->Volume() * (*eit)->Read( phi_key );
      }
      
    // fracture pore volume (of the lower-dimensional fractures) that is substracted from matrix volume
    const csmp::Index thi_key = sg.Database().StorageKey("thickness");
    double fracture_pore_volume(0.);
    for ( auto eit=fgref.CellsBegin(); eit!=fgref.CellsEnd(); ++eit ) {
         assert( (*eit)->IsSurface() );
         fracture_pore_volume += (*eit)->Volume() * (*eit)->Read( thi_key) * (*eit)->Read( phi_key );
      }

    const double fracture_fraction_of_void_space = fracture_pore_volume / (matrix_pore_volume - fracture_pore_volume);
    if ( verbose ) {
         cout <<"\nfractureMatrixVoidRatio: fraction of the total void space that is due to fractures: ";
         cout << fracture_fraction_of_void_space <<"\n\n";
      }
   
    return fracture_fraction_of_void_space;

} // end fractureMatrixVoidRatio

template double fractureMatrixVoidRatio( const Model<2U>&, const char*, const char*, bool );
template double fractureMatrixVoidRatio( const Model<3U>&, const char*, const char*, bool );


 
/** Determination of qf/qm ratio

    calculate how much of the total flow is in the matrix and how much is in the matrix, assuming that the
    matrix K has a single value.
    
    @attention it is assumed that the fractures make a negligible contribution to the flow cross-sectional area
    
    @attention this method only works for a single matrix domain with an isotropic permeability.
    
    @test refactored by SKM 2/11/2013: introduced thickness attribute to capture fracture aperture.
*/
template<uint32_t dim>
double fractureMatrixFluxRatio( const Model<dim>& sg,
                                  const char* matrix_region,
                                  double fluid_viscosity, 
                                  double grad_pf, double xsect_area, double total_flux )
 {
    const Region<dim>&  mgref(sg.Region(matrix_region));
    
    // finding the matrix permeability
    double   matrix_k1, matrix_k2; 
    mgref.MinMaxOf( "permeability", matrix_k1, matrix_k2 );
    if ( matrix_k1 != matrix_k2 ) {
         throw csmp::Exception( ERROR, "fractureMatrixFluxRatio", 
                                       "calculation not possible, 'rock matrix' has a range of permeability values." );
         return std::strtod("NAN",NULL);
      }
 
    // Contribution of matrix to total flux for the reduced flow cross-section and the given far-field
    // fluid pressure gradient
    const double matrix_flux =  xsect_area * (matrix_k1 / fluid_viscosity) * fabs(grad_pf);

    // now compute the fracture flux matrix flux ratio (if this is the case, there is no excess permeability)
    if ( total_flux - matrix_flux <= 0. ) return 0.;
    return (total_flux - matrix_flux) / matrix_flux;  

} // end fractureMatrixFluxRatio                           

template
double fractureMatrixFluxRatio( const Model<2U>&, const char*, double, double, double, double );
template
double fractureMatrixFluxRatio( const Model<3U>&, const char*, double, double, double, double );





/**
     Compute matrix permeability that the matrix would need to have to carry the same amount of flow
     as the fractures.
     
     @attention Method works only for a homogeneous isotropic matrix domain.
*/
double equivalentMatrixPermeability( const Model<3>& sg,
                                     const char* matrix_region, const char* fracture_region,
                                     double fluid_viscosity, 
                                     double grad_pf, double xsect_area, double total_flux )
 {
    // evaluating matrix permeability
    const Region<3>&  mgref(sg.Region(matrix_region));
    double  matrix_k1, matrix_k2;
    mgref.MinMaxOf( "permeability", matrix_k1, matrix_k2 );
    if ( matrix_k1 != matrix_k2 ) {
         throw csmp::Exception( ERROR, "equivalentMatrixPermeability", 
                               "calculation not possible, 'rock matrix' has a range of permeability values." );
         return std::strtod("NAN",NULL);
      }
 
    const double  fracture_volume = regionVolume(sg,fracture_region);
    // volume of lower-dimensional fractures must be subtracted because they occupy the same space as rock matrix
    const double  matrix_volume   = regionVolume(sg,matrix_region) - fracture_volume;
   
    double flux_multiplier = matrix_volume / (matrix_volume + fracture_volume);
    
    // Contribution of matrix to total flux for the reduced flow cross-section and the given far-field
    // fluid pressure gradient
    const double matrix_flux = flux_multiplier * xsect_area * (matrix_k1/fluid_viscosity) * fabs(grad_pf);
    if ( total_flux - matrix_flux <= 0. ) return matrix_k1;
   
    const double matrix_k = ((total_flux - matrix_flux) / (fabs(grad_pf) * flux_multiplier * xsect_area)) * fluid_viscosity;

    return matrix_k;

} // end equivalentMatrixPermeability                           





/// taking into account the thickness attribute correctly
static void computeVelocityAndVolumeFlux( Model<3U>& model, const char* region )
{
   Region<3U>&  ref(model.Region(region));
   const csmp::Index pf_key  = model.Database().StorageKey("fluid pressure");
   const csmp::Index K_key   = model.Database().StorageKey("hydraulic conductivity");
   const csmp::Index phi_key = model.Database().StorageKey("porosity");
   const csmp::Index vD_key  = model.Database().StorageKey("velocity");
   const csmp::Index vi_key  = model.Database().StorageKey("pore velocity");
   const csmp::Index vf_key  = model.Database().StorageKey("volume flux");
   VectorVariable<3U>   velo;
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
              velo(2)  += pf  * -DERIV(2,i) * K;
           }
        
         // storing the computed velocity (checking the variable status flags)
         (*it)->Store( vD_key, velo );
        
         // storing the computed velocity (checking the variable status flags)
         const double volume_flux = velo.Length();
         (*it)->Store( vf_key, makeScalar(PLAIN,volume_flux) );
        
         // pore velocity flux
         velo /= (*it)->Read( phi_key );
         (*it)->Store( vi_key, velo );
      }

} // end computeVelocityAndVolumeFlux




/**
    fracture - matrix equivalent permeability analysis
    
    @attention: regions must consist out of discontiguous fracture subregions that
    do not include multiple sets.
    
    @attention: for lower-dimensional elements fracture aperture must an attribure stored
                in the variable called 'thickness'

    @attention  If there are bumps in the boundaries, then the calculation will not be correct; this is reported.

*/
static int k_equiv_qfqm_Analysis3D( const char* model_name, bool with_vtk_output )
{ 
  double&  model_time( ModelTime::Instance().modelTime );
  const bool verbose(false);

  cout <<"\n\nCSMP++ FRACTURE - MATRIX EQUIVALENT PERMEABILITY & FRACTURE-MATRIX FLUX RATIO ANALYZER"<< endl;
  cout <<"copyright (c) Stephan K. Matthai (2013)."<< endl << endl;
  printSpecifications();

  Standard_IO_Handler    stdio( model_name );
  csmp::Model<3U>        model(string{model_name});
  Region<3U>&            sgref(model.Region("Model"));
  const double         model_volume(sgref.Volume()),
                         pore_volume(sgref.Volume(true));
  printModelDimensions( model );
  cout <<"\nk_equiv_qfqm_Analysis3D: model volume: "<< model_volume <<" and pore volume: "<< pore_volume << endl;
  cout <<"\nk_equiv_qfqm_Analysis3D: model contains the unique regions: ";
  for ( Model<3>::regionConstIterator rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); rit++ )
    cout << (*rit).first <<"  ";
  cout <<"\n\n";

  // screening for surface bumps greater than 1% of the model dimensions
  if ( detectOutOfPlaneNodesAtBoxBoudaries( model, model_name, 0.01 ) )
    throw csmp::Exception( FATAL_ERROR, "k_equiv_qfqm_Analysis3D:", "model boundaries contain irregularities; check model before continueing." );

// ------------------------------------------------------------------------------------
// 1. Establishing initial & essential conditions for testing
// ------------------------------------------------------------------------------------
  // fluid volume source
  model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );

// hydrostatic fluid pressure gradient (g x rho_water) 
  const double farfield_pf_gradient(9806.65);
  
  list<string> integral_properties;  integral_properties.push_back("volume flux");
  list<string> range_properties;     range_properties.push_back("fluid pressure");

  RegionMonitor<3U>  monitor( model, integral_properties, range_properties );

// --------------------------------------------------------------------------------------
// 2. Combine configuration-file defined fracture sets into the region 'fractures'
// --------------------------------------------------------------------------------------
// tested:
  bool         with_matrix(false), with_fractures(false);
  set<string>  fracture_regions, matrix_regions;
  double     fracture_matrix_interface_area(0.);
  // looping over all the unique regions that contain fractures and or matrix and combining thee
  for ( Model<3>::regionConstIterator rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); rit++ )
    {
       // if the region name contains 'FRACTURE' it is pooled into an non-unique region 'fractures'
       if ( (*rit).first.find("FRAC") != std::string::npos || (*rit).first.find("SET") != std::string::npos) {
            // recording the region name
            fracture_regions.insert((*rit).first);
            fracture_matrix_interface_area += (*rit).second.Volume();
            if ( with_fractures == false ) {
                 model.CopyRegion( (*rit).first.c_str(), "fractures" );
                 with_fractures = true;
              }
            else model.AssimilateRegion( (*rit).first.c_str(), "fractures" );
         }
       // same for matrix
       if ( (*rit).first.find("MATRIX") != std::string::npos ) {
            // recording the region name
            matrix_regions.insert((*rit).first);
            if ( with_matrix == false ) {
                 model.CopyRegion( (*rit).first.c_str(), "matrix" );
                 with_matrix = true;
              }
            else model.AssimilateRegion( (*rit).first.c_str(), "matrix" );
         }
    }
  
  if ( fracture_regions.empty() )
    throw csmp::Exception( FATAL_ERROR, "k_equiv_qfqm_Analysis3D", "there is no region tagged 'FRACTURE' or 'SET' in the model; quitting.");

  if ( matrix_regions.empty() )
    throw csmp::Exception( FATAL_ERROR, "k_equiv_qfqm_Analysis3D", "there is no region tagged 'MATRIX' in the model; quitting.");
  
  if ( matrix_regions.size() > 1 ) {
       cout <<"\nk_equiv_qfqm_Analysis3D: Model contains multiple matrix regions: ";
       for ( set<string>::const_iterator it=matrix_regions.begin(); it!=matrix_regions.end(); it++ )
         cout << (*it) <<" ";
       cout <<"\n\n";
       throw csmp::Exception( ERROR, "k_equiv_qfqm_Analysis3D", "unfortunately, this tool can only handle single-valued matrix permeability");
    }
  
  cout <<"\n\nk_equiv_qfqm_Analysis3D: Model has a fracture-matrix interface area, Af of: ";
  cout << fracture_matrix_interface_area * 2. <<" (m2), and has a P32 of: "<< fracture_matrix_interface_area / model_volume <<"\n\n";
  
  cout <<"\n\nk_equiv_qfqm_Analysis3D: fraction of void space due to fractures: ";
  cout << fractureMatrixVoidRatio( model, "matrix", "fractures" ) <<"\n\n";



// ------------------------------------------------------------------------------------
// 3. Calculating hydraulic conductivity from permeability
//    (thickness attribute is used to scale lower dimensional fractures)
// ------------------------------------------------------------------------------------
  VTK_Interface<3U>  vtk_output;
  if ( with_vtk_output ) vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
// tested:
  const csmp::Index k_key   = model.Database().StorageKey("permeability");
  const csmp::Index K_key   = model.Database().StorageKey("hydraulic conductivity");
  const csmp::Index thi_key = model.Database().StorageKey("thickness");
  printRangeOfVariable( model, stdio, "permeability" );
  printRangeOfVariable( model, stdio, "thickness" );
  cout <<"\nNote: 'thickness' here refers to 'fracture aperture'\n\n";
  const double  fluid_viscosity(1.6e-3); // Pa s-1
  for ( auto it=sgref.CellsBegin(); it!=sgref.CellsEnd(); it++ ) {
        const double k = ((*it)->Read(k_key) * (*it)->Read(thi_key)) / fluid_viscosity;
       (*it)->Store( K_key, makeScalar(PLAIN,k) );
    }
  printRangeOfVariable( model, stdio, "hydraulic conductivity" );




// ------------------------------------------------------------------------------------
// 4. Steady-state fluid pressure algorithm [K]{p} = {Q}
// ------------------------------------------------------------------------------------
  SAMG_Settings  settings;

  SteadyStateDiffusor<3U,Element>  fluid_pressure( model, "hydraulic conductivity", "fluid pressure", "fluid volume source" );

  // targeting SAMG DLL 1 for this pressure solver
  fluid_pressure.GetSolverSettings().SetSolverInstance(1);
  // iout
  fluid_pressure.GetSolverSettings().Set_iout1( 0 );
  fluid_pressure.GetSolverSettings().Set_iout2( 0 );
  if ( !verbose ) fluid_pressure.GetSolverSettings().Set_idmp( -1 );
  // one-time solver set-up: nothing is remembered for next try
  fluid_pressure.GetSolverSettings().Set_iswit(5); 
  // SAMG solution criteria
  fluid_pressure.GetSolverSettings().Set_eps(0.); // absolute criterion
  // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
  fluid_pressure.GetSolverSettings().Set_a_cmplx(2);
  // Pre-adjust SAMG mesh complexity, based on solver output
  fluid_pressure.GetSolverSettings().Set_g_cmplx(1.5);
  fluid_pressure.GetSolverSettings().Set_w_avrge(2);

#ifdef SAMG_OUTPUT_TO_FILE
  // trigger output to file
  fluid_pressure.GetSolverSettings().Set_idmp( 8 );
  // // define SAMG file output format for reduced file size, idmp > 1 is required
  fluid_pressure.GetSolverSettings().Set_ioform( "f" );
  // set filename for SAMG file output other than default "level", idmp > 1 is required
  fluid_pressure.GetSolverSettings().Set_filnam_dump( "ReservoirSimulator_steady_state_p" );
#endif



// ------------------------------------------------------------------------------------
// 5. Analysis of results (analyser is used in loop below)
// ------------------------------------------------------------------------------------
  StatisticalAnalyzer<3U>                                     flux_histogram( model );
  vector<pair<double,double> >                            bins;
  map<string,pair<vector<pair<double,double> >,size_t> >  results;
  flux_histogram.DefineBins( "flux_histogram.bins", bins );
            
  

// ------------------------------------------------------------------------------------
// 6. Loop over the three coordinate directions (0=x, 1=y, 2=z)
// ------------------------------------------------------------------------------------
  vector<double>  k_equivalent(3), flux_ratio(3), matrix_permeability(3);
  const double  pmin(0.);
  double        model_throughput;
  string    monitor_file(model_name);
  monitor_file += "-monitored-regions";

  // single-phase solute advection-ONLY constructor, 1st-order method
  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", model,
                                                   "porosity", "fluid pressure", "velocity",
                                                   "nodal fluid volume source", false, false, "thickness" );
                
  // cross-sectional areas normal to the direction of flow
  Point<3U>  min_xyx, max_xyz;
  sgref.MinMaxCoordinates( min_xyx, max_xyz );
  double xsec_area;
 
  // x, y, and z direction
  for ( unsigned int i=0U; i<3U; i++ ) {
    // releasing any pre-assigned boundary conditions
    sgref.ChangePropertyStatus( "fluid pressure", PLAIN, COMPLETE );
	  // calculate fluid pressure in x direction
	  if ( i == 0 ) {
	       cout << "\n\n\nk_equiv_qfqm_Analysis3D: Computing 'fluid pressure' and 'velocity' in x-direction" << endl;
	       const double pmax = farfield_pf_gradient * fabs(max_xyz[0] - min_xyx[0]);
	       xsec_area = (max_xyz[1] - min_xyx[1]) * (max_xyz[2] - min_xyx[2]);
	       model.InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH, pmax) );
	       model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH, pmin) );
	    }
	  // calculate fluid pressure in y direction
	  else if ( i == 1 ) {
	       cout << "\n\n\nk_equiv_qfqm_Analysis3D: Computing 'fluid pressure' and 'velocity' in y-direction" << endl;
	       const double pmax = farfield_pf_gradient * fabs(max_xyz[1] - min_xyx[1]);
	       xsec_area = (max_xyz[0] - min_xyx[0]) * (max_xyz[2] - min_xyx[2]);
	       model.InputBoundaryValue( BOTTOM,  "fluid pressure", makeScalar(DIRICH, pmax) );
	       model.InputBoundaryValue( TOP,     "fluid pressure", makeScalar(DIRICH, pmin) );
	    }
	  // calculate fluid pressure in z direction
	  else if ( i == 2 ) {
	       cout << "\n\n\nk_equiv_qfqm_Analysis3D: Computing 'fluid pressure' and 'velocity' in z-direction" << endl;
	       const double pmax = farfield_pf_gradient * fabs(max_xyz[2] - min_xyx[2]);
	       xsec_area = (max_xyz[0] - min_xyx[0]) * (max_xyz[1] - min_xyx[1]);
	       model.InputBoundaryValue( FRONT,  "fluid pressure", makeScalar(DIRICH, pmax) );
	       model.InputBoundaryValue( BACK,   "fluid pressure", makeScalar(DIRICH, pmin) );
	    }

	  // fluid pressure and velocity calculation
	  model.Apply( fluid_pressure );
    computeVelocityAndVolumeFlux( model, "Model" );
	  printRangeOfVariable( model, stdio, "fluid pressure" );
	  printRangeOfVariable( model, stdio, "velocity" );

      if ( with_vtk_output ) {
           // dependent variable ranges in groups
           for ( map<string,Region<3U> >::const_iterator 
                 grit=model.RegionsBegin(); grit!=model.RegionsEnd(); grit++ ) {
                  printRangeOfVariable( model, stdio,(*grit).first.c_str(),  "fluid pressure" );
                  printRangeOfVariable( model, stdio, (*grit).first.c_str(), "velocity" );
                  printRangeOfVariable( model, stdio, (*grit).first.c_str(), "volume flux" );
                // output results for visualization using VTK
                vtk_output.OutputDataToVTK( model, (*grit).first.c_str(), "fluid-pressure", "fluid pressure", static_cast<size_t>(i), true );
                vtk_output.OutputDataToVTK( model, (*grit).first.c_str(), "velocity",  "velocity", static_cast<size_t>(i), true );
                vtk_output.OutputDataToVTK( model, (*grit).first.c_str(), "volume-flux", "volume flux", static_cast<size_t>(i), true );
             }
           // output results for visualization using VTK
           vtk_output.OutputDataToVTK( model, "permeability",   "permeability",    static_cast<size_t>(i), true );
           vtk_output.OutputDataToVTK( model, "velocity",       "velocity",        static_cast<size_t>(i), true );
           vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure",  static_cast<size_t>(i), true );
        }
        
	  // monitor the group properties
	  monitor.ScalarPropertyIntegrals( model, model_time );
	  monitor.ScalarPropertyRanges( model, model_time );
	  if ( i == 0 )      monitor.Out( (monitor_file + "-x").c_str() );
	  else if ( i == 1 ) monitor.Out( (monitor_file + "-y").c_str() );
	  else               monitor.Out( (monitor_file + "-z").c_str() );
	  
	  // effective permeability
	  cout <<"\nk_equiv_qfqm_Analysis3D: total and specific cross-sectional flux through model: "<< (model_throughput = fabs(advector.ModelInflow()));
	  cout <<", "<< model_throughput / xsec_area << endl;
	  k_equivalent[i]  = model_throughput / farfield_pf_gradient;
	  k_equivalent[i] /= xsec_area;
	  k_equivalent[i] *= fluid_viscosity;
	  
	  flux_histogram.RegionPropertyHistograms( "volume flux", bins, results );
	  if      ( i == 0 ) flux_histogram.OutputRegionPropertyAbundancePolygonsMaple( "volume-flux-x", bins, results, true ); 
	  else if ( i == 1 ) flux_histogram.OutputRegionPropertyAbundancePolygonsMaple( "volume-flux-y", bins, results, true ); 
	  else               flux_histogram.OutputRegionPropertyAbundancePolygonsMaple( "volume-flux-z", bins, results, true ); 
	  

	// ------------------------------------------------------------------------------------
	// 7. calculate the relative proportions of the fracture and the matrix flux
	// ------------------------------------------------------------------------------------
	  flux_ratio[i] = fractureMatrixFluxRatio( model, "matrix", fluid_viscosity,
	                                           farfield_pf_gradient, xsec_area, model_throughput );

	  matrix_permeability[i] = equivalentMatrixPermeability( model, "matrix", "fractures", fluid_viscosity,
	                                                         farfield_pf_gradient, xsec_area, model_throughput );
  } // end of loop over three coordinate directions
  
  
  // display effective properties to screen
  cout <<"\nk_equiv_qfqm_Analysis3D: The permeability range of the model was:"<< endl;
  printRangeOfVariable( model, stdio, "permeability" );

  cout <<"\nk_equiv_qfqm_Analysis3D: The equivalent (flow-based upscaling derived) permeability [m2] of the model is: ";
  cout <<"\nin the x-direction: " << k_equivalent[0];
  cout <<"\nin the y-direction: " << k_equivalent[1];
  cout <<"\nin the z-direction: " << k_equivalent[2];
  cout << endl;
  
  cout <<"\nk_equiv_qfqm_Analysis3D: Contribution of the fractures relative to the matrix: ";
  cout <<"\nin the x-direction: " << flux_ratio[0];
  cout <<"\nin the y-direction: " << flux_ratio[1];
  cout <<"\nin the z-direction: " << flux_ratio[2];
  cout << endl;

  cout <<"\nk_equiv_qfqm_Analysis3D: The matrix would need to have a permeability of ";
  cout <<"\nin the x-direction: " << matrix_permeability[0];
  cout <<"\nin the y-direction: " << matrix_permeability[1];
  cout <<"\nin the z-direction: " << matrix_permeability[2];
  cout <<"\nm2 in order to carry the same amount of flow as the fractures."<< endl;
  
  // and write them to a file
  string effective_properties_file_name(model_name);
  effective_properties_file_name += "-equivalent-properties.text";
  ofstream  ofs;
  ofs.open( effective_properties_file_name.c_str(), ios::out|ios::trunc );

  ofs <<"\nk_equiv_qfqm_Analysis3D: The equivalent permeability [m2] of the model is: ";
  ofs <<"\nin the x-direction: " << k_equivalent[0];
  ofs <<"\nin the y-direction: " << k_equivalent[1];
  ofs <<"\nin the z-direction: " << k_equivalent[2];
  ofs << endl;
  
  ofs <<"\nk_equiv_qfqm_Analysis3D: Contribution of the fractures relative to the matrix: ";
  ofs <<"\nin the x-direction: " << flux_ratio[0];
  ofs <<"\nin the y-direction: " << flux_ratio[1];
  ofs <<"\nin the z-direction: " << flux_ratio[2];
  ofs << endl;

  ofs <<"\nk_equiv_qfqm_Analysis3D: The matrix would need to have a permeability of ";
  ofs <<"\nin the x-direction: " << matrix_permeability[0];
  ofs <<"\nin the y-direction: " << matrix_permeability[1];
  ofs <<"\nin the z-direction: " << matrix_permeability[2];
  ofs <<"\nm2 in order to carry the same amount of flow as the fractures."<< endl;
  
  ofs.close();
  
  cout <<"\nk_equiv_qfqm_Analysis3D: That's it..."<< endl;
  return 0;
  
} // end keff_qfqm_linear_FEM_3D




/**
    compares how much flow you actually get as compared to a prediction
    based on average permeabilty.
*/
template<uint32_t dim>
void  analyzeSensitivity( Model<dim>& sg, const char* group, Standard_IO_Handler& io, 
                          Interrelation<dim>& itr, PDE_Integrator<dim,Element>& algo )
 {
    for ( ;; ) {
          Region<dim>&  gref(sg.Region(group));
          double average_k = gref.Average( "permeability" ); 
          cout <<"\nanalyze_sensitivity: Current average permeability, k = "<< average_k;
          cout <<"; enter new k value (-1. to break loop): ";
          double perm;
          cin >> perm;
          if ( perm < 0. ) throw out_of_range("analyze_sensitivity: negative permeability cannot be entered.");
          gref.InputPropertyValue( "permeability", makeScalar(PLAIN,perm), COMPLETE );
          sg.Apply( itr );
          sg.Apply( algo );
          printRangeOfVariable( sg, io, group, "fluid pressure" );
          printRangeOfVariable( sg, io, group, "velocity" );
          printRangeOfVariable( sg, io, group, "volume flux" );
      }
 }
 
template
void  analyzeSensitivity( Model<2U>&, const char*, Standard_IO_Handler&, 
                          Interrelation<2U>&, PDE_Integrator<2U,Element>& );
template
void  analyzeSensitivity( Model<3U>&, const char*, Standard_IO_Handler&, 
                          Interrelation<3U>&, PDE_Integrator<3U,Element>& );
 
 
 



/** 
    WHO ???          works only for Surface element types
*/
void computeApertureDistribution( Model<3U>& sg, const char* groupname )
 {
    const ScalarVariable  bc(DIRICH,0.), perm(PLAIN, 1.), src(PLAIN, 0.001);
    Region<3U>&   gref(sg.Region(groupname));
    
    // 1. boundary conditions 
    gref.InputPropertyValue( "fluid pressure", bc, PERIMETER );
    gref.ChangePropertyStatus( "fluid pressure", DIRICH, PERIMETER );
    // an inflating source term
    gref.InputPropertyValue( "permeability",       perm, COMPLETE );
    gref.InputPropertyValue( "fluid volume source", src, COMPLETE );
    // 2. the computation
    SAMG_Solver solver;
    PDE_Integrator<3U,Element>  aperture( solver );

    NumIntegral_dNT_op_dN_dV<3U>  conductance( sg.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<3U>    source( sg.Database(), "fluid volume source", "fluid pressure" );
    aperture.Add( &conductance );
    aperture.Add( &source );
    aperture.IntegrateOver( gref );
 
 } // end computeApertureDistribution



} // end namespace csmp





























