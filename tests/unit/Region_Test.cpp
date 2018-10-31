#include "Region_Test.h"
#include "RegionInterface.h"

#include "Point.h"
#include "InputDataManager.h"
#include "Region.h"
#include "ANSYS_Model3D.h"
#include "PDE_Integrator.h"
#include "PropertyConstraints.h"
#include "PropertyHandle.h"
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"
#include "CSMP_highLevelUtilities.h"
#include "RegionMonitor.h"
#include "vsetMakers.h"

// algebraic multigrid solvers
#include "LinearSolver.h"

// Interrelations
#include "ConstantFactor.h"

// output
#include "VTK_Interface.h"

#include <iostream>

#define DIM 3U

using namespace std;

namespace csmp {

Region_Test::Region_Test( bool verbose )
 : verbose_(verbose)
{
}


Region_Test::~Region_Test()
{

}

/**
    SKM revised 8/3/2016 
    
    @todo test for Model as non-unique region is commented out because rest of code is broken
    
    @todo check again for code coverage: a number of important methods are not tested yet like:
         CreateBetween()
         CorrectLowDimRegionOrientation()
         VolumeIntegral_x_Thickness()
 
*/
void Region_Test::run()
{
    // 1. integrity check on regions written to disk
    // -----------------------------------------------------
    TestBoundaryFaceFunctionality(); // pyramid-hexa model generated with vsetMakers
    TestBoundaryFaceFunctionality("cube_flag");
	//JC: invalid data. check if it is a discontiguous model.
    //TestBoundaryFaceFunctionality("hex1_3");


    // 2. test of the functionality of regions
    // -----------------------------------------------------
    const char* model_name="prism_test";
    const string varFileName( "CSMP-1phase-variables.txt" );
    ANSYS_Model3D model( model_name, varFileName.data() );
    _test( consistencyCheckNeighborVersusPerimeterFaces( model ) );

    InputDataManager<DIM>  model_configuration;

    model_configuration.ConfigureFromFile( model, model_name,
                                           false, true, true, true, false );

    TestRegionFileInputOutput( model, "Model" );

    // we will solve a 3D steady state pressure problem (model contain 3 fractures)
    // calculate hydraulic conductivity, K = k / mu
    const double64 fluid_viscosity(1.0e-03);

    ConstantFactor<DIM,divides>  conductivity( model.Database(),
                                               "conductivity", "permeability",
                                               fluid_viscosity );
    model.Apply( conductivity );

    printRangeOfVariable( model, "conductivity" );

    VTK_Interface<DIM>  vtk_output;
    if ( verbose_ ) vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );

    // compute steady-state fluid pressure and Darcy velocity
    SteadyStateDiffusor<DIM,Region>  steady_state_pressure( model, "conductivity", "fluid pressure", "fluid volume source" );
    VelocityAndVolumeFlux<DIM,Element<DIM> >  postpro( model, "conductivity", "porosity", "fluid pressure" );
    steady_state_pressure.AddPostProcess( &postpro );

    steady_state_pressure.ComputeSteadyState( model.Region("Model") );

    // output result
    if ( verbose_ ) {
        printRangeOfVariable( model, "fluid pressure" );
        printRangeOfVariable( model, "velocity" );
        printRangeOfVariable( model, "pore velocity" );
        vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1, true );
        vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1, true );
        vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",    1, true );
        vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "FRAC_VOLUMES", "fluid pressure", 1, true );
      }
    // Testing start
    // 1. testing whether the model contains the right regions and test their uniqueness
    if ( verbose_ ) cout <<"\nRegion_Test::run: Does the model contain FRAC_VOLUMES and MATRIX regions? ";
    _test(  model.ContainsRegion("FRAC_VOLUMES") == true );
    _test(  model.ContainsRegion("MATRIX") == true );
    _test(  model.ContainsRegion("fracture") == false );

    if ( verbose_ ) cout <<"\nRegion_Test::run: Are they unique (have no overlap)? ";
    _test(  model.IsUnique("MATRIX") == true);
    _test(  model.IsUnique("FRAC_VOLUMES") == true );
    _test(  model.IsUnique("Model") == false ); // should be non-unique when other regions exist

    // creation of multiple regions from discrete property values
    // illustrated using the 'permeability'
    if ( verbose_ ) {
        printRangeOfVariable( model, "permeability" );
        vtk_output.OutputDataToVTK( model, "permeability",    "permeability", 1, true );
        vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "permeability_FRAC_VOLUMES", "permeability", 1, true );
      }

    // for each distinct value of the 'permeability' a region is created
    set<string>  group_names;

    model.FormRegionsFromPropertyValues( "permeability", group_names );

    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: FormRegionsFromPropertyValues()  creating 'permeability' regions: ";
    size_t count(0.);

    for ( set<string>::const_iterator it=group_names.begin(); it!=group_names.end(); it++ )
    {
        if ( verbose_ ) cout << (*it) <<" ";
        if ( verbose_ ) vtk_output.OutputDataToVTK( model, (*it).c_str(), "permeability", "permeability", 1, true );
        count++;
        size_t num_elements_perm_region(model.Region((*it).c_str()).Elements());
        _test((  num_elements_perm_region == model.Region("MATRIX").Elements()  )
                 or  ( num_elements_perm_region == model.Region("FRAC_VOLUMES").Elements() ));
    }
    if ( verbose_ ) cout << endl;
    cout.flush();
    _test( count == 2 ) ;


    // merging these regions into a single group 'permeability_regions'
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: MergeRegions()  merging these regions into new region of elements 'permeability_regions'."<< endl;
    model.MergeRegions( group_names, "permeability_regions" );

    if ( verbose_ ) vtk_output.OutputDataToVTK( model, "permeability_regions", "permeability_regions", "permeability", 1, true );

    _test( model.Region("permeability_regions").Elements() == model.Region("Model").Elements() );
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: RemoveRegion()  removing region 'permeability_regions' and other new regions."<< endl;

    model.RemoveRegion( "permeability_regions", false );

    //test of remove
    // removal of the new permeability-based unique regions:
    _test( model.ContainsRegion("permeability_regions") == false );

    for ( set<string>::const_iterator it=group_names.begin(); it!=group_names.end(); it++ )
    {
        model.RemoveRegion( (*it).c_str(), false );
       _test( model.ContainsRegion((*it).c_str()) == false );
    }

    // creation of a region from combined k & pf ranges using PropertyConstraints
    // --------------------------------------------------------------------------
    if ( verbose_ ) cout <<"\nRegion_Test::run: FormRegionFrom()  forming new region using permeability and fluid pressure constraints."<< endl;
    PropertyConstraints  rangeP( "permeability", 1.e-12, 1.e-11 );

    rangeP.AddConstraint("fluid pressure", 1.5e6, 1.8e6 );

    model.FormRegionFrom( "pressure_permeability_overlap1", rangeP );
    if ( verbose_ ) vtk_output.OutputDataToVTK( model, "pressure_permeability_overlap1",  "pressure_permeability_overlap1","fluid pressure", 1, true );

    vector<Node<DIM>*>::iterator nodes_end=model.Region("pressure_permeability_overlap1").NodesEnd();
    vector<Node<DIM>*>::iterator nodes_begin=model.Region( "pressure_permeability_overlap1" ).NodesBegin();

    vector<Element<DIM>*>::iterator Elements_end=model.Region("pressure_permeability_overlap1").ElementsEnd();
    vector<Element<DIM>*>::iterator Elements_begin=model.Region( "pressure_permeability_overlap1").ElementsBegin();

    Index PresKey( model.Database().StorageKey("fluid pressure" ) );
    Index PermKey( model.Database().StorageKey("permeability" ) );

    size_t tes(0);

    // loop over nodes for pressure test
    for ( vector<Node<DIM>*>::iterator npit= nodes_begin; npit!=nodes_end; npit++ )
      {
          double64 pres;
          pres = (*npit)->Read(PresKey);
          if ( (pres < 1.5e6 or pres > 1.8e6) ) tes++;
      }
    // test that the pressure region does not contain any nodes out of range
    _test( tes == 0 );

    // loop over elements for permeability test
    tes=0;
    for ( vector<Element<DIM>*>::iterator npit= Elements_begin; npit!=Elements_end; npit++ )
      {
          double64 perm;
          perm = (*npit)->Read(PermKey);
          if ( (perm > 1.e-11 or perm < 1.e-12 ) ) tes++;
      }
    _test( tes == 0 );

    // copy region
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: CopyRegion()  copying a region."<< endl;
    model.CopyRegion( "pressure_permeability_overlap1", "pressure_permeability_overlap2" );

    // visualisation
    if ( verbose_  ) {
        vtk_output.OutputDataToVTK( model, "pressure_permeability_overlap1", "region_kpf1", "fluid pressure", 1, true );
        vtk_output.OutputDataToVTK( model, "pressure_permeability_overlap2", "region_kpf2", "fluid pressure", 1, true );
      }
    _test( model.Region("pressure_permeability_overlap1").Elements() == model.Region("pressure_permeability_overlap1").Elements() );
    _test( model.ContainsRegion("pressure_permeability_overlap1") == true );
    _test( model.ContainsRegion("pressure_permeability_overlap2") == true );

    // removal
    model.RemoveRegion( "pressure_permeability_overlap1", false );
    model.RemoveRegion( "pressure_permeability_overlap2", false );
    // regions should no longer be there
    _test( model.ContainsRegion("pressure_permeability_overlap1") == false );
    _test( model.ContainsRegion("pressure_permeability_overlap2") == false );


    // creating a region from properties within a specified range
    // ----------------------------------------------------------
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: FormRegionFrom()  forming new region using a fluid pressure range."<< endl;
    model.FormRegionFrom( "pf_window", "fluid pressure", 1.5e6, 1.7e6 );
    _test( model.ContainsRegion("pf_window") == true );
    if ( verbose_ ) vtk_output.OutputDataToVTK( model, "pf_window", "pf_window", "fluid pressure", 1, true );
    model.RemoveRegion( "pf_window", false );


    // creating and testing rectangular regions
    // ----------------------------------------
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: FormRectangularRegion()  forming new box-shaped region inside the model."<< endl;
    Point<3U> pMin, pMax;
    model.MinMaxCoordinates( pMin, pMax );
    cout << "\npMin: " << pMin << endl;
    cout << "\npMax: " << pMax << endl;
    Point<3U> p1( pMin + pMax/4. );
    Point<3U> p2( pMax - pMax/4. );
    if ( verbose_ ) {
        cout << "\ncoordinate point 1: "<< p1 << endl;
        cout << "\ncoordinate point 2: "<< p2 << endl;
        cout << endl;
      }

    model.FormRectangularRegion( "rectangular region", p1, p2 );
    if ( verbose_ ) vtk_output.OutputDataToVTK( model, "rectangular region", "box", "fluid pressure", 1, true );
    Point<3U> p1_newRegion;
    Point<3U> p2_newRegion;
    model.Region("rectangular region").MinMaxCoordinates( p1_newRegion, p2_newRegion );
    if ( verbose_ ) {
        cout << "\ncoordinate point 1: "<< p1_newRegion << endl;
        cout << "\ncoordinate point 2: "<< p1_newRegion << endl;
        cout << endl;
      }

    // is the region contained in the bounding box with the corners p1 and and p2 that was used to create the region?
    _test( p1 == p1_newRegion or p1 < p1_newRegion );
    _test( p2 == p2_newRegion or p2 > p2_newRegion );

    model.RemoveRegion( "rectangular region", false );


    // Region Union Test
    if ( verbose_ ) cout <<"\n\n\nRegion_Test::run: RegionUnion()  combining 'FRAC_VOLUMES' and 'MATRIX' regions."<< endl;
    model.RegionUnion( "FRAC_VOLUMES", "MATRIX", "MODEL" );
    _test( model.Region("Model").Elements() == model.Region("MODEL").Elements() );

    if ( verbose_ )
      vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "volumetric_fractures", "fluid pressure", 1, true );

    // Region Include Test
    if ( verbose_ ) cout <<"\nRegion_Test::run: Is FRAC_VOLUMES a part of Model(yes):  "<< endl;
    _test( model.RegionIncludes( "MODEL","FRAC_VOLUMES" ) == true );

    if ( verbose_ ) cout <<"\nRegion_Test::run: Is Model a part of FRAC_VOLUMES(no):  "<< endl;
    _test( model.RegionIncludes( "FRAC_VOLUMES","MODEL" ) == false );

    model.RemoveRegion( "MODEL", false );

    // Region Intersection Test
    if ( verbose_ ) {
        cout <<"\n\n\nRegion_Test::run: RegionIntersection() between FRAC_VOLUMES and MATRIX."<< endl;
        cout <<"\nRegion_Test::run: Does region FRAC_VOLUMES overlap with MATRIX? (no): ";
      }
    _test( model.RegionIntersection( "FRAC_VOLUMES", "MATRIX", "empty_group" ) == false );
     // NB: empty_group will not be formed

    if ( verbose_ ) cout <<"\nRegion_Test::run: Does region FRAC_VOLUMES overlap with Model? (yes): ";
    _test( model.RegionIntersection( "FRAC_VOLUMES", "Model", "model and fractures" ) == true );

    // Region difference test (returns the data in A that is not in B)
    if ( verbose_ ) cout <<"\nRegion_Test::run: What is the difference between MATRIX and Model (FRAC_VOLUMES): ";
    _test( model.RegionDifference( "Model", "MATRIX", "difference" ) == true );
    _test( model.Region("difference").Elements() == model.Region("FRAC_VOLUMES").Elements() );

    if ( verbose_ )
      vtk_output.OutputDataToVTK( model, "difference", "fluid-pressure", "fluid pressure", 1, true );

    // Region symmetric difference Test (the symmetric difference is the union without the intersection)
    if ( verbose_ ) {
        cout <<"\n\n\nRegion_Test::run: RegionSymmetricDifference()"<< endl;
        cout <<"\nRegion_Test::run: What is the symmetric difference between MATRIX and Model (FRAC_VOLUMES): ";
      }
    _test( model.RegionSymmetricDifference( "MATRIX", "Model", "diff" ) == true );
    // should contain all volumetric elements from FRAC_VOLUMES
    _test( model.Region("diff").Elements() == model.Region("FRAC_VOLUMES").Elements() );
    model.RemoveRegion( "diff", false );

    // breaking a region into contiguous sub-regions
    // ---------------------------------------------
    if ( verbose_ ) cout << "\n\n\nRegion_Test::run: number of contiguous sub regions: ";
    _test(model.PartitionRegionIntoContiguousSubRegions( "FRAC_VOLUMES" ) == 0 );
    // removal of new partitions
	//JC: check it later due to the comment from SKM: logic of this method seems to be broken and it does not always work.Refactor!
    //_test(model.RemoveRegionPartitionsFor("FRAC_VOLUMES") == 0 );

} // end run




/**
    vsetMaker version of previous method
 
    uses SurfaceArea() to test whether the face connectivity is preserved
    when a model is saved to disk and reloaded afterwards.
 
    try input models "cube_flag" and a hexahedral one
 
    @attention the model must be box-shaped else, surface area computation will fail.
*/
bool Region_Test::TestBoundaryFaceFunctionality()
{
    // 1. building and saving test model to disk
    // -----------------------------------------
    bool intact_functionality(true);
    const string varFileName( "CSMP-1phase-variables.txt" );
    const string model_name("PyramidHexaPatch");
    const bool   skewed_elements(false); // otherwise model is not a box anymore
    VSet<3U>     vset;
  
    test_Create_Pyramid_Hexa_VSet( vset, skewed_elements );
  
    Model<3U>  model1( vset, varFileName.c_str(), true );
    _test( consistencyCheckNeighborVersusPerimeterFaces( model1 ) );

    // corner_points
    Point<3U> xyz_min, xyz_max;
    model1.MinMaxCoordinates( xyz_min, xyz_max );
    const double64 dx(xyz_max[0]-xyz_min[0]), dy(xyz_max[1]-xyz_min[1]), dz(xyz_max[2]-xyz_min[2]);
    const double64 surface_area = 2.*dx*dy + 2.*dx*dz + 2.*dy*dz;
  
    const Region<3U>& model1_domain(model1.Region("Model"));
    const double64 surface_area1 = model1_domain.SurfaceArea();
    const double64 volume1       = model1_domain.Volume();
  
    _equal( surface_area1, surface_area, numeric_limits<double64>::epsilon() * surface_area );
  
    // extracting perimeter elements to set for comparison
    size_t perimeter_elements(model1_domain.PerimeterElements());
    // extracting the perimeter face vector for comparison with re-read model2
    vector<vector<int8> > perimeter_faces;
    perimeter_faces.reserve(model1_domain.PerimeterElements());
    for ( size_t e=model1_domain.InteriorElements(); e<model1_domain.Elements(); ++e ) {
         vector<int8>  face_vec;
         for ( size_t i=0U; i<model1_domain.PerimeterFaces(e); ++i )
           face_vec.push_back( static_cast<int8>(model1_domain.PerimeterFace(e,i)) );
         perimeter_faces.push_back( move(face_vec) );
      }
    // counting the perimeter faces
    size_t n_perimeter_faces(0U);
    for ( auto it : perimeter_faces ) n_perimeter_faces += it.size();
  
    model1.OutputToBinaryFile("model1");
  
  
    // 2. reading the model back in and testing SurfaceArea again
    // ----------------------------------------------------------
    Model<3U>  model2(string("model1"));
    _test( consistencyCheckNeighborVersusPerimeterFaces( model2 ) );

    const Region<3U>& model2_domain(model2.Region("Model"));
 
    // test 0: same number of perimeter elementds
    size_t perimeter_elements2(model2_domain.PerimeterElements());
    _test( perimeter_elements = perimeter_elements2 );
  
    // test 1: re-read model2
    vector<vector<int8> > perimeter_faces2;
    perimeter_faces2.reserve(model2_domain.PerimeterElements());
    for ( size_t e=model2_domain.InteriorElements(); e<model2_domain.Elements(); ++e ) {
         vector<int8>  face_vec;
         for ( size_t i=0U; i<model2_domain.PerimeterFaces(e); ++i )
           face_vec.push_back( static_cast<int8>(model2_domain.PerimeterFace(e,i)) );
         perimeter_faces2.push_back( move(face_vec) );
      }
    // counting the perimeter faces
    size_t n_perimeter_faces2(0U);
    for ( auto it : perimeter_faces2 ) n_perimeter_faces2 += it.size();
  
    // test 3: same number of perimeter faces ?
    _test( n_perimeter_faces == n_perimeter_faces2 );

    // test 5: is the content of the perimeter face vectors actually the same ?
	_test( equal(perimeter_faces2.begin(), perimeter_faces2.end(), perimeter_faces.begin(), perimeter_faces.end() ) );

    // test 6: verifying that the outer surface area and volume of in the re-read CSMP native model is the same
    // region surface area
    const double64 surface_area2 = model2_domain.SurfaceArea();
    _equal( surface_area1, surface_area2, numeric_limits<double64>::epsilon() * surface_area1 );
    // region volume
    const double64 volume2 = model2_domain.Volume();
    _equal( volume1, volume2, numeric_limits<double64>::epsilon() * volume1 );

    // 3. getting extra diagnostics from the RegionMonitor
    // ----------------------------------------------------------
    model1.InputPropertyValue ( "fluid pressure", makeScalar(PLAIN,0.) );
    model1.InputPropertyValue ( "permeability", makeScalar(PLAIN,1.0e-12) );
    model1.InputPropertyValue ( "porosity", makeScalar(PLAIN,1.) );
    model1.InputPropertyValue ( "fluid volume source", makeScalar(PLAIN,0.) );
    model1.InputPropertyValue ( "nodal fluid volume source", makeScalar(PLAIN,0.) );
    model1.InputPropertyValue ( "concentration", makeScalar(PLAIN,0.) );
    // boundary conditions
    const double64 pressure (10*101325.);
    model1.InputBoundaryValue( LEFT,  "fluid pressure",   makeScalar(DIRICH, 0.) );
    model1.InputBoundaryValue( RIGHT, "fluid pressure",   makeScalar(DIRICH,pressure) ); // 1 bar

    const bool integrate_pore_volume_only(true);
    RegionMonitor<3U>  monitor( model2, "porosity", "fluid pressure",  integrate_pore_volume_only );
  
    return intact_functionality;

} // end TestBoundaryFaceFunctionality





bool Region_Test::TestBoundaryFaceFunctionality( const string& model_name )
{
    // 1. building and saving test model to disk
    bool intact_functionality(true);
    const string varFileName( "CSMP-1phase-variables.txt" );
  
    // ansys model
    ANSYS_Model3D model1( model_name.c_str(), varFileName.c_str() );
    _test( consistencyCheckNeighborVersusPerimeterFaces( model1 ) );
  
    // corner_points
    Point<3U> xyz_min, xyz_max;
    model1.MinMaxCoordinates( xyz_min, xyz_max );
    const double64 dx(xyz_max[0]-xyz_min[0]), dy(xyz_max[1]-xyz_min[1]), dz(xyz_max[2]-xyz_min[2]);
    const double64 surface_area = 2.*dx*dy + 2.*dx*dz + 2.*dy*dz;

    InputDataManager<DIM>  model_configuration;
    model_configuration.ConfigureFromFile( model1, model_name.c_str(),
                                           false, true, true, true, false  );
  
    const Region<3U>& model1_domain(model1.Region("Model"));
    double64 surface_area1 = model1_domain.SurfaceArea();
  
    _equal( surface_area1, surface_area, numeric_limits<double64>::epsilon() * surface_area * 100. );
  
    model1.OutputToBinaryFile("model1");
  
    // 2. reading the model back in and testing SurfaceArea again
    Model<3U>  model2(string("model1"));
    _test( consistencyCheckNeighborVersusPerimeterFaces( model2 ) );

    const Region<3U>& model2_domain(model2.Region("Model"));
    double64 surface_area2 = model2_domain.SurfaceArea();
  
    _equal( surface_area1, surface_area2, numeric_limits<double64>::epsilon() * surface_area1 * 100. );
  
    const bool integrate_pore_volume_only(true);
    RegionMonitor<3U>  monitor( model2, "porosity", "fluid pressure",  integrate_pore_volume_only);
  
    return intact_functionality;

} // end TestBoundaryFaceFunctionality





/**
    Compares whether the data of the region class are correctly recovered from binary file.
*/
bool Region_Test::TestRegionFileInputOutput( const Model<3U>& model, const char* region )
 {
     // 1. creating sets for interior and perimeter nodes, elements, boundary faces etc.
     // --------------------------------------------------------------------------------
     set<Node<3U>*>    interior_nodes, perimeter_nodes;
     set<Element<3U>*> interior_elements, perimeter_elements;
 
     const Region<3>&  domain(model.Region(region));
     for ( auto nit=domain.NodesBegin(); nit!=domain.PerimeterNodesBegin(); nit++ ) interior_nodes.insert( (*nit) );
     for ( auto nit=domain.PerimeterNodesBegin(); nit!=domain.NodesEnd(); nit++ )   perimeter_nodes.insert( (*nit) );
     for ( auto eit=domain.ElementsBegin(); eit!=domain.PerimeterElementsEnd(); eit++ ) interior_elements.insert( (*eit) );
     for ( auto eit=domain.PerimeterElementsBegin(); eit!=domain.ElementsEnd(); eit++ ) perimeter_elements.insert( (*eit) );
   
     size_t n_perimeter_nodes(domain.PerimeterNodes());
     size_t n_perimeter_elements(domain.PerimeterElements());
     // extracting the perimeter face vector for comparison with re-read model2
     vector<vector<int8> > perimeter_faces;
     perimeter_faces.reserve(domain.PerimeterElements());
     for ( size_t e=domain.InteriorElements(); e<domain.Elements(); ++e ) {
          vector<int8>  face_vec;
          for ( size_t i=0U; i<domain.PerimeterFaces(e); ++i )
            face_vec.push_back( static_cast<int8>(domain.PerimeterFace(e,i)) );
          perimeter_faces.push_back( move(face_vec) );
       }

     // saving and retrieving the model from file
     model.OutputToBinaryFile(model.Name());

  
     // 2. retrieving the model and getting the same diagnostics
     // --------------------------------------------------------------------------------
     Model<3U>  model2(string(model.Name()));
     const Region<3>&  domain2(model.Region(region));

     set<Node<3U>*>    interior_nodes2, perimeter_nodes2;
     set<Element<3U>*> interior_elements2, perimeter_elements2;

     for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); nit++ ) interior_nodes2.insert( (*nit) );
     for ( auto nit=domain2.PerimeterNodesBegin(); nit!=domain2.NodesEnd(); nit++ )   perimeter_nodes2.insert( (*nit) );
     for ( auto eit=domain2.ElementsBegin(); eit!=domain2.PerimeterElementsEnd(); eit++ ) interior_elements2.insert( (*eit) );
     for ( auto eit=domain2.PerimeterElementsBegin(); eit!=domain2.ElementsEnd(); eit++ ) perimeter_elements2.insert( (*eit) );
   
     size_t n_perimeter_nodes2(domain2.PerimeterNodes());
     size_t n_perimeter_elements2(domain2.PerimeterElements());
     // extracting the perimeter face vector for comparison with re-read model2
     vector<vector<int8> > perimeter_faces2;
     perimeter_faces2.reserve(domain2.PerimeterElements());
     for ( size_t e=domain2.InteriorElements(); e<domain2.Elements(); ++e ) {
          vector<int8>  face_vec;
          for ( size_t i=0U; i<domain2.PerimeterFaces(e); ++i )
            face_vec.push_back( static_cast<int8>(domain2.PerimeterFace(e,i)) );
          perimeter_faces2.push_back( move(face_vec) );
       }  
	 
     // 3. Testing
     // --------------------------------------------------------------------------------
     _test( interior_nodes == interior_nodes2 );
     _test( perimeter_nodes == perimeter_nodes2 );
     _test( interior_elements == interior_elements2 );
     _test( perimeter_elements == perimeter_elements2 );
     _test( n_perimeter_nodes == n_perimeter_nodes2 );
     _test( n_perimeter_elements == n_perimeter_elements2 );
	 _test( equal(perimeter_faces2.begin(), perimeter_faces2.end(), perimeter_faces.begin(), perimeter_faces.end()) );  
	 

    return true;

 } // end TestRegionFileInputOutput





/**
    Checks whether the neighbor information matches the boundary face info for region "Model"
    Two potential failures are detected and reported:
    1. the number of perimeter faces is not correct
    2. the ids of the perimeter faces are not correct
*/
bool consistencyCheckNeighborVersusPerimeterFaces( const Model<3U>& model )
 {
    size_t consistency_check_failures(0U);
   
    const Region<3U>& model_domain(model.Region("Model"));
   
    for ( size_t e=model_domain.InteriorElements(); e<model_domain.Elements(); ++e )
     {
         const size_t expected_perimeter_faces(model_domain.PerimeterFaces(e));
         size_t       perimeter_faces(0U);
         
         for ( size_t j=0U; j<model_domain.E(e)->Neighbors(); ++j )
           // if there is no neighbor, there should be a boundary face corresponding to this
           if ( model_domain.E(e)->Neighbor(j) == nullptr )
             {
                // checking the perimeter face information
                size_t face = model_domain.PerimeterFace( e, perimeter_faces );
                if ( j != face )
                  consistency_check_failures++;
                perimeter_faces++;
             }
           if ( expected_perimeter_faces != perimeter_faces )
             consistency_check_failures++;
      }
  
    return ( consistency_check_failures == 0U );
   
 } // end consistencyCheckNeighborVersusPerimeterFaces




} // end csmp
