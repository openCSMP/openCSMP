#include "Region_Example.h"

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

// Interrelations
#include "ConstantFactor.h"

// output
#include "VTK_Interface.h"
#include "FEM_Data.h"

#define DIM 3U

using namespace std;

namespace csmp{

void Region_Example::Specifications()
{
  SetTitle( "Operations involving csmp::Region objects" );
  SetDifficulty( 3 );
  SetCategory( "Software Functionality" );
  AddAuthor( "SKM" );
  AddDescription( "demonstrations of all methods of the RegionInterface" );
  AddDescription( "exercise: step through this example using the debugger to see what is happening" );
  AddDescription( "source in: Region_Example.cpp" );
  AddDescription( "application and features of csmp::Regions" );
  AddRequirement( "'fracs4' .asc, .dat, -regions. & -configuration.txt" );
  AddRequirement( "CSMP-1phase-variables.txt" );
}



void Region_Example::Run()
{
  //  Create a model from ANSYS mesh and configure it from file
   const char* model_name="fracs4";
   ANSYS_Model3D  model( model_name, "CSMP-1phase-variables.txt");
   printModelDimensions( model, true );

   // testing whether model contains the desired regions
   cout <<"\nmain: Does the model contain FRACS and MATRIX regions? ";
   cout << model.ContainsRegion("FRACS") <<" "<< model.ContainsRegion("MATRIX") << endl;
   cout <<"or a non-existing one called 'dummy'? - "<< model.ContainsRegion("dummy");
   cout <<"\nmain: Are they unique (have no overlap)? ";
   cout << model.IsUnique("FRACS") <<" "<< model.IsUnique("MATRIX") << endl;

   // configure model / regions from file 'fracs4-configuration.txt'
   InputDataManager<DIM>().ConfigureFromFile( model, model_name,
                                              false, true, true, true, false );

   // creation and destruction of a variable at runtime (to visualise permeability)
   PropertyHandle<DIM>   k( model, "permeability", SCALAR, ELEMENT );
   PropertyHandle<DIM>*  log_k = new PropertyHandle<3U>( model, "log10 permeability", SCALAR, ELEMENT );
   *log_k  = 0.;
   *log_k += k;
   log_k->Log10();

   // output in VTK format (visualisation toolkit, paraview, mayavi)
   VTK_Interface<DIM>  vtk_output;
   vtk_output.OutputDataToVTK( model, "log10-permeability", "log10 permeability", 0 );
   delete log_k;

   // calculate hydraulic conductivity, K = k / mu
   const double64 fluid_viscosity(1.0e-03);
   ConstantFactor<DIM,divides>  conductivity( model.Database(),
                                             "conductivity", "permeability",
                                              fluid_viscosity );
   model.Apply( conductivity );
   printRangeOfVariable( model, "conductivity" );

   vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );

   // compute steady-state fluid pressure and Darcy velocity
   SteadyStateDiffusor<DIM,Region>  steady_state_pressure( model, "conductivity", "fluid pressure", "fluid volume source" );
   VelocityAndVolumeFlux<DIM,Element<DIM> >  postpro( model, "conductivity", "porosity", "fluid pressure" );
   steady_state_pressure.AddPostProcess( &postpro );
   steady_state_pressure.IntegrateOver( model.Region("Model") );

   printRangeOfVariable( model, "fluid pressure" );
   printRangeOfVariable( model, "velocity" );
   printRangeOfVariable( model, "pore velocity" );

   vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
   vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1 );
   vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",    1 );


  // -----------------------------------------------------------------------
  //
  //  Demonstrate how to create model subregions
  /*
      - This is the task of RegionInterface policy of the Model
  */
  // -----------------------------------------------------------------------
   // VTK output of existing regions
   for ( map<string,Region<3> >::const_iterator
         it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); it++ )
     vtk_output.OutputDataToVTK (model, (*it).first.c_str(), "fluid-pressure", "fluid pressure", 1, true );

   // creation of multiple regions from discrete property values illustrated using the 'permeability'
   printRangeOfVariable( model, "permeability" );
   // for each distinct value of the 'permeability' a region is created
   set<string>  group_names;
   // test 1: O.K.
   model.FormRegionsFromPropertyValues( "permeability", group_names );
   //    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   cout <<"\n\n\nmain: FormRegionsFromPropertyValues()  creating 'permeability' regions: ";
   for ( set<string>::const_iterator it=group_names.begin(); it!=group_names.end(); it++ ) cout << (*it) <<" ";
   cout << endl;
   cout.flush();
   // merging these regions into a single region 'permeability_regions'
   // test 2: O.K.
   cout <<"\n\n\nmain: MergeRegions()  merging these regions into new regions of elements 'permeability_regions'."<< endl;
   model.MergeRegions( group_names, "permeability_regions" );
   //    ^^^^^^^^^^^^
   vtk_output.OutputDataToVTK( model, "permeability_regions", "permeability_regions", "permeability", 1, true );
   // test 3: O.K.
   cout <<"\n\n\nmain: RemoveRegion()  removing region 'permeability_regions' and other new regions."<< endl;
   model.RemoveRegion( "permeability_regions", false );
   //    ^^^^^^^^^^^^
   // removal of the new permeability-based unique regions:
   for ( set<string>::const_iterator it=group_names.begin(); it!=group_names.end(); it++ )
     model.RemoveRegion( (*it).c_str(), false );


   // creation of a region from combined property values using PropertyConstraints
   // test 4: O.K.
   cout <<"\nmain: FormRegionFrom()  forming new region using permeability and fluid pressure constraints."<< endl;
   PropertyConstraints  rangeP( "permeability", 1.e-11, 1.e-9 );
   rangeP.AddConstraint("fluid pressure", 1.5e6, 1.8e7 );
   model.FormRegionFrom( "pressure_permeability_overlap1", rangeP, false );
   //    ^^^^^^^^^^^^^^
   // test 5: O.K.
   cout <<"\n\n\nmain: CopyRegion()  copying a region."<< endl;
   model.CopyRegion( "pressure_permeability_overlap1", "pressure_permeability_overlap2" );
   //    ^^^^^^^^^^
   // visualisation
   vtk_output.OutputDataToVTK( model, "pressure_permeability_overlap1", "region_kpf1", "fluid pressure", 1, true );
   vtk_output.OutputDataToVTK( model, "pressure_permeability_overlap1", "region_kpf1", "fluid pressure", 1, true );
   // removal
   model.RemoveRegion( "pressure_permeability_overlap1", false );
   model.RemoveRegion( "pressure_permeability_overlap2", false );


   // creating a region from element-, node- or other properties within a specified range
   // test 6: O.K.
   cout <<"\n\n\nmain: FormRegionFrom()  forming new region using a fluid pressure range."<< endl;
   model.FormRegionFrom( "pf_window", "fluid pressure", 1.0e6, 1.0e7 );
   //    ^^^^^^^^^^^^^^
   // visualisation
   vtk_output.OutputDataToVTK( model, "pf_window", "pf_window", "fluid pressure", 1, true );
   // removal
   model.RemoveRegion( "pf_window", false );

   // breaking a region into contiguous sub-regions
   cout <<"\n\n\nmain: number of contiguous sub regions: ";
   // test 7:
   cout << model.PartitionRegionIntoContiguousSubRegions( "FRACS" );
   //            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
   // test 8:
   // merging 2 of the new regions into the region FRACS that was removed by the previous operation
   set<string>  merged_regions;
   merged_regions.insert("FRACS1");
   merged_regions.insert("FRACS2");
   merged_regions.insert("FRACS3");
   // now the 3 created subregions are removed
   model.MergeRegions( merged_regions, "FRACS" );
  
   // removal of new partitions
   cout <<"\n\n\nmain: removing the previously created subregions of FRACS:\n";
   cout <<"\n\tsubregions removed: "<<  model.RemoveRegionPartitionsFor("FRACS");
   //                                         ^^^^^^^^^^^^^^^^^^^^^^^^^

   cout <<"\n\n\nmain: FormRectangularRegion()  forming new box-shaped region inside the model."<< endl;
   Point<3U> pMin, pMax;
   model.MinMaxCoordinates( pMin, pMax );
   cout << "\npMin: " << pMin << endl;
   cout << "\npMax: " << pMax << endl;
   Point<3U> p1( pMin + pMax/4. );
   Point<3U> p2( pMax - pMax/4. );
   // test 9: O.K.
   model.FormRectangularRegion( "rectangular region", p1, p2 );
   //    ^^^^^^^^^^^^^^^^^^^^^^
   vtk_output.OutputDataToVTK( model, "rectangular region", "box", "fluid pressure", 1, true );
   model.RemoveRegion( "rectangular region", false );


   // building a region from a subset of element numbers
   // first, a range of element numbers 10 to 100 has to be generated
   cout <<"\n\n\nmain: FormRegionFrom()  forming a new region that contains the elements 10 to 100."<< endl;
   vector<size_t>  element_ids( 100U );
   size_t          ecount(10U);
   for ( size_t i=0U; i<element_ids.size(); i++ ) element_ids[i] = ecount++;
   // test 10: O.K.
   model.FormRegionFrom( "elements10to100", element_ids );
   //    ^^^^^^^^^^^^^^
   vtk_output.OutputDataToVTK( model, "elements10to100", "conductivity", "conductivity", 1, true );


   // assimilating this regions into a new region of elements 110 to 200, using a C++ lambda function
   transform( element_ids.begin(), element_ids.end(), element_ids.begin(), [&]( auto elmt ){ return elmt + 101U; } );
   model.FormRegionFrom( "elements110to200", element_ids );
   // test 11: O.K.
   cout <<"\n\n\nmain: AssimilateRegion()  assimilating this regions into FRACS region."<< endl;
   model.AssimilateRegion( "elements110to200", "elements10to100" );
   //    ^^^^^^^^^^^^^^^^
   cout <<"\nmain: number of elements before and after assimilation: "<< element_ids.size() <<" vs. ";
   cout << model.Region("elements10to100").Elements() << endl;

   vtk_output.OutputDataToVTK( model, "elements10to100", "fluid-pressure", "fluid pressure", 1, true );
   model.RemoveRegion( "elements110to200", false );



  // -----------------------------------------------------------------------
  //
  //  Demonstratation of boolean operations involving groups
  /*
      - This is again possible through the Model<>
  */
  // -----------------------------------------------------------------------
   // test: fails?
   cout <<"\nmain: Is FRACS a part of Model(yes):          "<<  model.RegionIncludes( "Model", "elements10to100" );
   cout <<"\nmain: How about the opposite? (should be no): "<<  model.RegionIncludes( "elements10to100", "Model" );
   //                                                                 ^^^^^^^^^^^^^^
   // here the union should be identical with "Model"
   cout <<"\n\n\nmain: RegionUnion()  combining 'FRACS' and 'MATRIX' regions."<< endl;
   // test :
   model.RegionUnion( "FRACS", "MATRIX", "MODEL" );
   //    ^^^^^^^^^^^
   assert( model.Region("Model").Elements() == model.Region("MODEL").Elements() );
   model.RemoveRegion( "MODEL", false );

   cout <<"\n\n\nmain: RegionIntersection() between FRACS and MATRIX."<< endl;
   // test: O.K.
   cout <<"\nmain: Does region FRACS overlap with MATRIX? (no): ";
   cout << model.RegionIntersection( "FRACS", "MATRIX", "empty_group" );
   //            ^^^^^^^^^^^^^^^^^^             empty_region will not be formed
   cout <<"\nmain: Does region FRACS overlap with Model? (yes): ";
   cout << model.RegionIntersection( "FRACS", "Model", "fractures" );

   cout <<"\n\n\nmain: RegionDifference() are there different parts?"<< endl;
   // test: O.K.
   cout <<"\nmain: What is the difference between MATRIX and Model: ";
   if ( model.RegionDifference( "MATRIX", "Model", "difference" ) )
   //         ^^^^^^^^^^^^^^^^
     vtk_output.OutputDataToVTK( model, "difference", "fluid-pressure", "fluid pressure", 1, true );

   cout <<"\n\n\nmain: RegionSymmetricDifference() see STL doc."<< endl;
   // test: (here an underscore should be appended to name since it already exists
   cout <<"\nmain: What is the symmetric difference between MATRIX and Model: ";
   size_t number_of_different_elements = model.RegionSymmetricDifference( "MATRIX", "Model", "difference" );
   //                                          ^^^^^^^^^^^^^^^^^^^^^^^^^
   cout << number_of_different_elements << endl;

   cout <<"\n\n\nmain: Add()  adding difference and FRACS."<< endl;
   model.Region("difference").Add( model.Region("FRACS") );
   //                         ^^^
   model.RemoveRegion( "difference", false );


  // -----------------------------------------------------------------------
  //
  //  Output of information and data from regions
  //
  // -----------------------------------------------------------------------
   Point<DIM>  xyz_min, xyz_max;
   cout <<"\n\n\nmain: Bounding box coordinates of region:"<< endl;
   // test: O.K.
   model.Region("FRACS").MinMaxCoordinates( xyz_min, xyz_max );
   // test: O.K.         ^^^^^^^^^^^^^^^^^
   xyz_min.Out();
   xyz_max.Out();

   cout <<"\n\n\nmain: fluid pressure min/max: ";
   double64  gmin, gmax;
   model.Region("FRACS").MinMaxOf( "fluid pressure", gmin, gmax );
   //  test: O.K.        ^^^^^^^^
   cout << gmin <<" - "<< gmax << endl;
   cout <<"\nmain: average 'fluid pressure'"<< model.Region("MATRIX").Average("fluid pressure") << endl;
   // test: O.K.                                                      ^^^^^^^
   cout <<"\nmain: surface area and volume of region: ";
   cout << model.Region("MATRIX").SurfaceArea() <<"  "<< model.Region("MATRIX").Volume() << endl;
   // test: O.K.                  ^^^^^^^^^^^                                   ^^^^^^
   cout <<"\n\n\nmain: porosity volume integral: "<< model.Region("MATRIX").VolumeIntegral( "porosity" ,false) << endl;
   // test: O.K.                                                            ^^^^^^^^^^^^^^

   // region_flag options are COMPLETE, INTERIOR or PERIMETER
   ScalarVariable phi(PLAIN,1.);
   model.Region("FRACS").InputPropertyValue( "porosity", phi, COMPLETE ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^
   model.Region("MATRIX").InputPropertyValue( "porosity", makeScalar(phi.Flag(),0.2), COMPLETE ); // test: O.K.
   printRangeOfVariable( model, "porosity" );

   FEM_Data<ScalarVariable>  input_data;
   input_data.Reset( model.Database().StorageKey("permeability"),
                     model.Region("Model").Elements(), makeScalar(PLAIN,1.0e-12) );

   model.Region("Model").InputVariableFrom( "permeability", input_data ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "permeability" );

   model.Region("Model").ChangePropertyStatus( "permeability", INIT_GUESS, INTERIOR );
   //                    ^^^^^^^^^^^^^^^^^^^^
   model.Region("Model").ChangePropertyStatusWhere( "permeability", INIT_COND, 1.0e-14, 1.0e-12 );
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^
   model.Region("Model").CopyReplace( "porosity", "permeability" ); // test: O.K.
   //                    ^^^^^^^^^^^
   printRangeOfVariable( model, "permeability" );

   PropertyHandle<DIM>  grad_pf( model, "fluid pressure gradient", VECTOR, ELEMENT );
   model.Region("Model").CopyGradientOfProperty_A_To_B( "fluid pressure", "fluid pressure gradient" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "fluid pressure gradient" );

   PropertyHandle<DIM>  nodal_perm( model, "nodal permeability", SCALAR, NODE );
   model.Region("Model").InputPropertyValue( "nodal permeability", makeScalar(PLAIN,1.0e-12), COMPLETE ); // test: O.K.
   model.Region("Model").InterpolateNodeToElementProperty( "nodal permeability", "permeability" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "permeability" );
   printRangeOfVariable( model, "nodal permeability" );

   PropertyHandle<DIM>  cpoint_k( model, "cpoint permeability", TENSOR, ELEMENT_INTEGRATION_POINT );
   TensorVariable<DIM> ts;
   ts      = 1.0e-12;
   ts(0,0) = ts(1,1) = 2.;
   ts(2,2) = 3.;
   ts.Out();
   cpoint_k = ts;
   printRangeOfVariable( model, "cpoint permeability" );
   PropertyHandle<DIM>  elmt_k( model, "element permeability", TENSOR, ELEMENT );
   model.Region("Model").InterpolateIntegrationPointToElementProperty( "cpoint permeability", "element permeability" );  // O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "element permeability" );

   PropertyHandle<DIM>  nodal_v( model, "nodal velocity", VECTOR, NODE );
   model.Region("Model").ExtrapolateElementToNodeProperty( "velocity","nodal velocity" ); // O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "velocity" );
   printRangeOfVariable( model, "nodal velocity" );  // interpolation by distance
   model.Region("Model").ExtrapolateElementToNodeProperty( "velocity","nodal velocity", false ); // O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "nodal velocity" );

   PropertyHandle<DIM>  nodal_k( model, "nodal tensor permeability", TENSOR, NODE );
   model.Region("Model").ExtrapolateIntegrationPointToNodeProperty( "cpoint permeability", "nodal tensor permeability" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "nodal tensor permeability" );

   PropertyHandle<DIM>  volume( model, "volume", SCALAR, ELEMENT );
   model.Region("Model").AssignElementCharacteristicsTo( "volume", "volume" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "volume" );

   model.Region("Model").MoveNodeCoordinatesBy( "nodal velocity" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^
   PropertyHandle<DIM>  nd_positions( model, "node positions", VECTOR, NODE );
   model.Region("Model").AssignNodeCoordinatesTo( "node positions" ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "node positions" );

   PropertyHandle<DIM>  depth( model, "depth", SCALAR, NODE );
   model.Region("Model").AssignNodeCoordinatesTo( "depth", 'y' ); // test: O.K.
   //                    ^^^^^^^^^^^^^^^^^^^^^^^
   printRangeOfVariable( model, "depth" );

   model.FormRegionFrom( "small", element_ids ); // see above
   Region<DIM>  gref(model.Region("small"));
   cout <<"\nmain: Region 'small' contains: "<< endl;
   cout <<"\n\tnodes:             "<< gref.Nodes(); // test: O.K.
   cout <<"\n\tconstraint points: "<< gref.IntegrationPoints(); // test: O.K.
   cout <<"\n\telements:          "<< gref.Elements(); // test: O.K.
   cout <<"\n\tinterior nodes:    "<< gref.InteriorNodes(); // test: O.K.
   cout <<"\n\tinterior elements: "<< gref.InteriorElements(); // test: O.K.
   cout <<"\n\tIs empty?          "<< gref.Empty() << endl; // test: O.K.

   const Element<DIM>&  e1(*(model.Mesh().RootElement(0)));
   cout <<"\nmain: Does the model contain a certain element? "<< model.Region("Model").Contains( &e1 ) << endl; // test: O.K.
   //                                                                                  ^^^^^^^^
   cout <<"\nmain: At its boundary? "<< model.Region("Model").IsPerimeterElement( &e1 ) << endl; // test: O.K.
   //                                                         ^^^^^^^^^^^^^^^^^^

   // output of group data
   // --------------------
   gref.Out();
   // establishing a contiguous node and element numbering for the region
   gref.UpdateMemberIndexes();
   gref.OutputVariableToScreen( "fluid pressure" ); // test: O.K.
   VSet<DIM>  vset;
   gref.OutputTo( vset, true ); // test: O.K. (true = including properties)
   vset.Out();
   FEM_Data<ScalarVariable>        sc_data;
   FEM_Data<VectorVariable<DIM> >  vc_data;
   gref.OutputVariableTo( "fluid pressure", sc_data ); // test: O.K.
   sc_data.Out();
   gref.OutputVariableTo( "velocity", vc_data ); // test: O.K.
   vc_data.Out();
} // Run()

} // csmp

