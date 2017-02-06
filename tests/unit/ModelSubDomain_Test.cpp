//
//  ModelSubDomain_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "ModelSubDomain_Test.hpp"
#include "Model.h"
#include "Boundary.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"
#include "VTK_Interface.h"
#include "vset_makers.h"
//#include "ModelComparator.h"
#include "Element.h"
#include "Box.h"

using namespace std;

namespace csmp {


/**
     Questions
     - are the flags correctly assigned
     - establish clear responsibility sharing between MeshManager and AllElements region !

     Observations
     - AllElements in single domain model is non-unique, but should be unique
     
     
     TODO: so far (4/1/17), this test only checks the writing and reading of the model from file;
     add tests for the other important parts of the functionality.
 
*/
void ModelSubDomain_Test::run()
  {
      std::ostream& os = getInfoStream();
     bool test_binary_file_recovery1(true),
          test_binary_file_recovery2(true), // with boundaries
          test_binary_file_recovery3(true), // complex model with multiple regions
          verbose(true);
    
     // Test 1: subdomain storage to file and recreation in a new model
     // ---------------------------------------------------------------
     if ( test_binary_file_recovery1 )
       {
         VSet<3U>   vset;
         const bool skewed(false), isoparametric(true);
         test_Create_Prism_Hexa_VSet( os, vset, skewed );
         Model<3U>   model1( vset, isoparametric );
        
         //model1.CreateProperty( "box flag", "none", SCALAR, NODE );
         //model1.CreateProperty( "box flag element", "none", SCALAR, ELEMENT );
         //boxFlagsToVariable( model1, "box flag", "box flag element" );
         //VTK_Interface<3U>  vtk_output;
         //vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
         //vtk_output.OutputDataToVTK( model1, "element_flag", "box flag element", 0 );
        
         cerr <<"\nModelSubDomain_Test::run: original model.";
         // TODO: numbering of boundary nodes does not seem to be correct
         vset.Out(os);
         //if ( verbose ) model1.Out(os);
         //model1.Region("All Elements").Out(os);
        
         model1.OutputToBinaryFile("ModelSubDomain_Test");
        
         // testing: model2.InputFromBinaryFile("ModelSubDomain_Test");
         Model<3U>  model2( string("ModelSubDomain_Test") );
         //vtk_output.OutputDataToVTK( model2, "node_flag", "box flag", 1 );
         //vtk_output.OutputDataToVTK( model2, "element_flag", "box flag element", 1 );
        
         cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
         //if ( verbose ) model2.Out(os);
        
         _test( CompareModelSubdomains( model1.Region("All Elements"), model2.Region("All Elements"), verbose ) );
         _test( CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), verbose ) );
       }
    
    
    // Test 2: box-shaped model with boundary information (Boundary->ModelSubDomain)
    // -----------------------------------------------------------------------------
    if ( test_binary_file_recovery2 )
      {
         ANSYS_Model3D model1( "cube_flag", "CSMP-variables.txt",
                                 false, /* irregular_mesh */
                                 true,  /* binary_file */
                                 true,  /* use_regions_file */
                                 true   /* create_boundaries */ );
        
         // loop over model boundary verifying consistency between AtBoundary() and region flags
         bool left_fail(false), right_fail(false), bottom_fail(false), top_fail(false), front_fail(false), back_fail(false);
	       for ( auto it = model1.BoundariesBegin(); it != model1.BoundariesEnd(); ++it ) {
               string bname((*it).first);
               for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.PerimeterNodesBegin(); ++nit ) {
                     if ( bname == "LEFT"   and (*nit)->AtBoundary() != LEFT )   left_fail = true;
                     if ( bname == "RIGHT"  and (*nit)->AtBoundary() != RIGHT )  right_fail = true;
                     if ( bname == "TOP"    and (*nit)->AtBoundary() != TOP )    top_fail = true;
                     if ( bname == "BOTTOM" and (*nit)->AtBoundary() != BOTTOM ) bottom_fail = true;
                     if ( bname == "FRONT"  and (*nit)->AtBoundary() != FRONT )  front_fail = true;
                     if ( bname == "BACK"   and (*nit)->AtBoundary() != BACK )   back_fail = true;
                 }
            }
         _test( left_fail   == false );
         _test( right_fail  == false );
         _test( bottom_fail == false );
         _test( top_fail    == false );
         _test( front_fail  == false );
         _test( back_fail   == false );

         // mapping flags to values to test assignments
// TODO: adding properties upsets indices for the access of the flagged array
// call  UpdateParametersAndDatabase(); or  UpdateIndexReferences(); but they are private?

/* CREATE PROPERTY UPSETS flagged array variable storage offset
         model1.CreateProperty( "box flag", "none", SCALAR, NODE );
         model1.CreateProperty( "box flag element", "none", SCALAR, ELEMENT );
         boxFlagsToVariable( model1, "box flag", "box flag element" );

         VTK_Interface<3U>  vtk_output;
         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
        
         // checking whether the side boundary interior and boundary flags are identified correctly
         //const csmp::Index  prop_key(model1.Database().StorageKey("box flag"));
         Boundary<3U>&      bref(model1.Boundary("FRONT"));
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double64>(FRONT)), INTERIOR );
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double64>(REGION_BOUNDARY)), PERIMETER );

         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 1 );

         Boundary<3U>&      brefl(model1.Boundary("LEFT"));
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double64>(LEFT)), INTERIOR );
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double64>(REGION_BOUNDARY)), PERIMETER );

         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 2 );
        
         // testing by comparison with BOX_BOUNDARY flags
*/
         model1.OutputToBinaryFile("ModelSubDomain_Test2");
         Model<3U>  model2( string("ModelSubDomain_Test2") );

         _test( CompareModelSubdomains( model1.Region("All Elements"), model2.Region("All Elements"), verbose ) );
         _test( CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), verbose ) );

      }

    if ( test_binary_file_recovery3 )
      {
         ANSYS_Model3D model1( "prism_test", "CSMP-variables.txt",
                                 false, /* irregular_mesh */
                                 true,  /* binary_file */
                                 true,  /* use_regions_file */
                                 true   /* create_boundaries */ );
        
         model1.OutputToBinaryFile("ModelSubDomain_Test3");
         Model<3U>  model2( string("ModelSubDomain_Test3") );

         _test( CompareModelSubdomains( model1.Region("FRAC_VOLUMES"), model2.Region("FRAC_VOLUMES"), verbose ) );
      }


    
  } // end run



//template bool compareModelSubdomains( const ModelSubDomain<3,Element>&, const ModelSubDomain<3,Element>&, false );



} // end csmp
