//
//  ModelSubDomain_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "ModelSubDomain_Test.h"
#include "Model.h"
#include "Boundary.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"
#include "VTK_Interface.h"
#include "vsetMakers.h"
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
     // Test 0: methods of subdomain in live model
     // ------------------------------------------
     _test( Test_EstablishNeighborConnectivity() );
  
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
         test_Create_Prism_Hexa_VSet( vset, skewed );
         Model<3U>   model1( vset, isoparametric );
        
         //model1.CreateProperty( "box flag", "none", SCALAR, NODE );
         //model1.CreateProperty( "box flag element", "none", SCALAR, ELEMENT );
         //boxFlagsToVariable( model1, "box flag", "box flag element" );
         //VTK_Interface<3U>  vtk_output;
         //vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
         //vtk_output.OutputDataToVTK( model1, "element_flag", "box flag element", 0 );
        
         cerr <<"\nModelSubDomain_Test::run: original model.";
         // TODO: numbering of boundary nodes does not seem to be correct
         //vset.Out();
         //if ( verbose ) model1.Out();
        
         model1.OutputToBinaryFile("ModelSubDomain_Test");
        
         // testing: model2.InputFromBinaryFile("ModelSubDomain_Test");
         Model<3U>  model2( string("ModelSubDomain_Test") );
         //vtk_output.OutputDataToVTK( model2, "node_flag", "box flag", 1 );
         //vtk_output.OutputDataToVTK( model2, "element_flag", "box flag element", 1 );
        
         cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
         //if ( verbose ) model2.Out();
        
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
               const string bname((*it).first);
               for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.PerimeterNodesBegin(); ++nit ) {
                     if ( bname == "LEFT"   and (*nit)->AtBoundary() != LEFT )   left_fail = true;
                     if ( bname == "RIGHT"  and (*nit)->AtBoundary() != RIGHT )  right_fail = true;
                     if ( bname == "TOP"    and (*nit)->AtBoundary() != TOP )    top_fail = true;
                     if ( bname == "BOTTOM" and (*nit)->AtBoundary() != BOTTOM ) bottom_fail = true;
                     if ( bname == "FRONT"  and (*nit)->AtBoundary() != FRONT )  front_fail = true;
                     if ( bname == "BACK"   and (*nit)->AtBoundary() != BACK )   back_fail = true;
                 }
            }
        
         // testing the side boundaries
         _test( left_fail   == false );
         _test( right_fail  == false );
         _test( bottom_fail == false );
         _test( top_fail    == false );
         _test( front_fail  == false );
         _test( back_fail   == false );
        
         // testing the edges and their perimeter nodes=model corners as well
         bool edge1_fail(false), edge2_fail(false), edge3_fail(false), edge4_fail(false),
              edge5_fail(false), edge6_fail(false), edge7_fail(false), edge8_fail(false),
              edge9_fail(false), edge10_fail(false), edge11_fail(false), edge12_fail(false);

	       for ( auto it = model1.BoundariesBegin(); it != model1.BoundariesEnd(); ++it ) {
               const string bname((*it).first);
               // the interior part of the edges comes first
               for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.PerimeterNodesBegin(); ++nit ) {
                     if ( bname == "EDGE1" and (*nit)->AtBoundary() != EDGE1 ) edge1_fail = true;
                     if ( bname == "EDGE2" and (*nit)->AtBoundary() != EDGE2 ) edge2_fail = true;
                     if ( bname == "EDGE3" and (*nit)->AtBoundary() != EDGE3 ) edge3_fail = true;
                     if ( bname == "EDGE4" and (*nit)->AtBoundary() != EDGE4 ) edge4_fail = true;
                     if ( bname == "EDGE5" and (*nit)->AtBoundary() != EDGE5 ) edge5_fail = true;
                     if ( bname == "EDGE6" and (*nit)->AtBoundary() != EDGE6 ) edge6_fail = true;
                     if ( bname == "EDGE7" and (*nit)->AtBoundary() != EDGE7 ) edge7_fail = true;
                     if ( bname == "EDGE8" and (*nit)->AtBoundary() != EDGE8 ) edge8_fail = true;
                     if ( bname == "EDGE9" and (*nit)->AtBoundary() != EDGE9 ) edge9_fail = true;
                     if ( bname == "EDGE10" and (*nit)->AtBoundary() != EDGE10 ) edge10_fail = true;
                     if ( bname == "EDGE11" and (*nit)->AtBoundary() != EDGE11 ) edge11_fail = true;
                     if ( bname == "EDGE12" and (*nit)->AtBoundary() != EDGE12 ) edge12_fail = true;
                 }
            }

         // testing the side boundaries
         _test( edge1_fail == false );
         _test( edge2_fail == false );
         _test( edge3_fail == false );
         _test( edge4_fail == false );
         _test( edge5_fail == false );
         _test( edge6_fail == false );
         _test( edge7_fail == false );
         _test( edge8_fail == false );
         _test( edge9_fail == false );
         _test( edge10_fail == false );
         _test( edge11_fail == false );
         _test( edge12_fail == false );
        
         // verifying that the perimeter nodes of the edges contain the correct corner nodes
         set<BOX_BOUNDARY> bnodes;
         // testing edge1
         const Boundary<3U> edge1(model1.Boundary("EDGE1"));
         _test( edge1.PerimeterNodes() == 2 );
         for ( auto nit=edge1.PerimeterNodesBegin(); nit!=edge1.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR2) > 0 );
         bnodes.clear();
         // testing edge2
         const Boundary<3U> edge2(model1.Boundary("EDGE2"));
         _test( edge2.PerimeterNodes() == 2 );
         for ( auto nit=edge2.PerimeterNodesBegin(); nit!=edge2.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR2) > 0 );
         bnodes.clear();
         // testing edge3
         const Boundary<3U> edge3(model1.Boundary("EDGE3"));
         _test( edge3.PerimeterNodes() == 2 );
         for ( auto nit=edge3.PerimeterNodesBegin(); nit!=edge3.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR4) > 0 );
         bnodes.clear();
         // testing edge4
         const Boundary<3U> edge4(model1.Boundary("EDGE4"));
         _test( edge4.PerimeterNodes() == 2 );
         for ( auto nit=edge4.PerimeterNodesBegin(); nit!=edge4.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR4) > 0 );
         bnodes.clear();
         // testing edge5
         const Boundary<3U> edge5(model1.Boundary("EDGE5"));
         _test( edge5.PerimeterNodes() == 2 );
         for ( auto nit=edge5.PerimeterNodesBegin(); nit!=edge5.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR5) > 0 );
         bnodes.clear();
         // testing edge6
         const Boundary<3U> edge6(model1.Boundary("EDGE6"));
         _test( edge6.PerimeterNodes() == 2 );
         for ( auto nit=edge6.PerimeterNodesBegin(); nit!=edge6.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR2) > 0 );
         _test( bnodes.count(CNR6) > 0 );
         bnodes.clear();
         // testing edge7
         const Boundary<3U> edge7(model1.Boundary("EDGE7"));
         _test( edge7.PerimeterNodes() == 2 );
         for ( auto nit=edge7.PerimeterNodesBegin(); nit!=edge7.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR7) > 0 );
         bnodes.clear();
         // testing edge1
         const Boundary<3U> edge8(model1.Boundary("EDGE8"));
         _test( edge8.PerimeterNodes() == 2 );
         for ( auto nit=edge8.PerimeterNodesBegin(); nit!=edge8.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR4) > 0 );
         _test( bnodes.count(CNR8) > 0 );
         bnodes.clear();
         // testing edge9
         const Boundary<3U> edge9(model1.Boundary("EDGE9"));
         _test( edge9.PerimeterNodes() == 2 );
         for ( auto nit=edge9.PerimeterNodesBegin(); nit!=edge9.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR5) > 0 );
         _test( bnodes.count(CNR6) > 0 );
         bnodes.clear();
         // testing edge10
         const Boundary<3U> edge10(model1.Boundary("EDGE10"));
         _test( edge10.PerimeterNodes() == 2 );
         for ( auto nit=edge10.PerimeterNodesBegin(); nit!=edge10.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR6) > 0 );
         _test( bnodes.count(CNR7) > 0 );
         bnodes.clear();
         // testing edge11
         const Boundary<3U> edge11(model1.Boundary("EDGE11"));
         _test( edge11.PerimeterNodes() == 2 );
         for ( auto nit=edge11.PerimeterNodesBegin(); nit!=edge11.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR7) > 0 );
         _test( bnodes.count(CNR8) > 0 );
         bnodes.clear();
         // testing edge12
         const Boundary<3U> edge12(model1.Boundary("EDGE12"));
         _test( edge1.PerimeterNodes() == 2 );
         for ( auto nit=edge12.PerimeterNodesBegin(); nit!=edge12.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR8) > 0 );
         _test( bnodes.count(CNR5) > 0 );
         bnodes.clear();

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

         // saving the model to csmp native binary file format
         model1.OutputToBinaryFile("ModelSubDomain_Test2");
        
         // bringing the model back from binary file
         Model<3U>  model2( string("ModelSubDomain_Test2") );

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

         _test( CompareModelSubdomains( model1.Region("FRAC_VOLUMES"),
                                        model2.Region("FRAC_VOLUMES"), verbose ) );
      }
    
  } // end run



/**
    Using vsets from 'vset_makers' as input data, this tests weither     ModelSubDoman::EstablishNeighborConnectivity()  recreates the correct neighbor connectivity
    
        @author SKM 5/12/20

*/
bool ModelSubDomain_Test::Test_EstablishNeighborConnectivity() 
 {
    // 0. creating the test model
    // --------------------------
    string varFileName("ModelSubDomain_Test-variables.txt");
    const bool   skewed_elements(false); // otherwise model is not a box anymore
    VSet<3U>     vset;
    test_Create_Pyramid_Hexa_VSet( vset, skewed_elements );
    // adding 'node number' as a variable
    PropertyData node_nums( NODE, SCALAR, 3U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );
    // adding 'element number' as a variable
    PropertyData elmt_nums( ELEMENT, SCALAR, 3U );
    node_nums.Reserve( vset.Elements() );
    for ( size_t i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );

    Model<3U> model( vset, varFileName.c_str(), true );

    // 1. recreating the neighbor connectivity and comparing
    // -----------------------------------------------------
    const csmp::Index eid_key(model.Database().StorageKey("element number"));
    Region<3U> domain = model.Region("Model");
    domain.EstablishNeighborConnectivity();
    
    bool no_mismatch(true);
    for ( vector<Element<3U>*>::const_iterator it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it )
      for ( size_t i=0U; i<(*it)->Neighbors(); ++i ) {
           if ( (*it)->Neighbor(i) != nullptr ) {
                const size_t elmt_id = static_cast<size_t>((*it)->Read( eid_key ));
                const size_t nbor_id = static_cast<size_t>((*it)->Neighbor(i)->Read( eid_key ));
               _test( nbor_id == vset.Pfvert( elmt_id, i ) );
                if ( nbor_id != vset.Pfvert( elmt_id, i ) ) {
                     cerr <<"\nelmt "<< elmt_id <<":"<< i <<": vset vs. reconstructed neighbor: ";
                     cerr << vset.Pfvert( elmt_id, i ) <<" vs. "<< nbor_id;
                     no_mismatch = false;
                  }
             }
           if ( no_mismatch == false ) cerr << endl;
        }

    return no_mismatch;
    
 } // end Test_EstablishNeighborConnectivity
     




} // end csmp
