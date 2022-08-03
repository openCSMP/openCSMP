//
//  ModelSubDomain_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "ModelSubDomain_Test.h"
#include "Model.h"
#include "Region.h"
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
         const bool skewed(false);
         test_Create_Prism_Hexa_VSet( vset, skewed );
         Model<3U>   model1( vset );
        
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
    
// TODO: the following test (using model 'cube_flag') should be moved to BoundaryInterface_Test
    // Test 2: box-shaped model with boundary information (Boundary->ModelSubDomain)
    // -----------------------------------------------------------------------------
    if ( test_binary_file_recovery2 )
      {
         ANSYS_Model3D model1( "cube_flag", "CSMP-variables.txt",
                                 false, /* irregular_mesh */
                                 true   /* binary_file */
                             );
        
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

         // E.P Need to create EDGE boundaries in the first place for this to Work!!! NEED TO CALL METHOD BELOW FIRST
         // model1.EstablishEdgeBoundariesOfBoxShapedModel();


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
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(FRONT)), INTERIOR );
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(REGION_BOUNDARY)), PERIMETER );

         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 1 );

         Boundary<3U>&      brefl(model1.Boundary("LEFT"));
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(LEFT)), INTERIOR );
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(REGION_BOUNDARY)), PERIMETER );

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
                                 true   /* binary_file */
                             );
        
         model1.OutputToBinaryFile("ModelSubDomain_Test3");
         Model<3U>  model2( string("ModelSubDomain_Test3") );

         _test( CompareModelSubdomains( model1.Region("FRAC_VOLUMES"),
                                        model2.Region("FRAC_VOLUMES"), verbose ) );
      }
    
  } // end run





/** compares node locations and connectivity
    - uses point locations for node comparison

*/
template<uint32_t dim,template<uint32_t> class simplicial_complex>
bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>& domain1,
                                                  const ModelSubDomain<dim,simplicial_complex>& domain2,
                                                  bool verbose )
 {
    // domain names
    if ( domain1.Name() != domain2.Name() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: names mismatch: "<< domain1.Name() <<" vs. "<< domain2.Name() <<"\n";
         //return false;
      }
    _test( domain1.Name() == domain2.Name() );
   
    // 1. nodes
    // --------
    // number of nodes
    if ( domain1.Nodes() != domain2.Nodes() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: node numbers do not match: "<< domain1.Nodes() <<" vs "<< domain2.Nodes() <<"\n";
         //return false;
      }
    _test( domain1.Nodes() == domain2.Nodes() );
   
    // interior nodes
    if ( domain1.InteriorNodes() != domain2.InteriorNodes() ) {
         if ( verbose ) {
              cerr <<"\ncompareModelSubdomains: number of interior nodes does not match: "<< domain1.InteriorNodes() <<" vs "<< domain2.InteriorNodes() <<"\n";
              cerr <<"\tinterior nodes domain 1 vs domain 2:\n";
              for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
                cerr << (*nit)->Idx() <<" ";
              cerr << endl;
              for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
                cerr << (*nit)->Idx() <<" ";
              cerr << endl;
           }
         //return false;
      }
    _test( domain1.InteriorNodes() == domain2.InteriorNodes() );
   
    // node coordinates (sorted ranges of points have to be created)
    // ranges for domain 1
    set<Point<dim> >  interior_points1;
    for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    set<Point<dim> >  perimeter_points1;
    for ( auto nit=domain1.PerimeterNodesBegin(); nit!=domain1.NodesEnd(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    // ranges for domain 2
    set<Point<dim> >  interior_points2;
    for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    set<Point<dim> >  perimeter_points2;
    for ( auto nit=domain2.PerimeterNodesBegin(); nit!=domain2.NodesEnd(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    // comparisons
    if ( interior_points1 != interior_points2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: locations of interior points do not match.\n";
         //return false;
      }
    _test( interior_points1 == interior_points2 );

    if ( perimeter_points1 != perimeter_points2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: locations of perimeter points do not match.\n";
         //return false;
      }
    _test( perimeter_points1 == perimeter_points2 );
   
    interior_points1.clear();
    interior_points2.clear();
    perimeter_points1.clear();
    perimeter_points2.clear();

    // node boundary flags
    // (creating maps where the nodes are ordered by their location because
    //  they do not necessarily have the same number)
    map<Point<dim>,BOX_BOUNDARY>  boundary_flags1;
    for ( auto nit=domain1.NodesBegin(); nit!=domain1.NodesEnd(); ++nit )
      boundary_flags1.insert( make_pair( (*nit)->Coordinate(), (*nit)->AtBoundary() ) );
    //
    map<Point<dim>,BOX_BOUNDARY>  boundary_flags2;
    for ( auto nit=domain2.NodesBegin(); nit!=domain2.NodesEnd(); ++nit )
      boundary_flags1.insert( make_pair( (*nit)->Coordinate(), (*nit)->AtBoundary() ) );
    //
    // may fail due to tolerance differences in point classification _test( boundary_flags1.size() == boundary_flags2.size() );
    int node_flag_mismatches(0);
    auto nit2=boundary_flags2.begin();
    for ( auto nit1=boundary_flags1.begin(); (nit1!=boundary_flags1.end() && nit2!=boundary_flags2.end()); ++nit1, ++nit2 )
      if ( (*nit1).second != (*nit2).second ) {
            cerr <<"\n\t"<< parseBoundary((*nit1).second) <<" vs "<< parseBoundary((*nit2).second);
            node_flag_mismatches++;
        }
    if ( node_flag_mismatches > 0 ) {
          cerr <<"\ncompareModelSubdomains: the BOX_BOUNDARY flags of the domains do not match.\n";
          //return false;
      }
    _test( node_flag_mismatches == 0 );
    boundary_flags1.clear();
    boundary_flags2.clear();

   
    // 2. elements, i.e. connectivity
    // ------------------------------
    // number of elements
    if ( domain1.Cells() != domain2.Cells() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: number of elements does not match.\n";
         //return false;
      }
    _test( domain1.Cells() == domain2.Cells() );
   
    // number of interior elements
    if ( domain1.PerimeterCells() != domain2.PerimeterCells() ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: mismatch in number of perimeter elements.\n";
         //return false;
      }
    _test( domain1.PerimeterCells() == domain2.PerimeterCells() );
   
    // are the interior elements the same ?
    set<size_t>  interior_elmts1, interior_elmts2;
    for ( auto it=domain1.CellsBegin(); it!=domain1.PerimeterCellsBegin(); ++it ) interior_elmts1.insert( (*it)->Idx() );
    for ( auto it=domain2.CellsBegin(); it!=domain2.PerimeterCellsBegin(); ++it ) interior_elmts2.insert( (*it)->Idx() );
    _test( interior_elmts1 == interior_elmts2 );
   
    // element connectivity (not assuming that elements are in same order)
    set<vector<size_t> > plist_entries1;
    for ( auto it=domain1.CellsBegin(); it!=domain1.CellsEnd(); ++it ) {
         vector<size_t> nodes( (*it)->Nodes() );
         for ( auto i{0U}; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
         plist_entries1.insert( move(nodes) );
      }
    set<vector<size_t> > plist_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<size_t> nodes( (*it)->Nodes() );
         for ( auto i{0U}; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
//cerr <<"\n"<< (*it)->Idx() <<": ";
//out(nodes);
         plist_entries2.insert( move(nodes) );
      }
    // comparing plists
    bool plists_are_the_same(true);
    if ( plist_entries1 != plist_entries2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: (plist) mismatch in element-to-node connectivity.\n";
         plists_are_the_same = false;
         //return false;
      }
    _test( plists_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // comparing the element neighbor connectivity
    set<vector<int64_t> > pfverts_entries1;
    for ( auto it=domain1.CellsBegin(); it!=domain1.CellsEnd(); ++it ) {
         vector<int64_t> nbors( (*it)->Neighbors(),0 );
         for ( auto i{0}; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = (*it)->Neighbor(i)->Idx();
             }
         pfverts_entries1.insert( move(nbors) );
      }
    set<vector<int64_t> > pfverts_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<int64_t> nbors( (*it)->Neighbors(),0 );
         for ( auto i{0}; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = (*it)->Neighbor(i)->Idx();
             }
         pfverts_entries2.insert( move(nbors) );
      }
    // comparing neighbor connectivity lists
    bool pfverts_are_the_same(true);
    if ( pfverts_entries1 != pfverts_entries2 ) {
         if ( verbose )
           cerr <<"\ncompareModelSubdomains: (pfverts) mismatch in element to neighbor-element connectivity.\n";
         pfverts_are_the_same = false;
         //return false;
      }
    _test( pfverts_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // all comparisons passed - domains have same nodes and elements.
    return true;
   
 } // end compareModelSubdomains

template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<3U,Element>&,
                                                           const ModelSubDomain<3U,Element>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<2U,Element>&,
                                                           const ModelSubDomain<2U,Element>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<1U,Element>&,
                                                           const ModelSubDomain<1U,Element>&, bool );

template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<3U,Face>&,
                                                           const ModelSubDomain<3U,Face>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<2U,Face>&,
                                                           const ModelSubDomain<2U,Face>&, bool );
template bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<1U,Face>&,
                                                           const ModelSubDomain<1U,Face>&, bool );


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
    for ( auto i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );
    // adding 'element number' as a variable
    PropertyData elmt_nums( ELEMENT, SCALAR, 3U );
    node_nums.Reserve( vset.Elements() );
    for ( auto i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );

    Model<3U> model( vset, varFileName.c_str() );

    // 1. recreating the neighbor connectivity and comparing
    // -----------------------------------------------------
    const csmp::Index eid_key(model.Database().StorageKey("element number"));
    Region<3U> domain = model.Region("Model");
    model.Mesh().UpdateConnectivity();
    
    bool no_mismatch(true);
    for ( vector<Element<3U>*>::const_iterator it=domain.CellsBegin(); it!=domain.CellsEnd(); ++it )
      for ( auto i{0}; i<(*it)->Neighbors(); ++i ) {
           if ( (*it)->Neighbor(i) != nullptr ) {
                const size_t elmt_id = static_cast<uint32_t>((*it)->Read( eid_key ));
                const size_t nbor_id = static_cast<uint32_t>((*it)->Neighbor(i)->Read( eid_key ));
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
