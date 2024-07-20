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
#include "Element.h"
#include "Box.h"

using namespace std;

namespace csmp {


/**
     Observations
     - AllElements in single domain model is non-unique, but should be unique
      
     TODO: test wide ranging functionality of ModelSubDomain: so far (4/1/17), this test only checks the writing and reading of the model from file;
     add tests for the other important parts of the functionality.
*/
void ModelSubDomain_Test::run()
  {
     // Test 0: methods of subdomain in live model
     // ------------------------------------------
     _test( Test_EstablishNeighborConnectivity() );
  
     bool test_binary_file_recovery1(true),
          verbose(true);
    
     // Test 1: subdomain storage to file and recreation in a new model
     // ---------------------------------------------------------------
     if ( test_binary_file_recovery1 )
       {
         VSet<3U>   vset;
         const bool skewed(false);
         create_Prism_Hexa_VSet( vset, skewed );
         vset.RemoveData("element number");
         Model<3U>   model1( vset );
        
         model1.CreateProperty( "box flag", "bFn", "none", SCALAR, NODE );
         model1.CreateProperty( "box flag element", "bFe", "none", SCALAR, ELEMENT );
         boxFlagsToVariable( model1, "box flag", "box flag element" );
         if ( verbose_ ) {
              VTK_Interface<3U>  vtk_output;
              vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
              vtk_output.OutputDataToVTK( model1, "element_flag", "box flag element", 0 );
           }
           
         cerr <<"\nModelSubDomain_Test::run: original model.";
         if ( verbose_ ) model1.Out();
        
         model1.OutputToBinaryFile("ModelSubDomain_Test");
        
         // testing: model2.InputFromBinaryFile("ModelSubDomain_Test");
         Model<3U>  model2( string("ModelSubDomain_Test") );
         if ( verbose_ ) {
              VTK_Interface<3U>  vtk_output;
              vtk_output.OutputDataToVTK( model2, "node_flag", "box flag", 1 );
              vtk_output.OutputDataToVTK( model2, "element_flag", "box flag element", 1 );
           }
         cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
         if ( verbose_ ) model2.Out();
        
        // TESTING
         _test( CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), verbose ) );
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
         plist_entries1.insert( std::move(nodes) );
      }
    set<vector<size_t> > plist_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<size_t> nodes( (*it)->Nodes() );
         for ( auto i{0U}; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
//cerr <<"\n"<< (*it)->Idx() <<": ";
//out(nodes);
         plist_entries2.insert( std::move(nodes) );
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
         pfverts_entries1.insert( std::move(nbors) );
      }
    set<vector<int64_t> > pfverts_entries2;
    for ( auto it=domain2.CellsBegin(); it!=domain2.CellsEnd(); ++it ) {
         vector<int64_t> nbors( (*it)->Neighbors(),0 );
         for ( auto i{0}; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = (*it)->Neighbor(i)->Idx();
             }
         pfverts_entries2.insert( std::move(nbors) );
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
    create_Pyramid_Hexa_VSet( vset, skewed_elements );
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
    Region<3U>& domain = model.Region("Model");
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
