//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include <ranges>
#include "CSMP_definitions.h"
#include "MeshManager_Test.h"
#include "MeshManagementUtilities.h"
#include "ErrorHandler.h"
#include "vsetMakers.h"
#include "VTU_Interface.h"
#include "Region.h"
#include "Element.h"
#include "compareFloats.h"

#include "ANSYS_Model3D.h"
#include "ANSYS_Interface.h"
#include "ModelTopology.h"

#include "IsoparametricLinearPyramid.h"
#include "VTK_Interface.h"
#include "ModelTopology.h"
#include "NodeManifold.h"

// #define CSMP_MESH_MANAGER_TEST_DEBUG

using namespace std;

namespace csmp {

/**
custom models
uses "CSMP-1phase-variables.txt"
*/
MeshManager_Test::MeshManager_Test()
{
}


/**
   Construct raw model using only the regions specified in -regions.txt file and eliminating lower-dim elements from these regions
*/
void MeshManager_Test::ConstructANSYS_ModelWithoutModification( const char* model_name )
 {
    const string mesh_file_set{ model_name };
    VSet<3U>  vset;
    bool isoparametric_elements( true );

    ModelTopology   mesh_topology( isoparametric_elements );
    ANSYS_Interface mesh_interface( isoparametric_elements );

    // 1. reading the mesh from ANSYS-CSMP-input files
    const bool binary_input_file{true};
    mesh_interface.Read_ANSYS_Mesh( mesh_file_set, vset, mesh_topology, binary_input_file );
    
    // 2. ATTENTION (comment from SKM): Since ANSYS does not output the neighbour connectivity correctly,
    //    the 'pfverts' neighbor container is zapped here so that VData does not think anymore that it has neighbor connectivity
    //    later on this connectivity will be recreated inside of the Model where suitable machinery exists.
    vset.RemovePfverts();

    // 3. keeping only the highest-dimensional elements in regions which have been specified in regions file
     {
        mesh_topology.RemoveLowDimCellsFromDomains( vset );
        mesh_topology.ReduceToDomains( model_name );
        //    reducing the element data to the desired elements specified in the topology object
        //    if the element numbers in the two are different.
        if ( mesh_topology.Cells() != vset.Elements() + vset.Faces() + vset.Interfaces() ) {
            map<size_t,size_t>  old_and_new_elmtids;
            mesh_topology.CreateNewCellNumbers( old_and_new_elmtids );
            vset.ReduceTo( old_and_new_elmtids );
            old_and_new_elmtids.clear();
          }
      }

   // constructor that ignores regions, but recreates box boundary flags for box-shaped tetrahedral models
   delete model3D_ptr_; // any earlier 3D model
   model3D_ptr_ = new Model<3U>( vset, "MeshManager_Test-variables.txt" );
   
 } // end



template<uint32_t dim>
bool TestNodeNeighborConnectivity( const Model<dim>& model )
 {
    const Region<dim>&  model_domain(model.Region("Model"));
    
    vector<set<size_t>> node_neighbors;
    nodeNeighbors( model_domain, node_neighbors );
    
    // determining typical number of node neighbors in mesh
    size_t n_neighbors{0};
    for ( auto nit : node_neighbors ) n_neighbors += nit.size();
    cout <<"\n\taverage number of neighbors per node: "<< n_neighbors / node_neighbors.size();
}




 
 
 
bool MeshManager_Test::CheckConnectivityOfModel3D( Model<3>& model )
 {
    set<Element<3>*> elements3;
    MeshManager<3>&  mesh(model.Mesh());
    cout << "\n\n\nMeshManager_Test::CheckConnectivityOfModel3D: '"<< model.Name() <<"'";
    cout <<"\nExamining the connectivity of  Nodes: " << mesh.Nodes() << "\n";
    set<Node<3>*>  contiguous_set_of_nodes;
    Node<3>*       nptr = &(*mesh.NodesBegin()); // from colony!
    cout << "\nInterconnected nodes: " << findInterconnectedNodeCluster( nptr, contiguous_set_of_nodes ) << "\n";
    cout << "\nExamining the connectivity of  Elements: " << mesh.Elements() << "\n";
    // checking the neighbor connectivity
#ifdef CSMP_MESH_MANAGER_TEST_DEBUG
    integrityCheck<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() );
#endif
    // checking whether the model is contiguous
    set<Element<3>*> contiguous_subset_of_cells;
    findContiguousMeshPatch( &(*mesh.ElementsBegin()), contiguous_subset_of_cells );
    if (  contiguous_subset_of_cells.size() == mesh.Elements() ) {
         cout <<"\nMeshManager_Test::CheckConnectivityOfModel3D: model is contiguous.";
      }
    else {
         map<string,vector<Element<3U>*> > elmt_map3;
         _test( findContiguousMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), elmt_map3 ) > 1 );
         cout <<"\nMeshManager_Test::CheckConnectivityOfModel3D: the model is discontiguous and consists of the mesh patches:";
         for ( auto& it : elmt_map3 ) {
              cout <<"\n\t\t'"<< it.first <<"': "<< it.second.size() <<" elements.";
           }
      }
    if ( mesh.Faces() > 0 ) {
        cout << "\n\nExamining the connectivity of  Faces: " << mesh.Faces() << "\n";
        map<string,vector<Face<3U>*> >  face_map3;
        _test( findContiguousMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) >= 1 );
        cout << "\nInterconnected faces: " << (*face_map3.begin()).second.size() << "\n";
      }
    if ( mesh.InterFaces() > 0 ) {
       cout << "\nExamining the connectivity of  Interfaces: " << mesh.InterFaces() << "\n";
        map<string,vector<InterFace<3U>*> >  iface_map3;
        _test( findContiguousMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) >= 1 );
        cout << "\nInterconnected Interfaces: " << (*iface_map3.begin()).second.size() << "\n";
      }
    
    // Junchul's test (goal is zero errors)
    _test( mesh.CheckElementConnectivity() == 0 );

    return true;
      
 } // end CheckConnectivityOfModel3D
 
 
 
 
 





/**
THE MASTER TEST FUNCTION
*/
void MeshManager_Test::run()
{
  // basics
   {
      cout << "\n----------------------------";
      cout << "\nMeshManager_Test::TestBasics";
      cout << "\n----------------------------";
      // ==========
      TestBasics();
      // ==========
   }

  // testing hex element consistency for meshes from ANSYS
  _test(TestNeigbourVersusFaceConsistency());
  
  _test( Test_BuildConnectivity() );
  _test( Test_UpdateConnectivity() );

  _test(Test_BuiltElementConnectivity2D()); // OK
  _test(Test_BuiltElementConnectivity3D()); // OK

	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::TestSavedElementFaceInterfaceModel2D()";
	cout << "\n(elements,faces,interfaces,regions,boundaries, split boundaries)";
	cout << "\n----------------------------------------------------";
  _test( TestSavedElementFaceInterfaceModel2D() );

  // tests whether traversal works for contiguous model
  _test( Test_MeshTraversal3D(/* Pyramid_Hexa_VSet */) );
  
  Test_DeleteNodesAndRepairNodeConnnectivity();
  
  // building more complex 'FracBox' model with Boundaries and lower-dimensional elements for further testing
  VSet<3U>      vset;
  ModelTopology topology;
  string        var_file("CSMP-variables.txt");
  const bool    regions_to_boundaries{true};
  create_FracBox( topology, vset );
  
  Model<3>      model( topology, vset, var_file.c_str(), regions_to_boundaries );
  
	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::CheckConnectivityOfModel3D";
	cout << "\n----------------------------------------------------";
  _test( CheckConnectivityOfModel3D(model) );

	cout << "\n------------------------------------------------";
	cout << "\nMeshManager_Test::TestEntityNumberingFunction";
	cout << "\n------------------------------------------------";
	_test(TestEntityNumberingFunction( model ) ); // originally using 'prism_test'

// MESH MODIFICATION AND REMESHING

	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::TestElementDeletionAndInsertion";
	cout << "\n----------------------------------------------------";
	_test(TestCellDeletionAndInsertion());

	cout << "\n-------------------------------------------------";
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion";
	cout << "\n-------------------------------------------------";
	_test(TestFaceDeletionAndInsertion());

	cout << "\n------------------------------------------------------";
	cout << "\nMeshManager_Test::TestInterFaceDeletionAndInsertion";
	cout << "\n------------------------------------------------------";
	_test(TestInterFaceDeletionAndInsertion());

	_test(TestEraseAllPrimitives());
 
  cout << endl;

} // end run




bool MeshManager_Test::Test_detachNeighborsFrom()
 {
    VSet<3>       vset;
    ModelTopology topology;
//    create_RubikCube( vset );
    create_FracBox( topology, vset );
    Model<3> model( vset, "MeshManager_Test-variables.txt" );
    model.Name("FracBox");
    
    // element 13 (14) sits in the center of the cube, having all neighbors
    Region<3U>& model_domain = model.Region("Model");
    size_t      errors{0ul};
 
    // finding element that has all its nbors, recording it and its neigbhors
    size_t               e_with_nbors{0ul};
    vector<Element<3U>*> e_nbors;
    for ( const auto& it : model_domain.CellVector() ) {
         if ( it->IsVolume() && it->ConnectedNeighbors() == it->Neighbors() ) {
              e_nbors.assign( it->NeighborsBegin(), it->NeighborsEnd() );
              break;
           }
         e_with_nbors++;
      }
    
    // recording the number of neighbors of all surrounding elements
    vector<int> n_nbors( model_domain.E(e_with_nbors)->Neighbors(), IRREGULAR_OUTSIDE );
    for ( uint32_t i{0u}; i<n_nbors.size(); i++ )
      n_nbors[i] = model_domain.E(e_with_nbors)->Neighbor(i)->ConnectedNeighbors();
    
    model.Mesh().DetachNeighborsFrom( &(*model_domain.E(e_with_nbors)) );
    //           ----------------------------------------------------

    // testing
    for ( uint32_t i{0u}; i<e_nbors.size(); i++ ) {
        _test( e_nbors[i]->ConnectedNeighbors() == n_nbors[i] - 1U );
        if ( !(e_nbors[i]->ConnectedNeighbors() == n_nbors[i] - 1U) ) errors++;
      }

    return ( errors == 0 );
    
 } // end Test_detachNeighborsFrom




// for a model with a splitboundary
// NOTE: 14/8/2024 - method not used anywhere yet
 void MeshManager_Test::Test_DeleteNodesAndRepairNodeConnnectivity()
  {
      VSet<2> vset;
      // model with boundaries and split boundary
      ModelTopology topology = create_BoundarySplitBoundaryPatch( vset );
      const bool    bds_from_regions{ false };
      Model<2> model( topology, vset, "MeshManager_Test-variables.txt", bds_from_regions );
      model.Name("BoundarySplitBoundaryPatch");
      Region<2> mdomain = model.Region("Model");
      auto n_original_nodes = vset.Vertices();
      
      // 1. getting some nodes to delete some of which are manifolds
      // -----------------------------------------------------------
      vector<Node<2>*> nodes;
      // nodes 13-17 from original mesh, containing split boundary and intersection point
      set<size_t> search_nodes{ 13, 14, 15, 16, 17 };
      const csmp::Index n_key = model.Database().StorageKey("node number");
      for ( const auto& node : mdomain.NodeVector() ) {
           size_t node_idx = static_cast<size_t>(node->Read(n_key));
           if ( search_nodes.find(node_idx) != search_nodes.end() )
             nodes.push_back( node );
        }
      _test( nodes.size() == search_nodes.size() );
      
      /*
         This should differentiate the treatment of nodes delelted versus their immediate neighbors
      */
      size_t deletions = model.Mesh().DeleteNodesAndRepairNodeConnnectivity( nodes.begin(), nodes.end() );
      _test( deletions == 0 ); // because they are all still interconnected
	    _test( model.Mesh().Nodes() == n_original_nodes );
  
  
      // 2. Now we chop away those elements that contain the nodes and then try deletion again
      //    Elements:  from original mesh, containing split boundary and intersection point
      // -----------------------------------------------------------
      set<size_t> search_elmts{ 6, 7, 8, 9, 10, 11 };
      const csmp::Index e_key = model.Database().StorageKey("element number");
      for ( auto& elmt : mdomain.CellVector() ) {
           size_t elmt_idx = static_cast<size_t>(elmt->Read(e_key));
           if ( search_elmts.find(elmt_idx) != search_elmts.end() )
             model.Mesh().Delete( elmt );
        }
        
      // and we try to delete the nodes again
      deletions = model.Mesh().DeleteNodesAndRepairNodeConnnectivity( nodes.begin(), nodes.end() );
      _test( deletions > 0 );
      _test( nodes.size() == search_nodes.size() );

  } // end Test_DeleteNodesAndRepairConnnectivity




/**
    Checks that numbers of elements etc. in mesh manager do indeed reflect those of input model
    Tests: Nodes(), Elements(), Faces(), InterFaces(), HybridElementMesh(), IsContiguous(), OutputMeshTo(vset) and node-nbor connectivity;
    detachNeighborsFrom()
*/
void MeshManager_Test::TestBasics()
{
	// vsetMakers: create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch");
  // bounding box
  Point<3> xyz_min, xyz_max;
  model.MinMaxCoordinates( xyz_min, xyz_max );

  const MeshManager<3U>& mesh = model.Mesh();

  // returns true if the mesh consists of multiple element types
  _test(mesh.HybridElementMesh() == true);
  
  // see whether manager correctly detects model that consists of disconnected mesh patches
  _test( Test_IsContiguous() );

  // counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
  set<size_t> test_set;
  _test( detectElementsWithAllNodesOnBoundary( mesh, test_set ) == 0U );

  // returns number of nodes=vertices in the current mesh
  _test(mesh.Nodes() == 65);

  // returns number of elements in the current mesh
  _test(mesh.Elements() == 32);

  // returns number of Faces=lower-dimensional elements in current mesh
  _test(mesh.Faces() == 0);

  // returns number of interfaces=faces with multiplicated nodes
  _test(mesh.InterFaces() == 0);
  
  _test( Test_detachNeighborsFrom() );
  
  // test model model volume (based on bounding box)
  const Region<3>&  model_domain(model.Region("Model"));
  const double volume = (xyz_max[0]-xyz_min[0]) * (xyz_max[1]-xyz_min[1]) *  (xyz_max[2]-xyz_min[2]);
  _test( approximatelyEqual(model_domain.Volume(),volume) );
  
  // does this also work correctly in 2D?
  VSet<2> vset2D;
  create_MeshPatchWithLineElements_VSet( vset2D );
  Model<2> model2D( vset2D, "CSMP-variables.txt" ); // needs "element number" and "node number"
  model2D.Name("MeshPatchWithLineElements");
  const Region<2>&  model_domain2D(model2D.Region("Model"));
  const double area{ 6. * 5. };
  _test( approximatelyEqual(model_domain2D.Volume(),area) );
  
  // do all nodes have the expected neighbors?
  vector<set<size_t>>  node_neighbors;
  nodeNeighbors( model_domain, node_neighbors );
  
  // testing OutputMeshTo(vset) and the node-neighbor connectivity from the VSet
  VSet<2> test_vset;
  const bool get_indices_from_stored_variables{true};
  model2D.OutputMeshTo( test_vset, get_indices_from_stored_variables );
  vector<set<size_t>> pnode1, pnode2;
  vset2D.EstablishNodeNeighborConnectivity( pnode1 );
  test_vset.EstablishNodeNeighborConnectivity( pnode2 );
  _test( pnode1.size() == pnode2.size() );
  for ( size_t i{0ul}; i<pnode1.size(); ++i )
    _test( pnode1[i].size() == pnode2[i].size() );
    
  // testing HasNodeManifolds()
  _test( !model2D.Mesh().HasNodeManifolds() );
  // example mesh with manifolds
  {
    VSet<2> vset_with_manifolds;
    ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset_with_manifolds );
    Model<2> split_model( topo, vset_with_manifolds, "MeshManager_Test-variables.txt", false );
    split_model.Name("BoundarySplitBoundaryPatch");
    // using vsetMaker's SPLIT22_BASIC for testing
    _test( split_model.Mesh().HasNodeManifolds() );
    _test( distance( split_model.Mesh().NodeManifoldsBegin(), split_model.Mesh().NodeManifoldsEnd() ) == 6 );
  }

} // end TestBasics






/**
      For discontiguous model 'DykePartiallySplit' and contiguous model '' tests whether contiguity is correctly diagonosed
*/
bool MeshManager_Test::Test_IsContiguous()
 {
     bool correct_diagnosis{true};
     
     // 3D discontiguous model DykeAllLayersSplit (1785 nodes) vs. DykeAllLayersSplit (1302 nodes)
     {
       const string model_name{"DykeAllLayersSplit"};
       // ANSYS_Model3D model( model_name.c_str(), model_name.c_str(), "MeshManager_Test-variables.txt", true );
       // IMPORTANT - model needs to be build differenty to avoid any modification before testing
       ConstructANSYS_ModelWithoutModification( model_name.c_str() );
        _test( model3D_ptr_->Mesh().IsContiguous() == false );
        if ( model3D_ptr_->Mesh().IsContiguous() ) correct_diagnosis = false;
     }
     // 2D discontiguous model
     {
        VSet<2> vset;
        create_Disconnected2D_VSet( vset );
        Model<2> model( vset, "CSMP-variables.txt" );
        _test( model.Mesh().IsContiguous() == false );
        // diagnostics for debugging
        if ( model.Mesh().IsContiguous() ) {
             Region<2> reg = model.Region("Model");
             for ( const auto& nit : reg.NodeVector() ) printNeighbors( nit );
             correct_diagnosis = false;
          }
     }
     // 2D model which should be contiguous
     {
        VSet<2> vset;
        create_MeshPatchWithLineElements_VSet( vset );
        Model<2> model( vset, "CSMP-variables.txt" );
        _test( model.Mesh().IsContiguous() == true );
        if ( !model.Mesh().IsContiguous() ) correct_diagnosis = false;
     }
     
     return correct_diagnosis;
 }



/**
       2D test case for rectangular model of 2 domains separated by SplitBoundary.
       Extra partitally penetrating split boundary crosses dividing SplitBoundary.
       
       @author SKM
       @date 19/3/22
*/
bool MeshManager_Test::TestSavedElementFaceInterfaceModel2D()
 {
    VSet<2U>       vset;
    ModelTopology  topo = create_BoundarySplitBoundaryPatch( vset );
    const bool     treat_all_domains_as_regions{false};
    Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", treat_all_domains_as_regions );
    
    // testing that the model has the right area (also checks element orientations)
    const Region<2>&  model_domain(model.Region("Model"));
    const double analytic_model_area{ 4.5 * 7. };
    const double computed_model_area{ model_domain.Volume() };
    const double tolerance{ analytic_model_area * numeric_limits<double>::epsilon() };
    _test( approximatelyEqual(analytic_model_area,computed_model_area,tolerance) );
    
    // can such a model be output to VTU?
    if ( verbose_ ) {
        VTU_Interface<2>  vtu_out( model );
        printRangeOfVariable( model, "permeability" );
        vtu_out.OutputDataToVTU( "SPLIT22_BASIC", "permeability", "Model", 0 );
      }
      
    // saving model to disk and bringing it back
    model.OutputToBinaryFile( model.Name() );
    set<string>  subset_variables; // all variables
    Model<2U>    restored_model( string{model.Name()}, subset_variables );
    double       pmin, pmax;
    restored_model.MinMaxOf( "permeability", pmin, pmax );
    _test( approximatelyEqual(pmin,1.0e-13) );
    _test( approximatelyEqual(pmax,1.0e-12) );
    restored_model.MinMaxOf( "node number", pmin, pmax );
    _test( approximatelyEqual(pmin,0) );
    _test( approximatelyEqual(pmax,vset.Vertices()-1U) );
    
    return true;
    
 } // end TestCompleteModel






// TODO: extend to Face and InterFace connectivity
bool MeshManager_Test::Test_BuiltElementConnectivity2D()
 {
    VSet<2U> vset, vset_orig;
    // using VSetMakers to create and compare input data
    create_TrianglePatch_VSet( vset );
    vset_orig = vset;
    Model<2> model( vset, "CSMP-variables.txt" );
    Region<2>& model_domain = model.Region("Model");
    
    // testing the reconstruction of element connectivity from face node pointers (old one gets removed)
    vector<vector<Element<2>*> > nbor_pointers;
    backupNeighborConnectivity( model_domain.CellsBegin(), model_domain.CellsEnd(), nbor_pointers );
    // rebuilding the connectivity
    model.Mesh().BuildConnectivity<Element>( model_domain.CellsBegin(), model_domain.CellsEnd() );
    // getting the connectivity that was recreated
    vector<vector<Element<2>*> > nbor_pointers2;
    backupNeighborConnectivity( model_domain.CellsBegin(), model_domain.CellsEnd(), nbor_pointers2 );
    // comparing the connectivity with the original VSet
    if ( nbor_pointers != nbor_pointers2 ) return false;
    return true;
    
 } // end Test_BuiltElementConnectivity2D




// TODO: extend to Face and InterFace connectivity
bool MeshManager_Test::Test_BuiltElementConnectivity3D()
 {
    VSet<3U> vset, vset_orig;
    create_Pyramid_Hexa_VSet( vset, false );
    vset_orig = vset;
    Model<3> model( vset, "CSMP-variables.txt" );
    Region<3>& model_domain = model.Region("Model");
    
    // testing the reconstruction of element connectivity from face node pointers (old one gets removed)
    vector<vector<Element<3>*> > nbor_pointers;
    backupNeighborConnectivity( model_domain.CellsBegin(), model_domain.CellsEnd(), nbor_pointers );
    // rebuilding the connectivity
    model.Mesh().BuildConnectivity<Element>( model_domain.CellsBegin(), model_domain.CellsEnd() );
    // getting the connectivity that was recreated
    vector<vector<Element<3>*> > nbor_pointers2;
    backupNeighborConnectivity( model_domain.CellsBegin(), model_domain.CellsEnd(), nbor_pointers2 );
    // comparing the connectivity with the original VSet
    if ( nbor_pointers != nbor_pointers2 ) return false;
    return true;
    
 } // end Test_BuiltElementConnectivity3D






// using VSetMakers to create and compare input data
bool MeshManager_Test::Test_MeshTraversal3D()
 {
    // building the test model
    VSet<3U> vset;
    create_Pyramid_Hexa_VSet( vset, false );
//    create_Hexahedra_VSet( vset, false );
//    testCreateTetra_VSet( vset );
    Model<3>   model( vset, "CSMP-variables.txt" );
    Region<3>& model_domain = model.Region("Model");
    
    /* tested: OK
    cout <<"\nNodes with their original indices:\n";
    for ( auto nit=model.Mesh().NodesBegin(); nit!=model.Mesh().NodesEnd(); ++nit )
      cout <<" "<< (*nit).Idx();
    cout << endl;
    cout <<"and in the order in the model subdomain\n";
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      cout <<" "<< (*nit)->Idx() <<":"<< parseBoundary( (*nit)->AtBoundary() );
    cout << endl;
    */
    
    // testing the different mesh traversal algorithms
    // -----------------------------------------------
    
    // are all elements of the contiguous model region discovered
    set<Element<3U>*>  discovered_elements;
    
    floodFill( model_domain.E(0), discovered_elements ); // OK
    _test( discovered_elements.size() == model_domain.Cells() );
    
    set<Element<3>*> elements;
    _test( findContiguousMeshPatch<3>( &(*model_domain.E(0)), elements ) == model_domain.Cells() ); // OK
    _test( elements.size() == model_domain.Cells() );
    
    map<string,vector<Element<3>*> > elmt_patches;
    size_t patches = findContiguousMeshPatches( model.Mesh().ElementsBegin(), model.Mesh().ElementsEnd(),
                                                elmt_patches );
    _test( patches == 1 );
    _test( (*elmt_patches.begin()).second.size() == model_domain.Cells() );
    for ( auto i : elmt_patches ) cout <<" "<< i.first;
    
    // checking that NodesOfSegment() and  CornerNodesPerSegmentForElementOfType() give the same answer
    for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it ) {
         for ( auto i{0}; i < (*it)->FE()->Segments(); ++i ) {
              vector<uint32_t> node_vec;
              (*it)->FE()->NodesOfSegment( i, node_vec );
              sort( node_vec.begin(), node_vec.end() );
              pair<size_t,size_t> node_pair =
                CSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType( (*it)->FE_Type(), i );
              if ( node_pair.first > node_pair.second ) swap( node_pair.first, node_pair.second );
              _test( node_vec[0] == node_pair.first );
              _test( node_vec[1] == node_pair.second );
              if ( node_vec[0] != node_pair.first || node_vec[1] != node_pair.second )
                cerr <<"\n"<< parseElementType( (*it)->FE_Type() ) <<": segm "<< i;
           }
      }
    
    
    // traversal of mesh via node-to-node connectivity
    // -----------------------------------------------
    
    // creating a node connectivity list to check the Node::Neighbor method
    vector<set<size_t>>  node_neighbor_ids;
    nodeNeighbors( model_domain, node_neighbor_ids );
    
//    cout <<"\nTest_MeshTraversal3D: storage requirements for the nodes (bytes):";  // OK - around 250 bytes per node
//    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
//      cout <<"\nNode "<< (*nit)->Idx() <<": "<< sizeOf( (*nit) );
//    cout << endl;

    // visually checking the node connectivity
//    cout <<"\nTest_MeshTraversal3D: node to parent connectivity:";
//    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
//      printNeighbors( (*nit) );
    
    // node-to-node, breadth first mesh traversal
    set<Node<3>*> set_of_interconnected_nodes;
    _test( findInterconnectedNodeCluster( &(*model_domain.N(0)), set_of_interconnected_nodes ) == model_domain.Nodes() );
    _test( set_of_interconnected_nodes.size() == model_domain.Nodes() );
    
    // rebuilding the node to node connectivity verifying that the same results are obtained
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         // recording the current connectivity
         vector<Node<3>*>  node_neighbors;
         node_neighbors.reserve( (*nit)->Neighbors() );
         for ( auto i{0}; i < (*nit)->Neighbors(); ++i )
           node_neighbors.push_back( (*nit)->Neighbor(i) );
         // rebuilding the connectivity
         _test( (*nit)->AssignNodeNeighbors() == node_neighbors.size() );
         // comparing the sorted node pointers with one another
         for ( auto i{0}; i < (*nit)->Neighbors(); ++i )
           _test( (*nit)->Neighbor(i) == node_neighbors[i] );
      }    

    return true;
    
 } // end Test_MeshTraversal3D



/**
        Mesh traversal as implemented by Junchul Kim in 2019
     uses 'prism_test' as test model
*/
bool MeshManager_Test::TestEntityNumberingFunction( Model<3>& model )
{
	// renumbering nodes via Model region
	Region<3U>& model_domain(model.Region("Model"));
	model_domain.UpdateMemberIndexes();
	vector<size_t>  node_numbers_Model;
	node_numbers_Model.reserve(model_domain.Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model.Name() << "': 'Model' numbered nodes:\n";
	for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); ++nit) {
      node_numbers_Model.push_back((*nit)->Idx());
      if ((*nit)->Idx() % 100 == 0) cout << (*nit)->Idx() << "...";
    }
	cout << "\n";

	// renumbering nodes
	const bool in_a_single_sequence(true);
  MeshManager<3>& mesh(model.Mesh());
	mesh.AssignUniqueNumbers(in_a_single_sequence);
	// checking the numbering
	vector<size_t>  nodes_renumbered;
	nodes_renumbered.reserve(model.Mesh().Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model.Name() << "': renumbered nodes:\n";
	{
		set<csmp::Node<3U>*>    discovered_nodes;
		deque<csmp::Node<3U>*>  current_nodes;
		discovered_nodes.insert( &(*mesh.NodesBegin()) );
		current_nodes.push_back( &(*mesh.NodesBegin()) );

		while (!current_nodes.empty()) {
			const csmp::Node<3U>*  n_ptr(*current_nodes.begin());
			for (auto i = 0U; i < n_ptr->Parents(); i++) {
				for (auto j = 0U; j < n_ptr->Parent(i)->Nodes(); j++) {
					if (j != n_ptr->ParentNodeNumber(i)) {
						pair<typename set<csmp::Node<3U>*>::iterator, bool>
							new_node = discovered_nodes.insert(n_ptr->Parent(i)->N(j));
						//cout << n_ptr->Parent(i)->N(j)->Idx() << " ";
						if (new_node.second) current_nodes.push_back(n_ptr->Parent(i)->N(j));
					}
				}
			}
			current_nodes.pop_front();
		}

		for (auto& nit : discovered_nodes )
			nodes_renumbered.push_back(nit->Idx());

		std::sort(nodes_renumbered.begin(), nodes_renumbered.end(), [](auto& lhs, auto& rhs) {return lhs < rhs; });

		for ( auto& nit : nodes_renumbered )
			if (nit % 100 == 0) cout << nit << "...";			
		cout << "\n";

		_test(equal(node_numbers_Model.begin(), node_numbers_Model.end(), nodes_renumbered.begin(), nodes_renumbered.end()));
	}

	return true;
}



/**
     Tests:
      - AddNodeAt()
      - AddElement(),
      - AddFace(),
      - ReplaceElementByFace(),
      - Delete( Element*)
      - connectNeighborsUsingNodeParents()
      - sharedNodes( eptr, eptr )
      - BuildConnectivity( iface)
      
*/
bool MeshManager_Test::TestCellDeletionAndInsertion()
{
  VSet<3U> vset;
  create_Hexahedra_VSet( vset, false );
  Model<3>   model( vset, "MeshManager_Test-variables.txt" );


  MeshManager<3U>& mesh(model.Mesh());
  const size_t     n_original_elmts(mesh.Elements());

  // TODO: prerequisite: test mesh must not be broken (element nbor connectivity issue)
  bool connectivity_is_broken = integrityCheck<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() );
  _test( !connectivity_is_broken );


	// 1. Checking that pointers to elements are not affected by element deletion
  // --------------------------------------------------------------------------
  // make a copy of first element
  Element<3>* e1ptr = &(*mesh.elements_.begin());
  auto elmt1 = (*mesh.elements_.begin());
  // get pointer to element 3
  Element<3>* const e3ptr = &(*next(mesh.elements_.begin(),3));
  const size_t e3idx = e3ptr->Idx();
  // erase element 1 from the plf_colony, repairing connectivity
	mesh.Delete( e1ptr );
  
  _test( e3idx == e3ptr->Idx() );
	_test(mesh.Elements() == n_original_elmts - 1);
  // --------------------------------------------
  
  
  // 2. recreate the missing the element and repair connectivity
  // -----------------------------------------------------------
	LocalVariables				     node_vars = model.Database().LocalVariablesAt(NODE);
	LocalVariables				     elmt_vars = model.Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 intp_vars = model.Database().IntegrationPointVariablesAt(ELEMENT);
  vector<Node<3>*>           nodes( elmt1.NodesBegin(), elmt1.NodesEnd() );

	e1ptr = mesh.AddElement( elmt1.FE_Type(), elmt_vars, intp_vars, nodes, elmt1.Material_ID() );
  //           -------------------------------------------------------------------------------
  mesh.ConnectNeighborsUsingNodeParents( e1ptr );
  _test( e1ptr->ConnectedNeighbors() == 3 );
  
  
  // 3. Create a new node and pyramid element and connected themm with mesh
  // ----------------------------------------------------------------------
  const size_t n_original_nodes(mesh.Nodes());
  // new node is placed at z=-1 (behind, in the center of backface of element 1,
  // getting its barycentre coordinates from it)
  Point<3U> bctr = e1ptr->BaryCenter();
  bctr[2]        = -1.;
	Node<3U>* n_ptr = mesh.AddNodeAt( bctr, node_vars );
  //                     ----------------------------
	_test( mesh.Nodes() == n_original_nodes + 1 );
  n_ptr->Idx( mesh.Nodes() ); // needs a number here

  // create new PYRAMID element on the backside of the mesh
  int32_t material_id(1); // new element's rock_tye
  vector<Node<3U>*>  nodes2 = { nodes[3], nodes[2], nodes[1], nodes[0], n_ptr };
	Element<3U>*	     py_ptr = mesh.AddElement( ISOPARAMETRIC_LINEAR_PYRAMID, elmt_vars, intp_vars, nodes2, material_id );
  //                               -------------------------------------------------------------------------------------
  // debugging diagnosts
  if ( verbose_ ) {
    cout <<"\n"<<"Nodes element 1:  ";
    for ( auto it=e1ptr->NodesBegin(); it!=e1ptr->NodesEnd(); ++it ) cout << (*it)->Idx() <<" ";
    cout <<"\n"<<"Nodes py element: ";
    for ( auto it=py_ptr->NodesBegin(); it!=py_ptr->NodesEnd(); ++it ) cout << (*it)->Idx() <<" ";
    cout << endl;
    size_t nbors_found = mesh.ConnectNeighborsUsingNodeParents( py_ptr );
    _test( nbors_found == 1 );
    _test( py_ptr->Neighbor(4) == e1ptr ); // opposite base-plane
   }
  // backing up the 4 nodes at pyramid base
  vector<Node<3U>*> py_base_nds( nodes2.begin(), next(nodes2.begin(),4) );
  // finding the shared nodes
  pair<vector<Node<3U>*>,bool> snodes = sharedNodes( e1ptr, py_ptr );
  _test( snodes.second == true ); // should be face nodes
  sort( nodes2.begin(), nodes2.end() );
  nodes2.erase( remove( nodes2.begin(), nodes2.end(), n_ptr ), nodes2.end() );
  sort( snodes.first.begin(), snodes.first.end() );
  _test( snodes.first == nodes2 ); // should be face nodes

  // getting rid of pyramid (and its connections so that face 0 of element 1 is on the outside of model again
  _test( mesh.Delete( py_ptr ) == true );
  _test( e1ptr->Neighbor(0) == nullptr );
  
  
  // 4. Creating an element that caps face of element 1 where pyramid was
  // --------------------------------------------------------------------
  //                                                                                           nodes in right order!
	Element<3U>* quad_ptr = mesh.AddElement( ISOPARAMETRIC_LINEAR_QUADRILATERAL, elmt_vars, intp_vars, py_base_nds, material_id );
 // no neighbor connectivity is required here because this is a stand-alone piece of surface mesh
  auto nrml = quad_ptr->UnitNormal();
  _test( nrml[2] < 1 ); // should be outward pointing


  // 4. Creating Face connected to hex at the back model boundary
  // ------------------------------------------------------------
	LocalVariables				     face_vars = model.Database().LocalVariablesAt(FACE);
	IntegrationPointVariables	 fitp_vars = model.Database().IntegrationPointVariablesAt(FACE);
  // testing pair<Element<dim>*>,uint32_t> parentElement<>()
  pair<Element<3U>*,uint32_t> adjacent_hex = parentElement<3U>( py_base_nds.begin(), py_base_nds.end() );
  _test( adjacent_hex.first == e1ptr );
  // creating duplicate faces
  mesh.AddBoundaryFace( adjacent_hex.first, adjacent_hex.second, face_vars, fitp_vars );
  mesh.AddBoundaryFace( adjacent_hex.first, adjacent_hex.second, face_vars, fitp_vars );
  _test( mesh.Faces() == 2 );
  // testing utility that identifies duplicated cells (verbose off)
  size_t duplicate_faces = detectDuplicateCells<3U,Face>( mesh.FacesBegin(), mesh.FacesEnd(), false );
  _test( duplicate_faces == 1 );
  // deleting quad
  mesh.Delete( quad_ptr );
  // after deleting the pyramid and the quad we should be back to original number of elements
  _test( mesh.Elements() == n_original_elmts );
  
	// checking the total number of nodes and elements (after insertion)
	cout << "\nNodes    of the model after insertion: " << mesh.Nodes();
	cout << "\nElements of the model after insertion: " << mesh.Elements();
  // one lost one gained
	_test( mesh.Elements() == n_original_elmts );
 
  // mesh should now be broken, but if it was already initially, this is of no interest
  if ( !connectivity_is_broken ) {
       connectivity_is_broken = integrityCheck<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() );
       _test( connectivity_is_broken );
    }
 
   // 5. Creating a Face between hex1 and a neighbor
  // ------------------------------------------------------------
  uint32_t valid_nbor_face1(UNSPECIFIED);
  for ( uint32_t i{0u}; i<e1ptr->Neighbors(); ++i ) if ( e1ptr->Neighbor(i) ) { valid_nbor_face1=i; break; }
  assert( valid_nbor_face1 < e1ptr->Neighbors() );
  // finding shared face
  pair<size_t,size_t> shared_faces_nbors = findAdjacentFacesFromNeighbors( e1ptr, e1ptr->Neighbor(valid_nbor_face1) );
  pair<size_t,size_t> shared_faces_nodes = findAdjacentFacesFromNodes( e1ptr, e1ptr->Neighbor(valid_nbor_face1) );
  _test( shared_faces_nbors == shared_faces_nodes );
  valid_nbor_face1          = (uint32_t) shared_faces_nodes.first;
  uint32_t valid_nbor_face2 = (uint32_t) shared_faces_nodes.second;
  // finding shared nodes
  auto shared_nodes = sharedNodes( e1ptr, e1ptr->Neighbor(valid_nbor_face1) );
  // putting a quadrilateral element between these elements
	quad_ptr = mesh.AddElement( ISOPARAMETRIC_LINEAR_QUADRILATERAL, elmt_vars, intp_vars, shared_nodes.first, material_id );
  // replacing this element by a face
  auto fptr = mesh.ReplaceElementByFace( quad_ptr, e1ptr, e1ptr->Neighbor(valid_nbor_face1), valid_nbor_face1, valid_nbor_face2, face_vars, fitp_vars );
  _test( fptr != nullptr );
  _test( mesh.Faces() == 3 );

	// 6. Node diagnostics and deletion
  // -----------------------------------------------------------------
  // create new nodes and return pointers to them (they will later be identified as orphan and deleted)
	Node<3U>*	ptr_n2 = mesh.AddNodeAt( Point<3U>(26., 27., 0.0), node_vars, NOT );
	Node<3U>*	ptr_n3 = mesh.AddNodeAt( Point<3U>(29., 30., 0.0), node_vars, NOT );
	Node<3U>*	ptr_n4 = mesh.AddNodeAt( Point<3U>(31., 32., 0.0), node_vars, NOT );
	Node<3U>*	ptr_n5 = mesh.AddNodeAt( Point<3U>(34., 35., 0.0), node_vars, NOT );
  _test( ptr_n2 != nullptr );
  _test( ptr_n3 != nullptr );
  _test( ptr_n4 != nullptr );
  _test( ptr_n5 != nullptr );
  //                   --------------
	_test( mesh.Nodes() == n_original_nodes + 5 );
  vector<Node<3U>*>  orphan_nodes1{ n_ptr, ptr_n2, ptr_n3, ptr_n4, ptr_n5 };
  sort( orphan_nodes1.begin(), orphan_nodes1.end() );
  // using this result to verify OrphanNodeVector()
  vector<const Node<3U>*> orphan_nodes2 = mesh.OrphanNodeVector();
  _test( orphan_nodes2.size() == orphan_nodes1.size() );

	mesh.DeleteNodesAndRepairNodeConnnectivity( orphan_nodes1.begin(), orphan_nodes1.end() );
	_test( mesh.Nodes() == n_original_nodes );
 
	return true;
  
} // end TestCellDeletionAndInsertion





/**
   create Faces at boundary and between prism and hexa elements and then deletes them again.
   Tests:
   
   - AddBoundaryFace()
   - AddFace()
   - ReplaceFaceByInterFace()
   - DeleteInterfacesAndRepairConnnectivity()
   - AssignUniqueNumbers()
   
*/
bool MeshManager_Test::TestFaceDeletionAndInsertion(/* "PyramidHexaPatch" */)
{
	// vsetMakers: create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch"); // only hexa and pyramid type of elements; pyramids=26..32

  MeshManager<3U>& mesh(model.Mesh());
  bool             boundary_face_constructed(false);
  bool             interior_face_constructed(false);
  LocalVariables				     fvars = model.Database().LocalVariablesAt(FACE);
	IntegrationPointVariables	 ivars = model.Database().IntegrationPointVariablesAt(FACE);
  const size_t n_original_nodes = mesh.Nodes();

  // 1. creating faces around the pyramid elements
  // ---------------------------------------------
  Region<3>&       model_domain(model.Region("Model"));
  vector<Face<3>*> ptrs_to_faces_created;
  for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); it++ )
    if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_PYRAMID )
      {
        for ( uint32_t i{0U}; i<(*it)->Faces(); i++ )
          {
             // if the face is at the boundary or has a hexahedral element neighbor, a Face is constructed there
             if ( (*it)->Neighbor(i) == nullptr ) {
                  // FACE CONSTRUCTION
                  ptrs_to_faces_created.push_back( mesh.AddBoundaryFace( (*it), i, fvars, ivars ) );
                  boundary_face_constructed = true; //  -------------------------------------------
               }
             // if the face is within model, we construct a normal face
             else if ( (*it)->Neighbor(i)->FE_Type() == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) {
                  uint32_t opposite_face = UNSPECIFIED;
                  for ( uint32_t j{0U}; j<(*it)->Neighbor(i)->Faces(); ++j ) {
                       if ( (*it)->Neighbor(i)->Neighbor(j) == (*it) ) {
                            opposite_face = j;
                            break;
                         }
                    }
                  // FACE_CONSTRUCTION
                  ptrs_to_faces_created.push_back( mesh.AddFace( (*it), i, (*it)->Neighbor(i), opposite_face, fvars, ivars ) );
                  interior_face_constructed = true; //  ----------------------------------------------------------------------
                  _test( ptrs_to_faces_created.back() != nullptr );
               }
          }
    
    } // end for all pyramids
    
  // interconnect the Face objects with their neighbors
  mesh.BuildConnectivity<Face>( ptrs_to_faces_created.begin(), ptrs_to_faces_created.end() );
    
  cout<<"\n\tcreated "<<  ptrs_to_faces_created.size() <<" faces."<< endl;
  _test( !boundary_face_constructed );
  _test( interior_face_constructed );
    
	_test( mesh.Faces() == ptrs_to_faces_created.size() );
 
  // 2. converting the Face into InterFace objects
  // ---------------------------------------------
  // (since the faces should include the pyramid region all of their nodes must be duplicated)
  // 2.1 testing that all inside elements of of the faces are pyramids
  for ( const auto& face : ptrs_to_faces_created ) {
       _test( face->InnerParent()->FE_Type() == ISOPARAMETRIC_LINEAR_PYRAMID );
       _test( face->OuterParent()->FE_Type() == ISOPARAMETRIC_LINEAR_HEXAHEDRON );
    }
  // 2.2 duplicating the face nodes for the creation on the interface
  LocalVariables				    nvars = model.Database().LocalVariablesAt(NODE);
  vector<vector<Node<3U>*>> outside_nodes( ptrs_to_faces_created.size(), vector<Node<3U>*>{nullptr} );
  unordered_set<Node<3U>*>  unique_outsiders;
  int face_counter{0};
  for ( const auto& face : ptrs_to_faces_created ) {
       outside_nodes[face_counter].resize( face->Nodes(), nullptr );
       for ( uint32_t j{0U}; j<face->Nodes(); ++j ) {
            auto nit = unique_outsiders.insert( face->N(j) );
            // if this is a new node that could be inserted it is duplicated and placed into the outside node container
            if ( nit.second ) { //               node duplication
                outside_nodes[face_counter][j] = mesh.Duplicate( face->N(j), nvars );
                //                                    -------------------------------
                cout <<"\n\t"<<"created node: "<< outside_nodes[face_counter][j]->Idx();
              }
            else // the node that already duplicated and simply has to be retrieved
              outside_nodes[face_counter][j] = (*nit.first);
         }
       face_counter++;
    }

  cout<<"\n\n\tcreated "<<  unique_outsiders.size() <<" nodes."<< endl;

  
  // 2.3 conversion of the last interior face to interface, but without duplication of nodes
  //     (using the inside nodes but in reverse order in order to create the outside nodes for the interface)
  // creating the interface
  LocalVariables				     ifvars = model.Database().LocalVariablesAt(INTER_FACE);
	IntegrationPointVariables	 iivars = model.Database().IntegrationPointVariablesAt(INTER_FACE);
   vector<InterFace<3>*>      ptrs_to_ifaces_created;
 
  // 'ranges' co-iteration of containers not yet: for ( const auto& [face,out_nodes] : zip( ptrs_to_faces_created, outside_nodes ) ) {
  for ( uint32_t i{0U}; i < ptrs_to_faces_created.size(); ++i ) {
        // getting the outside nodes in reverse order
        vector<Node<3U>*> iface_outside_nodes( outside_nodes[i].rbegin(), outside_nodes[i].rend() );
        
       // INTERFACE CONSTRUCTION
       ptrs_to_ifaces_created.push_back( mesh.ReplaceFaceByInterFace( ptrs_to_faces_created[i], ifvars, iivars, iface_outside_nodes ) );
       //                                     ---------------------------------------------------------------------------------------
       _test( ptrs_to_ifaces_created.back() != nullptr );
    }

  // connect interfaces to their neighbors
  mesh.BuildConnectivity<InterFace>( ptrs_to_ifaces_created.begin(), ptrs_to_ifaces_created.end() );
  cout<<"\n\tcreated "<<  ptrs_to_ifaces_created.size() <<" interfaces."<< endl;
    
  _test( mesh.Faces() == 0 );
  _test( mesh.InterFaces() == ptrs_to_faces_created.size() );
	
	// delete the new interface(s) again and repair the mesh connectivity
	mesh.DeleteInterfacesAndRepairConnnectivity( ptrs_to_ifaces_created.begin(), ptrs_to_ifaces_created.end() );
  _test( mesh.Nodes() == n_original_nodes );
 
  // Do all elements still have their nodes?
  mesh.AssignUniqueNumbers();
  for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it )
    for ( auto nit=(*it).NodesBegin(); nit!=(*it).NodesEnd(); ++nit )
      _test( (*nit) != nullptr && (*nit)->Idx() < mesh.Nodes() );
  
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion: model '" << model.Name() << "' (after deletion of faces):\n";
  cout << "\n\tFaces:      " << mesh.Faces() << "\n";
	cout << "\n\tInterFaces: " << mesh.InterFaces() << "\n";

	return true;
  
} // end TestFaceDeletionAndInsertion






/**
    Tests:
    - AddInterFace()
    - findAdjacentFacesFromNodes()
    - findAdjacentFacesFromNeighbors()
    - DeleteInterfacesAndRepairConnnectivity()
    - AddInterveningElement()
*/
bool MeshManager_Test::TestInterFaceDeletionAndInsertion(/* "PyramidHexaPatch" */)
{
	// vsetMakers: create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch");

  MeshManager<3U>&           mesh(model.Mesh());
  bool                       interface_constructed(false);
  LocalVariables				     nvars  = model.Database().LocalVariablesAt(NODE);
  LocalVariables             ifvars = model.Database().LocalVariablesAt(INTER_FACE);
  IntegrationPointVariables  iivars = model.Database().IntegrationPointVariablesAt(INTER_FACE);
  LocalVariables				     evars = model.Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 eivars = model.Database().IntegrationPointVariablesAt(ELEMENT);
  const size_t n_original_ifaces = mesh.InterFaces();

  // puts interfaces between the interior faces of Element # and Element #
  Element<3U>* const eptr( &(*next(mesh.ElementsBegin(),4)) );
  InterFace<3U>*     ifptr(nullptr);
  vector<InterFace<3U>*> iface_ptrs;
  iface_ptrs.reserve( eptr->Faces() );
  for ( uint32_t i{0}; i<eptr->Faces(); i++ )
    {
       if ( eptr->Neighbor(i) != nullptr && !interface_constructed ) {
            // getting the nodes for the inside of the future interface
            auto fnids = eptr->FE()->NodesOfFace( i );
            vector<Node<3U>*> inside_nodes;
            inside_nodes.reserve( fnids.size() );
            for ( uint32_t j=0U; j<fnids.size(); ++j )
              inside_nodes.push_back( eptr->N( fnids[j] ) );
            // duplicating these nodes to get nodes for the inside element and the other side
            vector<Node<3U>*> middle_nodes, outside_nodes;
            middle_nodes.reserve( fnids.size() );
            outside_nodes.reserve( fnids.size() );
            for ( auto& j : fnids ) {
                 Node<3U>* mnptr = mesh.Duplicate( eptr->N(j), nvars );
                 middle_nodes.push_back( mnptr );
                 Node<3U>* onptr = mesh.Duplicate( eptr->N(j), nvars );
                 outside_nodes.push_back( onptr );
              }
            // find matching faces via the shared nodes
            pair<uint32_t,uint32_t> face_ids1 = findAdjacentFacesFromNodes( eptr, eptr->Neighbor(i) );
            _test( i == face_ids1.first );
            // find matching faces via neighbor element pointers (faster)
            pair<uint32_t,uint32_t> face_ids2 = findAdjacentFacesFromNeighbors( eptr, eptr->Neighbor(i) );
            _test( face_ids1.first  == face_ids2.first );
            _test( face_ids1.second == face_ids2.second );
            //                        inner  outer
            ifptr = mesh.AddInterFace( eptr, face_ids1.first, eptr->Neighbor(i), face_ids1.second, ifvars, iivars );
            //           ------------------------------------------------------------------------------------------
            interface_constructed = true;
            // create an intervening element
            const int32_t material_id(5);
            Element<3U>*	ieptr = mesh.AddInterveningElement( ifptr, evars, eivars, middle_nodes, material_id );
            //                         -------------------------------------------------------------------------
            // is the InterFace connected to middle element
            _test( ifptr->InterveningElement() == ieptr );
            iface_ptrs.push_back( ifptr );
         }
    }
  if ( interface_constructed ) {
       cout<<"\n\tcreated interface: ";
       cout << ifptr->Idx() << endl;
       ifptr->Out();
    }
    
	_test( iface_ptrs.size() == 1 );
 
	// delete the new interface(s) again
	mesh.DeleteInterfacesAndRepairConnnectivity( iface_ptrs.begin(), iface_ptrs.end() );
	cout << "\nMeshManager_Test::TestInterFaceDeletionAndInsertion: model '" << model.Name() << "' (after deletion of interfaces):\n";
	cout << "\nFaces: " << mesh.InterFaces() << "\n";

	_test( mesh.InterFaces() == n_original_ifaces );

	return true;
}



  
/*
    Erasure of all elements, faces and interfaces from 2d model.
*/
bool MeshManager_Test::TestEraseAllPrimitives()
{
   // TODO: needs 2D VSet with SplitBoundary (throughgoing) and fracture type that ends inside of model
   /*
   MeshManager<2U>& mesh(model2d_->Mesh());
   const size_t     n_original_elmts(mesh.Elements());
   const size_t     n_original_faces(mesh.Faces());
   const size_t     n_orig_interfaces(mesh.InterFaces());
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
   _test( mesh.Erase( mesh.NodesBegin(),      mesh.NodesEnd() )      == n_original_elmts );
   _test( mesh.Erase( mesh.ElementsBegin(),   mesh.ElementsEnd() )   == n_original_elmts );
   _test( mesh.Erase( mesh.FacesBegin(),      mesh.FacesEnd() )      == n_original_faces );
   _test( mesh.Erase( mesh.InterFacesBegin(), mesh.InterFacesEnd() ) == n_orig_interfaces );
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
	 _test( mesh.Elements() + mesh.Faces() + mesh.InterFaces() + mesh.Nodes() == 0 );
  */
	return true;
}



bool MeshManager_Test::TestNeigbourVersusFaceConsistency()
 {
    // ErrorHandler& csmp_error{ ErrorHandler::Instance() };
    const string model_name("box1x1x1_hexa_struct");
    if ( verbose_ ) {
         cout <<"\n"<<"MeshManager_Test::TestNeigbourVersusFaceConsistency: "<<this->getName()<<endl<<endl;
         cout <<"Building Model: '"<< model_name <<"'"<< endl;
      }
    const string variablesFile("MeshManager_Test-variables.txt");
    constexpr uint32_t dim{3u};
    //                   fileset             regions-file
    ANSYS_Model3D model( model_name.c_str(), model_name.c_str(), variablesFile.c_str(), true );
    Region<dim>&  model_domain(model.Region("Model"));
    model_domain.UpdateMemberIndexes();
 //   VTK_Interface<3U>  vtk_output;
 //   VTU_Interface<3U>  vtu_output( model );

    size_t errors{0ul};
    
    // checking that the corner elements have the right neighbors on each model side
    // and that the node numbers match
    size_t n_corner_elmts_found{0ul};
    for ( const auto& nit : model_domain.NodeVector() )
      if ( isCorner(nit->AtBoundary()) )
        {
           // verify that this is a corner node
           _test( nit->Parents() == 1u );

           // getting the parent element of this node
           const Element<3u>* const eptr = nit->Parent(0u);
           if ( verbose_ ) {
                cout <<"\n"<<"Discovered corner element:";
                eptr->Out();
             }
           // corner hex should have 3 outside faces (6-3=3)
           _test( eptr->ConnectedNeighbors() == 3u );
           
           // for the element faces at model boundary get the face nodes and verify that they match the expected nodes of the face
           for ( uint32_t face_id{0u}; face_id<eptr->Neighbors(); ++face_id )
             // for each element face at the boundary
             if ( eptr->Neighbor(face_id) == nullptr ) {
                  vector<uint32_t> fnids = eptr->FE()->CornerNodesOfFace( face_id );
                  // All face nodes must:
                  //  - reside at model boundary
                  for ( const auto& local_node : fnids ) _test( eptr->N(local_node)->AtBoundary() != NOT );
                  //  - match the local node numbers returned for them from parent element
                  for ( const auto& local_node : fnids ) {
                       // for the face node, find parent element that is equal to corner element
                       Element<dim>* eptr_p{nullptr};
                       uint32_t      eparent{0u};
                       while ( eparent<eptr->N(local_node)->Parents() ) {
                            eptr_p = eptr->N(local_node)->Parent(eparent);
                            // once the parent element of interest has been found
                            if ( eptr_p == eptr ) break;
                            ++eparent;
                         }
                       assert( eptr_p != nullptr );
                       auto local_node_number_in_parent_elmt = eptr->N(local_node)->ParentNodeNumber(eparent);
                       _test( local_node_number_in_parent_elmt == local_node );
                       if ( local_node_number_in_parent_elmt != local_node ) errors++;
                    }
               }
            n_corner_elmts_found++;
        }
    _test( n_corner_elmts_found == 8u );
      
     return ( errors == 0U );
 }
 
 
/**
       Tests that a local connectivity update does not break the connectivity
       at the margin of the updated region.
*/
bool MeshManager_Test::Test_BuildConnectivity()
  {
     ModelTopology topology;
     VSet<3U>      vset;
     create_FracBox( topology, vset );
     const bool    create_boundaries_from_regions{ true };
     Model<3U>     model( topology, vset, "MeshManager_Test-variables.txt", create_boundaries_from_regions );
     int           errors(0ul);
     
     // 1. testing that a local connectivity update does not break neighbor connections else
     vector<Element<3U>*>  test_volumes;
     test_volumes.reserve( model.Mesh().Elements() );
     const Region<3U>&     model_domain = model.Region("Model");
     for ( auto& eit : model_domain.CellVector() )
       if ( eit->IsVolume() && eit->ConnectedNeighbors() == eit->Neighbors() )
         test_volumes.push_back( eit );
         
     // connectivity update
     model.Mesh().BuildVolumeConnectivity<Element>( next(test_volumes.begin(),30), next(test_volumes.begin(),60) );
         
     // testing
     for ( const auto& eit : test_volumes )
       if ( eit->ConnectedNeighbors() != eit->Neighbors() ) errors++;
     
     return !( errors > 0 );
  }
  
 
 
 
/**
    For a the box-shaped "FracBox" model with lower dimensional elements, the neighbor and parent  connectivity is removed and then rebuilt.
    The results are tested with integrityCheck.
    
    Also tests BuildVolumeConnectivty, BuidSurfaceConnectivity, BuildLineConnectivity
*/
bool MeshManager_Test::Test_UpdateConnectivity()
 {
     ModelTopology topology;
     VSet<3U>      vset;
     create_FracBox( topology, vset );
     const bool    create_boundaries_from_regions{ true };
     Model<3U>     model( topology, vset, "MeshManager_Test-variables.txt", create_boundaries_from_regions );
     int           errors(0ul);
     
     // 1. Removing existing connectivity
     // ---------------------------------
     // element
     Region<3U>&  model_domain = model.Region("Model");
     for ( auto& eit : model_domain.CellVector() )
       for ( uint32_t i{0u}; i<eit->Neighbors(); ++i )
         if ( eit->Neighbor(i) != nullptr )
           eit->Neighbor(i)->UnassignNeighbor(i);
     // node parents and node neighbors
     for ( auto& nit : model_domain.NodeVector() ) {
          nit->EraseParents();
          nit->EraseNeighbors();
          // manifolds
          if ( nit->IsManifold() )
            nit->Manifold()->Remove( nit );
       }
      // faces
      if ( model.Mesh().Faces() > 0 ) {
            for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit ) {
                for ( uint32_t i{0u}; i<fit->Neighbors(); ++i )
                  if ( fit->Neighbor(i) != nullptr )
                    fit->Neighbor(i)->UnassignNeighbor(i);
              }
         }
     
     // 2. Recreating connectivity
     // --------------------------
     // connectivity update
    auto t0 = chrono::high_resolution_clock::now();
    // _________________________________
    model.Mesh().UpdateConnectivity();
    auto t1 = chrono::high_resolution_clock::now();
    cout<<"\n"<<"MeshManager_Test::Test_UpdateConnectivity: time to reconnect the mesh: ";
    cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
    cout <<" milliseconds."<< endl;
         
     // testing
    const bool elements_ok = integrityCheck<3U,Element>( model.Mesh().ElementsBegin(), model.Mesh().ElementsEnd() );
    _test( elements_ok );
    if ( !elements_ok ) errors++;
    const bool faces_ok = integrityCheck<3U,Face>( model.Mesh().FacesBegin(), model.Mesh().FacesEnd() );
    _test( faces_ok );
    if ( !faces_ok ) errors++;
     
    return !( errors > 0 );
     
 } // end Test_UpdateConnectivity
 
 
 
 
 

} // end csmp
