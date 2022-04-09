//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "MeshManager_Test.h"
#include "MeshManagementUtilities.h"
#include "vsetMakers.h"
#include "VTU_Interface.h"
#include "Region.h"
#include "Element.h"
#include "compareFloats.h"

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
         _test( findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), elmt_map3 ) > 1 );
         cout <<"\nMeshManager_Test::CheckConnectivityOfModel3D: the model is discontiguous and consists of the mesh patches:";
         for ( auto& it : elmt_map3 ) {
              cout <<"\n\t\t'"<< it.first <<"': "<< it.second.size() <<" elements.";
           }
      }
    cout << "\n\nExamining the connectivity of  Faces: " << mesh.Faces() << "\n";
    if ( mesh.Faces() > 0 ) {
        map<string,vector<Face<3U>*> >  face_map3;
        _test( findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) >= 1 );
        cout << "\nInterconnected faces: " << (*face_map3.begin()).second.size() << "\n";
      }
    cout << "\nExamining the connectivity of  Interfaces: " << mesh.InterFaces() << "\n";
    if ( mesh.InterFaces() > 0 ) {
        map<string,vector<InterFace<3U>*> >  iface_map3;
        _test( findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) >= 1 );
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
   
  _test(Test_BuiltElementConnectivity2D()); // OK
  _test(Test_BuiltElementConnectivity3D()); // OK
  _test(Test_parentElementsSharedByFace()); // OK

  _test( Test_MeshTraversal3D(/* Pyramid_Hexa_VSet */) );
  
	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::TestCompleteModel2D";
	cout << "\n(elements,faces,interfaces,regions,boundaries, split boundaries)";
	cout << "\n----------------------------------------------------";
  _test( TestCompleteModel2D() );

  // building more complex 'FracBox' model with Boundaries and lower-dimensional elements for further testing
  VSet<3U>      vset;
  ModelTopology topology;
  string        var_file("CSMP-variables.txt");
  const bool    regions_to_boundaries{true};
  test_Create_FracBox( topology, vset );
  
  Model<3>      model( topology, vset, var_file.c_str(), regions_to_boundaries );
  
	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::CheckConnectivityOfModel3D";
	cout << "\n----------------------------------------------------";
  _test( CheckConnectivityOfModel3D(model) );

	cout << "\n----------------------------------------------------";
	cout << "\nMeshManager_Test::TestElementDeletionAndInsertion";
	cout << "\n----------------------------------------------------";
	_test(TestElementDeletionAndInsertion());

	cout << "\n-------------------------------------------------";
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion";
	cout << "\n-------------------------------------------------";
	_test(TestFaceDeletionAndInsertion());

	cout << "\n------------------------------------------------------";
	cout << "\nMeshManager_Test::TestInterFaceDeletionAndInsertion";
	cout << "\n------------------------------------------------------";
	_test(TestInterFaceDeletionAndInsertion());

	cout << "\n------------------------------------------------";
	cout << "\nMeshManager_Test::TestEntityNumberingFunction";
	cout << "\n------------------------------------------------";
	_test(TestEntityNumberingFunction( model ) ); // originally using 'prism_test'

	_test(TestEraseAllPrimitives());
 
  cout << endl;

} // end run





/**
Checks that numbers of elements etc. in mesh manager do indeed reflect those of input model
*/
void MeshManager_Test::TestBasics()
{
	// vsetMakers: test_Create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  test_Create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch");
  // bounding box
  Point<3> xyz_min, xyz_max;
  model.MinMaxCoordinates( xyz_min, xyz_max );

  const MeshManager<3U>& mesh = model.Mesh();

  // returns true if the mesh consists of multiple element types
  _test(mesh.HybridElementMesh() == true);

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
  
  // test model model volume (based on bounding box)
  const Region<3>&  model_domain(model.Region("Model"));
  const double volume = (xyz_max[0]-xyz_min[0]) * (xyz_max[1]-xyz_min[1]) *  (xyz_max[2]-xyz_min[2]);
  _test( approximatelyEqual(model_domain.Volume(),volume) );
  
  // does this also work correctly in 2D?
  VSet<2> vset2D;
  test_Create_MeshPatchWithLineElements_VSet( vset2D );
  Model<2> model2D( vset2D, "CSMP-variables.txt" ); // needs "element number" and "node number"
  model2D.Name("MeshPatchWithLineElements");
  const Region<2>&  model_domain2D(model2D.Region("Model"));
  const double area{ 6. * 5. };
  _test( approximatelyEqual(model_domain2D.Volume(),area) );
  
  // do all nodes have the expected neighbors?
  vector<set<size_t>>  node_neighbors;
  nodeNeighbors( model_domain, node_neighbors );
  // testing using the connectivity from the VSet

} // end TestBasics



/**
       2D test case for rectangular model of 2 domains separated by SplitBoundary.
       Extra partitally penetrating split boundary crosses dividing SplitBoundary.
       
       @author SKM
       @date 19/3/22
*/
bool MeshManager_Test::TestCompleteModel2D()
 {
    VSet<2U>       vset;
    ModelTopology  topo = test_Create_BoundarySplitBoundaryPatch( vset );
    const bool     treat_all_domains_as_regions{false};
    Model<2>       model( topo, vset, "CSMP-1phase-variables.txt", treat_all_domains_as_regions );
    
    // testing that the model has the right area (also checks element orientations)
    const Region<2>&  model_domain(model.Region("Model"));
    const double model_area{ 4.5 * 7. };
    _test( approximatelyEqual(model_domain.Volume(),model_area) );
    
    // can such a model be output to VTU?
    VTU_Interface<2>  vtu_out( model );
    vtu_out.OutputDataToVTU( "SPLIT22_BASIC", "permeability", "Model", 0 );
    
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



// tests method with the same name
bool MeshManager_Test::Test_parentElementsSharedByFace()
 {
    VSet<3U> vset;
    // testing with element 13 with face 4 on the LEFT outside
    test_Create_Prism_Hexa_VSet( vset, false );
    Model<3U>    model( vset, "CSMP-variables.txt" );
    const size_t ELMT{13}; // 13 in VSet
    Element<3>*  eptr = &(*next(model.Mesh().ElementsBegin(),ELMT));
    //eptr->Out();
    assert( eptr->Neighbor(4) == nullptr );
    
    // getting an inner face in the 3D model that is not on the boundary
    vector<Node<3>*> face_nodes;
    Element<3>*      inner_eptr(nullptr), *outer_eptr(nullptr);
    for ( auto i{0}; i<eptr->Neighbors(); ++i )
     if ( eptr->Neighbor(i) != nullptr ) {
          vector<uint32_t> fnids;
          eptr->FE()->NodesOfFace(i,fnids);
          face_nodes.reserve( fnids.size() );
          for (size_t j{0}; j<fnids.size(); ++j )
            face_nodes.push_back( eptr->N( fnids[j] ) );
          inner_eptr = eptr;
          outer_eptr = eptr->Neighbor(i);
          break;
       }
    
    // calling the function that is being tested
    pair<Element<3>*,Element<3>*> parents = parentElementsSharedByFace<3>( face_nodes.begin(), face_nodes.end() );
    
    // test that the correct neighbor elements were found (inner one should be first
    _test( parents.first  != nullptr );
    _test( parents.second != nullptr );
    _test( parents.first  == inner_eptr );
    _test( parents.second == outer_eptr );
    
    // now testing for face 4 that is on the left outside
    vector<uint32_t> fnids;
    eptr->FE()->NodesOfFace(4,fnids);
    face_nodes.resize( fnids.size() );
    for (size_t j{0}; j<fnids.size(); ++j )
      face_nodes[j] = eptr->N( fnids[j] );
    inner_eptr = eptr;
    outer_eptr = nullptr;

    // calling the function that is being tested
    parents = parentElementsSharedByFace<3>( face_nodes.begin(), face_nodes.end() );
    
    // test that the correct neighbor elements were found (inner one should be first
    _test( parents.first  != nullptr );
    _test( parents.second == nullptr );
    _test( parents.first  == eptr );

    return true;
    
 } // end Test_parentElementsSharedByFace




// using VSetMakers to create and compare input data
bool MeshManager_Test::Test_BuiltElementConnectivity2D()
 {
    // 2D functionality
    VSet<2U> vset, vset_orig;
    test_Create_TrianglePatch_VSet( vset );
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



// using VSetMakers to create and compare input data
bool MeshManager_Test::Test_BuiltElementConnectivity3D()
 {
    // 2D functionality
    VSet<3U> vset, vset_orig;
    test_Create_Pyramid_Hexa_VSet( vset, false );
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
    test_Create_Pyramid_Hexa_VSet( vset, false );
//    test_Create_Hexahedra_VSet( vset, false );
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
    size_t patches = findStandAloneMeshPatches( model.Mesh().ElementsBegin(), model.Mesh().ElementsEnd(),
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
         _test( (*nit)->ReassignNeighbors() == node_neighbors.size() );
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
     Starts with a hexahedral mini model
*/
bool MeshManager_Test::TestElementDeletionAndInsertion()
{
  VSet<3U> vset;
  test_Create_Hexahedra_VSet( vset, false );
  Model<3>   model( vset, "CSMP-variables.txt" );
  //Region<3>& model_domain = model.Region("Model");

  MeshManager<3U>& mesh(model.Mesh());
  const size_t     n_original_elmts(mesh.Elements());
  
	// 0. the first element is copy constructed and stored, and then deleted
	csmp::Element<3U> first_element(*mesh.ElementsBegin());
	// 1. delete elements 1 from the model
	mesh.Erase( mesh.ElementsBegin() );
  
	_test(mesh.Elements() == n_original_elmts - 1);

  const size_t n_original_nodes(mesh.Nodes());
	// 2. create a new element with its nodes
	IsoparametricLinearPyramid fe;
	LocalVariables				     node_vars = model.Database().LocalVariablesAt(NODE);
	LocalVariables				     elmt_vars = model.Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 intp_vars = model.Database().IntegrationPointVariablesAt(ELEMENT);
  // copy the first node
  Node<3U>      n1(*mesh.NodesBegin());
  const size_t  nearby_node(4);
	Node<3U>*		  ptr_n1 = mesh.AddNodeAtUniqueLocation( n1.Coordinate(), nearby_node, node_vars );
  //                     ----------------------------
  // method must return pointer to node 1 pointer
  _test( ptr_n1 == &(*mesh.NodesBegin()) );

  // create new nodes and return pointers to them
	Node<3U>*		ptr_n2 = mesh.AddNodeAt( Point<3U>(26., 27., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n3 = mesh.AddNodeAt( Point<3U>(29., 30., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n4 = mesh.AddNodeAt( Point<3U>(31., 32., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n5 = mesh.AddNodeAt( Point<3U>(34., 35., 0.0), node_vars, NOT );
  //                   --------------

	_test( mesh.Nodes() == n_original_nodes + 5 );

  // create new PYRAMID element
  int32_t material_id(1); // new element's rock_tye
  vector<Node<3U>*>  nodes = {ptr_n1,ptr_n2,ptr_n3,ptr_n4,ptr_n5};
	Element<3U>*	     ptr_e1 = mesh.AddElement( ISOPARAMETRIC_LINEAR_PYRAMID, elmt_vars, intp_vars, nodes, material_id );

	// assign new element as a parent to its nodes (TODO: should be done when nodes are connected
  ptr_n1->ResizeParentStorage( n1.Parents()+1 );
	ptr_n1->Assign( n1.Parents() - 1, ptr_e1);
  
  ptr_n2->ResizeParentStorage( n1.Parents()+1 );
	ptr_n2->Assign( n1.Parents() - 1, ptr_e1);
  
  ptr_n3->ResizeParentStorage( n1.Parents()+1 );
	ptr_n3->Assign( n1.Parents() - 1, ptr_e1);
  
  ptr_n4->ResizeParentStorage( n1.Parents()+1 );
	ptr_n4->Assign( n1.Parents() - 1, ptr_e1);
  
  ptr_n4->ResizeParentStorage( n1.Parents()+1 );
	ptr_n5->Assign( n1.Parents() - 1, ptr_e1);

	// checking the total number of nodes and elements (after insertion)
	std::cout << "\nNodes    of the model after insertion: " << mesh.Nodes();
	std::cout << "\nElements of the model after insertion: " << mesh.Elements();
  // one lost one gained
	_test( mesh.Elements() == n_original_elmts );

	// 2. deleting these nodes again
	mesh.Delete( nodes.begin(), nodes.end() );
	_test( mesh.Nodes() == n_original_nodes );

// TODO: test insertion / deletion / connection of Face and InterFace objects
	return true;
}




// create Faces at boundary and between elements and then delete them again.
bool MeshManager_Test::TestFaceDeletionAndInsertion(/* "PyramidHexaPatch" */)
{
	// vsetMakers: test_Create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  test_Create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch");

  MeshManager<3U>& mesh(model.Mesh());
  bool             boundary_face_constructed(false);
  bool             interior_face_constructed(false);
  LocalVariables				     fvars = model.Database().LocalVariablesAt(FACE);
	IntegrationPointVariables	 ivars = model.Database().IntegrationPointVariablesAt(FACE);
  const size_t n_original_faces = mesh.Faces();

  // element 4
  Element<3U>* const eptr( &(*next(mesh.ElementsBegin(),4)) );
  Face<3U>*    fptr1(nullptr), *fptr2(nullptr);
  for ( auto i{0U}; i<eptr->Faces(); i++ )
    {
       // if the face is at the boundary, we construct a boundary face
       if ( eptr->Neighbor(i) == nullptr && !boundary_face_constructed ) {
            fptr1 = mesh.AddBoundaryFace( eptr, i, fvars, ivars );
            boundary_face_constructed = true;
         }
       // if the face is within model, we construct a normal face
       if ( eptr->Neighbor(i) != nullptr && !interior_face_constructed ) {
            uint32_t opposite_face = UNSPECIFIED;
            for ( auto j{0}; j<eptr->Neighbor(i)->Faces(); ++j ) {
                 if ( eptr->Neighbor(i)->Neighbor(j) == eptr ) opposite_face = j;
                 break;
              }
            fptr2 = mesh.AddFace( eptr, i, eptr->Neighbor(i), opposite_face, fvars, ivars );
            interior_face_constructed = true;
         }
    }
  cout<<"\n\tcreated faces: ";
  if ( boundary_face_constructed ) fptr1->Out();
  if ( interior_face_constructed ) fptr2->Out();
    
	_test( mesh.Faces() == n_original_faces + boundary_face_constructed + interior_face_constructed );
 
  // conversion of interior face to interface, but without duplication of nodes
  vector<Node<3U>*>         outside_nodes( fptr2->Nodes(), nullptr );
  // using the same nodes but in reverse order
  for ( auto i{0U}; i<fptr2->Nodes(); i++ ) outside_nodes[ fptr2->Nodes()-i-1U ] = fptr2->N(i);
  LocalVariables				     ifvars = model.Database().LocalVariablesAt(INTER_FACE);
	IntegrationPointVariables	 iivars = model.Database().IntegrationPointVariablesAt(INTER_FACE);
  if ( interior_face_constructed ) {
       InterFace<3U>* ifptr = mesh.ReplaceFaceByInterFace( fptr2, ifvars, iivars, outside_nodes );
       ifptr->Out();
    }
	
	// delete the new face(s) again
  vector<Face<3U>*> face_ptrs{ fptr1, fptr2 };
	mesh.Delete( face_ptrs.begin(), face_ptrs.end() );
  _test( fptr1 == nullptr );
  _test( fptr2 == nullptr );
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion: model '" << model.Name() << "' (after deletion of faces):\n";
	cout << "\nFaces: " << mesh.Faces() << "\n";

	_test( mesh.Faces() == n_original_faces );

	return true;
}







bool MeshManager_Test::TestInterFaceDeletionAndInsertion(/* "PyramidHexaPatch" */)
{
	// vsetMakers: test_Create_Pyramid_Hexa_VSet
  VSet<3>    vset;
  const bool bSkewed{false};
  test_Create_Pyramid_Hexa_VSet( vset, bSkewed );
  Model<3> model( vset );
  model.Name("PyramidHexaPatch");

  MeshManager<3U>&           mesh(model.Mesh());
  bool                       interface_constructed(false);
  LocalVariables				     ifvars = model.Database().LocalVariablesAt(INTER_FACE);
	IntegrationPointVariables	 iivars = model.Database().IntegrationPointVariablesAt(INTER_FACE);
  const size_t n_original_ifaces = mesh.InterFaces();

  // puts interfaces between the interior faces of Element # and Element #
  Element<3U>* const eptr( &(*next(mesh.ElementsBegin(),4)) );
  InterFace<3U>*     ifptr(nullptr);
  vector<InterFace<3U>*> iface_ptrs;
  iface_ptrs.reserve( eptr->Faces() );
  for ( auto i{0}; i<eptr->Faces(); i++ )
    {
       if ( eptr->Neighbor(i) != nullptr && !interface_constructed ) {
            vector<uint32_t> fnids;
            // getting the nodes for the inside of the future interface
            eptr->FE()->NodesOfFace( i, fnids );
            vector<Node<3U>*> inside_nodes;
            inside_nodes.reserve( fnids.size() );
            for ( auto j=0U; j<fnids.size(); ++j )
              inside_nodes.push_back( eptr->N( fnids[j] ) );
            // duplicating these nodes to get nodes for the inside element and the other side
            vector<Node<3U>*> middle_nodes, outside_nodes;
            middle_nodes.reserve( fnids.size() );
            outside_nodes.reserve( fnids.size() );
            for ( auto j=0U; j<fnids.size(); ++j ) {
                 Node<3U>* mnptr = mesh.Duplicate( eptr->N( fnids[j] ), MIDDLE );
                 middle_nodes.push_back( mnptr );
                 Node<3U>* onptr = mesh.Duplicate( eptr->N( fnids[j] ), OUTSIDE );
                 outside_nodes.push_back( onptr );
              }
            // find matching faces via the shared nodes
            pair<uint32_t,uint32_t> face_ids1 = findAdjacentElementFaces( eptr, eptr->Neighbor(i) );
            _test( i == face_ids1.first );
            // find matching faces via neighbor element pointers (faster)
            pair<uint32_t,uint32_t> face_ids2 = findAdjacentFacesFromNeighbors( eptr, eptr->Neighbor(i) );
            _test( face_ids1.first == face_ids2.first );
            _test( face_ids1.first == face_ids2.second );
            //                        inner  outer
            ifptr = mesh.AddInterFace( eptr, face_ids1.first, eptr->Neighbor(i), face_ids1.second, ifvars, iivars );
            interface_constructed = true;
            // create an intervening element
            const int32_t material_id(5);
            Element<3U>*	ieptr = mesh.AddElement( ISOPARAMETRIC_LINEAR_TRIANGLE, ifvars, iivars, middle_nodes, material_id );
            
            // connecting the InterFace to the middle element
            ifptr->Assign( ieptr );
            iface_ptrs.push_back( ifptr );
         }
    }
  cout<<"\n\tcreated faces: ";
  if ( interface_constructed ) ifptr->Out();
    
	_test( mesh.Faces() == n_original_ifaces + interface_constructed );
 
	// delete the new interface(s) again
	mesh.Delete( iface_ptrs.begin(), iface_ptrs.end() );
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



/**
      finding the neighbors nodes of each node.
      
      This method is equivalent to creating a sparsity pattern for an accumulation.
      
      @test OK SKM 8/12/21
*/
template<uint32_t dim>
void nodeNeighbors( const Region<dim>& subdomain, vector<set<size_t>>& node_neighbors )
 {
    if ( !node_neighbors.empty() ) node_neighbors.clear();
    node_neighbors.resize( subdomain.Nodes() );
    
    // 1. get continuous indices to access the sets contained in the node_neighbor vectors
    subdomain.RenumberNodes();
 
    vector<uint32_t> segm_nodes;

    const auto elmtsEnd{ subdomain.CellsEnd() };
    for ( auto it=subdomain.CellsBegin(); it!=elmtsEnd; ++it ) {
         const auto n_segments{ (*it)->Segments() };
         for ( auto segm_id{0}; segm_id < n_segments; ++segm_id ) {
              (*it)->FE()->NodesOfSegment( segm_id, segm_nodes );
              // replacing local with global node ids
              for ( auto& sit : segm_nodes ) sit = static_cast<uint32_t>((*it)->N(sit)->Idx());
              // storing the node-to-node connections avoiding duplicates
              assert( segm_nodes.size() == 2 ); // only linear segments are considered by this function
              size_t segm_node1{ subdomain.N(*segm_nodes.begin())->Idx() };
              size_t segm_node2{ subdomain.N(*segm_nodes.rbegin())->Idx() };
              node_neighbors[ segm_node1 ].insert( segm_node2 );
              node_neighbors[ segm_node2 ].insert( segm_node1 );
           }
      }
      
    // printing the node-neighbor vector for testing
    /*
    cout <<"\n\nnodeNeighbors: connectivity created for "<< node_neighbors.size() <<" nodes:";
    size_t node{0};
    for ( auto nit : node_neighbors ) {
         cout <<"\n\t" << node <<": ";
         for (  auto i : nit )
           cout << subdomain.N(i)->Idx() <<" ";
         cout <<" ("<< parseBoundary( subdomain.N(node)->AtBoundary() ) <<")";
         node++;
      }
    */
      
 } // end nodeNeighbors

template void nodeNeighbors( const Region<3>&, vector<set<size_t>>& );
template void nodeNeighbors( const Region<2>&, vector<set<size_t>>& );


} // end csmp
