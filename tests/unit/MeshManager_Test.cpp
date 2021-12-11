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
#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "Region.h"
#include "Element.h"

#include "IsoparametricLinearPyramid.h"
#include "VTK_Interface.h"
#include "ModelTopology.h"
#include "NodeManifold.h"


using namespace std;

namespace csmp {

/**
custom models
uses "CSMP-1phase-variables.txt"
*/
MeshManager_Test::MeshManager_Test()
 : model2d_(nullptr), model3d_(nullptr)
{
}


void MeshManager_Test::Create_ANSYS2D_Model( bool reconstruct_from_file )
 {
    // ansys 2d model - contiguous
    cout << "\n------------------------------------------";
    cout << "\nMeshManager_Test: ANSYS model 'box2d_fault'";
    cout << "\n------------------------------------------";
    model2d_name_ = "box2d_fault";
    string varFileName = "CSMP-variables.txt";
    model2d_ = new ANSYS_Model2D(model2d_name_.c_str(), varFileName.c_str());
    MeshManager<2>& mesh(model2d_->Mesh());
    cout << "\nNodes: " << mesh.Nodes() << "\n";
    set<Element<2>*> elements;
    cout << "\nInterconnected elements: " << findContiguousMeshPatch<2,Element>( &(*mesh.ElementsBegin()), elements ) << "\n";
    cout << "\nElements: " << mesh.Elements() << "\n";
    std::map<std::string,std::vector<Element<2U>*> > patch_map;
    cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map ) << "\n";
    cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
    std::map<std::string,std::vector<Face<2U>*> >  face_map;
    cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map ) << "\n";
    cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
    std::map<std::string,std::vector<InterFace<2U>*> >  iface_map;
    cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map ) << "\n";
    
    if (reconstruct_from_file) {
        model2d_->OutputToBinaryFile(model2d_name_.c_str());
        delete model2d_;
        model2d_ = new Model<2U>(model2d_name_);
        MeshManager<2>& mesh(model2d_->Mesh());
        cout << "\nNodes: " << mesh.Nodes() << "\n";
        cout << "\nNode Groups: " << findContiguousMeshPatch<2,Element>( &(*mesh.ElementsBegin()), elements ) << "\n";
        cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
        cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map ) << "\n";
        cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
        cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), model2d_->Mesh().FacesEnd(), face_map ) << "\n";
        cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
        cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), model2d_->Mesh().InterFacesEnd(), iface_map ) << "\n";
      }
      
 } // end Create_ANSYS2D_Model
 
 
 
void MeshManager_Test::Create_ANSYS3D_Model( bool contiguous, bool reconstruct_from_file )
 {
    // ansys 3d model - discontiguous
   if ( !contiguous ) {
        cout << "\n-------------------------------------------------------";
        cout << "\nMeshManager_Test: ANSYS model 'ModelDykeAllLayersSplit'";
        cout << "\n-------------------------------------------------------";
        string varFileName = "ANSYS_SplitBoundaryMatch_Test-variables.txt";
        model3d_name_ = "ModelDykeAllLayersSplit";
        model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str(), true, true, true, true);
 
        MeshManager<3>& mesh(model3d_->Mesh());
        cout << "\nNodes: " << mesh.Nodes() << "\n";
        set<Element<3>*> elements3;
        cout << "\nNode Groups: " << findContiguousMeshPatch<3,Element>( &(*mesh.ElementsBegin()), elements3 ) << "\n";
        cout << "\nElements: " << mesh.Elements() << "\n";
        std::map<std::string,std::vector<Element<3U>*> > patch_map3;
        cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map3 ) << "\n";
        cout << "\nFaces: " << mesh.Faces() << "\n";
        std::map<std::string,std::vector<Face<3U>*> >  face_map3;
        cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) << "\n";
        cout << "\nInterfaces: " << mesh.InterFaces() << "\n";
        std::map<std::string,std::vector<InterFace<3U>*> >  iface_map3;
        cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) << "\n";
        
        //writing ansys model to file deleting it and then recreating a csmp native model from the file
        if ( reconstruct_from_file ) {
            model3d_->OutputToBinaryFile(model3d_name_.c_str());
            delete model3d_;
            model3d_ = new Model<3U>(model3d_name_);
            MeshManager<3>& mesh(model3d_->Mesh());
            cout << "\nNodes: " << mesh.Nodes() << "\n";
            cout << "\nNode Groups: " << findContiguousMeshPatch<3,Element>( &(*mesh.ElementsBegin()), elements3 ) << "\n";
            cout << "\nElements: " << mesh.Elements() << "\n";
            cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map3 ) << "\n";
            cout << "\nFaces: " << mesh.Faces() << "\n";
            cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) << "\n";
            cout << "\nInterfaces: " << mesh.InterFaces() << "\n";
            cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) << "\n";
          }
        delete model3d_;
        model3d_ = nullptr;
        return;
    }

	// ansys 3d model - contiguous
	cout << "\n-------------------------------------------------------";
	cout << "\nMeshManager_Test: ANSYS model 'prism_test'";
	cout << "\n-------------------------------------------------------";
  string varFileName = "CSMP-variables.txt";
	model3d_name_ = "prism_test";
	model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str());

  Region<3U>&  model_domain = model3d_->Region("Model");
  vector<set<size_t>> node_neighbors;
  nodeNeighbors( model_domain, node_neighbors );
  // determining typical number of node neighbors in mesh
  size_t n_neighbors{0};
  for ( auto nit : node_neighbors ) n_neighbors += nit.size();
  cout <<"\n\taverage number of neighbors per node: "<< n_neighbors / node_neighbors.size();

  set<Element<3>*> elements3;
  MeshManager<3>& mesh(model3d_->Mesh());
	cout << "\nNodes: " << mesh.Nodes() << "\n";
	cout << "\nNode Groups: " << findContiguousMeshPatch<3,Element>( &(*mesh.ElementsBegin()), elements3 ) << "\n";
	cout << "\nElements: " << mesh.Elements() << "\n";
  map<string,vector<Element<3U>*> > patch_map3;
	cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map3 ) << "\n";
	cout << "\nFaces: " << mesh.Faces() << "\n";
  if ( mesh.Faces() > 0 ) {
      map<string,vector<Face<3U>*> >  face_map3;
      cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) << "\n";
    }
	cout << "\nInterfaces: " << mesh.InterFaces() << "\n";
  if ( mesh.InterFaces() > 0 ) {
      map<string,vector<InterFace<3U>*> >  iface_map3;
      cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) << "\n";
    }
    
	if ( reconstruct_from_file ) {
      // writing ansys model to file deleting it and then recreating a csmp native model from the file
      model3d_->OutputToBinaryFile(model3d_name_.c_str());
      delete model3d_;
      model3d_ = new Model<3U>(model3d_name_);
      MeshManager<3>& mesh(model3d_->Mesh());
      Region<3U>&  model_domain2 = model3d_->Region("Model");
      cout << "\nNodes: " << mesh.Nodes() << "\n";
      cout << "\nNode Groups: " << findContiguousMeshPatch<3,Element>( (*model_domain2.ElementsBegin()), elements3 ) << "\n";
      cout << "\nElements: " << mesh.Elements() << "\n";
      cout << "\nElement Groups: " << findStandAloneMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map3 ) << "\n";
      cout << "\nFaces: " << mesh.Faces() << "\n";
      if ( mesh.Faces() > 0 ) {
          map<string,vector<Face<3U>*> >  face_map3;
          cout << "\nFace Groups: " << findStandAloneMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map3 ) << "\n";
        }
      cout << "\nInterfaces: " << mesh.InterFaces() << "\n";
      if ( mesh.InterFaces() > 0 ) {
          map<string,vector<InterFace<3U>*> >  iface_map3;
          cout << "\nInterface Groups: " << findStandAloneMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map3 ) << "\n";
        }
    }
  
 } // end Create_ANSYS3D_Model





/**
THE MASTER TEST FUNCTION
*/
void MeshManager_Test::run()
{
  // basics
  /*
   {
      cout << "\n----------------------------";
      cout << "\nMeshManager_Test::TestBasics";
      cout << "\n----------------------------";
      string varFileName("CSMP-1phase-variables.txt");
      model3d_name_ = "PyramidHexaPatch";
      const bool   skewed_elements(false); // otherwise model is not a box anymore
      VSet<3U>     vset;
      test_Create_Pyramid_Hexa_VSet(vset, skewed_elements);
      model3d_ = new Model<3U>(vset, varFileName.c_str(), true);
      TestBasics();
      delete model3d_;
      model3d_ = nullptr;
   }

  _test(Test_BuiltElementConnectivity2D()); // OK
  _test(Test_BuiltElementConnectivity3D()); // OK
  _test(Test_parentElementsSharedByFace()); // OK
*/

  _test( Test_MeshTraversal3D() );


	cout << "\n------------------------------------------------";
	cout << "\nMeshManager_Test::TestEntityNumberingFunction";
	cout << "\n------------------------------------------------";
	_test(TestEntityNumberingFunction());

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

	_test(TestEraseAllPrimitives());
 
  cout << endl;

} // end run





/**
Checks that numbers of elements etc. in mesh manager do indeed reflect those of input model
*/
void MeshManager_Test::TestBasics()
{
	// vsetMakers: test_Create_Pyramid_Hexa_VSet
	if (model3d_name_ == "PyramidHexaPatch") {
		const MeshManager<3U>& mesh = model3d_->Mesh();

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
	}
	else
		cerr << "\nMeshManager_Test::TestBasics; function only acts on 'test_Create_Pyramid_Hexa_VSet' model.\n";

} // end TestBasics




// method with the same name
bool MeshManager_Test::Test_parentElementsSharedByFace()
 {
    VSet<3U> vset;
    // testing with element 13 with face 4 on the LEFT outside
    test_Create_Prism_Hexa_VSet( vset, false );
    Model<3U>    model( vset, "CSMP-variables.txt", true );
    const size_t ELMT{13}; // 13 in VSet
    Element<3>*  eptr = &(*next(model.Mesh().ElementsBegin(),ELMT));
    //eptr->Out();
    assert( eptr->Neighbor(4) == nullptr );
    
    // getting an inner face in the 3D model that is not on the boundary
    vector<Node<3>*> face_nodes;
    Element<3>*      inner_eptr(nullptr), *outer_eptr(nullptr);
    for ( size_t i{0}; i<eptr->Neighbors(); ++i )
     if ( eptr->Neighbor(i) != nullptr ) {
          vector<size_t> fnids;
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
    vector<size_t> fnids;
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
    backupNeighborConnectivity( model_domain.ElementsBegin(), model_domain.ElementsEnd(), nbor_pointers );
    // rebuilding the connectivity
    model.Mesh().BuildConnectivity<Element>( model_domain.ElementsBegin(), model_domain.ElementsEnd() );
    // getting the connectivity that was recreated
    vector<vector<Element<2>*> > nbor_pointers2;
    backupNeighborConnectivity( model_domain.ElementsBegin(), model_domain.ElementsEnd(), nbor_pointers2 );
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
    backupNeighborConnectivity( model_domain.ElementsBegin(), model_domain.ElementsEnd(), nbor_pointers );
    // rebuilding the connectivity
    model.Mesh().BuildConnectivity<Element>( model_domain.ElementsBegin(), model_domain.ElementsEnd() );
    // getting the connectivity that was recreated
    vector<vector<Element<3>*> > nbor_pointers2;
    backupNeighborConnectivity( model_domain.ElementsBegin(), model_domain.ElementsEnd(), nbor_pointers2 );
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
    set<Element<3U>* const>  discovered_elements;
    
    floodFill( model_domain.E(0), discovered_elements ); // OK
    _test( discovered_elements.size() == model_domain.Elements() );
    
    set<Element<3>*> elements;
    _test( findContiguousMeshPatch<3>( &(*model_domain.E(0)), elements ) == model_domain.Elements() ); // OK
    _test( elements.size() == model_domain.Elements() );
    
    map<string,vector<Element<3>*> > elmt_patches;
    size_t patches = findStandAloneMeshPatches( model.Mesh().ElementsBegin(), model.Mesh().ElementsEnd(),
                                                elmt_patches );
    _test( patches == 1 );
    _test( (*elmt_patches.begin()).second.size() == model_domain.Elements() );
    for ( auto i : elmt_patches ) cout <<" "<< i.first;
    
    // checking that NodesOfSegment() and  CornerNodesPerSegmentForElementOfType() give the same answer
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         for ( size_t i{0}; i < (*it)->FE()->Segments(); ++i ) {
              vector<size_t> node_vec;
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
    vector<set<size_t>>  node_neighbors;
    nodeNeighbors( model_domain, node_neighbors );
    
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
         for ( size_t i{0}; i < (*nit)->Neighbors(); ++i )
           node_neighbors.push_back( (*nit)->Neighbor(i) );
         // rebuilding the connectivity
         _test( (*nit)->ReassignNeighbors() == node_neighbors.size() );
         // comparing the sorted node pointers with one another
         for ( size_t i{0}; i < (*nit)->Neighbors(); ++i )
           _test( (*nit)->Neighbor(i) == node_neighbors[i] );
      }    

    return true;
    
 } // end Test_MeshTraversal3D




bool MeshManager_Test::TestEntityNumberingFunction()
{
  const bool contiguous{true}, reconstruct_from_CSMP_binary_file{true};
  Create_ANSYS3D_Model( contiguous, reconstruct_from_CSMP_binary_file );
	// nodes numbered via Model region
	Region<3U>& model_domain(model3d_->Region("Model"));
 
	model_domain.UpdateMemberIndexes();
	vector<size_t>  node_numbers_Model;
	node_numbers_Model.reserve(model_domain.Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model3d_name_ << "': 'Model' numbered nodes:\n";
	for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); ++nit) {
      node_numbers_Model.push_back((*nit)->Idx());
      if ((*nit)->Idx() % 100 == 0) cout << (*nit)->Idx() << "...";
    }
	cout << "\n";

	// renumbering nodes
	const bool in_a_single_sequence(true);
  MeshManager<3>& mesh(model3d_->Mesh());
	mesh.AssignUniqueNumbers(in_a_single_sequence);
	// checking the numbering
	vector<size_t>  nodes_renumbered;
	nodes_renumbered.reserve(model3d_->Mesh().Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model3d_name_ << "': renumbered nodes:\n";

	{
		set<csmp::Node<3U>*>    discovered_nodes;
		deque<csmp::Node<3U>*>  current_nodes;
		discovered_nodes.insert( &(*mesh.NodesBegin()) );
		current_nodes.push_back( &(*mesh.NodesBegin()) );

		while (!current_nodes.empty()) {
			const csmp::Node<3U>*  n_ptr(*current_nodes.begin());
			for (size_t i = 0U; i < n_ptr->Parents(); i++) {
				for (size_t j = 0U; j < n_ptr->Parent(i)->Nodes(); j++) {
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



bool MeshManager_Test::TestElementDeletionAndInsertion()
{
  MeshManager<3U>& mesh(model3d_->Mesh());
  const size_t     n_original_elmts(mesh.Elements());
	// 0. the first element is copy constructed and stored, and then deleted

	csmp::Element<3U> first_element(*mesh.ElementsBegin());
	// 1. delete elements 1 from the model
	mesh.Erase( mesh.ElementsBegin() );
  
	_test(mesh.Elements() == n_original_elmts - 1);

  const size_t n_original_nodes(mesh.Nodes());
	// 2. create a new element with its nodes
	IsoparametricLinearPyramid fe;
	LocalVariables				     node_vars = model3d_->Database().LocalVariablesAt(NODE);
	LocalVariables				     elmt_vars = model3d_->Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 intp_vars = model3d_->Database().IntegrationPointVariablesAt(ELEMENT);
  // copy the first node
  Node<3U>      n1(*mesh.NodesBegin());
  const size_t  nearby_node(4);
	Node<3U>*		  ptr_n1 = mesh.AddNodeAtUniqueLocation( n1.Coordinate(), nearby_node, node_vars );
  // method must return pointer to node 1 pointer
  _test( ptr_n1 == &(*mesh.NodesBegin()) );

  // create new nodes and return pointers to them
	Node<3U>*		ptr_n2 = mesh.AddNodeAt( Point<3U>(26., 27., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n3 = mesh.AddNodeAt( Point<3U>(29., 30., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n4 = mesh.AddNodeAt( Point<3U>(31., 32., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n5 = mesh.AddNodeAt( Point<3U>(34., 35., 0.0), node_vars, NOT );

	_test( mesh.Nodes() == n_original_nodes + 5 );

  // create new PYRAMID element
  int32_t material_id(1); // new element's rock_tye
  vector<Node<3U>*>  nodes = {ptr_n1,ptr_n2,ptr_n3,ptr_n4,ptr_n5};
  Element<3U>*       neptr( &(*next(mesh.ElementsBegin(),4)) ); // just a neighbor to try
	Element<3U>*	ptr_e1 = mesh.AddElement( ISOPARAMETRIC_LINEAR_PYRAMID, elmt_vars, intp_vars, nodes, material_id );

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




// create Faces at boundary and between elements, converted to InterFace objects and then delete them.
bool MeshManager_Test::TestFaceDeletionAndInsertion()
{
  MeshManager<3U>& mesh(model3d_->Mesh());
  bool             boundary_face_constructed(false);
  bool             interior_face_constructed(false);
  LocalVariables				     fvars = model3d_->Database().LocalVariablesAt(FACE);
	IntegrationPointVariables	 ivars = model3d_->Database().IntegrationPointVariablesAt(FACE);
  const size_t n_original_faces = mesh.Faces();

  // element 4
  Element<3U>* const eptr( &(*next(mesh.ElementsBegin(),4)) );
  Face<3U>*    fptr1(nullptr), *fptr2(nullptr);
  for ( size_t i=0U; i<eptr->Faces(); i++ )
    {
       // if the face is at the boundary, we construct a boundary face
       if ( eptr->Neighbor(i) == nullptr && !boundary_face_constructed ) {
            fptr1 = mesh.AddBoundaryFace( eptr, i, fvars, ivars );
            boundary_face_constructed = true;
         }
       // if the face is within model, we construct a normal face
       if ( eptr->Neighbor(i) != nullptr && !interior_face_constructed ) {
            size_t opposite_face = UNSPECIFIED;
            for ( size_t j{0}; j<eptr->Neighbor(i)->Faces(); ++j ) {
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
 
 // conversion of interior face to interfaces
  LocalVariables				     ifvars = model3d_->Database().LocalVariablesAt(INTER_FACE);
	IntegrationPointVariables	 iivars = model3d_->Database().IntegrationPointVariablesAt(INTER_FACE);
  if ( interior_face_constructed ) {
       InterFace<3U>* ifptr = mesh.ReplaceFaceByInterFace( fptr2, ifvars, iivars );
       ifptr->Out();
    }
	
	// delete the new face(s) again
  vector<Face<3U>*> face_ptrs{ fptr1, fptr2 };
	mesh.Delete( face_ptrs.begin(), face_ptrs.end() );
  _test( fptr1 == nullptr );
  _test( fptr2 == nullptr );
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion: model '" << model3d_name_ << "' (after deletion of faces):\n";
	cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";

	_test( mesh.Faces() == n_original_faces );

	return true;
}







bool MeshManager_Test::TestInterFaceDeletionAndInsertion()
{
  MeshManager<3U>& mesh(model3d_->Mesh());
  bool                       interface_constructed(false);
  LocalVariables				     ifvars = model3d_->Database().LocalVariablesAt(INTER_FACE);
	IntegrationPointVariables	 iivars = model3d_->Database().IntegrationPointVariablesAt(INTER_FACE);
  const size_t n_original_ifaces = mesh.InterFaces();

  // puts interfaces between the interior faces of Element # and Element #
  Element<3U>* const eptr( &(*next(mesh.ElementsBegin(),4)) );
  InterFace<3U>*     ifptr(nullptr);
  vector<InterFace<3U>*> iface_ptrs;
  iface_ptrs.reserve( eptr->Faces() );
  for ( size_t i=0U; i<eptr->Faces(); i++ )
    {
       if ( eptr->Neighbor(i) != nullptr && !interface_constructed ) {
            vector<size_t> fnids;
            // getting the nodes for the inside of the future interface
            eptr->FE()->NodesOfFace( i, fnids );
            vector<Node<3U>*> inside_nodes;
            inside_nodes.reserve( fnids.size() );
            for ( size_t j=0U; j<fnids.size(); ++j )
              inside_nodes.push_back( eptr->N( fnids[j] ) );
            // duplicating these nodes to get nodes for the inside element and the other side
            vector<Node<3U>*> middle_nodes, outside_nodes;
            middle_nodes.reserve( fnids.size() );
            outside_nodes.reserve( fnids.size() );
            for ( size_t j=0U; j<fnids.size(); ++j ) {
                 Node<3U>* mnptr = mesh.Duplicate( eptr->N( fnids[j] ), MIDDLE, ManifoldType::INTERFACE );
                 middle_nodes.push_back( mnptr );
                 Node<3U>* onptr = mesh.Duplicate( eptr->N( fnids[j] ), OUTSIDE, ManifoldType::INTERFACE );
                 outside_nodes.push_back( onptr );
              }
            // find matching faces via the shared nodes
            pair<size_t,size_t> face_ids1 = findAdjacentElementFaces( eptr, eptr->Neighbor(i) );
            _test( i == face_ids1.first );
            // find matching faces via neighbor element pointers (faster)
            pair<size_t,size_t> face_ids2 = findAdjacentFacesFromNeighbors( eptr, eptr->Neighbor(i) );
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
	cout << "\nMeshManager_Test::TestInterFaceDeletionAndInsertion: model '" << model3d_name_ << "' (after deletion of interfaces):\n";
	cout << "\nFaces: " << model3d_->Mesh().InterFaces() << "\n";

	_test( mesh.InterFaces() == n_original_ifaces );

	return true;
}



  
/*
    Erasure of all elements, faces and interfaces from 2d model.
*/
bool MeshManager_Test::TestEraseAllPrimitives()
{
   bool reconstruct_from_CSMP_binary_file{true};
   Create_ANSYS2D_Model( reconstruct_from_CSMP_binary_file );
   MeshManager<2U>& mesh(model2d_->Mesh());
   const size_t     n_original_elmts(mesh.Elements());
   const size_t     n_original_faces(mesh.Faces());
   const size_t     n_orig_interfaces(mesh.InterFaces());
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
   _test( mesh.Erase( mesh.NodesBegin(),   mesh.NodesEnd() )   == n_original_elmts );
   _test( mesh.Erase( mesh.ElementsBegin(),   mesh.ElementsEnd() )   == n_original_elmts );
   _test( mesh.Erase( mesh.FacesBegin(),      mesh.FacesEnd() )      == n_original_faces );
   _test( mesh.Erase( mesh.InterFacesBegin(), mesh.InterFacesEnd() ) == n_orig_interfaces );
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
	 _test( mesh.Elements() + mesh.Faces() + mesh.InterFaces() + mesh.Nodes() == 0);

	return true;
}



/**
      finding the neighbors nodes of each node.
      
      This method is equivalent to creating a sparsity pattern for an accumulation.
      
      @test OK SKM 8/12/21
*/
template<size_t dim>
void nodeNeighbors( const Region<dim>& subdomain, vector<set<size_t>>& node_neighbors )
 {
    if ( !node_neighbors.empty() ) node_neighbors.clear();
    node_neighbors.resize( subdomain.Nodes() );
    
    // 1. get continuous indices to access the sets contained in the node_neighbor vectors
    subdomain.RenumberNodes();
 
    vector<size_t> segm_nodes;

    const auto elmtsEnd{ subdomain.ElementsEnd() };
    for ( auto it=subdomain.ElementsBegin(); it!=elmtsEnd; ++it ) {
         const size_t n_segments{ (*it)->Segments() };
         for ( size_t segm_id{0}; segm_id < n_segments; ++segm_id ) {
              (*it)->FE()->NodesOfSegment( segm_id, segm_nodes );
              // replacing local with global node ids
              for ( auto& sit : segm_nodes ) sit = (*it)->N(sit)->Idx();
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
