//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

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
MeshManager_Test::MeshManager_Test( bool reconstruct_model_from_CSMP_binary_file )
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

	// ansys 2d model - contiguous
	cout << "\n------------------------------------------";
	cout << "\nMeshManager_Test: ANSYS model 'box2d_fault'";
	cout << "\n------------------------------------------";
	model2d_name_ = "box2d_fault";
	model2d_ = new ANSYS_Model2D(model2d_name_.c_str(), varFileName.c_str());
	cout << "\nNodes: " << model2d_->Mesh().Nodes() << "\n";
  std::deque<Element<2>*> elements;
	cout << "\nInterconnected Nodes: " << findContiguousMeshPatch( model2d_->Mesh().N(0U), elements ) << "\n";
	cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
  std::map<std::string,std::deque<Element<2U>*> > patch_map;
	cout << "\nElement Groups: " << findStandAloneMeshPatches( model2d_->Mesh().ElementsBegin(), model2d_->Mesh().ElementsEnd(), patch_map ) << "\n";
	cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
  std::map<std::string,std::deque<Face<2U>*> >  face_map;
	cout << "\nFace Groups: " << findStandAloneMeshPatches( model2d_->Mesh().FacesBegin(), model2d_->Mesh().FacesEnd(), face_map ) << "\n";
	cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
  std::map<std::string,std::deque<InterFace<2U>*> >  iface_map;
	cout << "\nInterface Groups: " << findStandAloneMeshPatches( model2d_->Mesh().InterFacesBegin(), model2d_->Mesh().InterFacesEnd(), iface_map ) << "\n";
	if (reconstruct_model_from_CSMP_binary_file) {
		model2d_->OutputToBinaryFile(model2d_name_.c_str());
		delete model2d_;
		model2d_ = new Model<2U>(model2d_name_);
		cout << "\nNodes: " << model2d_->Mesh().Nodes() << "\n";
		cout << "\nNode Groups: " << findContiguousMeshPatch( model2d_->Mesh().N(0U), elements ) << "\n";
		cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
		cout << "\nElement Groups: " << findStandAloneMeshPatches( model2d_->Mesh().ElementsBegin(), model2d_->Mesh().ElementsEnd(), patch_map ) << "\n";
		cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
		cout << "\nFace Groups: " << findStandAloneMeshPatches( model2d_->Mesh().FacesBegin(), model2d_->Mesh().FacesEnd(), face_map ) << "\n";
		cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
		cout << "\nInterface Groups: " << findStandAloneMeshPatches( model2d_->Mesh().InterFacesBegin(), model2d_->Mesh().InterFacesEnd(), iface_map ) << "\n";
	}

	// ansys 3d model - discontiguous
	cout << "\n-------------------------------------------------------";
	cout << "\nMeshManager_Test: ANSYS model 'ModelDykeAllLayersSplit'";
	cout << "\n-------------------------------------------------------";
	varFileName = "ANSYS_SplitBoundaryMatch_Test-variables.txt";
	model3d_name_ = "ModelDykeAllLayersSplit";
	model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str(), true, true, true, true);
	cout << "\nNodes: " << model3d_->Mesh().Nodes() << "\n";
  std::deque<Element<3>*> elements3;
	cout << "\nNode Groups: " << findContiguousMeshPatch( model3d_->Mesh().N(0U), elements3 ) << "\n";
	cout << "\nElements: " << model3d_->Mesh().Elements() << "\n";
  std::map<std::string,std::deque<Element<3U>*> > patch_map3;
	cout << "\nElement Groups: " << findStandAloneMeshPatches( model3d_->Mesh().ElementsBegin(), model3d_->Mesh().ElementsEnd(), patch_map3 ) << "\n";
	cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";
  std::map<std::string,std::deque<Face<3U>*> >  face_map3;
	cout << "\nFace Groups: " << findStandAloneMeshPatches( model3d_->Mesh().FacesBegin(), model3d_->Mesh().FacesEnd(), face_map3 ) << "\n";
	cout << "\nInterfaces: " << model3d_->Mesh().InterFaces() << "\n";
  std::map<std::string,std::deque<InterFace<3U>*> >  iface_map3;
	cout << "\nInterface Groups: " << findStandAloneMeshPatches( model3d_->Mesh().InterFacesBegin(), model3d_->Mesh().InterFacesEnd(), iface_map3 ) << "\n";
	//writing ansys model to file deleting it and then recreating a csmp native model from the file
	if (reconstruct_model_from_CSMP_binary_file) {
		model3d_->OutputToBinaryFile(model3d_name_.c_str());
		delete model3d_;
		model3d_ = new Model<3U>(model3d_name_);
		cout << "\nNodes: " << model3d_->Mesh().Nodes() << "\n";
		cout << "\nNode Groups: " << findContiguousMeshPatch( model3d_->Mesh().N(0U), elements3 ) << "\n";
		cout << "\nElements: " << model3d_->Mesh().Elements() << "\n";
		cout << "\nElement Groups: " << findStandAloneMeshPatches( model3d_->Mesh().ElementsBegin(), model3d_->Mesh().ElementsEnd(), patch_map3 ) << "\n";
		cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";
		cout << "\nFace Groups: " << findStandAloneMeshPatches( model3d_->Mesh().FacesBegin(), model3d_->Mesh().FacesEnd(), face_map3 ) << "\n";
		cout << "\nInterfaces: " << model3d_->Mesh().InterFaces() << "\n";
		cout << "\nInterface Groups: " << findStandAloneMeshPatches( model3d_->Mesh().InterFacesBegin(), model3d_->Mesh().InterFacesEnd(), iface_map3 ) << "\n";
	}
	delete model3d_;

	// ansys 3d model - contiguous
	cout << "\n-------------------------------------------------------";
	cout << "\nMeshManager_Test: ANSYS model 'prism_test'";
	cout << "\n-------------------------------------------------------";
	model3d_name_ = "prism_test";
	model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str());
	cout << "\nNodes: " << model3d_->Mesh().Nodes() << "\n";
	cout << "\nNode Groups: " << findContiguousMeshPatch( model3d_->Mesh().N(0U), elements3 ) << "\n";
	cout << "\nElements: " << model3d_->Mesh().Elements() << "\n";
	cout << "\nElement Groups: " << findStandAloneMeshPatches( model3d_->Mesh().ElementsBegin(), model3d_->Mesh().ElementsEnd(), patch_map3 ) << "\n";
	cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";
	cout << "\nFace Groups: " << findStandAloneMeshPatches( model3d_->Mesh().FacesBegin(), model3d_->Mesh().FacesEnd(), face_map3 ) << "\n";
	cout << "\nInterfaces: " << model3d_->Mesh().InterFaces() << "\n";
	cout << "\nInterface Groups: " << findStandAloneMeshPatches( model3d_->Mesh().InterFacesBegin(), model3d_->Mesh().InterFacesEnd(), iface_map3 ) << "\n";

	//writing ansys model to file deleting it and then recreating a csmp native model from the file
	if (reconstruct_model_from_CSMP_binary_file) {
		model3d_->OutputToBinaryFile(model3d_name_.c_str());
		delete model3d_;
		model3d_ = new Model<3U>(model3d_name_);
		cout << "\nNodes: " << model3d_->Mesh().Nodes() << "\n";
		cout << "\nNode Groups: " << findContiguousMeshPatch( model3d_->Mesh().N(0U), elements3 ) << "\n";
		cout << "\nElements: " << model3d_->Mesh().Elements() << "\n";
		cout << "\nElement Groups: " << findStandAloneMeshPatches( model3d_->Mesh().ElementsBegin(), model3d_->Mesh().ElementsEnd(), patch_map3 ) << "\n";
		cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";
		cout << "\nFace Groups: " << findStandAloneMeshPatches( model3d_->Mesh().FacesBegin(), model3d_->Mesh().FacesEnd(), face_map3 ) << "\n";
		cout << "\nInterfaces: " << model3d_->Mesh().InterFaces() << "\n";
		cout << "\nInterface Groups: " << findStandAloneMeshPatches( model3d_->Mesh().InterFacesBegin(), model3d_->Mesh().InterFacesEnd(), iface_map3 ) << "\n";
	}
}




/**
THE MASTER TEST FUNCTION
*/
void MeshManager_Test::run()
{
  TestBasics();
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

	cout << "\n-------------------------------------------";
	cout << "\nMeshManager_Test::TestEraseAllPrimitives";
	cout << "\n-------------------------------------------";
	_test(TestEraseAllPrimitives());

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





// TODO: needs to operate on model with faces and interfaces
bool MeshManager_Test::TestEntityNumberingFunction()
{
	// nodes numbered via Model region
	Region<3U>& model_domain(model3d_->Region("Model"));
	model_domain.UpdateMemberIndexes();
	vector<size_t>  node_numbers_Model;
	node_numbers_Model.reserve(model_domain.Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model3d_name_ << "': 'Model' numbered nodes:\n";
	for (auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); ++nit) {
		node_numbers_Model.push_back((*nit)->Idx());
		if ((*nit)->Idx() % 100 == 0) cout << (*nit)->Idx() << "...";
	}
	cout << "\n";

	// renumbering nodes
	const bool in_a_single_sequence(true);
	model3d_->Mesh().AssignUniqueNumbers(in_a_single_sequence);
	// checking the numbering
	vector<size_t>  nodes_renumbered;
	nodes_renumbered.reserve(model3d_->Mesh().Nodes());
	cout << "\nMeshManager_Test::TestEntityNumberingFunction: model '" << model3d_name_ << "': renumbered nodes:\n";

	{
		set<csmp::Node<3U>*>    discovered_nodes;
		deque<csmp::Node<3U>*>  current_nodes;
		discovered_nodes.insert(model3d_->Mesh().N(0));
		current_nodes.push_back(model3d_->Mesh().N(0));

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

		for (auto nit : discovered_nodes)
			nodes_renumbered.push_back(nit->Idx());

		std::sort(nodes_renumbered.begin(), nodes_renumbered.end(), [](auto& lhs, auto& rhs) {return lhs < rhs; });

		for (auto nit : nodes_renumbered)
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

	csmp::Element<3U> first_element( (*mesh.E(0)) );
	// 1. delete elements 1 from the model
	mesh.Delete( mesh.ElementsBegin(), next(mesh.ElementsBegin(),1) );
  
	_test(mesh.Elements() == n_original_elmts - 1);

  const size_t n_original_nodes(mesh.Nodes());
	// 2. create a new element with its nodes
	IsoparametricLinearPyramid fe;
	LocalVariables				     node_vars = model3d_->Database().LocalVariablesAt(NODE);
	LocalVariables				     elmt_vars = model3d_->Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 intp_vars = model3d_->Database().IntegrationPointVariablesAt(ELEMENT);
  // copy the first node
  Node<3U>    n1( *mesh.N(0) );
  const size_t nearby_node(4);
	Node<3U>*		ptr_n1 = mesh.AddNodeAtUniqueLocation( n1.Coordinate(), nearby_node, node_vars );
  // method must return pointer to node 1 pointer
  _test( ptr_n1 == mesh.N(0) );

  // create new nodes and return pointers to them
	Node<3U>*		ptr_n2 = mesh.AddNodeAt( Point<3U>(26., 27., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n3 = mesh.AddNodeAt( Point<3U>(29., 30., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n4 = mesh.AddNodeAt( Point<3U>(31., 32., 0.0), node_vars, NOT );
	Node<3U>*		ptr_n5 = mesh.AddNodeAt( Point<3U>(34., 35., 0.0), node_vars, NOT );

	_test( mesh.Nodes() == n_original_nodes + 5 );

  // create new PYRAMID element
  int32 material_id(1); // new element's rock_tye
  vector<Node<3U>*>    nodes = {ptr_n1,ptr_n2,ptr_n3,ptr_n4,ptr_n5};
  Element<3U>*         neptr( mesh.E(4) ); // just a neighbor to try
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
	mesh.Delete<Node>( nodes.begin(), nodes.end() );
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
  Element<3U>* const eptr( mesh.E(4) );
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
	mesh.Delete( next(mesh.FacesBegin(),n_original_faces), mesh.FacesEnd() );
  mesh.EraseNullPointerCells();
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
  Element<3U>* const eptr( mesh.E(4) );
  InterFace<3U>*     ifptr(nullptr);
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
            const int32 material_id(5);
            Element<3U>*	ieptr = mesh.AddElement( ISOPARAMETRIC_LINEAR_TRIANGLE, ifvars, iivars, middle_nodes, material_id );
            
            // connecting the InterFace to the middle element
            ifptr->Assign( ieptr );
         }
    }
  cout<<"\n\tcreated faces: ";
  if ( interface_constructed ) ifptr->Out();
    
	_test( mesh.Faces() == n_original_ifaces + interface_constructed );
 
	// delete the new interface(s) again
	mesh.Delete( next(mesh.InterFacesBegin(),n_original_ifaces), mesh.InterFacesEnd() );
  mesh.EraseNullPointerCells();
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
   MeshManager<2U>& mesh(model2d_->Mesh());
   const size_t     n_original_elmts(mesh.Elements());
   const size_t     n_original_faces(mesh.Faces());
   const size_t     n_orig_interfaces(mesh.InterFaces());
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
   _test( mesh.Delete( mesh.ElementsBegin(),   mesh.ElementsEnd() )   == n_original_elmts );
   _test( mesh.Delete( mesh.FacesBegin(),      mesh.FacesEnd() )      == n_original_faces );
   _test( mesh.Delete( mesh.InterFacesBegin(), mesh.InterFacesEnd() ) == n_orig_interfaces );
   _test( mesh.EraseNullPointerCells() == n_original_elmts );
	 cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
	 _test( mesh.Elements() + mesh.Faces() + mesh.InterFaces() + mesh.Nodes() == 0);

	return true;
}

} // end csmp
