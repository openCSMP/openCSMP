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
#include "Element.h"

#include "IsoparametricLinearPyramid.h"
#include "VTK_Interface.h"
#include "ModelTopology.h"
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
	// variables file
	string varFileName("CSMP-1phase-variables.txt");
	model3d_name_ = "PyramidHexaPatch";
	const bool   skewed_elements(false); // otherwise model is not a box anymore
	VSet<3U>     vset;
	test_Create_Pyramid_Hexa_VSet(vset, skewed_elements);
	vset.Out();
	model3d_ = new Model<3U>(vset, varFileName.c_str(), true);
	TestBasics();	
	delete model3d_;

	// ansys 2d model - contiguous
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
	// check the numbering of the nodes and elements from the mesh
	set<csmp::Element<3U>*> discovered_elmts;
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
					if (new_node.second) current_nodes.push_back(n_ptr->Parent(i)->N(j));
				}
			}
			discovered_elmts.insert(n_ptr->Parent(i));
		}
		current_nodes.pop_front();
	}
	// 1. delete an element from the model
	// the first element from the model is stored into a new object, and then delete it
	auto it = discovered_elmts.begin();
	csmp::Element<3U> first_element(*(*it));

	model3d_->Mesh().Erase(*it);
	_test(model3d_->Mesh().Elements() == discovered_elmts.size() - 1);

	model3d_->Mesh().Add( csmp::Element<3U>(first_element) );
	_test(model3d_->Mesh().Elements() == discovered_elmts.size());

	// 2. create a new element with its nodes
	IsoparametricLinearPyramid fe;
	LocalVariables				     node_vars = model3d_->Database().LocalVariablesAt(NODE);
	LocalVariables				     elmt_vars = model3d_->Database().LocalVariablesAt(ELEMENT);
	IntegrationPointVariables	 intp_vars = model3d_->Database().IntegrationPointVariablesAt(ELEMENT);

	Node<3U>&	n1 = *(model3d_->Mesh().RootNode(0));
	Node<3U>*		ptr_n1 = model3d_->Mesh().AddIfUnique(n1); // if it is already in the mesh, returns the existing node's pointer
	Node<3U>*		ptr_n2 = model3d_->Mesh().Add( Node<3U>(model3d_->Mesh().Nodes() + 1, Point<3U>(26., 27., 0.0), node_vars, NOT) ); // otherwise, create new node and return its pointer
	Node<3U>*		ptr_n3 = model3d_->Mesh().Add( Node<3U>(model3d_->Mesh().Nodes() + 1, Point<3U>(29., 30., 0.0), node_vars, NOT) );
	Node<3U>*		ptr_n4 = model3d_->Mesh().Add( Node<3U>(model3d_->Mesh().Nodes() + 1, Point<3U>(31., 32., 0.0), node_vars, NOT) );
	Node<3U>*		ptr_n5 = model3d_->Mesh().Add( Node<3U>(model3d_->Mesh().Nodes() + 1, Point<3U>(34., 35., 0.0), node_vars, NOT) );

  size_t elmt_id = model3d_->Mesh().Elements(); // new element's id is equal to the number of the existing elements in the model
	Element<3U>*	ptr_e1 = model3d_->Mesh().Add( csmp::Element<3U>(elmt_id, &fe, elmt_vars, intp_vars, NOT) );

	// assign node connectivity where it is connected one of the last element's nodes
	ptr_n1->ResizeParentStorage(ptr_n1->Parents() + 1);
	ptr_n2->ResizeParentStorage(1);
	ptr_n3->ResizeParentStorage(1);
	ptr_n4->ResizeParentStorage(1);
	ptr_n5->ResizeParentStorage(1);

	// assign node's parent element
	ptr_n1->Assign(n1.Parents() - 1, ptr_e1);
	ptr_n2->Assign(0, ptr_e1);
	ptr_n3->Assign(0, ptr_e1);
	ptr_n4->Assign(0, ptr_e1);
	ptr_n5->Assign(0, ptr_e1);

	// assign new element's neighbors
	ptr_e1->Assign(0, ptr_n1);
	ptr_e1->Assign(1, ptr_n2);
	ptr_e1->Assign(2, ptr_n3);
	ptr_e1->Assign(3, ptr_n4);
	ptr_e1->Assign(4, ptr_n5);
	ptr_e1->Assign(0, model3d_->Mesh().RootNode(0)->Parent(0));

	// checking the total number of nodes and elements (after insertion)
	std::cout << "\nNodes    discovered   after insertion: " << discovered_nodes.size();
	std::cout << "\nElements discovered   after insertion: " << discovered_elmts.size();
	std::cout << "\nNodes    of the model after insertion: " << model3d_->Mesh().Nodes();
	std::cout << "\nElements of the model after insertion: " << model3d_->Mesh().Elements();

	_test(model3d_->Mesh().Nodes() == (discovered_nodes.size() + 4));
	_test(model3d_->Mesh().Elements() == (discovered_elmts.size() + 1));

	// 2. deleting these nodes again
	model3d_->Mesh().Erase(ptr_e1);

	// checking the total number of nodes and elements (after deletion)
	std::cout << "\nNodes    discovered   after deletion: " << discovered_nodes.size();
	std::cout << "\nElements discovered   after deletion: " << discovered_elmts.size();
	std::cout << "\nNodes    of the model after deletion: " << model3d_->Mesh().Nodes();
	std::cout << "\nElements of the model after deletion: " << model3d_->Mesh().Elements();

	_test(model3d_->Mesh().Elements() == discovered_elmts.size());

	discovered_elmts.clear();
	discovered_nodes.clear();
	current_nodes.clear();

	return true;
}



bool MeshManager_Test::TestFaceDeletionAndInsertion()
{
	// delete a face 
	csmp::Face<3U>*			target_face(NULL); // face to be deleted
	set<csmp::Face<3U>*>	discovered_faces;
	deque<csmp::Face<3U>*>	current_faces;
	for (size_t i = 0U; i < model3d_->Mesh().FaceGroups(); i++) {
		auto root_face = model3d_->Mesh().RootFace(i);
		discovered_faces.insert(root_face);
		current_faces.push_back(root_face);
		while (!current_faces.empty()) {
			csmp::Face<3U>*  n_ptr(*current_faces.begin());
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;
				if (target_face == NULL) target_face = n_ptr->Neighbor(i);
				auto new_face = discovered_faces.insert(n_ptr->Neighbor(i));
				if (new_face.second) current_faces.push_back(n_ptr->Neighbor(i));
			}
			current_faces.pop_front();
		}
	}
	
	_test(discovered_faces.size() == model3d_->Mesh().Faces());
	
	// delete the target face
	model3d_->Mesh().Erase(target_face);
	cout << "\nMeshManager_Test::TestFaceDeletionAndInsertion: model '" << model3d_name_ << "' (after deletion of faces):\n";
	cout << "\nFaces: " << model3d_->Mesh().Faces() << "\n";
	cout << "\nFace Groups: " << model3d_->Mesh().FaceGroups() << "\n";

	_test(discovered_faces.size() - 1 == model3d_->Mesh().Faces());

	return true;
}


bool MeshManager_Test::TestInterFaceDeletionAndInsertion()
{
	if (model3d_->Mesh().InterFaces() == 0) return true;

	// delete a interface 
	csmp::InterFace<3U>*		target_interface(NULL); // interface to be deleted
	set<csmp::InterFace<3U>*>	discovered_interfaces;
	deque<csmp::InterFace<3U>*>	current_interfaces;
	for (size_t i = 0U; i < model3d_->Mesh().InterFaceGroups(); i++) {
		auto root_interface = model3d_->Mesh().RootInterFace(i);
		discovered_interfaces.insert(root_interface);
		current_interfaces.push_back(root_interface);
		while (!current_interfaces.empty()) {
			csmp::InterFace<3U>*  n_ptr(*current_interfaces.begin());
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;
				if (target_interface == NULL) target_interface = n_ptr->Neighbor(i);
				auto new_interface = discovered_interfaces.insert(n_ptr->Neighbor(i));
				if (new_interface.second) current_interfaces.push_back(n_ptr->Neighbor(i));
			}
			current_interfaces.pop_front();
		}
	}
	csmp::InterFace<3U> new_interface(*target_interface);

	_test(discovered_interfaces.size() == model3d_->Mesh().InterFaces());

	model3d_->Mesh().Erase(target_interface);
	cout << "\nMeshManager_Test::TestInterFaceDeletionAndInsertion: model '" << model3d_name_ << "' (after deletion of interfaces):\n";
	cout << "\nInterfaces: " << model3d_->Mesh().InterFaces() << "\n";
	cout << "\nInterface Groups: " << model3d_->Mesh().InterFaceGroups() << "\n";

	_test(discovered_interfaces.size() - 1 == model3d_->Mesh().InterFaces());

	return true;
}




bool MeshManager_Test::TestEraseAllPrimitives()
{
	bool ret = false;
	ret = model3d_->Mesh().EraseElements();
	cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
	_test(model3d_->Mesh().Elements() == 0);

	return true;
}

} // end csmp
