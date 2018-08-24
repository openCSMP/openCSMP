//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "MeshManager_Test.h"
#include "vsetMakers.h"
#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp {

/**
     creates test for a vsetMakers model
     uses "CSMP-1phase-variables.txt"
*/
MeshManager_Test::MeshManager_Test()
 : model_name_("PyramidHexaPatch")
 {
    const string varFileName( "CSMP-1phase-variables.txt" );
    const bool   skewed_elements(false); // otherwise model is not a box anymore
    VSet<3U>     vset;
  
    test_Create_Pyramid_Hexa_VSet( vset, skewed_elements );
   
    vset.Out();
  
    model_ = new Model<3U>( vset, varFileName.c_str(), true );
   
    TestBasics();
 }


/**
    custom models 
    uses "CSMP-1phase-variables.txt"
*/
MeshManager_Test::MeshManager_Test( const char* ansys_test_model, bool reconstruct_model_from_CSMP_binary_file )
 : model_name_(ansys_test_model)
 {
    // ansys model
    const string varFileName( "CSMP-1phase-variables.txt" );
    model_ = new ANSYS_Model3D( model_name_.c_str(), varFileName.c_str() );
   
    // writing ansys model to file deleting it and then recreating a csmp native model from the file
    if ( reconstruct_model_from_CSMP_binary_file ) {
         model_->OutputToBinaryFile(model_name_.c_str());
         delete model_;
         model_ = new Model<3U>(model_name_);
      }
 }




/**
     THE MASTER TEST FUNCTION
*/
void MeshManager_Test::run()
 {
    // TestBasics() already called in model constructor
 
    _test( TestEntityNumberingFunction() );
  
//    _test( TestElementDeletionAndInsertion() );
  
//    _test( TestFaceDeletionAndInsertion() );
  
 } // end run





/**
    Checks that numbers of elements etc. in mesh manager do indeed reflect those of input model
*/
void MeshManager_Test::TestBasics()
  {
     // vsetMakers: test_Create_Pyramid_Hexa_VSet
     if ( model_name_ == "PyramidHexaPatch" ) {
          const MeshManager<3U>& mesh = model_->Mesh();
       
         // returns true if the mesh consists of multiple element types
         _test( mesh.HybridElementMesh() == true );
  
         // counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
         set<size_t> test_set;
         _test( mesh.DetectElementsWithAllNodesOnBoundary(test_set) == 0U );

         // returns number of nodes=vertices in the current mesh
         _test( mesh.Nodes() == 65 );
  
         // returns number of elements in the current mesh
         _test( mesh.Elements() == 32 );
  
         // returns number of Faces=lower-dimensional elements in current mesh
         _test( mesh.Faces() == 0 );
  
         // returns number of interfaces=faces with multiplicated nodes
         _test( mesh.InterFaces() == 0 );
       
         // the valdity of mesh connectivity tree
         // _test( mesh.RootNode() )
         _test( mesh.RootElement().FE()   != nullptr );
         if ( mesh.Faces() > 0 )      _test( mesh.RootFace().FE()      != nullptr );
         if ( mesh.InterFaces() > 0 ) _test( mesh.RootInterFace().FE() != nullptr );
       }
     else
     cerr <<"\nMeshManager_Test::TestBasics; function only acts on 'test_Create_Pyramid_Hexa_VSet' model.\n";
    
  } // end TestBasics







 
/**
    Tests whether the elements are correctly numbered consecutively by
    the mesh manager.
    
@code
    /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
    void AssignUniqueNumbers( bool in_a_single_sequence ) const;
 @endcode

*/
bool MeshManager_Test::TestEntityNumberingFunction()
  {
     // nodes numbered via Model region
     Region<3U>& model_domain(model_->Region("Model"));
     model_domain.UpdateMemberIndexes();
     vector<size_t>  node_numbers_Model;
     node_numbers_Model.reserve( model_domain.Nodes() );
     cout <<"\nMeshManager_Test::TestEntityNumberingFunction: model '"<< model_name_ <<"': 'Model' numbered nodes:\n";
     for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
          node_numbers_Model.push_back( (*nit)->Idx() );
          cout << (*nit)->Idx() <<" ";
       }
     cout <<"\n";
  
     // renumbering nodes
     const bool in_a_single_sequence(true);
     model_->Mesh().AssignUniqueNumbers( in_a_single_sequence );
     // checking the numbering
     vector<size_t>  nodes_renumbered;
     nodes_renumbered.reserve( model_->Mesh().Nodes() );
     cout <<"\nMeshManager_Test::TestEntityNumberingFunction: model '"<< model_name_ <<"': renumbered nodes:\n";
     for ( auto nit=model_->Mesh().NodesBegin(); nit!=model_->Mesh().NodesEnd(); ++nit ) {
          nodes_renumbered.push_back( (*nit).Idx() );
          cout << (*nit).Idx() <<" ";
       }
     cout <<"\n";
    
     _test( equal( node_numbers_Model.begin(), node_numbers_Model.end(), nodes_renumbered.begin(), nodes_renumbered.end() ) );
    
     // additing some nodes (at the end) and deleting the same number again
     // Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY=NOT );
     LocalVariables  node_vars = model_->Database().LocalVariablesAt( NODE );
     Node<3U>  extra_node1(66,Point<3>(23.,24.,25.),node_vars,NOT),
               extra_node2(67,Point<3>(26.,27.,28.),node_vars,NOT);
    
     Node<3U>* nd_ptr1 = model_->Mesh().PushBack( Node<3U>(66,Point<3>(23.,24.,25.),node_vars,NOT) );
     Node<3U>* nd_ptr2 = model_->Mesh().PushBack( Node<3U>(67,Point<3>(26.,27.,28.),node_vars,NOT) );
     // fails to compile: model_->Mesh().EmplaceNode( Node<3U>(68,Point<3>(26.,27.,28.),node_vars,NOT), Node<3U>(69,Point<3>(26.,27.,28.),node_vars,NOT) );
    
     // deleting these nodes again
     model_->Mesh().Erase( *nd_ptr1 );
     model_->Mesh().Erase( *nd_ptr2 );
    
     // checking the numbering (after insertion and deletion)
     vector<size_t>  nodes_resized;
     nodes_resized.reserve( model_->Mesh().Nodes() );
     cout <<"\nMeshManager_Test::TestEntityNumberingFunction: model '"<< model_name_ <<"': renumbered nodes (after insertion and deletion of nodes):\n";
     for ( auto nit=model_->Mesh().NodesBegin(); nit!=model_->Mesh().NodesEnd(); ++nit ) {
          nodes_resized.push_back( (*nit).Idx() );
          cout << (*nit).Idx() <<" ";
       }
     cout <<"\n";
 
     _test( equal( node_numbers_Model.begin(), node_numbers_Model.end(), nodes_resized.begin(), nodes_resized.end() ) );

     return true;
  }

} // end csmp
