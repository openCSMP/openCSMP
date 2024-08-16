//
//  MeshManager_Test.hpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MESH_MANAGER_TEST_H
#define CSMP_MESH_MANAGER_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

/**
   Testing checklist
   - buid mesh from VSet
   - verify node connectivity (nodes connected to eachother)
   - verify node numbering through basic calculations
   - verify neighbor connectivity (all volume elements with eachother, all surface ...)
   - creation and deletion of Face objects
   


  Test Dependencies
    @todo FiniteElementManager_Test
    @todo FiniteVolumeManager_Test

    @author SKM
    @date 23/8/2018
*/
class MeshManager_Test : public Test {
  public:
    MeshManager_Test();
    virtual ~MeshManager_Test() { delete model3D_ptr_; delete model2D_ptr_; }
    
    virtual void run();
  
  private:
  
    const bool verbose_ = true; // turn on or off for detailed reporting
 
    /// constructs and returns 3D model "as is" from ANSYS input files not trying to make any corrections (except for ignoring broken pfverts record) or construct boundaries
    void ConstructANSYS_ModelWithoutModification( const char* model_name );

    // TESTING PUBLIC INTERFACES OF MESH MANAGER
    
    /// Tests: Nodes(), Elements(), Faces(), InterFaces(), HybridElementMesh(), IsContiguous(), OutputMeshTo(vset),  HasNodeManifolds(),
    /// and node-nbor connectivity, checking that mesh matches what was in the VSet
    void TestBasics();
    
    // not tested but used: FiniteElements(), FiniteVolumes(), InitializeFiniteVolumeStencils(-tested elsewhere) iterators
    
    /// uses ANSYS 'DykePartiallySplit' (discontiguous) and vsetMaker models in the test; returns true if tests pass
    bool Test_IsContiguous();
    
    // SUPPORTING FUNCTIONALITY FROM MESH_MANAGEMENT_UTILITIES
    bool Test_connectNeighborsUsingNodeParents(); // TODO: make separate test for all these functions
    
    // hex 27 VSet - disconnected Neighbors from central element
    bool Test_detachNeighborsFrom();
    
    // TODO: test the uncommented ones...
    //bool Test_AddNodeAt();                  // <- TestCellDeletionAndInsertion() tests this
    //bool Test_Duplicate(); // node
    //bool Test_AddElement();                 // <- TestCellDeletionAndInsertion() tests this
    //bool Test_AddInterveningElement();
    //bool Test_ReplaceElementByFace();       // <- TestCellDeletionAndInsertion() tests this
    bool Test_ReplaceElementByInterFace();
    //bool Test_AddFace();
    bool Test_AddEdgeFace();
    //bool Test_AddBoundaryFace();            // <- TestCellDeletionAndInsertion() tests this
    bool Test_ReplaceFaceByInterFace();
    //bool Test_AddInterFace(); // two versions
    
    // MESH MODIFICATION (elements to faces etc.)

    bool TestCellDeletionAndInsertion();      // done
    //bool TestElementDeletionAndInsertion();
    bool TestFaceDeletionAndInsertion();      // done
    bool TestInterFaceDeletionAndInsertion(); // done

    bool Test_ReplaceInteriorElementsByFaces();
    bool Test_ReplaceBoundaryElementsByFaces();
    bool Test_CreateFacesBetweenNodeSharingElements();
    bool Test_ReplaceFacesByInterFaces();
    bool Test_ReplaceElementsByInterFaces();
    bool Test_CreateInterfacesBetweenNodeSharingElements();
    bool Test_CreateInterfacesBetweenNodeMatchingElements();

    bool Test_BuildConnectivity(); // test UpdateConnectivity() done which calls these
    //bool Test_BuildVolumeConnectivity();
    //bool Test_BuildSurfaceConnectivity();
    //bool Test_BuildLineConnectivity();
    //bool Test_BuildInterFaceConnectivity();
    bool Test_ConnectNodesToParentsAndNeighbors();
    
    void Test_DeleteNodesAndRepairNodeConnnectivity(); // done
    
    bool Test_DeleteCellsAndRepairConnnectivity(); // TODO: refactor: 4 versions for Node, Element, Face, InterFace - split ino RepairConnectivity()
    
    bool Test_CheckElementConnectivity();
    bool Test_OutputStoredVariablesTo();   // covered by binary file saving and reading tests
    bool Test_InputStoredVariablesFrom();  // covered by binary file saving and reading tests
    
    // TODO: test RepairConnectivity( input from SplitB generation, nodes, nbors with new nodes etc.
    
    // indirectly tested underlying functionality
    /*
       plf_colony - there is a separate unit test for it
       AssignUniqueNumbers();
       ReorderObjectsByIndexes();
       OutputMeshTo(vset); // done many times over in next test section
    */


    // COMPOSITE FUNCTIONALITY
    
    /// check that the nodes of an element face are the ones shared with the corresponding element neighbor
    bool TestNeigbourVersusFaceConsistency();

    /// tests BuildConnectivityMethod() that creates neighbor and node connectivity
    bool Test_BuiltElementConnectivity2D(); // INCOMPLETE
    bool Test_BuiltElementConnectivity3D(); // INCOMPLETE
    
    // TODO: add missing test:  DeleteCellsAndRepairConnnectivity()

    // TODO: add missing test:
    template<uint32_t dim>
    bool TestNodeNeighborConnectivity( const Model<dim>& );

    /// tests floodfill algorithm (in MeshManagementUtilities) etc.
    bool Test_MeshTraversal3D();

    /// JCK traversal of mesh tree, numbering Node, Element etc.
    bool TestEntityNumberingFunction( Model<3>& );

    /// tests building of CONTIGUOUS vs DISCONTIGUOUS models (consisting of isolated mesh patches) and whether all nodes, elements etc can be reached
    bool CheckConnectivityOfModel3D( Model<3>& ); // via mesh traversal, includes call to CheckElementConnectivity()
    
    /// tests model saving with SPLIT22_BASIC (rectangle with through-going subhorizontal split boundary and oblique s.b. cutting it)
    bool TestSavedElementFaceInterfaceModel2D();
 
    Model<3U>* model3D_ptr_ = nullptr;
    Model<2U>* model2D_ptr_ = nullptr;
 
 // JCK STUFF - DEPRECATE?
 
// TODO: meant to test erasure of elements (JCK)
    bool TestEraseAllPrimitives();

};

} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
