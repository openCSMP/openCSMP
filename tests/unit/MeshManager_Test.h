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
    virtual ~MeshManager_Test() {}
    
    virtual void run();
  
  private:
  
    void TestBasics();
    
    bool TestEntityNumberingFunction();
    bool TestElementDeletionAndInsertion();
    bool TestFaceDeletionAndInsertion();
    bool TestInterFaceDeletionAndInsertion();
    bool TestEraseAllPrimitives();
    
    // method with the same name
    bool Test_parentElementsSharedByFace();

    // create test models that are subsequently used for the testing
    void Create_ANSYS2D_Model( bool reconstruct_from_CSMP_binary_file );
    void Create_ANSYS3D_Model( bool contiguous, bool reconstruct_from_CSMP_binary_file );

    // checks whether all nodes, elements etc can be reached
    bool CheckConnectivityOfModel3D( Model<3>& );
    
    template<uint32_t dim>
    bool TestNodeNeighborConnectivity( const Model<dim>& );
    
    // using VSetMakers to create and compare input data
    bool Test_BuiltElementConnectivity2D();
    bool Test_BuiltElementConnectivity3D();
    
    // floodfill etc.
    bool Test_MeshTraversal3D();

};

// build sparsity pattern for testing the connectivity among nodes
template<uint32_t dim> void nodeNeighbors( const Region<dim>&, std::vector<std::set<size_t>>& node_neighbors );

} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
