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
  
    const bool verbose_ = true; // turn on or off for detailed reporting
  
    void TestBasics();
    
    /// check that the nodes of an element face are the ones shared with the corresponding element neighbor
    bool TestNeigbourVersusFaceConsistency();
 
    bool TestElementDeletionAndInsertion();
    bool TestFaceDeletionAndInsertion();
    bool TestInterFaceDeletionAndInsertion();

    /// rectangle with through-going subhorizontal split boundary and oblique s.b. cutting it
    bool TestCompleteModel2D();
    
    bool TestEntityNumberingFunction( Model<3>& );
    
    bool TestEraseAllPrimitives();

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

} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
