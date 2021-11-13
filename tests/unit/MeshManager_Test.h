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
    @author SKM
    @date 23/8/2018
*/
class MeshManager_Test : public Test {
  public:
    MeshManager_Test();
    virtual ~MeshManager_Test() { delete model2d_; delete model3d_; }
    
    virtual void run();
  
  private:
    // create test models that are subsequently used for the testing
    void Create_ANSYS2D_Model( bool reconstruct_from_CSMP_binary_file );
    void Create_ANSYS3D_Model( bool contiguous, bool reconstruct_from_CSMP_binary_file );
    
    void TestBasics();
    bool TestEntityNumberingFunction();
    bool TestElementDeletionAndInsertion();
    bool TestFaceDeletionAndInsertion();
    bool TestInterFaceDeletionAndInsertion();
    bool TestEraseAllPrimitives();
    
    // method with the same name
    bool Test_parentElementsSharedByFace();
    
    // using VSetMakers to create and compare input data
    bool Test_BuiltElementConnectivity2D();
    bool Test_BuiltElementConnectivity3D();

  private:
    std::string model2d_name_;
    std::string model3d_name_;

    Model<2U>*  model2d_ = nullptr;
    Model<3U>*  model3d_ = nullptr;
};


} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
