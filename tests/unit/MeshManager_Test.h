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
    MeshManager_Test( const char* ansys_test_model, bool reconstruct_model_from_CSMP_binary_file );
    MeshManager_Test( /* vset maker based model */ );
    ~MeshManager_Test() { delete model_; }
    virtual void run();
  
  private:
    void TestBasics();
    bool TestEntityNumberingFunction();
    bool TestElementDeletionAndInsertion();
    bool TestFaceDeletionAndInsertion();
    bool TestMemoryManagement();
  
  private:
    std::string model_name_;
    Model<3U>*  model_ = nullptr;
};


} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
