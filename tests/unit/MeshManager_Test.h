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
    MeshManager_Test( bool reconstruct_model_from_CSMP_binary_file );
    MeshManager_Test( /* vset maker based model */ );
    ~MeshManager_Test() { delete model3d_; }
    virtual void run();
  
  private:
    void TestBasics();
    
	bool TestEntityNumberingFunction_2D();
	bool TestElementDeletionAndInsertion_2D();
	bool TestFaceDeletionAndInsertion_2D();
	bool TestInterFaceDeletionAndInsertion_2D();
	bool TestEraseAllPrimitives_2D();
  
	bool TestEntityNumberingFunction_3D();
	bool TestElementDeletionAndInsertion_3D();
	bool TestFaceDeletionAndInsertion_3D();
	bool TestInterFaceDeletionAndInsertion_3D();
	bool TestEraseAllPrimitives_3D();

  private:
	std::string model2d_name_;
	std::string model3d_name_;

    Model<2U>*  model2d_ = nullptr;
	Model<3U>*  model3d_ = nullptr;
};


} // end csmp

#endif /* CSMP_MESH_MANAGER_TEST_H */
