//
//  MeshManagementUtilities_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 16/1/2026.
//

#ifndef CSMP_MESH_MANAGEMENT_UTILITIES_TEST_H
#define CSMP_MESH_MANAGEMENT_UTILITIES_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

/**
   Testing checklist

    @author SKM
    @date 16/1/2026
*/
class MeshManagementUtilities_Test : public Test {
  public:
    MeshManagementUtilities_Test() {}
    
    virtual void run();
  
  private:
  
    void Test_transformMeshIntoRefinedLinearAndQuadraticMeshes();
  
    static const bool verbose_ = true;
  
};

} // end csmp

#endif
