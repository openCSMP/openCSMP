#ifndef BOX_TEST_H
#define BOX_TEST_H

#include "Test.h"
#include "Box.h"

namespace csmp{

class Box_Test : public Test
{
  public:
    Box_Test( const char* name="cube_flag" )
      : model_name_(name) {}

    virtual void run();
  
    bool TestBoundaryFlagAssigment2D();
    
    bool TestWhetherSideBoundaryFlagsArePresent();
    
    bool TestWhetherAllBoxFlagsArePresent();
    
    bool TestBoundaryFlagging();
  
    bool TestBoundaryFlagRecreation();
    
    bool TestWhetherBoundaryFlagsArePreservedInBinaryFile();
  
    // new format of binary IO
    bool TestWhetherBoundaryFlagsArePreservedInBinaryFile1();

    bool TestConsistencyOfBoxFlaggingWithBoundaryIdentification();
  
    /// testing the normals of the volumetric elements and the boundary faces
    void TestWhetherSimplexNormalsAreOutwardPointing();

  private:
    std::string model_name_;
};

} // csmp

#endif // BOX_TEST_H
