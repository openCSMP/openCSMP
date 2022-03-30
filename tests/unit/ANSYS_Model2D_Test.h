#ifndef ANSYS_MODEL2D_TEST_H
#define ANSYS_MODEL2D_TEST_H

#include "Test.h"

namespace csmp{

// P. Lang 2011
class ANSYS_Model2D_Test : public Test
  {
    const bool verbose_ = true;
    
    public:
      virtual void run();
      
      // using line element VSet
      void Test_printLineElementRegion();
       
      // same test as in VData_Test, but with a more realistic model created in ANSYS
      void Test_CreateConsistentLineElementOrientations2D();
      
      // tests the creation of an internal boundary from the line-element regions
      void Test_CreatInternalBoundary();
  };

} // csmp

#endif // ANSYS_MODEL2D_TEST_H
