#ifndef ANSYS_MODEL2D_TEST_H
#define ANSYS_MODEL2D_TEST_H

#include "Test.h"

namespace csmp{

// P. Lang 2011
class ANSYS_Model2D_Test : public Test
  {
    public:
      virtual void run();
      
    // same test as in VData_Test, but with a more realistic model created in ANSYS
    void Test_CreateConsistentLineElementOrientations2D();
  };

} // csmp

#endif // ANSYS_MODEL2D_TEST_H
