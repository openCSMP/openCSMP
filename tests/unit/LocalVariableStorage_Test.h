#ifndef LOCAL_VARIABLE_STORAGE_TEST_H
#define LOCAL_VARIABLE_STORAGE_TEST_H

#include "Test.h"

namespace csmp
  {

// TODO: add test of the move constructor
  class LocalVariableStorage_Test : public Test
    {
    public:
      virtual void run();
    private:
      template<size_t dim>
      void runTest();

      void run3D();
    };

  } // csmp

#endif
