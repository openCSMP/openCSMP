#ifndef LOCAL_VARIABLE_STORAGE_TEST_H
#define LOCAL_VARIABLE_STORAGE_TEST_H

#include "Test.h"

namespace csmp
  {

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