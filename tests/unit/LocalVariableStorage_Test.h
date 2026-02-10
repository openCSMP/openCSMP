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
      template<uint32_t dim>
      void runTest();

      void run3D();
      
      /// read/store for INDEX class eliminating dynamic dispatching on type and placement (2020)
      void run3D_with_templatized_INDEX();
      
      /// tests function that distinguishes different CSMP_variables
      void Test_variableType();
      
      /// ReadVector, ReadTensor, ReadArray directly returning csmp variables (skm 12/24)
      void runExtensionsForDeclarativeProgramming();
      void runExtensionsForDeclarativeProgramming_INDEX();
    };

  } // csmp

#endif
