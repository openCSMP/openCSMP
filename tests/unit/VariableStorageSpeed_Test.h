#ifndef CSMP_VARIABLE_STORAGE_SPEED_TEST_H
#define CSMP_VARIABLE_STORAGE_SPEED_TEST_H

#include "Test.h"

namespace csmp
  {
  
  template<uint32_t> class Model;

  class VariableStorageSpeed_Test : public Test
    {
    public:
      virtual void run();

      void ScalarReadWriteWithSmallDataset();
      void CompareIndexWithINDEX( Model<3U>& );
      void TensorReadWithINDEXvsIndex( Model<3U>& );
      void TestReadingArrayVariableVersusVector( Model<3U>& );
   //    void TestReadingArrayNew_vector_vs_reused_vector( const Model<3U>& );
      
      static const bool verbose = true;
    };
    
  } // csmp

#endif
