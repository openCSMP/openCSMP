#ifndef DYNAMIC_ARRAY_TEST_H
#define DYNAMIC_ARRAY_TEST_H

#include "Test.h"

namespace csmp
{
    class DynamicArray_Test : public Test
      {
        public:
          virtual void run();

          void Test_DynamicArray2D();
          void Test_DynamicArray3D();
   
        private:
          const bool verbose_ = true;
     };
}
#endif // DYNAMIC_ARRAY_3D_TEST_H
