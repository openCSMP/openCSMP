/**
@file CompressedRowMatrix_Test.h
@author J. E. Mindel

Tests Compressed Row Matrix storage (CSR) for compatibility
with external solvers.
Methods OffsetNumbering,SetRowToZero,RemoveRowAndColumn,RemoveRowsFromMatrix are
not used and therefore not tested. They probably exist in the TODO list.

*/
#ifndef COMPRESSED_ROW_MATRIX_TEST_H
#define COMPRESSED_ROW_MATRIX_TEST_H

#include "CompressedRowMatrix.h"
#include "Test.h"

namespace csmp
{
  /**
  @addtogroup CSMPUnitTests
  @{
  */
  class CompressedRowMatrix_Test : public Test
    {
       public:
         explicit CompressedRowMatrix_Test( bool verbose = true,
                                            double tolerance = std::numeric_limits<double>::epsilon() * 10. )
           : tolerance_(tolerance),
              verbose_(verbose) {}
              
         virtual void run();
         
         void TestCRM_Multiplication();
         
         void Test_generateSparsityPattern();

    private:
        double     tolerance_;
        const bool verbose_;
  };
  /**
  @}
  */
} // end namespace csmp

#endif // COMPRESSED_ROW_MATRIX_TEST_H
