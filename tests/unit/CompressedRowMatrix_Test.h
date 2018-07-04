/**
@file CompressedRowMatrix_Test.h
@author J. E. Mindel

Tests Compressed Row Matrix storage (CSR) for compatibility
with external solvers.
Methods OffsetNumbering,SetRowToZero,RemoveRowAndColumn,RemoveRowsFromMatrix are
not used and therefore not tested. They probably exist in the TODO list.

*/
#ifndef COMPRESSEDROWMATRIX_TEST_H
#define COMPRESSEDROWMATRIX_TEST_H

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
         explicit CompressedRowMatrix_Test( bool verbose = false,
                                            double64 tolerance = std::numeric_limits<double64>::epsilon() * 10. )
           : tolerance_(tolerance),
              verbose_(verbose) {}
              
         virtual void run();

    private:
        const double64 tolerance_;
        const bool verbose_;
  };
  /**
  @}
  */
} // end namespace csmp

#endif // COMPRESSEDROWMATRIX_TEST_H
