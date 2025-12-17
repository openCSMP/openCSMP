#ifndef LU_SOLVER_H
#define LU_SOLVER_H

#include "DenseMatrix.h"

namespace csmp {

/** Trivial LU decompositiono Solver

as in Numerical Recipes, Cambridge 3rd Edition. Adapted to comply
with csmp::SparseMatrix, STL and according functionality

@author Philipp Lang
@modified by Shaho for DenseMatrix
@date 2011
*/
class LU_Solver
{
    public:
      LU_Solver();
      ~LU_Solver();
      void Solve(DenseMatrix<DM_MIN>& A,
                 std::vector<double>& b,
                 std::vector<double>& x,
                 uint32_t no_unknowns);
    private:
      void ludcmp(DenseMatrix<DM_MIN>& a,
                  uint32_t n,
                  std::vector<uint32_t>& indx);
      void lubksb(DenseMatrix<DM_MIN>& a,
                  uint32_t n,
                  std::vector<uint32_t>& indx,
                  std::vector<double>& x,
                  std::vector<double>& b );

      double                 tiny_; /**< represents chosen numerical limit */
      std::vector<uint32_t>  index_; /**< cache vector for indices */
};

} // end namespace csmp



#endif // LU_SOLVER_H
