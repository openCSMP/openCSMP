#ifndef LUDCMP_SOLVER_H
#define LUDCMP_SOLVER_H


#include "Solver.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

namespace csmp {

/** Trivial LU decomposition Solver

This solver was implemented during the 2010 CSMP retreat by MUL
in Starigrad, Croatia, and used to solve FEM matrices with up to 300 dof.
As in Numerical Recipes C++, Cambridge 3rd Edition. Adapted to comply
with csmp::SparseMatrix, STL and according functionality

@attention This is a direct solver and therefore works only for small systems.

@author P. Lang
@date 2010
*/
  class LUdcmp_Solver : public Solver {
    public:
      LUdcmp_Solver();
      virtual ~LUdcmp_Solver();

    protected:
      virtual void SolveMatrixEquation( SparseMatrix& A,
                                        std::vector<double>& b,
                                        std::vector<double>& x,
                                        size_t no_unknowns );
                                        
      virtual void SolveMatrixEquation( CompressedRowMatrix& A,
                                        std::vector<double>& b,
                                        std::vector<double>& x,
                                        size_t no_unknowns );
    private:
      void ludcmp( SparseMatrix& a,
                   long n,
                   std::vector<size_t>& indx);
                   
      void lubksb( SparseMatrix& a,
                   long n,
                   std::vector<size_t>& indx,
                   std::vector<double>& x,
                   std::vector<double>& b );
      void luout();

      double             tiny_; /**< represents chosen numerical limit */
      long                 n_;
      std::vector<size_t>  index_; /**< cache vector for indices */
};

} // end namespace csmp



#endif // LUDCMP_SOLVER_H
