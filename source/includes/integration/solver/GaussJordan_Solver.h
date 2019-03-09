#ifndef GAUSS_JORDAN_SOLVER_H
#define GAUSS_JORDAN_SOLVER_H

#include "Solver.h"
#include "SparseMatrix.h"
#include "CompressedSparseRowMatrix.h"

namespace csmp {

/** Trivial Gauss Solver

Basic Gauss-Jordan solver adapted to comply
with csmp::SparseMatrix, STL and according functionality

@author P. Lang
@date 2010

*/
class GaussJordan_Solver : public Solver {
public:
  GaussJordan_Solver();
  virtual ~GaussJordan_Solver();

protected:
  virtual void SolveMatrixEquation( SparseMatrix& A,
                                    std::vector<double64>& b,
                                    std::vector<double64>& x,
                                    size_t no_unknowns );

  virtual void SolveMatrixEquation( CompressedSparseRowMatrix& A,
                                    std::vector<double64>& b,
                                    std::vector<double64>& x,
                                    size_t no_unknowns );

private:
  void GaussJordan( SparseMatrix& A, std::vector<double64>& b );

  void SwapSparseMatrixElements( SparseMatrix& M,
                                 long row1,
                                 long col1,
                                 long row2,
                                 long col2 ) const;
};

} // end namespace csmp


#endif // GAUSSJORDAN_SOLVER_H
