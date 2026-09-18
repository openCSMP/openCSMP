// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GAUSS_JORDAN_SOLVER_H
#define GAUSS_JORDAN_SOLVER_H

#include "Solver.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

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
  
  std::string Name() const override { return "GaussJordan_Solver"; }

protected:
  virtual void SolveMatrixEquation( SparseMatrix& A,
                                    std::vector<double>& b,
                                    std::vector<double>& x,
                                    size_t no_unknowns ) override final;

  virtual void SolveMatrixEquation( CompressedRowMatrix& A,
                                    std::vector<double>& b,
                                    std::vector<double>& x,
                                    size_t no_unknowns ) override final;

private:
  void GaussJordan( SparseMatrix& A, std::vector<double>& b );

  void SwapSparseMatrixElements( SparseMatrix& M,
                                 long row1,
                                 long col1,
                                 long row2,
                                 long col2 ) const;
};

} // end namespace csmp


#endif // GAUSSJORDAN_SOLVER_H
