// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_EIGEN_SOLVER_H
#define CSMP_EIGEN_SOLVER_H

#include "Solver.h"

/**
  @brief a csmp-eigen interface solver class
  @author Luat Khoa Tran
 */
namespace csmp {

class SparseMatrix;
class CompressedRowMatrix;

class EigenSolver : public Solver
  {
  public:
    EigenSolver() = default;
    virtual ~EigenSolver() = default;

  std::string Name() const override { return "EigenSolver"; }

  protected:

   /**
     @brief solves matrix equation A x = b using direct Eigen Matrix solver.
     
     @param A csmp sparse matrix
     @param b RHS vector
     @param x vector of unknowns
     @param no_unknowns degrees of freedom of the problem
   */
    virtual void SolveMatrixEquation(SparseMatrix& A,
                                     std::vector<double>& b,
                                     std::vector<double>& x,
                                     size_t no_unknowns) override final;

    virtual void SolveMatrixEquation(CompressedRowMatrix& A,
                                     std::vector<double>& b,
                                     std::vector<double>& x,
                                     size_t no_unknowns) override final;
  };

} // end namespace csmp

#endif
