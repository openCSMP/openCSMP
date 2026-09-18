// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LINEAR_SOLVER_EXAMPLE_H
#define LINEAR_SOLVER_EXAMPLE_H

#include "Example.h"

#include "CSMP_definitions.h"
#include "LinearSolver.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

namespace csmp {

class  LinearSolver_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();

  void SetupMatrix( csmp::SparseMatrix& A );
  void SetupVectors( csmp::SparseMatrix& A,
                     std::vector<double>&,
                     std::vector<double>&);

  #if defined(CSMP_WITH_DUNE_ISTL)
  void SetupMatrix( Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >&);
  void SetupVectors( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >&,
                     Dune::BlockVector<Dune::FieldVector<double,1> >&,
                     Dune::BlockVector<Dune::FieldVector<double,1> >&);
  #endif
};

} // csmp

#endif // LINEAR_SOLVER_EXAMPLE_H
