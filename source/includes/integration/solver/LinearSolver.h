#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

/*
=======================
Linear Solvers Interface
=======================
*/
// define CSMP_DEFAULT_LINEAR_SOLVER by the priority order
// TODO: relieve this from the pre-processor to get better control on what happens

#define CSMP_DEFAULT_LINEAR_SOLVER csmp::Solver
#define CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS csmp::SolverSettings
#include "Solver.h"

#ifdef CSMP_WITH_MESCHACH
#undef CSMP_DEFAULT_LINEAR_SOLVER
#undef CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS
#define CSMP_DEFAULT_LINEAR_SOLVER csmp::Gauss_Solver
#define CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS csmp::Gauss_SolverSettings
#include "Gauss_Solver.h"
#endif

#ifdef CSMP_WITH_DUNE_ISTL
#undef CSMP_DEFAULT_LINEAR_SOLVER
#undef CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS
#define DUNE_VECTOR_BLOCK_FORMAT Dune::FieldVector<double,1>
#define DUNE_MATRIX_BLOCK_FORMAT Dune::FieldMatrix<double,1,1>
#define DUNE_VECTOR_FORMAT Dune::BlockVector< DUNE_VECTOR_BLOCK_FORMAT >
#define DUNE_MATRIX_FORMAT Dune::BCRSMatrix< DUNE_MATRIX_BLOCK_FORMAT >
#define DUNE_SOLVER_INTERFACE csmp::DuneISTL_SEQ_BCGS_AMG_SOR< 3, DUNE_MATRIX_FORMAT, DUNE_VECTOR_FORMAT >
#define CSMP_DEFAULT_LINEAR_SOLVER csmp::DuneISTL_Solver< DUNE_MATRIX_FORMAT, DUNE_VECTOR_FORMAT, DUNE_VECTOR_FORMAT, DUNE_SOLVER_INTERFACE >
#define CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS csmp::DuneISTL_Settings
#include "DuneISTL_Solver.h"
#endif

#ifdef CSMP_WITH_SAMG_SOLVER
#undef CSMP_DEFAULT_LINEAR_SOLVER
#undef CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS
#define CSMP_DEFAULT_LINEAR_SOLVER csmp::SAMG_Solver
#define CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS csmp::SAMG_Settings
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif

#endif
