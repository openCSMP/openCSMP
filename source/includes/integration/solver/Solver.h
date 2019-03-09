#ifndef CSMP_SOLVER_H
#define CSMP_SOLVER_H

#include "CSMP_number_types.h"
#include <cstdlib>
#include <vector>
#include "SolverSettings.h"

namespace csmp {

class SparseMatrix;
class CompressedSparseRowMatrix;

/**

@brief Inverts global solution matrices matrices of the finite-element
computations using a variety of solution procedures sparse matrices (class SparseMatrix).
Solver base class in strategy pattern.

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999

@section design Design Intent

To encapsulate numerical algebraic methods in an object which services
the PDE_Integrator and other class objects.

@section applicability Applicability

Currently the solver offers an algebraic multigrid method (AMG), a
conjugate gradient method (CG), and a Gaussian elimination method (GE)
as solution schemes. Gaussian elimination can be applied to non-symmetric
matrices as they arise if advection-dispersion matrices are accumulated
into the global solution matrix. The algebraic multigrid and conjugate
gradient methods require symmetric positive-definite (no negative
elements in the diagonal) solution matrices in order to be applicable.

Dependent on the size of the computation different solution methods are
appropriate: If you have a very small test problem, like a 30 x 30
matrix, CG ist as fast or faster than AMG and it will require less
storage. Above a size of 50 x 50, AMG becomes much faster than the
other methods.


@section participants Participants

Each inversion involves a sparse global solution matrix of type SpMat,
a righthand, and a solution vector of type vector<double64>.


@section implementation Implementation

The Solver uses the Meschach++ and the AMG libraries to carry out the
matrix operations. The user sets the private variable solution method,
whose default is algebraic multigrid.


@section examples Application Examples

The Solver is called by the Algorithm and Algorithm objects, when
they are submitted to the Model through its Apply() interfaces.
In this process the Algorithm hands the solver the global solution
matrix 'G', the solution vector 'x' and the righthand vector 'rh'.:

@code
SolveMatrixEquation( G, rh, x );
@endcode

The Solver then returns the solution into 'x' and it is tested further by
the Algorithm before it is mapped back to the Model variable
storage.

*/
class Solver {
public:
  Solver();
  Solver( SolverSettings* settings );
  virtual ~Solver();

  void  Verbose( bool verb );
  bool  Verbose() const;
  
  void  Solve( SparseMatrix& G,
               std::vector<double64>& rh,
               std::vector<double64>& x,
               size_t no_unknowns = 1U );

  void  Solve( CompressedSparseRowMatrix& G,
               std::vector<double64>& rh,
               std::vector<double64>& x,
               size_t no_unknowns = 1U );

  void  Out( const SparseMatrix& mat,
             const char* fname = "SparseMatrix" ) const;

  void  Out( const std::vector<double64>& vec,
             const char* fname = "CSP_Vec" ) const;

  double64  CalculateResidual( const SparseMatrix& A,
                               const std::vector<double64>& b,
                               const std::vector<double64>& x ) const;

  virtual SolverSettings* GetSolverSettings();

protected:
  virtual void  SolveMatrixEquation( SparseMatrix& A,
                                     std::vector<double64>& b,
                                     std::vector<double64>& x,
                                     size_t no_unknowns ) = 0;

  virtual void  SolveMatrixEquation( CompressedSparseRowMatrix& A,
                                     std::vector<double64>& b,
                                     std::vector<double64>& x,
                                     size_t no_unknowns ) = 0;

  SolverSettings* solver_settings_;

private:
  bool  verbose_;
};

} // csmp
#endif
