#ifndef CSMP_PARDISO_SOLVER_H
#define CSMP_PARDISO_SOLVER_H

#include "SolverSettings.h"
#include "Solver.h"

namespace csmp {

class SparseMatrix;

struct PARDISO_SolverSettings : public SolverSettings {
};


/**
    Interface to the direct PARDISO solver from the MKL library.
    Derived from Solver class, following CSMP's strategy design pattern for the solution of linear algebraic systems of equations.
 
    @note works only if there is access to the MKL library.
 
    @author draft SKM
    @date 30/7/2019
*/
class PARDISO_Solver : public Solver {
  public:
      explicit PARDISO_Solver( PARDISO_SolverSettings& );
      virtual ~PARDISO_Solver();
  
  protected:
      // TODO: this does not override the base class!
      virtual void SolveMatrixEquation( csmp::SparseMatrix& A,
                                        std::vector<double>& b,
                                        std::vector<double>& x,
                                        double solver_tolerance );

  private:
      PARDISO_SolverSettings& settings_;
};

} // end namespace csmp

#endif // CSMP_PARDISO_SOLVER_H
