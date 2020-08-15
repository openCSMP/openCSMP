//
//  Integrator.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 19/7/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_INTEGRATOR_H
#define CSMP_INTEGRATOR_H

#include "IntegralEquation.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"

namespace csmp {

class SAMG_Solver;

/**
    Integrator - to solve the linear algebraic system arising from the upwind implicit form of the transport equations.
    
    Operations 
    + IntegrateOver( equation, domain, time )
    + AssignInitialConditions( domain )
    + AssignEssentialConditions( domain )
    + InsertInMatrix()
    + InsertInVector()
    + Solve( solver );
    + PostProcess( domain );
    + OutputResults( domain );
    Attributes
    + SparseMatrix  A_;
    + Vector             rhs_;
    + Vector             x_;
*/
template<size_t dim, template<size_t> class USER>
class Integrator {
  public:
    /// m=n since solution matrix is square, limits should include permitted tolerances
    Integrator( size_t m_x_n, double64 lower_limit, double64 upper_limit );

    /// uses IntegralEquation and Accumulator to build system and Solver to solve it
    void IntegrateOver( double64 time_increment );
   
    /// transform solver settings from initial use to reuse of settings and previous solution 
    void ReconfigureSolverForRepeatedUse();

  protected:
    void AssignInitialConditions();
    void AssignEssentialConditions();
    
    void SolveLinearAlgebraicSystem();
    
    /// transfer results updating concentration, zeroing out 'new concentration' values, and performing range checks; returns error
    double64 VerifyAndAssignResults( bool show_range, bool do_range_check );

  private:
    /// shorthand for accessing the class that this is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
    
    SAMG_Settings   settings_;
    SAMG_Solver     solver_;
    const double64  upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
};

} // end csmp

#endif /* CSMP_INTEGRATOR_H */


//   void InsertInMatrix( size_t i, siz_t j );
//   void InsertInVector( size_t i );
