#ifndef Upwind_Integral_dNT_rhsop_dN_dV_h
#define Upwind_Integral_dNT_rhsop_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<size_t dim,class SIMPLEX=Element<dim> >
class Upwind_Integral_dNT_rhsop_dN_dV : public MathOperatorRHS<dim> {
  public:
    Upwind_Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& p,
                                     const char* oper, 
                                     const char* basic,
                                     const char* test,
                                     const char* upwind,
                                     const char* trigger,
                                     const double prefactor = 1. );
    
    virtual void GetOperands(  SIMPLEX& e );
    virtual void ComputeContribution( SIMPLEX& e );
  
  private:
    DenseMatrix<DM_MIN>  DN, DNT;
    Parameter basic_;
    Parameter upwind_;
    Parameter trigger_;
    std::vector<ScalarVariable > basic_var_;
    std::vector<ScalarVariable > upwind_var_;
    std::vector<ScalarVariable > trigger_var_;
    const double prefactor_;
};



/**
 
@class Upwind_Integral_dNT_rhsop_dN_dV Upwind_Integral_dNT_rhsop_dN_dV "pde_operators/Upwind_Integral_dNT_rhsop_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999

PDE-operator for off-diagonal terms on the right hand side of a matrix
equation.
 
@section motivation Motivation
 
In transient coupled problems, the following matrix formulation is seen
often:

	1/dt(M u^n+1 - M u^n) = A u^n

The term

	1/dt M u^n

can be collapsed into the right hand side vector, since u^n is known.
In the case of a coupled problem, where u = [u1, u2], the matrix M can
contain off-diagonal blocks which are multiplied with u1 (u2), but the
contribution must be assembled into the u2 (u1) part of the global rhs
vector! Since Integral_NT_op_N_dV cannot deal with that, I wrote this
animal here...
 

@section applicability Applicability
 
Only in transient problems 
 
 
@section structure Structure
 
Additional variable "basic" compared to other classes derived from
MathOperatorRHS. Corresponds to the u1 variable mentioned in the example
above. The variable "test" is needed to specify where the entry is 
assembled in the global rhs vector. Additional variables "upwind" and
"trigger" for upwinding
 */


} // csp

#endif
















