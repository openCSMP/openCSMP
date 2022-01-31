#ifndef INTEGRAL_NT_RHSOP_N_DV_H
#define INTEGRAL_NT_RHSOP_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_NT_rhsop_N_dV : public MathOperatorRHS<dim> {
  public:
    Integral_NT_rhsop_N_dV( const PropertyDatabase<dim>& p,
                            const char* oper,
                            const char* basic,
                            const char* test );
    
    virtual void GetOperands( const SIMPLEX& );
    virtual void ComputeContribution( const SIMPLEX& );
  
  private:
    DenseMatrix<DM_MIN>           INN;
    Parameter                     basic_;
    std::vector<ScalarVariable >  basic_var_;
    ScalarVariable                sc;
};



/**
 
@class Integral_NT_rhsop_N_dV Integral_NT_rhsop_N_dV "pde_operators/Integral_NT_rhsop_N_dV.h"
 
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
 
Only in transient problems (actually the problem with Integral_NT_op_N_dV
only arises there...), using the late accumulate option. The late
accumulate is switched on automatically and you'll be damned if you
switch it off again.
 
 
@section structure Structure
 
Additional variable "basic" compared to other classes derived from
MathOperatorRHS. Corresponds to the u1 variable mentioned in the example
above. The variable "test" is needed to specify where the entry is 
assembled in the global rhs vector.
*/

} // csmp

#endif
















