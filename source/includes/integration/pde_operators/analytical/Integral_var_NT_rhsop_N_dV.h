#ifndef CSMP_INTEGRAL_VAR_NT_RHSOP_N_DV_H
#define CSMP_INTEGRAL_VAR_NT_RHSOP_N_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief Integral_var_NT_rhsop_N_dV Integral_var_NT_rhsop_N_dV "pde_operators/Integral_var_NT_rhsop_N_dV.h"
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
alternative version.
 
@section applicability Applicability
 
Only in transient problems
 
@section structure Structure
 
Additional variable "basic" compared to other classes derived from
MathOperatorRHS. Corresponds to the u1 variable mentioned in the example
above. The variable "test" is needed to specify where the entry is
assembled in the global rhs vector. Additional variables "upwind" and
"trigger" for upwinding

@attention Implementation works only for linear triangles

@author A. Burri
@date 2004

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_var_NT_rhsop_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_var_NT_rhsop_N_dV( const PropertyDatabase<dim>&,
                                const char* oper,
                                const char* basic,
                                const char* test,
                                const char* var,
                                const double prefactor = 1.);
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_var_NT_rhsop_N_dV<dim,CELL>* clone() const override final
       { return new Integral_var_NT_rhsop_N_dV<dim,CELL> (*this); }

  private:
    void ComputeIntegral( const CELL<dim>& );
    
    ScalarVariable op_;
    
    Parameter basic_;
    Parameter var_;
    
    std::vector<ScalarVariable > basic_var_;
    std::vector<ScalarVariable > vvar_;
    const double prefactor_;
};


} // csmp

#endif
















