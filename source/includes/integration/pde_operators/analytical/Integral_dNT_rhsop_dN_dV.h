#ifndef INTEGRAL_DNT_RHSOP_DN_DV_H
#define INTEGRAL_DNT_RHSOP_DN_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// Known as: streaming potential source term: grad . [L grad pf]
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_rhsop_dN_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>&,
                              const char* oper, const char* test, const char* grad_var );

    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

  private:
    DenseMatrix<DM_MIN>  DN, DNT, VAR; 
    csmp::Index              grad_key;
};


/**
@class Integral_dNT_rhsop_dN_dV Integral_dNT_rhsop_dN_dV "pde_operators/Integral_dNT_rhsop_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999
 
@section motivation Motivation

PDE operator representing the divergence squared of the dependent 
variable.
*/

} // csmp

#endif
