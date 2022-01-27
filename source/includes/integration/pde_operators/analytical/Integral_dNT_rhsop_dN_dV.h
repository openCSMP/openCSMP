#ifndef INTEGRAL_DNT_RHSOP_DN_DV_H
#define INTEGRAL_DNT_RHSOP_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/// Known as: streaming potential source term: grad . [L grad pf]
template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_dNT_rhsop_dN_dV : public MathOperatorRHS<dim> {
  public:
    Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& pref, 
                               const char* oper, const char* test, const char* grad_var );

    virtual void GetOperands( const SIMPLEX& e );
    virtual void ComputeContribution( const SIMPLEX& e );

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
