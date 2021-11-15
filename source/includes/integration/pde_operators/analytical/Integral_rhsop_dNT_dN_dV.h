#ifndef INTEGRAL_RHSOP_DNT_DN_DV_H
#define INTEGRAL_RHSOP_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/// Known as: streaming potential source term
template<size_t dim, class SIMPLEX=Element<dim>, typename var=ScalarVariable >
class Integral_rhsop_dNT_dN_dV : public MathOperatorRHS<dim> {
  public:
    Integral_rhsop_dNT_dN_dV( const PropertyDatabase<dim>& pref, 
                              const char* oper,
                              const char* basis,
                              const char* test,
                              double prefactor=1. );

    virtual void GetOperands( SIMPLEX& e );
    
    virtual void ComputeContribution( SIMPLEX& e );

  private:
    DenseMatrix<DM_MIN> DN, DNT;
    Parameter           basic_;
    std::vector<var>    basic_var_;
    const double      prefactor_;
};


/**
@class Integral_rhsop_dNT_dN_dV Integral_rhsop_dNT_dN_dV "pde_operators/Integral_rhsop_dNT_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999
 
@section motivation Motivation

PDE operator representing the divergence squared of the dependent 
variable.
*/

} // csmp

#endif
