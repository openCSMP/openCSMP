#ifndef INTEGRAL_DNT_DN_DV_H
#define INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_dNT_dN_dV : public MathOperatorRHS<dim> {
  public:
    Integral_dNT_dN_dV( const PropertyDatabase<dim>& pref, 
                              const char* oper, const char* test );

    virtual void ComputeContribution( SIMPLEX& e );

  private:
    DenseMatrix<DM_MIN>  DN, DNT, UNITY; 
};


/**
 
@class Integral_dNT_dN_dV Integral_dNT_dN_dV "pde_operators/Integral_dNT_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999

@section motivation Motivation
 
PDE operator representing the divergence squared of the dependent 
variable.
*/

} // csmp

#endif
