#ifndef INTEGRAL_DNT_RHSOP_DN_DV_H
#define INTEGRAL_DNT_RHSOP_DN_DV_H

#include "MathOperatorRHS.h"
#include "ScalarVariable.h"

namespace csmp {

template<uint32_t> class Element;

/// integral over product of interpolation function derivatives
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_rhsop_dN_dV : public MathOperatorRHS<dim,CELL> {
  public:
    /// the gradient variable is an extra scalar placed on the nodes that gets multiplied with
    Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>&,
                              const char* oper, const char* test );

    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

     virtual Integral_dNT_rhsop_dN_dV<dim,CELL>* clone() const override final
        { return new Integral_dNT_rhsop_dN_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>          DN_, DNT_, TEMP_;
    std::vector<ScalarVariable>  u_; ///< test function variable
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
