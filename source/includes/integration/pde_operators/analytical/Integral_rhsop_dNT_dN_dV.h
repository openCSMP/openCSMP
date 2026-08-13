#ifndef INTEGRAL_RHSOP_DNT_DN_DV_H
#define INTEGRAL_RHSOP_DNT_DN_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// Known as: streaming potential source term
template<uint32_t dim, template<uint32_t> class CELL=Element, typename var=ScalarVariable>
class Integral_rhsop_dNT_dN_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_rhsop_dNT_dN_dV( const PropertyDatabase<dim>&,
                              const char* oper,
                              const char* basis,
                              const char* test,
                              double prefactor=1. );

    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

    Integral_rhsop_dNT_dN_dV<dim,CELL>* clone() const override final
      { return new Integral_rhsop_dNT_dN_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN> DN, DNT;
    Parameter           basic_;
    std::vector<var>    basic_var_;
    const double        prefactor_;
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
