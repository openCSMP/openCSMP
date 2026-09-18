// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INTEGRAL_DNT_RHSOP_DN_V_DV_H
#define INTEGRAL_DNT_RHSOP_DN_V_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// Known as: streaming potential source term: grad . [L grad pf] -> RHS[j]=∫ Ω e (∇N j ) T [σ]DN⋅{u}dV
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_rhsop_dN_v_dV : public MathOperatorRHS<dim,CELL> {
  public:
    /// the gradient variable is an extra scalar placed on the nodes that gets multiplied with
    Integral_dNT_rhsop_dN_v_dV( const PropertyDatabase<dim>&,
                                const char* oper, const char* test, const char* grad_var );

    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

     virtual Integral_dNT_rhsop_dN_v_dV<dim,CELL>* clone() const override final
        { return new Integral_dNT_rhsop_dN_v_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  DN, DNT;
    DenseMatrix<dim>     VAR;
    csmp::Index          grad_key;
};


/**
@class Integral_dNT_rhsop_dN_v_dV Integral_dNT_rhsop_dN_v_dV "pde_operators/Integral_dNT_rhsop_dN_v_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999
 
@section motivation Motivation

PDE operator representing the divergence squared of the dependent 
variable.
*/

} // csmp

#endif
