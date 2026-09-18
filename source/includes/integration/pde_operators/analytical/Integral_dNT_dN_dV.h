// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INTEGRAL_DNT_DN_DV_H
#define INTEGRAL_DNT_DN_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    Interpolation function derivatives squared.
    No material operand required.
    
    @author S.K. Matthai
    @author S. Roberts
    @date 1999

    RHS PDE operator representing the divergence squared of the dependent
    variable.

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_dN_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_dNT_dN_dV( const PropertyDatabase<dim>&, const char* test );

    void ComputeContribution( const CELL<dim>& )  override final;
  
    virtual Integral_dNT_dN_dV<dim,CELL>* clone() const override final
       { return new Integral_dNT_dN_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  DN, DNT, UNITY; 
};


/**
 
@class Integral_dNT_dN_dV Integral_dNT_dN_dV "pde_operators/Integral_dNT_dN_dV.h"
*/

} // csmp

#endif
