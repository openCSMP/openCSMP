// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_INTEGRAL_SET_RHS_TO_ZERO_H
#define CSMP_INTEGRAL_SET_RHS_TO_ZERO_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

  /**
  @author S.K. Matthaei
  @author S. Roberts
  @date 1999 */

template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_SetRHS_to_Zero : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_SetRHS_to_Zero( const PropertyDatabase<dim>&, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final { /* no data need to be collected */ }
    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_SetRHS_to_Zero<dim,CELL>* clone() const override final { return new Integral_SetRHS_to_Zero<dim,CELL> (*this); }
};

} // csmp

#endif // CSMP_INTEGRAL_SET_RHS_TO_ZERO_H
















