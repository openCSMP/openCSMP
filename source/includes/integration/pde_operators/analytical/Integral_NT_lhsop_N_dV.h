// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INTEGRAL_NT_LHSOP_N_DV_H
#define INTEGRAL_NT_LHSOP_N_DV_H

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
  @author S.K. Matthaei
  @author S. Roberts
  @date 1999 */

/// known as mass or capacitance matrix
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_NT_lhsop_N_dV : public MathOperatorLHS<dim,CELL> {
  public:
    Integral_NT_lhsop_N_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* basic, const char* test );

    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_NT_lhsop_N_dV<dim,CELL>* clone() const override final
      { return new Integral_NT_lhsop_N_dV<dim,CELL> (*this); }

  private:
    double scalar_value_;
};



/**

copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















