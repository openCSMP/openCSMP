// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DNT_DN_DV_H
#define NUM_INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// Laplacian squared PDE operator:  div^2 N = interpolation function derivate matrix squared. No material operand.
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_dN_dV final : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_dN_dV( const PropertyDatabase<dim>&, 
                           const char* basic, 
                           const char* test );
  
    /// no operands need to be fetched from computational domain
    void GetOperands( const CELL<dim>& ) override final {}  
    void ComputeContribution( const CELL<dim>& ) override final;
  
    NumIntegral_dNT_dN_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNT_dN_dV<dim,CELL >(*this); }
    
  private:
    DenseMatrix<DM_MIN>  B, BT; 
};

} // csmp

#endif
