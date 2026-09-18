// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NumIntegral_dNT_rhsop_dN_dV_h
#define NumIntegral_dNT_rhsop_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    @brief RHS = vector operator - interpolation function derivative matrix collapsed
    into righthand side vector.

    Known as Diffusive source / flux divergence RHS

    f_j = ∫_Ω (∇N_j)ᵀ [σ] ∇φ dV

    where φ is a known scalar field (e.g. a reference pressure or gravity head).

    Operand: Scalar or tensor — element-placed.
    Test variable: Scalar, node-placed.
    Application: Gravity-driven flow (hydrostatic body force in pressure formulation), reference state subtraction.
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_rhsop_dN_dV final : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& p, 
                                 const char* oper,          
                                 const char* test );       

    NumIntegral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& p, 
                                 const char* integral_multiplier,
                                 const char* oper,          
                                 const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
  
    NumIntegral_dNT_rhsop_dN_dV<dim,CELL>* clone() const override final
      { return new NumIntegral_dNT_rhsop_dN_dV<dim,CELL> (*this); }
      
  private:
    DenseMatrix<DM_MIN>           DN, DNT, OPMAT, TEMP;
    std::vector<ScalarVariable >  noperand;
    ScalarVariable                eoperand, multiplier;
    bool                          has_multiplier_;
    csmp::Index                   mult_key;
};


} // csmp

#endif
















