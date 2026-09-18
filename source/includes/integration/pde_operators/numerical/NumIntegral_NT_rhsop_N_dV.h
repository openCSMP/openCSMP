// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_NUM_INTEGRAL_NT_RHSOP_N_DV_H
#define CSMP_NUM_INTEGRAL_NT_RHSOP_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief Known as: "mass matrix", "fluid sources or sinks", or "capacitance matrix"

f_j = ∫_Ω N_j [σ] N_k φ_k dV

where φ_k are known nodal values of a scalar field.

Operand: Scalar — element or integration point-placed.

Test variable: Scalar, node-placed.

Application: Volumetric source or sink terms (fluid injection/production, heat generation).
This operator collapses the matrix into the right-hand vector.

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_rhsop_N_dV final : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_rhsop_N_dV( const PropertyDatabase<dim>&,
                               const char* oper, const char* test );

    void GetOperands(        const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

    NumIntegral_NT_rhsop_N_dV<dim,CELL>* clone() const override final
      { return new NumIntegral_NT_rhsop_N_dV<dim,CELL>(*this); }

  private:
    std::vector<double> op_;   ///< scalar material operand values —
                               ///< size 1 for element/face/region placement,
                               ///< size n_integration_points otherwise
    DenseMatrix<DM_MIN> NT_;   ///< interpolation function column vector (nodes x 1)
    DenseMatrix<DM_MIN> N_;    ///< interpolation function row vector    (1 x nodes)
    DenseMatrix<DM_MIN> RHS_TEMP_; ///< accumulated mass matrix          (nodes x nodes)
};


} // csmp

#endif
















