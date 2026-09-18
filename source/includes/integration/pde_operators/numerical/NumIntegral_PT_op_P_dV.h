// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_PT_OP_P_DV_H
#define NUM_INTEGRAL_PT_OP_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**
@brief Vector capacitance / storage capacity RHS

f_k = ∫_Ω Pₖᵀ [σ] Pⱼ dV · 1  (consistent form)

where P is the vector interpolation matrix, [σ] is the scalar or diagonal tensor operand, and the result is condensed to a vector by contraction with a unity vector — equivalent to row-summing the consistent mass matrix:

f_k = [σ] · V / n_nodes (lumped form)

where V is the physical element volume and n_nodes is the number of element nodes. The operand-weighted volume is distributed equally among all nodes.

Operand: Scalar or diagonal tensor — element-placed (lumped and consistent) or node-placed (consistent only). Only the diagonal entries MTRL(i,i) are used, making this operator suitable for isotropic or diagonally anisotropic material properties.

Test variable: Vector (displacement or velocity), node-placed.

Key distinction from NumIntegral_NT_lhsop_N_dV: That operator assembles the scalar mass matrix (LHS) for scalar DOFs. This operator assembles the right-hand side vector for vector DOFs by projecting a known property field onto the nodal basis — it does not involve an unknown solution variable.

Key distinction from NumIntegral_op_NT_N_dV: That operator operates on scalar DOFs using the scalar interpolation functions N and Nᵀ. This operator operates on vector DOFs using the vector interpolation matrix P and Pᵀ, making it appropriate for mechanics problems where the solution variable is a displacement or velocity vector.

Key distinction from NumIntegral_PT_op_dV: That operator computes ∫ Pᵀ [f] dV where [f] is a known body force vector — a single contraction of P with the operand. This operator computes ∫ Pᵀ [σ] P dV condensed to a vector — a double contraction of P with the operand, equivalent to the diagonal of the consistent mass matrix projected onto the RHS.

Physical meaning: Distributes an element-level scalar or tensor property (e.g. density, compressibility, storage coefficient) to the nodal degrees of freedom weighted by the element volume. When the operand is unity, the result is the nodal volume fractions summing to the total element volume.

Application:

- Capacitance matrix assembly for vector field problems
- Storage capacity distribution in poromechanics with vector primary variable
- Inertia term assembly in dynamic mechanics (consistent or lumped mass)
- Distributing element-averaged density to nodal DOFs

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_PT_op_P_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_PT_op_P_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* test );

    void ComputeContribution( const CELL<dim>& ) override final;

  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















