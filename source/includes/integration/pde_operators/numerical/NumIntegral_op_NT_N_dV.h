// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_OP_NT_N_DV_H
#define NUM_INTEGRAL_OP_NT_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**
    @brief Lumped property projection RHS (not implemented here)

    f_k = ∫_Ω Nₖ [σ] Nⱼ φⱼ dV

    but in its lumped form (as implemented)

    f_k = [σ] · V / n_nodes

    where V is the physical element volume and n_nodes is the number of element nodes.
    The element contribution is distributed equally among all nodes,
    equivalent to row-sum lumping of the consistent integral.

    For nodal operands the quadrature weights are included:

    f_k = Σᵢ wᵢ [σ(ξᵢ)] · V / n_nodes

    Operand: Scalar, vector or tensor — element-placed (uniform over the element) or node/integration-point-placed (spatially varying). The diagonal entries of MTRL are used for vector operands; all entries for tensor operands.

    Test variable: Scalar or vector, node-placed.

    Lumped formulation: Always used in practice — the consistent formulation (LumpedFormulation(false)) is not yet implemented and will print a warning.

    Key distinction from NumIntegral_NT_lhsop_N_dV: That operator assembles the consistent or lumped mass matrix (LHS). This operator assembles the right-hand side vector by projecting a known property field onto the nodal basis in lumped form
    — it does not involve an unknown solution variable.

    Key distinction from NumIntegral_NT_rhsop_N_dV: That operator multiplies the operand by known nodal values φⱼ of a scalar field. This operator distributes the operand itself directly to the nodes without multiplying by a separate scalar field.

    Application:

    - Distributing element-averaged material properties to nodes for post-processing or initialisation
    - Lumped projection of integration-point quantities onto the nodal basis
    - Assembling nodal equivalent forces from element-level stress or strain fields in linear elasticity

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_op_NT_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_op_NT_N_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* test );

    void GetOperands( const CELL<dim>& ) override final { /* do not load any data */ }
    void ComputeContribution( const CELL<dim>& ) override final;
    
  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger
} // csmp

#endif
















