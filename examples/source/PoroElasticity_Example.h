/*
 *  PoroElasticity_Example.h
 *  opencsmp
 *
 *  Demonstrates coupled poromechanics (Biot consolidation) using the
 *  block-structured PDE_Integrator framework.
 */

#ifndef CSMP_PORO_ELASTICITY_EXAMPLE_H
#define CSMP_PORO_ELASTICITY_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Biot poromechanics example — coupled displacement and fluid pressure.

    Demonstrates the assembly and solution of a block-structured coupled
    finite element system for quasi-static Biot consolidation in 2D plane
    strain, using the openCSMP PDE_Integrator framework.

    @par Governing equations

    The coupled Biot system in weak form leads to the block matrix problem:

        | K    Q  | { u }   { f_u }
        |         |       =
        | Qᵀ   S  | { p }   { f_p }

    where:
    - K  is the mechanical stiffness matrix (displacement-displacement)
    - Q  is the Biot coupling matrix (displacement-pressure)
    - Qᵀ is the transpose coupling matrix (pressure-displacement)
    - S  is the storage / compressibility matrix (pressure-pressure)
    - u  is the nodal displacement vector
    - p  is the nodal fluid pressure vector
    - f_u is the mechanical force vector (body forces, tractions)
    - f_p is the fluid source vector

    @par CSMP operators used

    | Operator | Block | Integral |
    |----------|-------|---------|
    | `NumIntegral_BT_D_B_dV` | K | ∫ Bᵀ D(E,ν) B dV |
    | `NumIntegral_NT_lhsop_N_dV` | S | ∫ Nᵀ (1/M) N dV |
    | `NumJacobianIntegral_BT_m_N_dV` | Q | -α ∫ Bᵀ m N dV |
    | `NumJacobianIntegral_N_mT_B_dV` | Qᵀ | -α ∫ Nᵀ mᵀ B dV |
    | `NumIntegral_PT_op_dV` | f_u | ∫ Pᵀ f dV |
    | `NumIntegral_NT_rhsop_N_dV` | f_p | ∫ Nᵀ q N dV |

    @par Model setup

    The example uses a small 2D mesh (TINY) with 6 nodes, 2 triangles,
    1 quadrilateral and 1 line element. Material properties:

    - Young's modulus E = 10000 Pa
    - Poisson's ratio ν = 0.25 (plane strain)
    - Biot coefficient α = 1.0
    - Biot modulus 1/M = 1.0e-5 Pa⁻¹

    Boundary conditions:
    - Left boundary:   u = 0 (fixed), p = 0 (drained)
    - Right boundary:  u_x = 0.01 m (prescribed), p = 500 Pa
    - Bottom boundary: u_y = 0 (roller)

    @par Dirichlet elimination

    Essential boundary conditions are applied by full elimination of
    constrained DOFs from the global system. The reduced system has
    5 free DOFs:

        [0] u_x(1)   [1] u_y(1)   [2] u_y(2)   [3] u_x(4)   [4] p(4)

    @par Output

    Results are written to VTU format via VTU_Interface. Output variables:
    - displacement vector [m]
    - fluid pressure [Pa]
    - stress tensor [Pa]
    - strain tensor [-]
    - mean stress [Pa]

    @par References

    - Biot, M.A. (1941). General theory of three-dimensional consolidation.
      Journal of Applied Physics, 12(2), 155-164.
    - Zienkiewicz & Taylor, The Finite Element Method, Vol. 1, Ch. 11.
    - Smith & Griffiths, Programming the Finite Element Method, Ch. 9.
*/
class PoroElasticity_Example : public Example {
  public:
    void Specifications();
    void Run();
};

} // namespace csmp

#endif /* CSMP_PORO_ELASTICITY_EXAMPLE_H */

