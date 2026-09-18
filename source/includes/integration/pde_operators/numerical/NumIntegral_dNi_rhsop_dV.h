// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DNI_RHSOP_DV_H
#define NUM_INTEGRAL_DNI_RHSOP_DV_H

#include "MathOperatorRHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

std::string parse( SPATIAL_DERIVATIVE );

template<uint32_t> class Element;

/**
    @brief Directional gradient projection RHS
    
    f_k = ∫_Ω N_k (∂φ/∂xᵢ) dV

    where φ is a scalar field interpolated from nodal values φⱼ via the shape functions:
    
    ∂φ/∂xᵢ = Σⱼ (∂N_j/∂xᵢ) φⱼ

    and xᵢ is the spatial direction selected by the xyz_ member (x, y or z).

    Operand: Scalar — node-placed. Nodal values are collected into op_vec_ via NodePropertyVector
    and used to reconstruct the gradient of the scalar field at each integration point.

    Test variable: Scalar, node-placed.

    Key distinction from NumIntegral_dNT_op_dV: That operator contracts the transposed gradient of the
    test function (∇N)ᵀ with a pre-computed body force vector stored as an element property.
    This operator instead computes the gradient of the operand itself from its nodal values and projects it
    onto the scalar test function N — the roles of the gradient and the interpolation are swapped.

    Application:

    Directional derivative of a nodal scalar field (e.g. pressure gradient contribution to a transport equation)
    - Advective flux terms where the transported quantity is node-placed
    - Gradient recovery operators
    
    Volume integral over the gradient of the Operand in the direction i (i=x,y,z).
    Accumulation into the right-hand side of the linear algebraic system Ax=b
    
    @note created to compute a grad P right-handside for a 2-step Stokes
    lubrication solver.
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNi_rhsop_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNi_rhsop_dV( const PropertyDatabase<dim>&,
                              const char* oper,
                              const char* test );
    
    ~NumIntegral_dNi_rhsop_dV() = default;
  
    /// reads the operand values from the nodes and stores them in a vector
    void GetOperands( const CELL<dim>& ) override final;
  
    void ComputeContribution( const CELL<dim>& ) override final;
  
    /// to chose the spatial derivate direction of interest; default is Y-axis
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    NumIntegral_dNi_rhsop_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNi_rhsop_dV<dim,CELL> (*this); }
  
  private:
    NumIntegral_dNi_rhsop_dV();
    
    SPATIAL_DERIVATIVE           xyz_;     ///< direction of partial derivative of interest
    std::vector<ScalarVariable>  op_vec_;  ///< nodal operand values
    DenseMatrix<DM_MIN>          DERIV_;   ///< shape function derivative matrix
};

} // csmp

#endif
