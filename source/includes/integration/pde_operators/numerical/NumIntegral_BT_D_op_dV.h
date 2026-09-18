// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NumIntegral_BT_D_op_dV_h
#define NumIntegral_BT_D_op_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    Volume strain operator.

    f_j = ∫_Ω Bⱼᵀ D(E, ν) [ε₀] dV

    where D is the isotropic constitutive matrix assembled from Young's modulus E and Poisson's ratio ν,
    B is the strain-displacement matrix, and [ε₀] is a known volumetric strain field
    (e.g. thermal expansion, swelling, or dilatation induced by fluid pressure changes).

    Operands: Young's modulus,
    Poisson's ratio (scalar, element-placed) and the
    volumetric strain operand [ε₀] (scalar or tensor, element or integration point-placed).

    Test variable: Vector (displacement), node-placed.

    Application:

    - Thermal loading (thermoelastic coupling)
    - Swelling or shrinkage strains in geomechanics
    - Dilatation-induced body forces in poromechanics
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_BT_D_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_BT_D_op_dV( const PropertyDatabase<dim>&,
                            const char* oper,     ///< volume strain
                            const char* youngs,   ///< Young's modulus
                            const char* poissons, ///< Poisson's ratio
                            const char* test, bool plane_strain=true ); ///< displacement

    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
    
    /// switch from default plane strain to plane stress
    void PlaneStress() noexcept { plane_strain_ = false; }
    
    NumIntegral_BT_D_op_dV<dim,CELL>* clone() const override final { return new NumIntegral_BT_D_op_dV<dim,CELL> (*this); }
    
  private:  
    csmp::Index         nu_key_; ///< Poisson's ratio
    csmp::Index         E_key_;  ///< Young's modulus
    std::vector<double> E_;      ///< Young's modulus
    std::vector<double> nu_;     ///< variable in which Poisson's ratio will be stored
    bool                plane_strain_;
    DenseMatrix<DM_MIN> B_;
    static constexpr auto MATDIM = (dim == 3) ? DM6 : DM3;
    std::vector<std::array<double,MATDIM>> initial_strains_;
    DenseMatrix<MATDIM> D_;
    ScalarVariable      sc_;
    VectorVariable<dim> vc_;
    TensorVariable<dim> ts_;
};



} // csmp

#endif
