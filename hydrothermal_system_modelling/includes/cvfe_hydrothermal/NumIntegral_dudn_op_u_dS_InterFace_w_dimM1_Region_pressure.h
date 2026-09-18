// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DUDN_OP_U_DS_P_H
#define NUM_INTEGRAL_DUDN_OP_U_DS_P_H

#include "MathOperatorLHS.h"
#include "Index.h"
#include "InterFace.h"

namespace csmp {

class FiniteElement;
template<uint32_t> class Model;
template<uint32_t> class InterFace;

/**
    Assembles the implicit pressure coupling matrix and the explicit gravity RHS
    for a split boundary with an intervening lower-dimensional region.

    The transfer coefficient for each phase φ ∈ {l, v, [a]} and side s ∈ {in, out} is:

        TC_s_φ = ml/mv/ma_upstream * (k * kr_φ / μ_φ) / dx_s * area3

    The gravity RHS contribution (scattered into pivotVector in AssignToGlobal) is:

        grav_s_φ = rho_mid_φ * g * (dy/dl)_s * (k * kr_φ / μ_φ) * ml/mv/ma_upstream * area3

    Both terms use identical upstream weighting, ensuring consistency with
    SplitRegionLeakage::ComputeGravity_PressureSourceTerm.  Assembling gravity here
    replaces the former "split region gravity mass source" explicit injection step.

    Variable names mirror SplitRegionLeakage throughout.

    @param with_air   Must match the with_air flag of the paired SplitRegionLeakage instance.
    @param open_space If true, use host-rock permeability instead of fault permeability.

    @author SKM / Benoit LC
*/
template<uint32_t dim>
class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure
    : public MathOperatorLHS<dim, InterFace>
{
public:
    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure(
        const Model<dim>&,
        const char* oper,
        const char* basic,
        const char* test,
        bool   open_space           = false,
        double fault_perm_anisotropy = 1.,
        bool   with_air             = false);

    virtual void GetOperands(const InterFace<dim>&) override;
    virtual void ComputeContribution(const InterFace<dim>&) override;
    virtual void AssignToGlobal(const InterFace<dim>&,
                                SparseMatrix&,
                                std::vector<double>&,
                                const std::vector<size_t>&) override;

    void ActivateAir(const Model<dim>& model);

    void SetSplitRegionParameters(bool open_space, double fault_perm_anisotropy)
    {
        open_space_            = open_space;
        fault_perm_anisotropy_ = fault_perm_anisotropy;
    }

    void UpdateTimeIncrement(double dt) { delta_t_ = dt; }

    virtual NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>* clone() const override
    { return new NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>(*this); }

    void ResetMinMaxTransferTerms() {
        transfer_min_ =  1e30;
        transfer_max_ = -1e30;
    }

private:
    double Diffusion1D(double val_farfield, double diffusivity, double x, double t);
    double Grad_Var_AtX0(double val_farfield, double diffusivity, double t);
    double Flux1DAtX0(double val_farfield, double diffusivity, double t, double conductivity);

private:
    // All keys are plain Index (same as SplitRegionLeakage) so that conditional
    // body-assignment in the constructor compiles without issue.

    csmp::Index  thickness_key;
    csmp::Index  p_key;

    csmp::Index sat_halite_key;

    // --- Liquid ---
    csmp::Index  density_l_key;
    csmp::Index  ml_key;
    csmp::Index  viscosity_l_key;
    csmp::Index  sat_l_key;

    // --- Vapor ---
    csmp::Index  density_v_key;
    csmp::Index  mv_key;
    csmp::Index  viscosity_v_key;
    csmp::Index sat_v_key;

    // --- Air (valid only when with_air_ == true) ---
    csmp::Index  density_a_key;
    csmp::Index  ma_key;
    csmp::Index  viscosity_a_key;
    csmp::Index sat_a_key;

    // --- Permeability key (same string "permeability" as SplitRegionLeakage k_key) ---
    // Retrieved via MaterialOperandKey() from the base class; no separate member needed.

    /// Transfer coefficients for the implicit pressure LHS.
    std::vector<double> transfer_coefficients_;

    /// Gravity RHS contributions, scattered into pivotVector in AssignToGlobal.
    std::vector<double> gravity_rhs_;

    double delta_t_             = 0.;
    double transfer_min_        =  1e30;
    double transfer_max_        = -1e30;
    const bool order1_FD_approximation_ = true;
    bool   open_space_;
    double fault_perm_anisotropy_;
    bool   with_air_;
};

} // namespace csmp

#endif
