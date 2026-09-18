// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure.h"
#include "Index.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Model.h"
#include "Exception.h"

#include <cmath>

using namespace std;
namespace csmp {

// =============================================================================
// Constructor
// =============================================================================

template<uint32_t dim>
NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::
    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure(
        const Model<dim>& model,
        const char*       oper,
        const char*       basic,
        const char*       test,
        bool              open_space,
        double            fault_perm_anisotropy,
        bool              with_air)
    : MathOperatorLHS<dim, InterFace>(model.Database(), oper, basic, test)
    // keys — all plain Index, matching SplitRegionLeakage naming exactly
    , thickness_key      (model.Database().StorageKey("thickness"))
    , p_key              (model.Database().StorageKey("fluid pressure"))
    , sat_halite_key     (model.Database().StorageKey("saturation halite"))
    , density_l_key      (model.Database().StorageKey("density liquid"))
    , ml_key             (model.Database().StorageKey("fluid mass liquid"))
    , viscosity_l_key    (model.Database().StorageKey("viscosity liquid"))
    , sat_l_key          (model.Database().StorageKey("saturation liquid"))
    , density_v_key      (model.Database().StorageKey("density vapor"))
    , mv_key             (model.Database().StorageKey("fluid mass vapor"))
    , viscosity_v_key    (model.Database().StorageKey("viscosity vapor"))
    , sat_v_key          (model.Database().StorageKey("saturation vapor"))

    , open_space_            (open_space)
    , fault_perm_anisotropy_ (fault_perm_anisotropy)
    , with_air_              (with_air)
{
    MathOperatorLHS<dim, InterFace>::Name(
        "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure", oper, basic, test);

    // Air keys assigned in body — plain Index assignment, no template mismatch
    if (with_air_)
    {
        density_a_key       = model.Database().StorageKey("density air");
        ma_key              = model.Database().StorageKey("fluid mass air");
        viscosity_a_key     = model.Database().StorageKey("viscosity air");
        sat_a_key            = model.Database().StorageKey("saturation air");
    }

    if (MathOperatorLHS<dim, InterFace>::BasicOperandPlacement() != NODE ||
        MathOperatorLHS<dim, InterFace>::BasicOperandType() != SCALAR)
        throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::(constructor)",
                              basic, "Operand (basic) must be a scalar property placed on the nodes.");

    if (MathOperatorLHS<dim, InterFace>::TestOperandPlacement() != NODE ||
        MathOperatorLHS<dim, InterFace>::TestOperandType() != SCALAR)
        throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::(constructor)",
                              test, "Operand (test) must be a scalar property placed on the nodes.");

    if (MathOperatorLHS<dim, InterFace>::MaterialOperandPlacement() != ELEMENT ||
        MathOperatorLHS<dim, InterFace>::MaterialOperandType() != SCALAR)
        throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim,CELL>::GetOperands",
                              "material operand placements other than on the CELL are not handled yet.");
}


// =============================================================================
// GetOperands
// =============================================================================

/**
 * Variable naming mirrors SplitRegionLeakage throughout:
 *   permeability_m / _i / _o   — fault / inner-parent / outer-parent permeability
 *   mul_in/mid/out             — liquid viscosity at in/mid/out node
 *   muv_, mua_                 — vapor / air viscosity
 *   rhol_mid, rhov_mid, rhoa_mid — phase densities at mid (used for gravity, same as SRL)
 *   ml_in/mid/out etc.         — mass concentrations (kg/m³ pore space)
 *   grad_pot_in_l etc.         — computed inline from current conditions (mirrors SplitRegionLeakage<dim>::ComputeFlowPotential())
 *   relative_perm_*            — computed inline from current conditions (mirrors SplitRegionLeakage<dim>::ComputeRelPerm())
 *   hydraulic_conductivity_*   — k * kr / mu, upstream-selected
 *   grav_component_in/out_*    — rho_mid * g * dy/dl, same formula as SRL::ComputeFlowPotential
 *
 * Results:
 *   transfer_coefficients_[i]             in-side  TC  = ml/mv/ma_up * hc / dx * area3
 *   transfer_coefficients_[i + n_nodes]   out-side TC
 *   transfer_coefficients_[i + 2*n_nodes] mid      TC  = TC_in + TC_out  (conservation)
 *   gravity_rhs_[i]                       in-side  gravity source
 *   gravity_rhs_[i + n_nodes]             out-side gravity source
 *   gravity_rhs_[i + 2*n_nodes]           mid      = -(in + out)         (conservation)
 */

template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::GetOperands(
    const InterFace<dim>& iface)
{
    assert(iface.HasInterveningElement());

    const uint32_t n_nodes = iface.InterveningElement()->Nodes();
    transfer_coefficients_.assign(n_nodes * 3U, 0.);
    gravity_rhs_.assign(n_nodes * 3U, 0.);

    if (!order1_FD_approximation_) return;

    // ----- Cell-level quantities (naming matches SplitRegionLeakage) -----
    const double thickness      = iface.InterveningElement()->Read(thickness_key);
    const double area3          = iface.InterveningElement()->Volume()
                         / static_cast<double>(n_nodes);

    assert(iface.FE()->MidSideNodes() == 0U);

    double permeability_m = iface.InterveningElement()->Read(
        MathOperatorLHS<dim, InterFace>::MaterialOperandKey());
    const double permeability_i = iface.InnerParent()->Read(
        MathOperatorLHS<dim, InterFace>::MaterialOperandKey());
    const double permeability_o = iface.OuterParent()->Read(
        MathOperatorLHS<dim, InterFace>::MaterialOperandKey());
    permeability_m *= fault_perm_anisotropy_;

    // Unit normal inside → outside, same as SplitRegionLeakage: unrml = (*sb_c)->UnitNormal(MIDDLE)
    const Point<dim> unrml = iface.UnitNormal(MIDDLE);

    // dy/dl along in_mid  =  unrml * distance_in_mid  → dy_in_mid  = dotProduct(unrml, vertical) * distance_in_mid
    // dy/dl along out_mid = -unrml * distance_out_mid → dy_out_mid = -dotProduct(unrml, vertical) * distance_out_mid
    // Gravity components divide by distance, leaving just: rho * g * dotProduct(unrml, vertical) (±)
    double dy_over_dl = 0.;
    if constexpr (dim == 3U) dy_over_dl = dotProduct(unrml, Point<3U>(0., 1., 0.));
    if constexpr (dim == 2U) dy_over_dl = dotProduct(unrml, Point<2U>(0., 1.));

    // ---- Inline relperm lambda — mirrors SplitRegionLeakage::ComputeRelPerm exactly ----
    //
    // Hard-coded residual liquid saturation of 30%, same as SRL.
    // Kr_l linear; Kr_v linear (no-air) or gas partitioned by effective saturation (with-air).
    // Sh >= 1 → all phases blocked.
    //
    // Sa is ignored (and may be passed as 0.) when !with_air_.

    auto computeRelPerm = [this](
                              double  Sl, double  Sv, double  Sh, double  Sa, double  p,
                              double& kr_l, double& kr_v, double& kr_a)
    {
        kr_l = kr_v = kr_a = 0.;

        if (Sh >= 1.) return;

        const double eff_l = Sl / (1. - Sh);
        const double eff_v = Sv / (1. - Sh);

        kr_l = (eff_l < 0.3) ? 0. : (eff_l - 0.3) / 0.7;

        if (!with_air_)
        {
            kr_v = (eff_l < 0.3) ? 1. : eff_v / 0.7;
        }
        else
        {
            const double eff_a  = Sa / (1. - Sh);
            const double kr_gas = 1. - kr_l;
            const double s_gas  = eff_v + eff_a;
            kr_v = (s_gas > 0.) ? kr_gas * eff_v / s_gas : 0.;
            kr_a = (s_gas > 0.) ? kr_gas * eff_a / s_gas : 0.;
        }
    };

    // ----- Per-node loop -----
    for (auto n = 0U; n < n_nodes; ++n)
    {
        auto* in_node  = iface.MatchingN(n, INSIDE);
        auto* mid_node = iface.MatchingN(n, MIDDLE);
        auto* out_node = iface.MatchingN(n, OUTSIDE);

        const bool at_periphery = (in_node == out_node);

        // ---- Distances (same logic as SplitRegionLeakage::ComputeFlowPotential) ----
        double distance_in_mid, distance_out_mid;
        if (at_periphery)
        {
            distance_in_mid  = thickness / 2.;
            distance_out_mid = thickness / 2.;
        }
        else
        {
            distance_in_mid  = (in_node->Coordinate()  - mid_node->Coordinate()).Length();
            distance_out_mid = (out_node->Coordinate() - mid_node->Coordinate()).Length();
        }

        //distance_in_mid = distance_out_mid = 0.5;//FOR BENCHMARKS ONLY!!!!!!!!!!!

        // ---- Saturations — read from current node state ----
        const double Sl_in  = in_node ->Read(sat_l_key);
        const double Sl_mid = mid_node->Read(sat_l_key);
        const double Sl_out = out_node->Read(sat_l_key);
        const double Sv_in  = in_node ->Read(sat_v_key);
        const double Sv_mid = mid_node->Read(sat_v_key);
        const double Sv_out = out_node->Read(sat_v_key);
        const double Sh_in  = in_node ->Read(sat_halite_key);
        const double Sh_mid = mid_node->Read(sat_halite_key);
        const double Sh_out = out_node->Read(sat_halite_key);

        // ---- Pressure — read from current node state ----
        const double p_in  = in_node ->Read(p_key);
        const double p_mid = mid_node->Read(p_key);
        const double p_out = out_node->Read(p_key);

        // ---- Densities-mid — read from current node state ----
        const double rhol_mid = mid_node->Read(density_l_key);
        const double rhov_mid = mid_node->Read(density_v_key);
        //24-04-2026 test fix
        const double rhol_in  = in_node->Read(density_l_key);
        const double rhov_in  = in_node->Read(density_v_key);
        const double rhol_out = out_node->Read(density_l_key);
        const double rhov_out = out_node->Read(density_v_key);
        //---

        double Sa_in = 0., Sa_mid = 0., Sa_out = 0.;
        if (with_air_)
        {
            Sa_in  = in_node ->Read(sat_a_key);
            Sa_mid = mid_node->Read(sat_a_key);
            Sa_out = out_node->Read(sat_a_key);
        }

        // ---- Inline relative permeabilities — mirrors SplitRegionLeakage::ComputeRelPerm exactly ----
        double relative_perm_in_l  = 0., relative_perm_in_v  = 0., relative_perm_in_a  = 0.;
        double relative_perm_mid_l = 0., relative_perm_mid_v = 0., relative_perm_mid_a = 0.;
        double relative_perm_out_l = 0., relative_perm_out_v = 0., relative_perm_out_a = 0.;
        computeRelPerm(Sl_in,  Sv_in,  Sh_in,  Sa_in,  p_in,  relative_perm_in_l,  relative_perm_in_v,  relative_perm_in_a);
        computeRelPerm(Sl_mid, Sv_mid, Sh_mid, Sa_mid, p_mid, relative_perm_mid_l, relative_perm_mid_v, relative_perm_mid_a);
        computeRelPerm(Sl_out, Sv_out, Sh_out, Sa_out, p_out, relative_perm_out_l, relative_perm_out_v, relative_perm_out_a);

        // ---- Viscosities ----
        const double mul_in  = in_node ->Read(viscosity_l_key);
        const double mul_mid = mid_node->Read(viscosity_l_key);
        const double mul_out = out_node->Read(viscosity_l_key);
        const double muv_in  = in_node ->Read(viscosity_v_key);
        const double muv_mid = mid_node->Read(viscosity_v_key);
        const double muv_out = out_node->Read(viscosity_v_key);

        // ---- Mass concentrations (kg/m³ pore space) ----
        const double ml_in  = in_node ->Read(ml_key);
        const double ml_mid = mid_node->Read(ml_key);
        const double ml_out = out_node->Read(ml_key);
        const double mv_in  = in_node ->Read(mv_key);
        const double mv_mid = mid_node->Read(mv_key);
        const double mv_out = out_node->Read(mv_key);

        // ---- Inline potential gradients and gravity source terms ----
        // Mirrors SplitRegionLeakage::ComputeFlowPotential exactly:
        //   in_mid  =  unrml * distance_in_mid  → dy_in_mid  =  dy_over_dl * distance_in_mid
        //   out_mid = -unrml * distance_out_mid → dy_out_mid = -dy_over_dl * distance_out_mid
        //   grav_component = rho_mid * g * dy/dl  (distances cancel)
        //   gradP = (p_mid - p_node) / distance
        //   grad_pot = gradP + grav_component
        //
        // Periphery guard matches SRL:
        //   gradP_in and grav_component_in_*: only when distance_in_mid > 0 && !at_periphery
        //   grav_component_out_*:             only when !at_periphery
        //   gradP_out:                        always when distance_out_mid > 0
        //   (in-side hydraulic block is already skipped at periphery below)

        double grad_pot_in_l  = 0., grad_pot_in_v  = 0., grad_pot_in_a  = 0.;
        double grad_pot_out_l = 0., grad_pot_out_v = 0., grad_pot_out_a = 0.;
        double grav_component_in_l  = 0., grav_component_in_v  = 0., grav_component_in_a  = 0.;
        double grav_component_out_l = 0., grav_component_out_v = 0., grav_component_out_a = 0.;

        if (distance_in_mid > 0. && !at_periphery)
        {
            const double gradP_in   = (p_mid - p_in) / distance_in_mid;

            //grav_component_in_l     =  9.80665 * rhol_mid * dy_over_dl;
            //grav_component_in_v     =  9.80665 * rhov_mid * dy_over_dl;

            // //24-04-2026 test improve
            double avg_rhol_in_mid = (rhol_in + rhol_mid) / 2.0;
            double avg_rhov_in_mid = (rhov_in + rhov_mid) / 2.0;

            grav_component_in_l     =  9.80665 * avg_rhol_in_mid * dy_over_dl;
            grav_component_in_v     =  9.80665 * avg_rhov_in_mid * dy_over_dl;

            grad_pot_in_l           = gradP_in + grav_component_in_l;
            grad_pot_in_v           = gradP_in + grav_component_in_v;

            if (with_air_)
            {
                //const double rhoa_mid = mid_node->Read(density_a_key);
                //grav_component_in_a  =  9.80665 * rhoa_mid * dy_over_dl;

                //24-04-2026 test improve
                const double rhoa_mid = mid_node->Read(density_a_key);
                const double rhoa_in  = in_node->Read(density_a_key);

                double avg_rhoa_in_mid = (rhoa_in + rhoa_mid) / 2.0;
                grav_component_in_a  =  9.80665 * avg_rhoa_in_mid * dy_over_dl;

                grad_pot_in_a        = gradP_in + grav_component_in_a;
            }
        }

        if (distance_out_mid > 0.)
        {
            const double gradP_out  = (p_mid - p_out) / distance_out_mid;
            if (!at_periphery)
            {
                //grav_component_out_l =  -9.80665 * rhol_mid * dy_over_dl;
                //grav_component_out_v =  -9.80665 * rhov_mid * dy_over_dl;

                //24-04-2026 test improve
                double avg_rhol_out_mid = (rhol_out + rhol_mid) / 2.0;
                double avg_rhov_out_mid = (rhov_out + rhov_mid) / 2.0;

                grav_component_out_l =  -9.80665 * avg_rhol_out_mid * dy_over_dl;
                grav_component_out_v =  -9.80665 * avg_rhov_out_mid * dy_over_dl;

                if (with_air_)
                {
                    //const double rhoa_mid  = mid_node->Read(density_a_key);
                    //grav_component_out_a   = -9.80665 * rhoa_mid * dy_over_dl;

                    //24-04-2026 test improve
                    const double rhoa_mid  = mid_node->Read(density_a_key);
                    const double rhoa_out  = out_node->Read(density_a_key);

                    double avg_rhoa_out_mid = (rhoa_out + rhoa_mid) / 2.0;
                    grav_component_out_a   = -9.80665 * avg_rhoa_out_mid * dy_over_dl;
                }
            }

            grad_pot_out_l = gradP_out + grav_component_out_l;
            grad_pot_out_v = gradP_out + grav_component_out_v;
            if (with_air_) grad_pot_out_a = gradP_out + grav_component_out_a;
        }

        // =========================================================
        // IN-SIDE  (skip at periphery — matches SplitRegionLeakage)
        // =========================================================
        if (distance_in_mid > 0. && !at_periphery)
        {
            // --- Liquid ---
            {
                double hydraulic_conductivity_in_l = 0.;
                double mC_in_l = 0.;
                const double perm = open_space_ ? permeability_i : permeability_m;
                if (grad_pot_in_l <= 0.)
                {
                    // flow l in → mid: upstream = in
                    if (mul_in > 0.) hydraulic_conductivity_in_l = perm * relative_perm_in_l / mul_in;
                    mC_in_l = ml_in;
                }
                else
                {
                    // flow l mid → in: upstream = mid
                    if (mul_mid > 0.) hydraulic_conductivity_in_l = perm * relative_perm_mid_l / mul_mid;
                    mC_in_l = ml_mid;
                }
                transfer_coefficients_[n] += mC_in_l * hydraulic_conductivity_in_l / distance_in_mid * area3;
                gravity_rhs_[n]           += grav_component_in_l * hydraulic_conductivity_in_l * mC_in_l * area3;
            }

            // --- Vapor ---
            {
                double hydraulic_conductivity_in_v = 0.;
                double mC_in_v = 0.;
                const double perm = open_space_ ? permeability_i : permeability_m;
                if (grad_pot_in_v <= 0.)
                {
                    if (muv_in > 0.) hydraulic_conductivity_in_v = perm * relative_perm_in_v / muv_in;
                    mC_in_v = mv_in;
                }
                else
                {
                    if (muv_mid > 0.) hydraulic_conductivity_in_v = perm * relative_perm_mid_v / muv_mid;
                    mC_in_v = mv_mid;
                }
                transfer_coefficients_[n] += mC_in_v * hydraulic_conductivity_in_v / distance_in_mid * area3;
                gravity_rhs_[n]           += grav_component_in_v * hydraulic_conductivity_in_v * mC_in_v * area3;
            }

            // --- Air ---
            if (with_air_)
            {
                // grad_pot_in_a and grav_component_in_a computed inline above
                const double mua_in  = in_node ->Read(viscosity_a_key);
                const double mua_mid = mid_node->Read(viscosity_a_key);
                const double ma_in   = in_node ->Read(ma_key);
                const double ma_mid  = mid_node->Read(ma_key);

                double hydraulic_conductivity_in_a = 0.;
                double mC_in_a = 0.;
                const double perm = open_space_ ? permeability_i : permeability_m;
                if (grad_pot_in_a <= 0.)
                {
                    if (mua_in  > 0.) hydraulic_conductivity_in_a = perm * relative_perm_in_a  / mua_in;
                    mC_in_a = ma_in;
                }
                else
                {
                    if (mua_mid > 0.) hydraulic_conductivity_in_a = perm * relative_perm_mid_a / mua_mid;
                    mC_in_a = ma_mid;
                }
                transfer_coefficients_[n] += mC_in_a * hydraulic_conductivity_in_a / distance_in_mid * area3;
                gravity_rhs_[n]           += grav_component_in_a * hydraulic_conductivity_in_a * mC_in_a * area3;
            }
        }

        // =========================================================
        // OUT-SIDE
        // =========================================================
        if (distance_out_mid > 0.)
        {
            // --- Liquid ---
            {
                double hydraulic_conductivity_out_l = 0.;
                double mC_out_l = 0.;
                const double perm = open_space_ ? permeability_o : permeability_m;
                if (grad_pot_out_l <= 0.)
                {
                    // flow l out → mid: upstream = out
                    if (mul_out > 0.) hydraulic_conductivity_out_l = perm * relative_perm_out_l / mul_out;
                    mC_out_l = ml_out;
                }
                else
                {
                    // flow l mid → out: upstream = mid
                    if (mul_mid > 0.) hydraulic_conductivity_out_l = perm * relative_perm_mid_l / mul_mid;
                    mC_out_l = ml_mid;
                }
                transfer_coefficients_[n + n_nodes] += mC_out_l * hydraulic_conductivity_out_l / distance_out_mid * area3;
                gravity_rhs_[n + n_nodes]           += grav_component_out_l * hydraulic_conductivity_out_l * mC_out_l * area3;
            }

            // --- Vapor ---
            {
                double hydraulic_conductivity_out_v = 0.;
                double mC_out_v = 0.;
                const double perm = open_space_ ? permeability_o : permeability_m;
                if (grad_pot_out_v <= 0.)
                {
                    if (muv_out > 0.) hydraulic_conductivity_out_v = perm * relative_perm_out_v / muv_out;
                    mC_out_v = mv_out;
                }
                else
                {
                    if (muv_mid > 0.) hydraulic_conductivity_out_v = perm * relative_perm_mid_v / muv_mid;
                    mC_out_v = mv_mid;
                }
                transfer_coefficients_[n + n_nodes] += mC_out_v * hydraulic_conductivity_out_v / distance_out_mid * area3;
                gravity_rhs_[n + n_nodes]           += grav_component_out_v * hydraulic_conductivity_out_v * mC_out_v * area3;
            }

            // --- Air ---
            if (with_air_)
            {
                // grad_pot_out_a and grav_component_out_a computed inline above
                const double mua_out  = out_node->Read(viscosity_a_key);
                const double mua_mid  = mid_node->Read(viscosity_a_key);
                const double ma_out   = out_node->Read(ma_key);
                const double ma_mid   = mid_node->Read(ma_key);

                double hydraulic_conductivity_out_a = 0.;
                double mC_out_a = 0.;
                const double perm = open_space_ ? permeability_o : permeability_m;
                if (grad_pot_out_a <= 0.)
                {
                    if (mua_out > 0.) hydraulic_conductivity_out_a = perm * relative_perm_out_a  / mua_out;
                    mC_out_a = ma_out;
                }
                else
                {
                    if (mua_mid > 0.) hydraulic_conductivity_out_a = perm * relative_perm_mid_a / mua_mid;
                    mC_out_a = ma_mid;
                }
                transfer_coefficients_[n + n_nodes] += mC_out_a * hydraulic_conductivity_out_a / distance_out_mid * area3;
                gravity_rhs_[n + n_nodes]           += grav_component_out_a * hydraulic_conductivity_out_a * mC_out_a * area3;
            }
        }

        // =========================================================
        // MID-NODE: sum of both sides (conservation)
        // =========================================================
        transfer_coefficients_[n + 2*n_nodes] =
            transfer_coefficients_[n] + transfer_coefficients_[n + n_nodes];

        gravity_rhs_[n + 2*n_nodes] =
            -(gravity_rhs_[n] + gravity_rhs_[n + n_nodes]);
    }
}


// =============================================================================
// ComputeContribution
// =============================================================================

template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::ComputeContribution(
    const InterFace<dim>& iface)
{
    assert(iface.UsesLocalCoordinates());
    assert(iface.HasInterveningElement());

    const uint32_t n_elmt_nodes  = iface.InterveningElement()->Nodes();
    const uint32_t n_total_nodes = iface.Nodes() + n_elmt_nodes;

    MathOperatorLHS<dim, InterFace>::LHS.Resize(n_total_nodes, n_total_nodes);
    MathOperatorLHS<dim, InterFace>::LHS.Zero();

    // Diagonal
    for (auto i = 0U; i < n_total_nodes; ++i)
        MathOperatorLHS<dim, InterFace>::LHS(i, i) = transfer_coefficients_[i];

    // Off-diagonal: in-nodes ↔ mid-nodes
    for (auto i = 0U; i < n_elmt_nodes; ++i)
    {
        MathOperatorLHS<dim, InterFace>::LHS(i,                  i + 2*n_elmt_nodes) = -transfer_coefficients_[i];
        MathOperatorLHS<dim, InterFace>::LHS(i + 2*n_elmt_nodes, i)                  = -transfer_coefficients_[i];
    }

    // Off-diagonal: out-nodes ↔ mid-nodes
    for (auto i = 0U; i < n_elmt_nodes; ++i)
    {
        MathOperatorLHS<dim, InterFace>::LHS(i + n_elmt_nodes,   i + 2*n_elmt_nodes) = -transfer_coefficients_[i + n_elmt_nodes];
        MathOperatorLHS<dim, InterFace>::LHS(i + 2*n_elmt_nodes, i + n_elmt_nodes)   = -transfer_coefficients_[i + n_elmt_nodes];
    }
}


// =============================================================================
// AssignToGlobal
// =============================================================================

/**
 * Scatters the LHS pressure coupling matrix into G and scatters gravity_rhs_
 * into pivotVector, replacing the former "split region gravity mass source"
 * explicit injection step.
 */
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::AssignToGlobal(
    const InterFace<dim>& iface,
    SparseMatrix&         G,
    vector<double>&       pivotVector,
    const vector<size_t>& DOF_indexes)
{
    const uint32_t n_elmt_nodes  = iface.InterveningElement()->Nodes();
    const uint32_t n_total_nodes = iface.Nodes() + n_elmt_nodes;
    const uint32_t n_elmt_nodes2 = n_elmt_nodes * 2U;

    vector<size_t> IDT(n_total_nodes);
    vector<size_t> IDB(n_total_nodes);

    for (auto i = 0U; i < n_elmt_nodes; ++i)
    {
        IDT[i]                = iface.MatchingN(i, INSIDE )->Idx();
        IDT[i + n_elmt_nodes] = iface.MatchingN(i, OUTSIDE)->Idx();
        IDT[i + n_elmt_nodes2]= iface.MatchingN(i, MIDDLE )->Idx();
        IDB[i]                = IDT[i];
        IDB[i + n_elmt_nodes] = IDT[i + n_elmt_nodes];
        IDB[i + n_elmt_nodes2]= IDT[i + n_elmt_nodes2];
    }

    assert(this->BasicOperandType() == SCALAR);
    assert(this->TestOperandType()  == SCALAR);

    for (auto i = 0U; i < IDT.size(); ++i)
    {
        IDT[i] += this->TestOperandOffset();
        IDT[i]  = DOF_indexes[IDT[i]];
    }
    for (auto i = 0U; i < IDB.size(); ++i)
    {
        IDB[i] += this->BasicOperandOffset();
        IDB[i]  = DOF_indexes[IDB[i]];
    }

    vector<double> nodal_values(IDB.size(), 1.);
    for (auto nIdx = 0U; nIdx < n_elmt_nodes; ++nIdx)
    {
        nodal_values[nIdx]                = iface.MatchingN(nIdx, INSIDE )->Read(this->TestOperandKey());
        nodal_values[nIdx + n_elmt_nodes] = iface.MatchingN(nIdx, OUTSIDE)->Read(this->TestOperandKey());
        nodal_values[nIdx + n_elmt_nodes2]= iface.MatchingN(nIdx, MIDDLE )->Read(this->TestOperandKey());
    }

    if (this->multiply_accumulate_)
    {
        throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::AssignToGlobal",
                              "Multiply Accumulate not supported.");
    }
    else if (this->add_accumulate_ || this->add_accumulate_later_)
    {
        for (auto i = 0U; i < this->LHS.Rows(); ++i)
        {
            if (IDT[i] == NULL_IDX) continue;
            for (auto j = 0U; j < this->LHS.Cols(); ++j)
            {
                if (IDB[j] == NULL_IDX)
                    pivotVector[IDT[i]] -= this->LHS(i, j) * nodal_values[j] * this->factor_;
                else
                    G.Add(IDT[i], IDB[j], this->LHS(i, j) * this->factor_);
            }
        }
        //24-04-2026 commented, this was an attempt to integrate gravity terms directly, does not work!
        for (auto i = 0U; i < n_total_nodes; ++i)
        {
            if (IDT[i] == NULL_IDX) continue;
            pivotVector[IDT[i]] += gravity_rhs_[i] * this->factor_ * delta_t_;
        }
    }
    else if (this->subtract_accumulate_ || this->subtract_accumulate_later_)
    {
        for (auto i = 0U; i < this->LHS.Rows(); ++i)
        {
            if (IDT[i] == NULL_IDX) continue;
            for (auto j = 0U; j < this->LHS.Cols(); ++j)
            {
                if (IDB[j] == NULL_IDX)
                    pivotVector[IDT[i]] += this->LHS(i, j) * nodal_values[j] * this->factor_;
                else
                    G.Add(IDT[i], IDB[j], -this->LHS(i, j) * this->factor_);
            }
        }
        //24-04-2026 commented, this was an attempt to integrate gravity terms directly, does not work!
        for (auto i = 0U; i < n_total_nodes; ++i)
        {
            if (IDT[i] == NULL_IDX) continue;
            pivotVector[IDT[i]] -= gravity_rhs_[i] * this->factor_ * delta_t_;
        }
    }
    else
    {
        throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::AssignToGlobal",
                              "accumulation instructions could not be parsed.");
    }
}

template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::ActivateAir(
    const Model<dim>& model)
{
    with_air_     = true;
    density_a_key   = model.Database().StorageKey("density air");
    ma_key          = model.Database().StorageKey("fluid mass air");
    viscosity_a_key = model.Database().StorageKey("viscosity air");
    sat_a_key       = model.Database().StorageKey("saturation air");
}

// =============================================================================
// Analytic diffusion helpers (unchanged)
// =============================================================================

template<uint32_t dim>
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::Diffusion1D(
    double val_farfield, double diffusivity, double x, double t)
{
    return val_farfield * erf(x / (2. * sqrt(diffusivity * t)));
}

template<uint32_t dim>
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::Grad_Var_AtX0(
    double val_farfield, double diffusivity, double t)
{
    constexpr double pi     = 3.14159265358979323846;
    const double     sqrtPi = sqrt(pi);
    double dvaldx = val_farfield / (diffusivity * t);
    dvaldx /= sqrtPi * sqrt(diffusivity * t);
    assert(fabs(dvaldx) < 1e20);
    return dvaldx;
}

template<uint32_t dim>
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>::Flux1DAtX0(
    double val_farfield, double diffusivity, double t, double conductivity)
{
    return Grad_Var_AtX0(val_farfield, diffusivity, t) * conductivity;
}


// =============================================================================
// Explicit template instantiations
// =============================================================================

template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<1U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<2U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<3U>;

} // namespace csmp
