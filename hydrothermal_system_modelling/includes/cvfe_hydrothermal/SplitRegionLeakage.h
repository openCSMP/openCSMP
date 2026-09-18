// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SPLIT_REGION_LEAKAGE_H
#define SPLIT_REGION_LEAKAGE_H


#include "Model.h"
#include "SplitBoundary.h"
#include <map>
#include <vector>
#include <set>
#include <string>

namespace csmp {

template<uint32_t> class Model;

/**
        @class SplitRegionLeakage SplitRegionLeakage.h

        Calculating leakage of mass and energy between fault (split boundary) and host.

        @author Benoit LC, ETH Zuerich
        @section contact Contact


        @changes changes Latest Changes

        @section motivation Motivation

        @section usage Usage

        @code

        @endcode

        @section dependencies Dependencies

        @section issues Known issues

        @section testing Testing
        testing was done in 2026

    */

template<uint32_t dim>
class SplitRegionLeakage {

public:
    SplitRegionLeakage(Model<dim> &model,
                       bool open_space = false,
                       double fault_perm_anisotropy = 1.);

    virtual ~SplitRegionLeakage();

    void WithSpeciesTransport(
        Model<dim> &model,
        bool with_air = false,
        bool with_magmatic_fluids = false,
        bool with_tracer=false,
        bool with_gold=false,
        bool with_lithium=false);

    void SetTimeStep(double dt_ext);
    bool ComputeGravity_PressureSourceTerm();

    // ── Leakage-in-the-loop API (frozen-rate split) ──────────────────────
    // PrepareStageRates(): run the full pass ONCE per transport stage with
    // dt=1 in capture mode — geometry, upwinding, viz stores as usual, but
    // per-triple leg fluxes are recorded as RATES (per second) and nothing
    // is applied. ExchangeSubstep(sub_dt): apply rate*sub_dt per triple with
    // fresh node state and per-BEAT drain budgets. Same frozen-flux contract
    // as the matrix transport stage (rates from stage-start state; per-beat
    // positivity comes from the budgets, replacing the per-leg caps).
    void PrepareStageRates();
    void ExchangeSubstep(double sub_dt);

    // Per-node, per-phase drain rates in PER-PORE-VOLUME units (1/s scale of
    // the transport LHS), captured at Prepare. Lets CheckDry add the frozen
    // leakage demand to the transport outflow when sizing sub-steps, with no
    // pore-volume plumbing on the transport side.
    struct PhaseDrainRates { double l = 0., v = 0., a = 0.; };
    const std::map<Node<dim>*, PhaseDrainRates>& NodeDrainRatesPerPV() const
    { return node_drain_rates_perpv_; }
    size_t StageRateCount() const { return stage_rates_.size(); }

    // ── interface clamp ledger (read by the transport stage ledger) ──────
    // ExchangeSubstep accumulates the flux WITHHELD by the per-beat budget
    // clamp (requested removal minus applied), and the mass actually applied,
    // per phase, in absolute (kg) units. Transport resets these once per stage
    // and reads them into its combined ledger, split out as the interface's
    // contribution vs the matrix transport's.
    // Interface standing-queue measurement (parallels transport's): each
    // beat's total budget-withheld IS the queue standing at the fault right
    // then. parked_now = the LAST beat's withheld (stage carryover);
    // parked_max = the largest per-beat withheld (worst transient queue);
    // applied = total drained over the stage (for the ratio). Reading:
    // max >> carryover -> the fault queue formed and drained (good);
    // carryover ~ max -> it persisted to stage end (the real error proxy).
    struct InterfaceLedger { double parked_now = 0., parked_max = 0., applied = 0.;
        double worst = 0.; size_t worst_idx = 0; };
    void   ResetInterfaceLedger() { iface_ledger_ = InterfaceLedger(); }
    const InterfaceLedger& GetInterfaceLedger() const { return iface_ledger_; }

    // Reason for the most recent dt-cut request (mode, node, fraction), for logging.
    const std::string& GetCutReason() const { return sb_cut_reason_; }

private:

    void ReadVariables();
    void ComputeRelPerm();
    void ComputeHydraulicConductivity();

    void TotalNetFluxes_for_visualization();
    void UpdateTransportVariables();
    void ComputeFlowPotential();

    Model< dim >                  &model_ref_;
    const PropertyDatabase<dim>   &prop_ref_;

    bool cut_dt, overdraining, open_space_;
    bool
        with_air_     = false,
        with_tracer_  = false,
        with_gold_    = false,
        with_lithium_ = false,
        with_magmatic_fluids_ = false;

    double fault_perm_anisotropy_;

    Index
        p_key,pp_key, nfvs_key,
        qM_key, qE_key, qS_key,
        k_key,
        qSO2_key, qAuHS0_key,

        viscosity_l_key, viscosity_v_key, viscosity_a_key,
        sat_l_key, sat_v_key, sat_a_key, sat_halite_key,
        density_l_key, density_v_key, density_a_key,

        pore_volume_key,

        ml_key, mv_key, ma_key,
        hCl_key, hCv_key, hCa_key,
        hv_key,   // specific vapor enthalpy (diagnostic: hCv=hv*rv*sv)
        xCl_key, xCv_key,
        TrCl_key, TrCv_key,
        LiCl_key, LiCv_key,
        so2Cl_key, so2Cv_key,
        auhs0Cl_key, auhs0Cv_key,

        thickness_key,

        grav_mass_source_key,
        split_grad_pot_l_key, split_grad_pot_v_key, split_grad_pot_a_key,
        split_relperm_l_key, split_relperm_v_key, split_relperm_a_key;


    ScalarVariable
        qM_in, qM_out, qM_mid, qE_in, qE_out, qE_mid,
        qS_in, qS_out, qS_mid, qT_in, qT_out, qT_mid,
	qSO2_in, qSO2_mid, qSO2_out,
	qAuHS0_in, qAuHS0_mid, qAuHS0_out;

    ScalarVariable
        mul_in, muv_in, mua_in,
        Sl_in, Sv_in, Sa_in, Sh_in,
        rhol_in, rhov_in, rhoa_in;

    ScalarVariable
        mul_out, muv_out, mua_out,
        Sl_out, Sv_out, Sa_out, Sh_out,
        rhol_out, rhov_out, rhoa_out;

    ScalarVariable
        mul_mid, muv_mid, mua_mid,
        Sl_mid, Sv_mid, Sa_mid, Sh_mid,
        rhol_mid, rhov_mid, rhoa_mid;

    ScalarVariable pore_volume_in, pore_volume_out, pore_volume_mid;

    ScalarVariable ml_in, ml_out, ml_mid;
    ScalarVariable mv_in, mv_out, mv_mid;
    ScalarVariable ma_in, ma_out, ma_mid;

    ScalarVariable hCl_in, hCl_out, hCl_mid;
    ScalarVariable hCv_in, hCv_out, hCv_mid;
    ScalarVariable hCa_in, hCa_out, hCa_mid;

    ScalarVariable xCl_in, xCl_out, xCl_mid;
    ScalarVariable xCv_in, xCv_out, xCv_mid;
    ScalarVariable TrCl_in, TrCl_out, TrCl_mid;
    ScalarVariable TrCv_in, TrCv_out, TrCv_mid;
    ScalarVariable LiCl_in, LiCl_out, LiCl_mid;
    ScalarVariable LiCv_in, LiCv_out, LiCv_mid;
    ScalarVariable so2Cl_in, so2Cl_out, so2Cl_mid;
    ScalarVariable so2Cv_in, so2Cv_out, so2Cv_mid;
    ScalarVariable auhs0Cl_in, auhs0Cl_out, auhs0Cl_mid;
    ScalarVariable auhs0Cv_in, auhs0Cv_out, auhs0Cv_mid;

    ScalarVariable thickness;

    ScalarVariable pressure_in, pressure_out, pressure_mid;
    ScalarVariable ppressure_in, ppressure_out, ppressure_mid;
    ScalarVariable nfvs_in, nfvs_out, nfvs_mid;
    ScalarVariable permeability_m, permeability_i, permeability_o;

    ScalarVariable
        grav_mass_source_in_tot,
        grav_mass_source_out_tot,
        grav_mass_source_mid_tot;

    ScalarVariable
        grad_pot_in_l,
        grad_pot_in_v,
        grad_pot_in_a,
        grad_pot_out_l,
        grad_pot_out_v,
        grad_pot_out_a;

    ScalarVariable
        relative_perm_in_l,
        relative_perm_in_v,
        relative_perm_in_a,
        relative_perm_mid_l,
        relative_perm_mid_v,
        relative_perm_mid_a,
        relative_perm_out_l,
        relative_perm_out_v,
        relative_perm_out_a;

    Node<dim> *in_node;
    Node<dim> *out_node;
    Node<dim> *mid_node;

    Point<dim> unrml;

    bool with_SB_leakage_report = false;

    double
        grav_mass_source_in_l,
        grav_mass_source_out_l,
        grav_mass_source_mid_l;

    double
        grav_mass_source_in_v,
        grav_mass_source_out_v,
        grav_mass_source_mid_v;

    double
        grav_mass_source_in_a,
        grav_mass_source_out_a,
        grav_mass_source_mid_a;

    double
        dt,
        distance_in_mid,
        distance_out_mid,
        gradP_in, gradP_out,
        area,

        hydraulic_conductivity_in_l, hydraulic_conductivity_in_v, hydraulic_conductivity_in_a,
        hydraulic_conductivity_out_l, hydraulic_conductivity_out_v, hydraulic_conductivity_out_a,

        effective_sat_in_l, effective_sat_in_v, effective_sat_in_a,
        effective_sat_mid_l, effective_sat_mid_v, effective_sat_mid_a,
        effective_sat_out_l, effective_sat_out_v,  effective_sat_out_a,

        net_Mflux_in_l, net_Mflux_in_v,  net_Mflux_in_a,
        net_Mflux_out_l, net_Mflux_out_v, net_Mflux_out_a,

        CFL_outflow_l_mid, CFL_outflow_v_mid, CFL_outflow_a_mid,
        CFL_outflow_l_in, CFL_outflow_v_in, CFL_outflow_a_in,
        CFL_outflow_l_out, CFL_outflow_v_out, CFL_outflow_a_out,

        CFL_inflow_l_mid, CFL_inflow_v_mid, CFL_inflow_a_mid,
        CFL_inflow_l_in, CFL_inflow_v_in, CFL_inflow_a_in,
        CFL_inflow_l_out, CFL_inflow_v_out, CFL_inflow_a_out,

        CFL_drain_max, CFL_stuff_max,

        net_Eflux_in_l, net_Eflux_in_v, net_Eflux_in_a,
        net_Eflux_out_l, net_Eflux_out_v, net_Eflux_out_a,

        net_Xflux_in_l, net_Xflux_in_v,
        net_Xflux_out_l, net_Xflux_out_v,

        net_Trflux_in_l, net_Trflux_out_l, net_Trflux_in_v, net_Trflux_out_v,
        net_Liflux_in_l, net_Liflux_out_l, net_Liflux_in_v, net_Liflux_out_v,
        net_so2flux_in_l, net_so2flux_out_l, net_so2flux_in_v, net_so2flux_out_v,
        net_auhs0flux_in_l, net_auhs0flux_out_l, net_auhs0flux_in_v, net_auhs0flux_out_v,

        mass_exchange_term_in_l, mass_exchange_term_out_l, mass_exchange_term_mid_l,
        mass_exchange_term_in_v, mass_exchange_term_out_v, mass_exchange_term_mid_v,
        mass_exchange_term_in_a, mass_exchange_term_out_a, mass_exchange_term_mid_a,

        energy_exchange_term_in_l, energy_exchange_term_out_l, energy_exchange_term_mid_l,
        energy_exchange_term_in_v, energy_exchange_term_out_v, energy_exchange_term_mid_v,
        energy_exchange_term_in_a, energy_exchange_term_out_a, energy_exchange_term_mid_a,

        salt_mass_exchange_term_in_l, salt_mass_exchange_term_out_l, salt_mass_exchange_term_mid_l,
        salt_mass_exchange_term_in_v, salt_mass_exchange_term_out_v, salt_mass_exchange_term_mid_v,

        tracer_mass_exchange_term_in_l, tracer_mass_exchange_term_out_l, tracer_mass_exchange_term_mid_l,
        tracer_mass_exchange_term_in_v, tracer_mass_exchange_term_out_v, tracer_mass_exchange_term_mid_v,
        
	lithium_mass_exchange_term_in_l, lithium_mass_exchange_term_out_l, lithium_mass_exchange_term_mid_l,
        lithium_mass_exchange_term_in_v, lithium_mass_exchange_term_out_v, lithium_mass_exchange_term_mid_v,
        
	so2_mass_exchange_term_in_l, so2_mass_exchange_term_out_l, so2_mass_exchange_term_mid_l,
        so2_mass_exchange_term_in_v, so2_mass_exchange_term_out_v, so2_mass_exchange_term_mid_v,
        
	auhs0_mass_exchange_term_in_l, auhs0_mass_exchange_term_out_l, auhs0_mass_exchange_term_mid_l,
        auhs0_mass_exchange_term_in_v, auhs0_mass_exchange_term_out_v, auhs0_mass_exchange_term_mid_v,

        grav_component_in_l, grav_component_in_v, grav_component_in_a,
        grav_component_out_l, grav_component_out_v, grav_component_out_a,

        volume_l_removed_from_mid, volume_v_removed_from_mid, volume_a_removed_from_mid,
        volume_l_removed_from_in, volume_v_removed_from_in, volume_a_removed_from_in,
        volume_l_removed_from_out, volume_v_removed_from_out, volume_a_removed_from_out,

        volume_l_added_to_mid, volume_v_added_to_mid, volume_a_added_to_mid,
        volume_l_added_to_in, volume_v_added_to_in, volume_a_added_to_in,
        volume_l_added_to_out, volume_v_added_to_out, volume_a_added_to_out,

        mass_l_removed_from_mid, mass_v_removed_from_mid, mass_a_removed_from_mid,
        mass_l_removed_from_in, mass_v_removed_from_in, mass_a_removed_from_in,
        mass_l_removed_from_out, mass_v_removed_from_out, mass_a_removed_from_out,

        available_liquid_mid, available_vapor_mid, available_air_mid,
        available_liquid_in, available_vapor_in, available_air_in,
        available_liquid_out, available_vapor_out, available_air_out;

    std::set<decltype(mid_node)> reported_nodes;

    // Per-pass drain budgets (beat_mode_ clamp): Node* -> remaining allowed
    // removal (mass) per phase. Set on a node's FIRST visit of the leakage pass
    // to (1-sb_drain_margin)*start-of-pass available; every visit draws from it.
    // Keying by Node* makes all sharing automatic: parent-cell revisits, the
    // in/out sides of a periphery doublet, and matrix nodes shared between
    // adjacent triples all draw from ONE budget, so every node lands at exactly
    // margin*(start-of-pass holdup) regardless of topology.
    struct DrainBudget { double l = 0., v = 0., a = 0.; bool init = false; };
    std::map<decltype(mid_node), DrainBudget> drain_budgets_;

    // ── frozen per-triple leg VELOCITIES (rebuilt every stage by
    //    PrepareStageRates) ─────────────────────────────────────────────────
    // FREEZE BOUNDARY = the velocity field, mirroring the matrix transport
    // sub-cycle. The Darcy mass flux is  grad_pot * K * area * m_donor * dt.
    // grad_pot depends on the frozen pressure field, K and area are geometry —
    // so  vel = grad_pot * K * area  is the genuinely frozen (per-second)
    // coefficient. The donor mass m_donor is NOT frozen: it is re-read fresh
    // every beat in ExchangeSubstep, so the flux shrinks as a donor drains
    // (exactly as RescaleFacetFluxToTimestep does for facet fluxes). This is
    // what keeps E/M consistent: enthalpy leg = net_Mflux * hCv_donor/mv_donor
    // = vel * hCv_donor * sub_dt — the frozen mass in the denominator that
    // caused the mode-2 desync is gone.
    //
    // Per leg we store the velocity coefficient and the donor SIDE (which node
    // supplies mass and intensive content), fixed by the frozen grad_pot sign
    // at capture. 'mid' donates when grad_pot > 0, matrix side when <= 0.
    enum DonorSide { DONOR_MATRIX, DONOR_MID };
    struct TripleRates {
        Node<dim> *in_n, *out_n, *mid_n;
        // velocity coefficients (grad_pot * K * area), per leg, per phase
        double vel_in_l, vel_out_l, vel_in_v, vel_out_v, vel_in_a, vel_out_a;
        // donor side per leg (frozen from grad_pot sign at capture)
        DonorSide don_in_l, don_out_l, don_in_v, don_out_v, don_in_a, don_out_a;
        // frozen grad_pot (kept for the budget clamp's donor selection, which
        // reads grad_pot_* members — restored from these each beat)
        double gp_in_l, gp_out_l, gp_in_v, gp_out_v, gp_in_a, gp_out_a;
        // stage-start donor holdups (drift diagnostic only)
        double mv_in0, mv_out0, mv_mid0, ml_in0, ml_out0, ml_mid0;
    };
    std::vector<TripleRates> stage_rates_;
    std::map<decltype(mid_node), PhaseDrainRates> node_drain_rates_perpv_;
    bool capture_rates_ = false;

    // Cut-vs-clamp is decided by HOW leakage is entered, not by a file-scope
    // knob (which could be set inconsistently with the run mode). ExchangeSubstep
    // sets beat_mode_=true  -> per-beat BUDGET CLAMP, never cuts.
    // The one-shot pass (legacy) leaves beat_mode_=false -> DETECT + request
    // dt CUT on overdrain (and optional CFL cuts). A beat physically cannot
    // request a cut; a one-shot physically cannot clamp.
    bool beat_mode_ = false;
    InterfaceLedger iface_ledger_;
    bool ledger_active_ = false;   // true only inside ExchangeSubstep beats
    double beat_withheld_ = 0.;    // per-beat standing-queue accumulator
    std::set<Node<dim>*> stale_reported_;   // STALE-ENTRY one-shot-per-node

    void ApplyDrainBudgets();   // extracted budget clamp (shared: apply path + ExchangeSubstep)

    std::string sb_cut_reason_;   // why the last cut fired; exposed via GetCutReason()
    std::string overdrain_detail_; // phase-side list that tripped overdraining, e.g. "liquid-mid,vapor-in"

};
} // end namespace csmp


#endif
