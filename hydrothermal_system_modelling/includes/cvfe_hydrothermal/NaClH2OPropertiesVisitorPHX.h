// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NACLH2OPROPERTIESVISITORPHX_H
#define NACLH2OPROPERTIESVISITORPHX_H

#include "CSMP_definitions.h"
#include "Model.h"
#include "Boundary.h"
#include "Visitor.h"
#include "ErrorHandler.h"

#include "H2ONaClThermalEquilibrator.h"
#include "Rock.h"
#include "Air.h"

namespace csmp
{

template<uint32_t dim>
class NaClH2OPropertiesVisitorPHX : public Visitor<dim>
{
public:

    NaClH2OPropertiesVisitorPHX( Model<dim>& );
    ~NaClH2OPropertiesVisitorPHX();

    // ── Visitor interface ─────────────────────────────────────────────────────
    virtual void Visit( Node<dim>* n );
    virtual void Visit( Model<dim>* n );
    virtual void Visit( Region<dim>* n );

    // ── Setup ─────────────────────────────────────────────────────────────────
    void SetVerbose(bool verbose) { verbose_output = verbose; }
    void SetTimeIncrement( double time_increment );
    void InitialPropertiesFromPTX( Model<dim>& );

    // ── Open boundaries ───────────────────────────────────────────────────────
    void OpenBoundaries( double ref_wt );                              // variable p-t top

    // ── Air phase ─────────────────────────────────────────────────────────────
    void ActivateAir(
        double humidity                       = 0.,     // relative humidity [0–1]
        bool   treat_air_and_vapor_as_mixture = false,
        bool   treat_air_and_vapor_as_hydraulic_mixture = false,
        double rain_mm_per_year               = 0. );   // liquid recharge at TOP [mm/yr]

    // ── Brooks-Corey relative permeability ────────────────────────────────
    // Replaces linear kr_liquid and kr_gas with a self-consistent Brooks-
    // Corey pair derived from a pore-size distribution model. Strongly
    // suppresses gas mobility in wet conditions (preventing air percolation
    // through saturated zones) and suppresses liquid mobility in dry
    // conditions (slowing drainage front propagation). This is the
    // standard macroscopic capillary retention model.
    //
    // Parameters:
    //   enable : true → Brooks-Corey, false → legacy linear
    //   lambda : Brooks-Corey pore-size distribution index [-]
    //            ~0.5  fine-grained, poorly sorted (shale, clay)
    //            ~1.0  fractured crystalline rock (reasonable default)
    //            ~2.0  coarse, well-sorted (sand, gravel)
    //
    // Both curves use the same lambda and the same residual saturations
    // (res_sl, res_sv, res_sa) already defined as member variables.
    void SetBrooksCoreyKr( bool b, double lambda = 2. );

    // ── Physics options ───────────────────────────────────────────────────────
    void SetThermodynamicDensityInGravityTermTo( bool b );
    void SetThermodynamicDensityInCompressibilityTermTo( bool b );
    void SetAvoidPressureOscillationsAtVLHTo( bool b );
    void SetConvergenceSpeedUpTo( bool b );
    void AddFluidContributionToHeatCapacity( bool b );
    void AttemptToSurviveFluidPropertiesError( bool b );

    // ── Rock model ────────────────────────────────────────────────────────────
    void WithRockLiquidusSolidus( bool b );
    void SetRockHeatCapacity( double mini_cp );
    void SetRockCrystallizationCurve(
        double nu, double sigma1, double latent_heat, double b_coeff, std::string crystallization_curve_ );


    //----
    /** One tracked extreme: the worst value of a monitored quantity this pass,
        plus the state at the node that produced it. `p` and `dp` are what the
        scheme needs to invert the criterion into a timestep scale. */
    struct DpExtreme {
        double value = 0.;   // the monitored quantity at its worst
        double p     = 0.;   // fluid pressure at that node
        double dp     = 0.;  // expected_dp at that node
        size_t node  = 0;
    };

private:

    // =====================================================================
    //  INFRASTRUCTURE
    // =====================================================================

    PropertyDatabase<dim>&         pref;
    csmp::MeshManager<dim>&        pmesh;
    ErrorHandler&                  csmp_error;

    // =====================================================================
    //  CONFIGURATION FLAGS
    // =====================================================================

    bool verbose_output                               = false;
    bool bogus_variables                              = false;
    bool pure_halite                                  = false;
    bool fixed_temperature                            = false;

    bool open_boundaries                              = false;

    bool thermodymamic_density_in_gravity_term         = true;
    bool thermodynamic_density_in_compressibility_term = true;

    bool avoid_pressure_oscillations_at_VLH            = true;

    bool with_rock_liquidus_solidus                   = false;
    bool add_fluid_contribution_to_heat_capacity      = false;
    bool thermal_eq_succeeded                         = false;
    bool verbose_equilibration                        = false;

    double expected_dp     = 0.;

    // ── Air configuration ─────────────────────────────────────────────────
    bool   with_air                                   = false;
    bool   full_gas_mixture                           = false;
    bool   hydraulic_gas_mixture                      = false;

    double ref_humidity                               = 0.;   // relative humidity [0–1]
    double rain_mm_per_year_top                       = 0.;   // liquid recharge [mm/yr]

    // Minimum H2O-NaCl pore fraction — always reserved for the Equilibrator.
    // Prevents phi_ → 0 and hdummy blowing up when air dominates the pore space.
    // Used in ComputeAirSaturationAndDensity (sa_max cap) and BoundaryFlow
    // (direct excess air removal). Must be consistent between both.
    double min_h2o_nacl_fraction = 0.02;
    //double sa_maximum = 0.;

    // ── Brooks-Corey relative permeability ────────────────────────────────
    // When enabled, both kr_liquid and kr_gas are computed from a Brooks-
    // Corey pore-size distribution model instead of the legacy linear forms.
    // The two curves form a self-consistent pair derived from the same
    // pore-network assumption:
    //
    //   kr_l_BC(Se)   = Se^((2 + 3λ)/λ)
    //   kr_gas_BC(Se) = (1 - Se)^2 · (1 - Se^((2+λ)/λ))
    //
    //   Se = (sl - res_sl·(1-sh)) / ((1-sh)·(1 - res_sl - res_sv - res_sa))
    //
    // Captures the dominant macroscopic effect of capillary retention at
    // the continuum scale: the gas phase cannot efficiently invade wet
    // pores, and the liquid phase loses mobility rapidly as it approaches
    // residual saturation. This is the standard model for unsaturated flow
    // in simulators that do not explicitly solve pc(sl), and is what every
    // major groundwater/geothermal code uses (TOUGH, PFLOTRAN, HYDRUS).
    //
    // The gas kr is applied to the combined gas phase (vapor + air) and
    // then split between the two components by volume fraction in the
    // existing kr_gas → krv, kra allocation — consistent with the fact
    // that vapor and air share the gas pore network at the pore scale.
    //
    // Set via SetBrooksCoreyKr(true, lambda).
    bool   use_brooks_corey_kr                        = false;  // opt-in
    double lambda_bc                                  = 2.;    // pore-size index default value

    // ── Boundary inflow reference values ──────────────────────────────────
    double ref_frac_air                               = 0.;   // [kg_air/kg_total] diagnostic
    double ref_frac_vapor_in_air                      = 0.;   // [kg_vap/kg_gas]
    double ref_spec_h_liquid                          = 0.;   // [J/kg]
    double ref_spec_h_vapor_in_air                    = 0.;   // [J/kg]
    double ref_spec_h_air                             = 0.;   // [J/kg]
    double ref_sal                                    = 0.;   // boundary salinity [mass frac]
    int boundary_flow_trial                           = 0;
    double rain_mass_deposited                        = 0.;

    // =====================================================================
    //  EQUILIBRATOR WORKING VARIABLES  (passed by reference)
    // =====================================================================

    double m_rock_       = 0.;
    double cp_rock_      = 0.;
    double rho_rock_     = 0.;
    double m_air_        = 0.;
    double cp_air_       = 0.;
    double phi_          = 0.;
    double m_fluid_      = 0.;
    double wt_           = 0.;
    double tp_           = 0.;
    double p_current_    = 0.;
    double H_current_    = 0.;
    double H_previous_   = 0.;
    double t_fixed       = 0.;
    double t_diffusion_  = 0.;

    // =====================================================================
    //  PER-NODE WORKING VARIABLES
    // =====================================================================

    // ── Increments ────────────────────────────────────────────────────────
    double dh_rock_diff_  = 0.;
    double dh_fluid_diff_ = 0.;
    double dh_air_diff_   = 0.;
    double dhCl_          = 0.;
    double dhCv_          = 0.;
    double dhCa_          = 0.;
    double dP             = 0.;
    double dxCl_          = 0.;
    double dxCv_          = 0.;
    double dx_diff_       = 0.;
    double dml_           = 0.;
    double dmv_           = 0.;
    double dma_           = 0.;

    // ── Scalars ───────────────────────────────────────────────────────────
    double x_              = 0.;
    double dt_             = 0.;
    double h_fluid_        = 2.086e6;  // init guess for EOS
    double t_              = 0.;
    double tl_             = -1.;      // rock liquidus
    double ts_             = -1.;      // rock solidus
    double p_              = 0.;
    double previous_state  = 0.;
    //double pore_volume     = 0.;
    //double rock_volume     = 0.;
    double volume_factor_LHS = 0.;
    double volume_factor_RHS = 0.;

    double hrock_prev      = 0.;
    double hrock_curr      = 0.;

    static constexpr double kelvin = 273.15;

    // ── Residual saturations ──────────────────────────────────────────────HARDCODED HERE, COULD BE MOVED TO SETTER FUNCTION
    double res_sl = 0.3;   // liquid (hardcoded)
    double res_sv = 0.;    // vapor
    double res_sa = 0.;    // air (fully mobile)

    double sat_min = 1.e-10;  // below this a phase is numerically empty

    // ── Relative permeabilities ───────────────────────────────────────────
    double krl = 0.;
    double krv = 0.;
    double kra = 0.;

    // ── Phase state ───────────────────────────────────────────────────────
    int old_state     = 0;
    int current_state = 0;

    // =====================================================================
    //  SCALAR VARIABLES  (read from / stored to CSMP mesh)
    // =====================================================================

    // ── Rock liquidus/solidus ─────────────────────────────────────────────
    ScalarVariable tl, ts;

    // ── Primary state ─────────────────────────────────────────────────────
    ScalarVariable t, tp, p, pp;

    // ── Saturations ───────────────────────────────────────────────────────
    ScalarVariable sl, sv, sa, sh;
    ScalarVariable slp, svp, sap;

    // ── Mass concentrations [kg/m³ pore] ──────────────────────────────────
    ScalarVariable mt, mtp;
    ScalarVariable ml, mv, ma;
    ScalarVariable mlp, mvp, map;
    ScalarVariable mh;

    // ── Phase densities [kg/m³] ───────────────────────────────────────────
    ScalarVariable rl, rv, ra, rh;
    ScalarVariable rlp, rvp, rap;

    // ── Viscosities [Pa·s] ────────────────────────────────────────────────
    ScalarVariable mul, muv, mua;
    ScalarVariable mulp, muvp, muap;

    // ── Specific enthalpies [J/kg] ────────────────────────────────────────
    ScalarVariable hf;
    ScalarVariable hl, hv, ha, hh;
    ScalarVariable hlp, hvp, hap;
    ScalarVariable cpf;//heat capacity fluid

    // ── Volumetric enthalpies [J/m³] ──────────────────────────────────────
    ScalarVariable hVl, hVv, hVa, hVh;
    ScalarVariable hCl, hCv, hCa, hCh;
    ScalarVariable hClp, hCvp, hCap;

    // ── Previous total enthalpy [J/m³] ──────────────────────────────────────
    ScalarVariable Htp;

    // ── Compressibility ───────────────────────────────────────────────────
    ScalarVariable beta, CT;

    // ── Rock / pore ───────────────────────────────────────────────────────
    ScalarVariable cpr, phi, rr, cpv;
    ScalarVariable beta_rock, ncp;

    // ── Air thermodynamic ─────────────────────────────────────────────────
    ScalarVariable cpa, beta_air;

    // ── Fluid sources ─────────────────────────────────────────────────────
    ScalarVariable src_h, src_rate, src_wt, nQ;

    // ── Transport densities ───────────────────────────────────────────────
    ScalarVariable rl_transport, rv_transport, ra_transport;
    ScalarVariable vol_fac, rho_bulk;

    // ── Salt ──────────────────────────────────────────────────────────────
    ScalarVariable wt;
    ScalarVariable xf, xl, xv, xh;
    ScalarVariable xVl, xVv, xVh;
    ScalarVariable xCl, xCv, xCh;
    ScalarVariable xClp, xCvp, xCf;
    ScalarVariable ms, msp;
    ScalarVariable xlp, xvp, hhp;

    // ── State tracking ────────────────────────────────────────────────────
    ScalarVariable state, state_p;
    ScalarVariable dangerous_phase_change;
    ScalarVariable after_phasechange_counter;

    // ── Phase mobilities ──────────────────────────────────────────────────
    ScalarVariable mml,  mmv,  mma;
    ScalarVariable mmld, mmvd, mmad;
    ScalarVariable eml,  emv,  ema;
    ScalarVariable emld, emvd, emad;
    ScalarVariable xml, xmv;
    ScalarVariable rvl, rvv, rva;

    // ── Boundary flow ─────────────────────────────────────────────────────
    ScalarVariable bfm, bfe, bfs, bfa;

    // ── Boundary reference (per-node) ─────────────────────────────────────
    ScalarVariable ref_enthalpy_liquid_top;
    ScalarVariable ref_enthalpy_vapor_in_air_top;
    ScalarVariable ref_enthalpy_air_top;
    ScalarVariable ref_frac_vapor_in_air_top;
    ScalarVariable ref_frac_air_top;
    ScalarVariable rain_recharge_top;

    // ── Gas mixture ───────────────────────────────────────────
    ScalarVariable xa_mix, xv_mix;
    ScalarVariable p_partial_air, p_partial_vapor;
    ScalarVariable rho_gas_mix, mu_gas_mix;
    ScalarVariable h_gas_mix, beta_gas_mix;

    // =====================================================================
    //  INDEX KEYS
    // =====================================================================

    csmp::Index tl_key, ts_key;
    csmp::Index cpr_key, rr_key, phi_key, cpv_key, beta_rock_key, ncp_key;
    csmp::Index t_key, tp_key, p_key, pp_key;
    csmp::Index state_key, state_p_key, apc_key, dpc_key;
    csmp::Index sl_key, sv_key, sa_key, sh_key;
    csmp::Index mt_key, mtp_key;
    csmp::Index ml_key, mv_key, ma_key, mlp_key, mvp_key, map_key, mh_key;
    csmp::Index rl_key, rv_key, ra_key, rh_key;
    csmp::Index rl_transport_key, rv_transport_key, ra_transport_key, rho_bulk_key;
    csmp::Index mul_key, muv_key, mua_key;
    csmp::Index hf_key, hl_key, hv_key, ha_key, hh_key;
    csmp::Index hCl_key, hClp_key, hCv_key, hCvp_key, hCh_key;
    csmp::Index hCa_key, hCap_key;
    csmp::Index hVl_key, hVv_key, hVa_key, hVh_key;
    csmp::Index cpf_key, Htp_key;
    csmp::Index beta_key, CT_key;
    csmp::Index beta_air_key, cpa_key;
    csmp::Index nQ_key, vol_fac_key;
    csmp::Index wt_key;
    csmp::Index xf_key, xl_key, xv_key, xh_key;
    csmp::Index xVl_key, xVv_key, xVh_key;
    csmp::Index xCl_key, xCv_key, xCh_key, xClp_key, xCvp_key, xCf_key;
    csmp::Index msp_key;
    csmp::Index src_h_key, src_rate_key, src_wt_key;
    csmp::Index mml_key, mmv_key, mma_key, mmld_key, mmvd_key, mmad_key;
    csmp::Index eml_key, emv_key, ema_key, emld_key, emvd_key, emad_key;
    csmp::Index xml_key, xmv_key;
    csmp::Index rvl_key, rvv_key, rva_key;
    csmp::Index bfm_key, bfe_key, bfs_key, bfa_key;
    csmp::Index ref_enthalpy_liquid_top_key, ref_enthalpy_vapor_in_air_top_key, ref_enthalpy_air_top_key;
    csmp::Index ref_frac_vapor_in_air_top_key, ref_frac_air_top_key;
    csmp::Index rain_recharge_top_key;
    csmp::Index xa_mix_key, xv_mix_key;
    csmp::Index p_partial_air_key, p_partial_vapor_key;
    csmp::Index rho_gas_mix_key, mu_gas_mix_key, h_gas_mix_key, beta_gas_mix_key;

    // =====================================================================
    //  PHYSICS / EOS OBJECTS
    // =====================================================================

    Rock                        rock;
    Air                         air;
    H2ONaClThermalEquilibrator  equilibrator;
    H2ONaClFluidProperties      fluid;
    Fluidproperties             Liquid, Vapor, Bulk, Salt;

    // =====================================================================
    //  PRIVATE MEMBER FUNCTIONS
    // =====================================================================

    double EffectiveLiquidSaturationHalitePresent(
        double sat_liquid, double sat_vapor, double sat_air = 0. ) const;
    double RelativePermeabilityLiquid(
        double sat_liquid, double sat_vapor, double sat_air = 0. ) const;
    double RelativePermeabilityGas(
        double sat_liquid, double sat_vapor, double sat_air = 0. ) const;

    void ReadAllVariables( Node<dim>* n );
    void CheckForOutOfRange( Node<dim>* n );
    void CheckBoundaryFlags( Node<dim>* n );
    void CalculateAbsoluteVariables( Node<dim>* n );
    void UpdateEquilibratorVariables( Node<dim>* n );
    void Equilibrate( Node<dim>* n );
    void UpdateCSMPVariables( Node<dim>* n );
    void VolumeFactorComputations( Node<dim>* n );
    void PrepareVariablesForStorage();
    void StorePropertiesAndFlags( Node<dim>* n );
    void StoreInitialPropertiesAndFlags( Node<dim>* n );
    void SetPrecipitationRecharge( Model<dim>& model );

    double BoundaryFlow();
    void   BoundaryIteration();

    double ComputeMaxAirSaturation() const;
    double ComputeTotalCompressibility();
    void ComputeAirSaturationAndDensity(
        double ma_, double t_, double p_, double p_air_partial);

    //   Authoritative mixture-property computation called by both
    //   Pass 1 (pre-Equil: muv_in=beta_v_in=hv_in=0) and
    //   Pass 2 (post-Equil: full Vapor properties supplied).
    void ComputeGasMixtureProperties(
        double t,               double p,
        double sv_in,           double sa_in,
        double rv_in,
        double muv_in,          double beta_v_in,    double hv_in,
        double p_partial_air_in, double p_partial_vapor_in,
        double& xa_out,         double& xv_out,
        double& ra_out,         double& mu_gas_mix_out,
        double& rho_gas_mix_out, double& h_gas_mix_out,
        double& beta_gas_mix_out);

    void ComputeGasMobilities(
        double& mob_v, double& mob_a,
        double krv, double kra, double rvp, double rap,
        double mu_v, double mu_a,
        double sl, double sv, double sa) const;

    double TwoPhasePureWaterCompressibility( double cpl, double cpv );

    void DumpNodeState(const char* call_site);
    void PAUSE();

    // ── Physical-state validator ───────────────────────────────────────────
    // Call at pipeline checkpoints with a stage label.  All violations are
    // written to cerr; if any are found the full equilibrator state is dumped.
    // Stage labels used in Visit():
    //   "ReadAllVariables"   — raw CSMP inputs, before any computation
    //   "PostEquilibrate"    — after EOS / BoundaryIteration + UpdateCSMPVariables
    //   "PostVolumeFactors"  — after VolumeFactorComputations
    //   "PostPrepare"        — after PrepareVariablesForStorage (final state)
    void CheckPhysicalState( Node<dim>* n, const char* stage );
};


// =========================================================================
//  INLINE IMPLEMENTATIONS
// =========================================================================

template<uint32_t dim>
inline double NaClH2OPropertiesVisitorPHX<dim>::EffectiveLiquidSaturationHalitePresent(
    double sat_liquid, double sat_vapor, double sat_air ) const
{
    const double sat_total = sat_liquid + sat_vapor + sat_air;
    const double my_res_sl = res_sl * sat_total;
    const double my_res_sv = res_sv * sat_total;
    const double my_res_sa = res_sa * sat_total;

    const double denom = sat_total - my_res_sl - my_res_sv - my_res_sa;
    if (denom <= 0.) return 0.;
    const double seff = (sat_liquid - my_res_sl) / denom;
    if (seff > 1.) return 1.;
    if (seff < 0.) return 0.;
    return seff;
}

// template<uint32_t dim>
// inline double NaClH2OPropertiesVisitorPHX<dim>::RelativePermeabilityLiquid(
//     double sat_liquid, double sat_vapor, double sat_air ) const
// {
//     const double sat_total = sat_liquid + sat_vapor + sat_air;
//     const double my_res_sl = res_sl * sat_total;
//     const double my_res_sv = res_sv * sat_total;
//     const double my_res_sa = res_sa * sat_total;

//     const double denom = sat_total - my_res_sl - my_res_sv - my_res_sa;
//     if (denom <= 0.) return 0.;
//     const double seff = sat_total * (sat_liquid - my_res_sl) / denom;
//     if (seff > sat_total) return sat_total;
//     if (seff < 0.) return 0.;
//     return seff;
// }

template<uint32_t dim>
inline double NaClH2OPropertiesVisitorPHX<dim>::RelativePermeabilityLiquid(
    double sat_liquid, double sat_vapor, double sat_air ) const
{
    // =========================================================================
    // Liquid-phase relative permeability.
    //
    // Returns kr_liquid scaled by sat_total = sat_liquid + sat_vapor + sat_air
    // = (1 - sh), consistent with the convention of RelativePermeabilityGas.
    // Callers use the returned value directly without additional (1-sh) scaling.
    //
    // LEGACY LINEAR FORM (use_brooks_corey_kr = false):
    //   kr_l = sat_total * Se
    //
    // BROOKS-COREY FORM (use_brooks_corey_kr = true):
    //   kr_l_BC = Se^((2 + 3λ)/λ)
    //   returned value = sat_total * kr_l_BC
    //
    //   Limiting behaviour:
    //     Se = 1 (sl at full pore):   kr_l = sat_total (full mobility)
    //     Se = 0 (sl at residual):    kr_l = 0         (no flow)
    //
    //   Much steeper than linear near residual — captures the physical fact
    //   that as liquid drains, the remaining water occupies smaller and more
    //   tortuous pores with rapidly decreasing connectivity.
    //
    //   At λ = 1:  kr_l(Se=0.8) ≈ 0.33   vs linear 0.80
    //              kr_l(Se=0.5) ≈ 0.03   vs linear 0.50
    //              kr_l(Se=0.3) ≈ 0.002  vs linear 0.30
    // =========================================================================
    const double sat_total = sat_liquid + sat_vapor + sat_air;
    if (sat_total <= 0.) return 0.;

    const double my_res_sl = res_sl * sat_total;
    const double my_res_sv = res_sv * sat_total;
    const double my_res_sa = res_sa * sat_total;
    const double denom     = sat_total - my_res_sl - my_res_sv - my_res_sa;
    if (denom <= 0.) return 0.;

    const double Se = std::max(0., std::min(1.,
                                            (sat_liquid - my_res_sl) / denom));

    if (!use_brooks_corey_kr)
    {
        // Legacy linear form: kr_l = sat_total * Se
        if (Se >= 1.) return sat_total;
        if (Se <= 0.) return 0.;
        return sat_total * Se;
    }

    // Brooks-Corey liquid: kr_l = Se^((2 + 3*lambda)/lambda)
    if (Se >= 1.) return sat_total;
    if (Se <= 0.) return 0.;

    const double lambda_safe = std::max(lambda_bc, 1.e-6);
    const double exponent    = (2.0 + 3.0 * lambda_safe) / lambda_safe;
    return sat_total * std::pow(Se, exponent);
}

template<uint32_t dim>
inline double NaClH2OPropertiesVisitorPHX<dim>::RelativePermeabilityGas(
    double sat_liquid, double sat_vapor, double sat_air ) const
{
    // =========================================================================
    // Gas-phase (vapor + air combined) relative permeability.
    //
    // Returns kr_gas scaled by sat_total = (1 - sh), matching the convention
    // of RelativePermeabilityLiquid. The returned value represents the total
    // gas-phase mobility; the caller splits it between vapor and air by
    // volume fraction (krv = kr_gas·sv/s_gas, kra = kr_gas·sa/s_gas) in the
    // existing kr split machinery. This is consistent with vapor and air
    // sharing the gas pore network at the pore scale.
    //
    // LEGACY LINEAR FORM (use_brooks_corey_kr = false):
    //   kr_gas = sat_total * (1 - Se) = 1 - sh - kr_l_linear
    //   (bit-identical to the original code when kr_l is also linear)
    //
    // BROOKS-COREY FORM (use_brooks_corey_kr = true):
    //   kr_gas_BC = (1 - Se)^2 * (1 - Se^((2+λ)/λ))
    //   returned value = sat_total * kr_gas_BC
    //
    //   Limiting behaviour:
    //     Se = 1 (sl at full pore):  kr_gas = 0         (no invasion)
    //     Se = 0 (sl at residual):   kr_gas = sat_total (full mobility)
    //
    //   At λ = 1:  kr_gas(Se=0.5) ≈ 0.22  vs linear 0.50
    //              kr_gas(Se=0.7) ≈ 0.065 vs linear 0.30
    //              kr_gas(Se=0.9) ≈ 0.002 vs linear 0.10
    //
    //   The factor-of-50-to-100 reduction in wet conditions is the dominant
    //   fix for air percolating too freely through saturated interior zones.
    //
    // KR SUM NOTE:
    //   Linear:       kr_l + kr_gas = sat_total exactly.
    //   Brooks-Corey: kr_l + kr_gas < sat_total in the intermediate regime.
    //   The "deficit" represents capillary trapping where neither phase flows
    //   efficiently — physical, not a bug. Downstream code does not rely on
    //   the sum equalling sat_total: each phase's mobility is computed
    //   independently and fractional allocation is done by normalisation.
    // =========================================================================
    const double sat_total = sat_liquid + sat_vapor + sat_air;
    if (sat_total <= 0.) return 0.;

    const double my_res_sl = res_sl * sat_total;
    const double my_res_sv = res_sv * sat_total;
    const double my_res_sa = res_sa * sat_total;
    const double denom     = sat_total - my_res_sl - my_res_sv - my_res_sa;
    if (denom <= 0.) return sat_total;   // degenerate: treat all space as gas

    const double Se = std::max(0., std::min(1.,
                                            (sat_liquid - my_res_sl) / denom));

    if (!use_brooks_corey_kr)
    {
        // Legacy linear form: kr_gas = sat_total * (1 - Se)
        return sat_total * (1.0 - Se);
    }

    // Brooks-Corey gas: kr_gas = (1 - Se)^2 * (1 - Se^((2+λ)/λ))
    if (Se >= 1.) return 0.;
    if (Se <= 0.) return sat_total;

    const double one_minus_Se = 1.0 - Se;
    const double lambda_safe  = std::max(lambda_bc, 1.e-6);
    const double exponent     = (2.0 + lambda_safe) / lambda_safe;
    const double kr_gas_bc    = one_minus_Se * one_minus_Se
                             * (1.0 - std::pow(Se, exponent));

    return sat_total * kr_gas_bc;
}

} // namespace csmp

#endif