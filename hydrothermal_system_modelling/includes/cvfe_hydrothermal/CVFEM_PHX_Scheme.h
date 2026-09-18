// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_PHX_Scheme_H
#define CVFEM_PHX_Scheme_H

#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <map>

#include "CVFEM_PHX_VariableNames.h"

#include "ThreePhaseTransportPHX.h"

#include "CVFEM_PressureGradientVisitor.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
//#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h" // JK26
#include "CVFEM_NumIntegral_dNT_op_dN_dV.h"
#include "CVFEM_PointSource_rhsop.h"
#include "NumIntegral_NT_lhs_nodal_op_N_dV.h"
#include "NumIntegral_NT_rhs_nodal_op_N_dV.h"
#include "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h"
#include "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h"
#include "NaClH2OPropertiesVisitorPHX.h"
#include "LimitVisitor.h"
#include "PDE_Integrator.h"
#include "WellModelPrototype.h" //Benoit

#include "SplitRegionLeakage.h"

//Benoit add for split boundary:
// LHS version
#include "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure.h"
#include "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature.h"

#include "CVFEM_SolverChoice.h"

namespace csmp {

template<uint32_t dim>
class CVFEM_PHX_Scheme {
public:

    explicit CVFEM_PHX_Scheme(Model<dim> &model_ref,
                                                     bool with_gravity = true,
                                                     bool thermodynamic_density_in_gravity_term = true,
                                                     const std::vector<std::string> wells = {},
                                                     bool with_zinc = false,
                                                     bool with_lithium = false,
                                                     SolverKind solver_kind = SolverKind::PETSc,
                                                     SolverKind well_solver_kind = SolverKind::Eigen,
                                                     const std::map<std::string, WellConfiguration> &well_configs = {});

    ~CVFEM_PHX_Scheme();

    void SetLargestTimeStep(double timestep);   // modifying maximum size of time step at the start of the simulation
    void SetCutLogLabel(const std::string& label) { dt_cut_log_label_ = label; } // names the dt-cut CSV (public: called from the driver)
    void ResetLargestTimeStep(double timestep);  // modifying in again later during the simulation
    void ChangeTimeStepTo(double timestep);   // adjusting timestep
    void InitialFluidPropertiesFromPTX(); // initialize fluid properties from current PTX conditions
    void PrepareTransientCalculations(); // preparation before transient calculations, calculating pressure gradient and updwind nodes from current status
    double Apply(); // main function to apply CVFEM scheme in transeint calculations, returns time step used for calculations

    // modify calculations of temperature-dependent heat capacity of the rock
    void TemperatureDependentHeatCapacityRock(double cpr_min_ext, double t_min_ext,
                                              double cpr_max_ext, double t_max_ext);

    // switch on open top, specifying salinity of inflowing fluid
    void OpenBoundaries(double wt_top);

    void SetBrickWallLimiterTo(bool limit);   // switch brick wall limiter for fluid pressure on or off

    void AddAdvectionVariable(const char *balanced_variable,
                              const char *new_lhs_liquid, const char *new_rhs_liquid,
                              const char *new_lhs_vapor,  const char *new_rhs_vapor);   // add further variables for FV calculations

    void Adjust_CFL_Criterion(double scale_factor, bool take_pore_velocity);  // modifying cfl criterion

    // access function to transient fluxes
    void GetFacetFluxFromInsideNodeToOutsideNode(Element<dim> &e,
                                                 unsigned int facet_idx,
                                                 double &flux_liquid,
                                                 double &flux_vapor,
                                                 double &flux_air);// dormant until ActivateAirPhase()

    // PW May 2016
    void SetEquilibratorConvergenceSpeedUpTo    (bool boost);

    // PW September 2018
    void SetAvoidPressureOscillationsAtVLHTo    (bool avoid_p_VLH);

    //Benoit add
    void WithRockLiquidusSolidus                (bool with_rock_liquidus_solidus = false);
    void AddFluidContributionToHeatCapacity     (bool add_fluid_contribution_to_heat_capacity);
    void AttemptToSurviveFluidPropertiesError   (bool attempt_to_survive_fluid_properties_error_);
    void SetRockHeatCapacity                    (double mini_cp);
    void SetRockCrystallizationCurve            (double nu_coefficient,
                                     double sigma1_coefficient,
                                     double latent_heat_of_fusion,
                                     double b_coefficient,
                                     std::string crystallization_curve_);

    void ResetAndCutTimestep(double new_dt, const char* reason,
                             bool smooth_decrease, bool reset_counter,
                             bool scale_nfvs = true, bool pause = false);

    //Benoit 2020, adding well simulator options
    /**
                @brief Initialize well, positive "delta_T" overrides "init_T", "follow_reservoir_T" overrides "delta_T" and "init_T"
            */
    /** Initialises one well. Geometry, completions and the initial temperature
        come from the WellConfiguration given to this scheme's constructor;
        init_bulk_smf is per well. */
    /** Initialises one well entirely from its WellConfiguration, given to this
        scheme's constructor keyed by well name. */
    void Initialize_well(const std::string &well_name);

    /** Requests a wellhead pressure and puts the well under wellhead-pressure
        control. Applied immediately unless the well's configuration sets a ramp
        rate. (The old `forced` flag is gone: it selected a hardcoded
        0.001 bar/min ramp that was invisible from the call site.) */
    void Set_well_top_pressure(const std::string &well_name,
                               const double &pressure);

    void Set_well_top_temperature(const std::string &well_name,
                                  const double &temperature);

    void Set_well_water_table(const std::string &well_name,
                              const bool &toggle_water_table_displacement,
                              const double &initial_water_table_depth);

    /** Holds a target MASS rate on this well: POSITIVE produces, NEGATIVE
        injects. Call again at any time to change it. `trigger == false` returns
        the well to wellhead-pressure control.

        (A `volumetric` argument used to be accepted and silently ignored — there
        is no volumetric control — so it has been removed.) */
    void Set_target_rate(const std::string &well_name,
                         const bool &trigger,
                         const double &ext_target_rate);

    /** Puts a well into top-injection mode at the given MASS rate [kg/s], which
        must be POSITIVE. The well is then filled from the top every step, with
        fluid at the temperature and salinity in its WellConfiguration.

        There is no "off" value: to stop top injection, put the well under another
        control mode with Set_well_top_pressure() or Set_target_rate().
        A non-positive rate throws. (It used to mean "off", which made a disabled
        feature indistinguishable from a typo — and, once modes existed, silently
        reset whichever mode the driver had chosen.) */
    void Set_injection_rate(const std::string &well_name, const double &ext_injection_rate);
    /** Switches well calculations on or off for ALL wells at once.

        Per-well control would mean moving the flag into WellModelPrototype and
        having WellReservoirExplicitCoupling skip a well whose flag is false —
        straightforward, but nothing needs it yet, so it stays global. */
    void IncludeWellCalculations(const bool &ext_with_well);
    void OutputWellResults(const std::string &well_name, const std::string &filename);
    void PAUSE();
    void SetModelTime(double time);

    //Benoit add air phase
    void ActivateAir(
        double humidity_air_top = 0.,    // Constant relative humidity [0–1]
        bool treat_air_and_vapor_as_mixture = false,
        bool treat_air_and_vapor_as_hydraulic_mixture = false,
        double rain_mm_per_year = 0.);  //Rain infiltration rate

    //Benoit add split boundary options
    void ActivateSplitRegionCoupling(
        bool open_space,
        double fault_perm_anisotropy,
        bool with_air = false,
        bool with_tracer = false,
        bool with_lithium = false,
        bool with_magmatic_fluids = false,
        bool with_gold = false);

    std::vector<double> injection_rates;


    // Used for 2D Benchmarks
    void ActivatePointInjectionBenchmarks(const char* val);

private:

    // sequence matters
    Model<dim> &model;
    const PropertyDatabase<dim> &p_ref;

    CVFEM_PHX_VariableNames
        names;

    bool
        with_gravity_,
        thermodynamic_density_in_gravity_term_,
        verbose,
        open_top,
        brick_wall_limiter,
        well_coupling_rejected,
        with_well_calculations,
        with_lithium_,
        with_air_,
        smooth_timestep_decrease,
        with_split_region_,
        open_space_,
        output_CVFEM_timing;


    double fault_perm_anisotropy_;

    uint32_t counter_smooth_timestep_decrease;

    CVFEM_PressureGradientVisitor<dim>
        pres_grad;

    ExplicitFiniteVolumeTransportPHX<dim>
        fv_transport_vapor,
        fv_transport_liquid,
        fv_transport_air;// dormant until ActivateAirPhase()

    UpwindControlVisitor<dim>
        upwind_control;

    ThreePhaseTransportPHX<dim>
        transport;

    // One bundle per system; backend chosen at construction (see CVFEM_SolverChoice.h).
    SolverBundle P_bundle, T_bundle;


    PDE_Integrator<dim, Element, SparseMatrix> // Benoit add SparseMatrix
        P_FE,
        T_FE;

    NumIntegral_NT_lhs_nodal_op_N_dV<dim, Element>
        capacitance_lhs;

    CVFEM_NumIntegral_dNT_op_dN_dV<dim, Element>
        conductance;

    NumIntegral_NT_rhs_nodal_op_N_dV<dim, Element>
        capacitance_rhs;

    CVFEM_PointSource_rhsop<dim, Element>
        heat_bottom;

    NumIntegral_NT_lhs_nodal_op_N_dV<dim, Element>
        capacitance_lhs_p;

    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element>
        conductance_p_upwind_liquid,
        conductance_p_upwind_vapor,
        conductance_p_upwind_air;// dormant until ActivateAirPhase()

    NumIntegral_NT_rhs_nodal_op_N_dV<dim, Element>
        capacitance_rhs_p;

    CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>
        grav_liq,
        grav_vap,
        grav_air;// dormant until ActivateAirPhase()

    CVFEM_PointSource_rhsop<dim, Element>
        source_p,
        source_p2;

    //Benoit add for splitboundary
    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_pressure<dim>
        iface_transferLHS_p;// no air until ActivateAirPhase()

    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>
        iface_transferLHS_t;

    CVFEM_PointSource_rhsop<dim, Element> iface_grav;

    NaClH2OPropertiesVisitorPHX<dim>
        equilibrator_properties;// no air until ActivateAirPhase()


    std::vector<PropertyHandle<dim>* >
        reset_properties;

    double
        minimum_pressure;

    LimitVisitor<dim>
        pressure_limiter_transport,
        pressure_limiter_fluid;

    std::vector<std::unique_ptr<WellModelPrototype<dim>>>
        well_models; //vector of well models

    std::vector<std::unique_ptr<WellModelPrototype<dim>>>InstantiateWellModels(Model<dim> &model, const std::vector<std::string> &wells, SolverKind well_solver_kind, const std::map<std::string, WellConfiguration> &well_configs);

    SplitRegionLeakage<dim> *Split_region_leakage = nullptr;

    double
        cfl_dt,
        largest_timestep,
        current_dt,
        timestep_increment_factor,
        control_dt,
        old_dt,
        min_value,
        max_value,

        cfl_min_l, cfl_max_l,
        cfl_min_v, cfl_max_v,
        cfl_min_a, cfl_max_a,

        well_source,
        total_tracer_mass_injected,
        injection_rate,
        model_time_;

    void AdvanceTransientVariables(); // book keeping of variables for transient calculations
    void ResetVariables(); // reset variables for transient pressure calculations
    void AdvectionDiffusionLoops(); // outer loop including advection and pressure diffion until mass-based time step criterion is met
    void PressureLoop(); // inner loop including pressure diffion until cfl-based time step criterion is met
    void FluidRockEquilibration(); // thermal quilibration between fluid and rock

    /** One-way (explicit) well -> reservoir exchange for this step. Solves each
        well, turns the result into reservoir source terms */
    void WellReservoirExplicitCoupling();

    /** Source terms for ONE well, in two phases.
        final_pass == false : compute them, and reject the step if any is
                              impossibly large.
        final_pass == true  : store them and accumulate the exchange totals.
        (Formerly Well_Reservoir_Iterative_Coupling(bool pressure_converged) —
        there is no iteration, and the flag was never a convergence result.) */
    void ComputeWellSourceTerms(bool final_pass, const std::string well_name);
    void PointInjection();
    void PointInjectionBenchmarks(std::string test);
    void nfvs_scaling(double scaling_factor);

    // ── Upwind freezing within a timestep (see the block above Apply) ──────
    // Latches the DIRECTION pass only: it runs once per timestep and its donor
    // map is reused by every later solve of that step (dt-cut retries,
    // split-boundary retries). Cleared at the top of Apply. Needed because that
    // pass runs BEFORE the solve, so after a ResetVariables it would decide
    // donors from STENCIL_DATA velocities left by the REJECTED solve — internal
    // state no property reset can rewind — and a donor flip changes the matrix
    // discretely. The VELOCITY pass is deliberately NOT latched: it runs after
    // UpdateProjection, so its inputs are always current, and freezing it made
    // cfl_dt stale on retries.
    bool upwind_done_this_step_ = false;

    // ── Timestep-cut logging ──────────────────────────────────────────
    std::ofstream dt_cut_log_;
    std::string   dt_cut_log_label_;
    std::string   dt_cut_log_stamp_;
    void LogTimestepCut(const char* reason, double new_dt);

    // Used for 2D Benchmarks
    std::string benchmark_;
    bool with_point_injection_ = false;

};

} // end namespace csmp

#endif
