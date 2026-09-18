// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef WELL_MODEL_PROT_H
#define WELL_MODEL_PROT_H

#include "CSMP_definitions.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "H2ONaClThermalEquilibrator_Fluid_Only.h"
#include "Model.h"
#include "PropertyDatabase.h"
#include "Region.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
// #include "Brine.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "EigenSolver.h"
#include "CVFEM_SolverChoice.h"
#include "VTU_Interface.h"
#include "compareFloats.h"
#include <string>
#include <vector>
#include <optional>

namespace csmp {

/**
        @author BLC
        @date 2020
    */

// 1D well N-R solver
// Only works properly within a 3D mesh


/** One perforated interval of the well, in absolute depth (negative downward).

    A well can have any number of these; a node is a completion if its depth lies
    inside any interval. Only completion nodes exchange with the reservoir
    (VO_well_completion gates the exchange), so this list IS the well's connection
    topology, not a convenience.
*/
struct CompletionInterval {
    double top_depth       = 0.;    ///< m, upper (shallower) limit
    double bottom_depth    = 0.;    ///< m, lower (deeper) limit

    /** Skin for this interval only. Left unset (the default), the interval uses
        WellConfiguration::default_skin — so writing `{ -1900., -10000. }` does
        NOT silently zero the well's skin. Set it explicitly, including to 0, to
        override the well default for this interval. */
    std::optional<double> skin;
};

/** How a well is driven. Exactly one at a time — the three used to be independent
    flags (a wellhead pressure, a target-rate trigger, an injection rate with -1
    meaning "off"), so nothing stopped two being active and quietly fighting over
    top_pressure. As an enum they are mutually exclusive by construction.

    Switch mode at any time during a run; the well reports each change. */
enum class WellControl {
    /** Run on the wellhead pressure. The default, and what a well does after
        Initialize_Well until told otherwise. */
    WellheadPressure,

    /** Hold a target mass rate. A controller adjusts the wellhead pressure each
        step until the rate is matched (to within ~0.3 kg/s). Sign convention:
        POSITIVE produces, NEGATIVE injects. Change the target whenever you like;
        that is the normal way to drive a doublet through a schedule. */
    TargetRate,

    /** Fill the well from the top at a given mass rate, with fluid at a
        configured temperature and salinity. Intended for filling a well and
        letting it find a steady state, typically together with water-table
        displacement. NOT compatible with TargetRate — both would drive the top
        of the well at once. */
    TopInjection
};

inline const char *WellControlName(WellControl c) noexcept {
    switch (c) {
    case WellControl::WellheadPressure: return "wellhead pressure";
    case WellControl::TargetRate:       return "target rate";
    case WellControl::TopInjection:     return "top injection";
    }
    return "unknown";
}

/** How the well temperature is initialised.

    Replaces the old (delta_T, init_T, follow_reservoir_T) triple, which was
    applied in sequence so each silently overrode the previous — and whose
    `init_T` was only honoured above a magic threshold of 5 degrees, while
    `delta_T` was not a delta at all but a lower bound. */
/** How conductive heat exchange between the well and the formation is modelled. */
enum class RadialHeatModel {
    /** Steady shape factor to the first reservoir control volume — the thermal
        analogue of the Peaceman well index already used for mass:

            Q_loss = 2 pi kappa (T_well - T_node) / ln(r_e / r_w)     [W/m]

        with the same r_e and r_w as the well index. Recommended, and the default.

        Rationale: Ramey's solution exists to stand in for a formation whose
        temperature is NOT solved — it assumes an infinite medium at a fixed
        far-field temperature, heating from t = 0, and its f(t) encodes how far the
        thermal front has diffused. Here the reservoir mesh already resolves and
        evolves the near-well temperature, so the transient is otherwise solved
        twice. That is also why neither the timestep nor the elapsed time is the
        right "t" to feed Ramey: the question should not arise.

        Consequences: no time variable, no validity limit at short timesteps, no
        negative denominator (r_e > r_w always), and the coupling tightens as the
        mesh is refined, exactly as it does for flow. It is mesh-dependent in the
        same way the well index is, which is correct — a coarse control volume
        means the near-well gradient is not resolved and the exchange should be
        weaker. */
    ShapeFactor,

    /** Ramey (1962) transient line-source solution, as in Eq. (6) of
        Lamy-Chappuis et al. (2022):

            Q_loss = 2 pi kappa (T_well - T_inf) / (ln(2 sqrt(alpha t) / r_w) - 0.29)

        Appropriate when the mesh does NOT resolve the near-well formation, which
        is the case Ramey was written for. Uses `radial_heat_time` below for t.

        NOTE it is a long-time approximation: f(t) is negative until roughly
        alpha t / r_w^2 > 1, i.e. hours to days for rock. Below that the model
        returns zero heat loss rather than a negative one. */
    Ramey
};

enum class InitialWellTemperatureMode {
    FollowReservoir,       ///< Twell = Treservoir             (was follow_reservoir_T)
    Uniform,               ///< Twell = value everywhere       (was init_T)
    ReservoirWithMinimum   ///< Twell = max(Treservoir, value) (was delta_T)
};

struct InitialWellTemperature {
    InitialWellTemperatureMode mode = InitialWellTemperatureMode::FollowReservoir;
    double value = 0.;     ///< oC; ignored for FollowReservoir
};

/** How the completion length of a node is obtained. It scales the well index
    linearly, so it matters more than the radius, which enters logarithmically. */
enum class CompletionLengthRule {
    FromElementLengths,  ///< half of each adjoining well element (default)
    FromHostThickness    ///< where the node's host is a lower-dimensional element
    ///< (a fault mid-region), use that element's `thickness`,
    ///< i.e. the aperture. Falls back to FromElementLengths
    ///< for nodes whose host is a volume element.
};

/** Per-well configuration. Everything here describes the CASE, not the numerics —
    values that a user could reasonably want to change without reading the solver.
    Numerical machinery (DAMP, n_of_eq, gravity, the residual initialisers) stays
    inside WellModelPrototype.

    Populate it in the driver (or, later, from a config file) and pass it to the
    constructor. Every member has the value the model used before this struct
    existed, so a default-constructed WellConfiguration reproduces the old
    behaviour — except the radius, which was previously 0.15 m for the well index
    and a hardcoded 0.17/0.12/0.11 m depth profile for the flow equations. Those
    two disagreed; see radius_depths below.

    Benoit DD/MM/YYYY
*/
struct WellConfiguration {

    // ── Geometry ────────────────────────────────────────────────────────────
    /** Radius as a function of depth, as parallel vectors.

        radius_depths  = { z1, z2, ... }   strictly DECREASING (negative downward)
        radius_values  = { r1, r2, r3 ... }

        r1 applies from the well top down to z1, r2 from z1 to z2, and the LAST
        radius continues to the toe — so radius_values.size() must be
        radius_depths.size() + 1.

        For a uniform well leave radius_depths empty and give a single radius.

        Depths are ABSOLUTE (same datum as the mesh, negative downward), not
        relative to the wellhead. Initialize_Well reports the well's own top and
        toe and which interval each segment fell into, so a profile that does not
        reach the well is visible rather than silent. */
    std::vector<double> radius_depths;
    std::vector<double> radius_values{ 0.15 };          // m

    /** Radius of the virtual segments above the mesh top, when
        number_of_virtual_top_segments > 0. Defaults to the topmost radius. */
    double virtual_top_radius = -1.;                    // m; <0 = use radius_values.front()

    double pipe_roughness = 3.e-3;                      // m

    /** Skin applied to completion nodes that do not carry their own. A
        CompletionInterval::skin overrides it for that interval. */
    /** Starting skin for completion nodes; CompletionInterval::skin overrides it
        for that interval. This SEEDS the "well skin" node property — the property
        is what the model reads every step, so writing it during a simulation
        changes the skin from that step on. Changing this struct mid-run does
        nothing. */
    double default_skin   = 0.;

    /** Fraction of the full 2*PI radial flow angle available to this well.

        LEAVE THIS AT 1.0 unless you know you need it. It is 1.0 for any well in
        the interior of the model, which is the normal case. It is only wrong
        there if the well lies ON a model boundary plane — a symmetry model of
        half or a quarter of a pattern — where the well drains through PI (0.5),
        PI/2 on an edge (0.25), or PI/4 in a corner (0.125).

        Not auto-detected: whether a boundary is a symmetry plane is modelling
        intent, and a well that merely touches a boundary is not the same thing.
        Detecting it from the node's boundary membership would be possible and is
        worth doing if symmetry models become common.

        Separate from the control-volume truncation, which `re` already accounts
        for through bulk_volume: correcting one does not correct the other. */
    double angle_fraction = 1.;

    // ── Wellhead pressure control ───────────────────────────────────────────
    /** Rate at which the wellhead pressure moves toward a newly requested value
        [Pa/s]. 0 (the default) applies it immediately. This replaces a hardcoded
        ramp of 0.001 bar/min that used to apply whenever the old `forced` flag
        was false, invisibly from the call site. */
    double wellhead_pressure_ramp = 0.;

    // ── Top injection (WellControl::TopInjection only) ──────────────────────
    /** Temperature and salt mass fraction of the fluid injected at the top.
        Previously hardcoded to 10 oC and 0 inside Top_Injection. */
    double injection_temperature = 10.;    ///< oC
    double injection_smf         = 0.;     ///< salt mass fraction [kg NaCl/kg]

    // ── Conductive heat exchange with the formation ─────────────────────────
    /** See RadialHeatModel. ShapeFactor is consistent with how the mass exchange
        is computed and needs no time variable; Ramey is for meshes that do not
        resolve the near-well formation. */
    RadialHeatModel radial_heat_model = RadialHeatModel::ShapeFactor;

    /** Formation thermal conductivity used by the heat-loss model [W/(m.K)].

        Set it to match the reservoir's own value — 2.0 is the default, which is
        what Lamy-Chappuis et al. (2022) used and what the lithium case runs with.
        It is NOT read from the reservoir automatically: "thermal conductivity" is
        an ELEMENT property and this model only has nodal values to hand. Creating
        a "nodal thermal conductivity" (interpolated from elements, as "nodal
        permeability" already is) would remove this second copy. */
    double formation_thermal_conductivity = 2.0;

    /** Formation thermal diffusivity [m2/s], alpha = kappa / (rho c). Used ONLY by
        RadialHeatModel::Ramey. Around 9.3e-7 for the rock properties in
        Lamy-Chappuis et al. (2022): kappa 2 W/(m.K), rho 2500 kg/m3, c 860 J/(kg.K).

        (The code previously used 1.0 here, about a million times too large, under
        the name "thermal_dispersivity". That inflated value is what kept the Ramey
        function positive at short timesteps; it was compensating for the wrong
        time being passed, not a property of the rock.) */
    double formation_thermal_diffusivity = 9.3e-7;

    /** Initial bulk SALT mass fraction in the well [kg NaCl / kg bulk] — one of
        the model's four primary variables, NOT a steam fraction. Differs between
        an injector, which carries the brine it re-injects, and a producer, which
        is why the configuration is per well rather than shared. */
    double init_bulk_smf  = 0.;

    // ── Completions ─────────────────────────────────────────────────────────
    /** Perforated intervals. An empty list is rejected by Validate(): it would
        give a well with no connection to the reservoir, which was the silent
        outcome of the old -1e8 default arguments. */
    std::vector<CompletionInterval> completions;

    CompletionLengthRule completion_length_rule = CompletionLengthRule::FromElementLengths;

    // ── Discretisation ──────────────────────────────────────────────────────
    double   target_segment_length          = 50.;      // m
    unsigned number_of_virtual_top_segments = 0;        // UNTESTED, see docs

    /// How the well temperature field is set up; see InitialWellTemperature.
    InitialWellTemperature initial_temperature;

    // ── Boundary defaults (run-time changes stay as setters) ────────────────
    double top_pressure    = 1.01325e5;                 // Pa
    double top_temperature = 15.;                       // oC

    // ── Physics ─────────────────────────────────────────────────────────────
    bool with_friction                   = true;
    bool with_grav_pot_energy            = true;
    bool with_kinetic_energy             = true;
    bool with_inertial_terms_in_momentum = true;
    bool with_radial_heat                = true;
    // NOTE: water-table displacement is deliberately NOT here. It is switched on
    // together with the depth it starts from, through Set_well_water_table(),
    // which is also how it changes during a run — a config field would be a
    // second source for the same state. It defaults to OFF.

    /** Where a node's reservoir exchange is placed within the sub-segments of its
        well element. false (default): the whole rate goes on the LAST sub-segment
        of the element. true: the element rate — the mean of its two end nodes — is
        spread evenly over all of that element's sub-segments, which represents a
        distributed completion better. Left off by default because that is what
        every result to date was produced with. Hidden from the config file on
        purpose; change it deliberately and re-check global mass balance, since the
        toe node is a special case that adds rather than averages. */
    bool spread_source = false;

    /** Throws with a readable message if the radius profile is inconsistent:
        wrong sizes, depths not strictly decreasing, or non-positive radii.
        Called by Initialize_Well before anything reads the profile. */
    void Validate() const;

    /** Radius at a given absolute depth, following the interval rules above. */
    double RadiusAtDepth(double depth) const;

    /** The completion interval containing this depth, or nullptr. Overlapping
        intervals are permitted; the first match wins. */
    const CompletionInterval *CompletionAtDepth(double depth) const;
};

template <uint32_t dim> class WellModelPrototype {

public:
    /** @param well_solver_kind  linear solver for the Newton Jacobian.
                    Defaults to Eigen (direct SparseLU): the per-well system is
                    small and dense-ish, which is where a direct solve wins and
                    where PETSc's per-call setup cost would dominate. */
    /** @param config  case description — geometry, discretisation and the physics
                         switches that have a working alternative branch. A
                         default-constructed WellConfiguration reproduces the
                         previous hardcoded behaviour except for the radius; see
                         WellConfiguration::radius_depths. */
    WellModelPrototype(Model<dim> &model, const std::string &well_name,
                       SolverKind well_solver_kind = SolverKind::Eigen,
                       const WellConfiguration &config = WellConfiguration());

    virtual ~WellModelPrototype();

    bool Apply(int count);

    /** Builds the well discretisation, geometry and completions entirely from the
        WellConfiguration supplied at construction — no arguments.

        (The old arguments went into the configuration: min/max_depth_feedzone ->
        completions; delta_T / init_T / follow_reservoir_T -> initial_temperature;
        init_bulk_smf -> init_bulk_smf.) */
    void Initialize_Well();

    void Write_results();
    void Display();
    void Display_short();
    void Display_extended();
    void Advance_well_variables();
    void NR_Advance_well_variables();

    void Reset();

    void WriteSolutionToFile(const std::string &filename);
    void SetTimeStep(const double &time_step);
    // ── Operational control ─────────────────────────────────────────────────
    // Exactly one WellControl mode is active at a time; each setter below
    // switches to its own mode and reports the change. See WellControl.

    /** Requests a wellhead pressure and switches to WellControl::WellheadPressure.
        Applied immediately unless WellConfiguration::wellhead_pressure_ramp is
        non-zero, in which case the pressure moves toward it at that rate — which
        replaces a hardcoded 0.001 bar/min ramp that used to apply whenever the
        old `forced` argument was false. */
    void Set_well_top_pressure(const double &pressure);

    void Set_well_top_temperature(const double &temperature);

    /** Water-table displacement, and the depth it starts from. This is the ONLY
        place the toggle is set — deliberately not a configuration field, so the
        state has one source. Off by default.

        Most useful with WellControl::TopInjection: the table drops, the well is
        refilled, and the pair settles. It is not expected to behave sensibly
        together with WellControl::TargetRate, which drives the top of the well
        at the same time; that combination is warned about. */
    void Set_water_table(const bool &toggle_water_table_displacement,
                         const double &initial_water_table_depth);

    void Initialize_Water_Table();

    /** Fills the well from the top at `mass_rate_injected` [kg/s], with fluid at
        WellConfiguration::injection_temperature and injection_smf, and switches
        the well to WellControl::TopInjection.

        The rate must be POSITIVE; a non-positive value throws. To stop top
        injection, put the well under another control mode. Called every step by
        the scheme for as long as the well stays in this mode, since it delivers
        rate * dt per step. */
    void Top_Injection(double mass_rate_injected);

    /** Holds a target mass rate: POSITIVE produces, NEGATIVE injects. A
        controller adjusts the wellhead pressure each step until the rate is
        matched. `trigger == false` returns the well to
        WellControl::WellheadPressure at whatever pressure it had reached.

        Call it again at any time to change the target — that is the intended way
        to drive a well through a schedule. */
    void Set_target_rate(const bool &trigger,
                         const double &ext_target_rate);

    /// Which mode this well is currently in.
    WellControl Control() const noexcept { return control_; }

    /// One line describing how this well is being driven right now.
    void ReportControl() const;

private:
    /// Switches mode, keeps the legacy flags in step, reports and warns on
    /// combinations that are not expected to work.
    void SwitchControl(WellControl to);

public:
    void PAUSE();

    bool Check_too_large_change_in_total_rate();

    const std::string well_name;
    const std::string getName() const {
        return well_name;
    }

private:
    H2ONaClThermalEquilibrator_Fluid_Only fluid_equilibrator;

    H2ONaClFluidProperties fluid;

    Fluidproperties Liquid;
    Fluidproperties Vapor;
    Fluidproperties Bulk;
    Fluidproperties Salt;

    // Brine    brine;

    // Linear solver for the Newton Jacobian; backend chosen at construction.
    SolverBundle well_bundle_;

    /// Case description supplied at construction; see WellConfiguration.
    WellConfiguration config_;

    /// Active operational mode; see WellControl. The old independent flags
    /// (target_rate_active, injection_mode) are kept in step with it.
    WellControl control_ = WellControl::WellheadPressure;

    /// Per completion node, from the interval it falls in. See CompletionInterval.
    std::vector<double> VO_angle_fraction;

    /// Drainage radius r_e of the reservoir control volume at each node, and its
    /// interpolation to segment centres. Same quantity the well index uses.
    std::vector<double> VO_drainage_radius, VO_D_drainage_radius;

    /// Skin as configured, per node. Seeds the "well skin" property in
    /// Initialize_Well; VO_skin itself is filled from that property every step.
    std::vector<double> VO_configured_skin;

    /// 1 where the node's host reservoir element is lower-dimensional (a fault
    /// mid-region), 0 otherwise; and that element's `thickness` (the aperture).
    /// Filled once in Initialize_Well, used by CompletionLengthRule.
    std::vector<double> V_host_is_lower_dimensional, VO_host_is_lower_dimensional;
    std::vector<double> V_host_thickness, VO_host_thickness;

    /// Well-region node pointers in V_ order, kept so per-node properties (skin)
    /// can be written from VO_-ordered data via node_ordering_vector.
    std::vector<Node<dim> *> V_nodes;

    void T_equilibration();
    void Initial_T_equilibration(const InitialWellTemperature &t0);
    void Initial_P();
    void Interpolate_pressure();
    void Velocity_model();
    void ComputeDensityDerivatives();

    void NR_iteration();
    void SolveWell();
    void NR_update_solution();

    void SetupResidualVectorAndJacobian();
    void SetupResidualVector();
    void SetupJacobianMatrix();

    void Compute_Momentum_Residuals(unsigned n, bool verbose);
    void Compute_Continuity_Residuals(unsigned n, bool verbose);
    void Compute_Energy_Residuals(unsigned n, bool verbose);
    void Compute_Salt_Mass_Residuals(unsigned n, bool verbose);

    void Fill_momentum_Jacobian_part_salt(unsigned n);
    void Fill_continuity_Jacobian_part_salt(unsigned n);
    void Fill_energy_Jacobian_part_salt(unsigned n);
    void Fill_salt_mass_Jacobian_part(unsigned n);

    void ComputeFrictionFactor();
    void Compute_radial_heat();

    void Find_Water_Table(bool move);
    void Element_average();
    void Weighted_Next_Guess(double w1, double w2);
    void Read_reservoir_variables(bool at_initialization);
    void Compute_Equivalent_Radius();
    void Compute_Inflow_Outflow();

    void Check_mass_balance();
    bool Check_density_and_rate_convergence();

    Model<dim> &model_ref_;
    const PropertyDatabase<dim> &prop_ref_;

    csmp::Index
        // Primary variables
        well_pressure_key_, // well fluid pressure
        wellhead_pressure_key_,
        well_velocity_liquid_key_,   // well liquid velocity
        well_velocity_vapor_key_,    // well vapor velocity
        well_velocity_fluid_key_,    // well mixture velocity
        well_s_enthalpy_liquid_key_, // well liquid s_enthalpy
        well_s_enthalpy_vapor_key_,  // well vapor s_enthalpy
        well_s_enthalpy_fluid_key_,  // well bulk s_enthalpy
        well_smf_liquid_key_,        // well liquid smf
        well_smf_vapor_key_,         // well vapor smf
        well_smf_fluid_key_,         // well mixture smf
        well_fluid_state_key_,       // well fluid state

        pressure_key_,          // reservoir fluid pressure
        s_enthalpy_liquid_key_, // reservoir liquid s_enthalpy
        s_enthalpy_vapor_key_,  // reservoir vapor s_enthalpy
        smf_liquid_key_,        // reservoir liquid smf
        smf_vapor_key_,         // reservoir vapor smf
        viscosity_liquid_key_,  // reservoir liquid dyn. viscosity
        viscosity_vapor_key_,   // reservoir vapor dyn. viscosity
        sat_liquid_key_,        // reservoir liquid saturation
        sat_vapor_key_,         // reservoir vapor saturation
        sat_halite_key_,        // reservoir halite saturation
        fluid_state_key_,       // reservoir fluid state

        // Secondary variables
        density_liquid_key_,       // reservoir liquid density
        density_vapor_key_,        // reservoir vapor density
        density_bulk_key_,         // reservoir vapor density
        well_density_liquid_key_,  // well liquid density
        well_density_vapor_key_,   // well vapor density
        well_density_fluid_key_,   // well mixture density
        well_viscosity_fluid_key_, // well mixture viscosity
        well_Re_key_,              // well Reynolds number
        well_sat_liquid_key_,      // well liquid saturation
        well_sat_vapor_key_,       // well vapor saturation
        well_sat_halite_key_,      // well halite saturation
        well_friction_factor_key_,
        well_temperature_key_, // well temperature
        temperature_key_,      // reservoir temperature
        well_sat_pressure_key_, well_sat_temperature_key_,

        // Well geometry
        node_depth_key_, element_depth_key_, well_cross_area_key_,
        well_completion_key_, well_completion_length_key_,
        well_equivalent_radius_key_,

        // Exchange
        well_injectivity_key_, // well injectivity in kg/s/Pa (Pa diferrential
        // pressure well-reservoir at rest)
        mass_transfer_rate_key_, energy_transfer_rate_l_key_,
        energy_transfer_rate_v_key_, energy_transfer_rate_key_,
        salt_mass_transfer_rate_key_, well_index_key_, radial_heat_key_,
        rel_perm_visc_liquid_key_, // relative permeability divided by viscosity
        rel_perm_visc_vapor_key_,  // relative permeability divided by viscosity
        formation_permeability_key_, permeability_key_, bulk_volume_key_,
        well_skin_key_;

    ScalarVariable well_pressure, wellhead_pressure, well_sat_pressure,
        well_velocity_liquid, well_velocity_vapor, well_velocity_fluid,
        well_s_enthalpy_liquid, well_s_enthalpy_vapor, well_s_enthalpy_bulk,
        well_s_enthalpy_fluid, well_smf_liquid, well_smf_vapor, well_fluid_state,
        well_smf_fluid, well_density_liquid, well_density_vapor,
        well_density_fluid, well_viscosity_liquid, well_viscosity_vapor,
        well_viscosity_fluid, well_Re, well_sat_liquid, well_sat_vapor,
        well_sat_halite, well_friction_factor, well_temperature,
        well_sat_temperature, temperature, radial_heat, well_index,

        well_cross_area, well_completion, well_completion_length,
        well_equivalent_radius,

        // transfer_coeff,
        well_injectivity, mass_transfer_rate, energy_transfer_rate_l,
        energy_transfer_rate_v, energy_transfer_rate, salt_mass_transfer_rate;

    double worst_relative_residual, previous_worst_relative_residual,
        pressure_correction_damp, DAMP, Damp, DAMP_JACOBIAN,
        lowest_ever_time_step, test_p, test_v, test_h, test_smf,
        target_segment_length, WT_VELOCITY, velocity_wt, real_velocity_wt,
        enthalpy_below_wt, smf_below_wt, temp_below_wt,

        TOTAL_MASS, Previous_TOTAL_MASS,
        WORST_MASS_MISMATCH = 0.,

        target_rate, total_rate, true_total_rate, previous_total_rate,
        NR_previous_total_rate, abs_p_min, abs_p_max, inflexion_pressure, radius,
        dt,
        // max_value,
        bottom_pressure, pipe_roughness, water_table_depth,
        previous_water_table_depth, NR_previous_water_table_depth, top_of_well,
        time_factor,

        seg_air_saturation_scaling, prev_seg_air_saturation_scaling,

        p, p_prev_seg, p_next_seg, p_prev_time, vm, vm_prev_seg, vm_next_seg,
        vm_prev_time, vm_prev_time_prev_seg, velocity_sign,
        velocity_sign_prev_seg, s_kinetic_energy, s_kinetic_energy_prev_seg,
        s_kinetic_energy_next_seg, s_kinetic_energy_prev_time,

        hf, hf_prev_seg, hf_next_seg, hf_prev_time, hh, hh_prev_seg, hh_next_seg,
        hh_prev_time, hb, hb_prev_seg, hb_next_seg, hb_prev_time,

        rhof, rhof_prev_seg, rhof_next_seg, rhof_prev_time,
        rhof_prev_time_prev_seg, rhoh, rhoh_prev_seg, rhoh_next_seg,
        rhoh_prev_time, rhob, rhob_prev_seg, rhob_next_seg, rhob_prev_time,

        smf_f, smf_f_prev_seg, smf_f_next_seg, smf_f_prev_time, smf_b,
        smf_b_prev_seg, smf_b_next_seg, smf_b_prev_time,

        // satf, satf_prev_seg, satf_next_seg, satf_prev_time,
        sath, sath_prev_seg, sath_next_seg, sath_prev_time,
        sath_prev_time_prev_seg, sath_prev_time_next_seg,

        drhofdp, drhofdhb, drhofdsmfb, drhobdp, drhobdhb, drhobdsmfb, dhdp,
        drhofdp_prev_seg, drhofdhb_prev_seg, drhofdsmfb_prev_seg,
        drhobdp_prev_seg, drhobdhb_prev_seg, drhobdsmfb_prev_seg, dhdp_prev_seg,
        drhofdp_next_seg, drhofdhb_next_seg, drhofdsmfb_next_seg,
        drhobdp_next_seg, drhobdhb_next_seg, drhobdsmfb_next_seg, dhdp_next_seg,

        drhofdp_, drhofdhb_, drhofdsmfb_, drhobdp_, drhobdhb_, drhobdsmfb_, dhdp_,
        drhofdp_prev_seg_, drhofdhb_prev_seg_, drhofdsmfb_prev_seg_,
        drhobdp_prev_seg_, drhobdhb_prev_seg_, drhobdsmfb_prev_seg_,
        dhdp_prev_seg_, drhofdp_next_seg_, drhofdhb_next_seg_,
        drhofdsmfb_next_seg_, drhobdp_next_seg_, drhobdhb_next_seg_,
        drhobdsmfb_next_seg_, dhdp_next_seg_,

        // old:
        drhodp, drhodh, drhodsmf, drhodp_prev_seg, drhodh_prev_seg,
        drhodsmf_prev_seg, drhodp_next_seg, drhodh_next_seg, drhodsmf_next_seg,

        drhodp_, drhodh_, drhodsmf_, drhodp_prev_seg_, drhodh_prev_seg_,
        drhodsmf_prev_seg_, drhodp_next_seg_, drhodh_next_seg_,
        drhodsmf_next_seg_,

        // not used:
        dhdsmf, dhdsmf_prev_seg, dhdsmf_next_seg, dsmfdp, dsmfdp_prev_seg,
        dsmfdp_next_seg, dsmfdh, dsmfdh_prev_seg, dsmfdh_next_seg,

        seg_length, next_seg_length, prev_seg_length, seg_crossarea,
        next_seg_crossarea, prev_seg_crossarea, prev_time_seg_crossarea,
        seg_inclination, prev_seg_inclination, next_seg_inclination,
        friction_factor, friction_factor_prev_seg,

        Momentum_residual, Continuity_residual, Energy_residual,
        Salt_mass_residual,

        // Variables for the "H2ONaClThermalEquilibrator" and
        // "H2ONaClFluidProperties" objects
        m_fluid_, h_fluid_,

        wt_,  // weight percent NaCl []
        smf_, // mass fraction NaCl, i.e. wt/100 []
        x_,   // mole fraction NaCl
        t_, tp_, p_, H_current_, H_previous_, prev_rel_res,

        top_pressure, previous_top_pressure, NR_previous_top_pressure,
        top_velocity, top_s_enthalpy, top_density, top_smf, top_temperature,

        // grav_acc,
        volume_injected, volume_injected_below, volume_injected_above,
        volume_injected_at_wt, enthalpy_injected, enthalpy_injected_below,
        enthalpy_injected_above, enthalpy_injected_at_wt, mass_injected,

        worst_value_p, worst_value_v, worst_value_h, worst_value_smf,
        worst_value_p_, worst_value_v_, worst_value_h_, worst_value_smf_,
        largest_increment;

    const double grav_acc;

    unsigned COUNTER, COUNTER_GOOD_CONVERGENCE, total_iterations,
        highest_ever_number_of_iterations_to_converge, n_of_eq,
        number_of_iterations, number_of_target_prod_rate_iterations,
        number_of_well_elements, number_of_well_nodes, number_of_well_segments,
        number_of_virtual_top_segments, converged_count;

    int wt_index;

    std::string message, message_at_end;

    std::vector<unsigned> element_ordering_vector, node_ordering_vector,
        element_deordering_vector, node_deordering_vector,
        number_of_segment_divisions, index_of_last_segment, VO_D_well_water_table,
        VO_D_previous_well_water_table, VO_D_NR_previous_well_water_table,
        VO_D_well_fluid_state, VO_D_previous_well_fluid_state,
        VO_D_NR_previous_well_fluid_state, VO_well_segment_activated,
        VO_previous_well_segment_activated;

    std::vector<double>
        // Well variables
        //// On depth-ordered mesh elements
        VO_interp_well_pressure, VO_well_velocity_liquid, VO_well_velocity_vapor,
        VO_well_velocity_fluid, VO_well_s_enthalpy_liquid,
        VO_well_s_enthalpy_vapor, VO_well_s_enthalpy_halite,
        VO_well_s_enthalpy_bulk, VO_well_s_enthalpy_fluid, VO_well_smf_liquid,
        VO_well_smf_vapor, VO_well_fluid_state, VO_well_smf_fluid,
        VO_well_smf_bulk,

        VO_well_temperature, VO_well_density_liquid, VO_well_density_vapor,
        VO_well_density_halite, VO_well_density_bulk, VO_well_density_fluid,
        VO_well_sat_liquid, VO_well_sat_vapor, VO_well_sat_halite,
        VO_previous_well_sat_halite, VO_well_mf_liquid, VO_well_mf_vapor,
        VO_well_mf_halite, VO_well_friction_factor, VO_well_cross_area,

        VO_well_viscosity_liquid, VO_well_viscosity_vapor,
        VO_well_viscosity_fluid,

        VO_well_Re,

        VO_Psat, VO_Tsat,

        //// Ordered by depth, downscaled to sub-segments, + virtual top segment
        VO_D_well_pressure, VO_D_interp_well_pressure, VO_D_well_velocity_liquid,
        VO_D_well_velocity_vapor, VO_D_well_velocity_fluid,
        VO_D_well_s_enthalpy_liquid, VO_D_well_s_enthalpy_vapor,
        VO_D_well_s_enthalpy_halite, VO_D_well_s_enthalpy_bulk,
        VO_D_well_s_enthalpy_fluid, VO_D_well_smf_liquid, VO_D_well_smf_vapor,
        VO_D_well_smf_fluid, VO_D_well_smf_bulk,

        VO_D_well_temperature, VO_D_well_density_liquid, VO_D_well_density_vapor,
        VO_D_well_density_halite, VO_D_well_density_bulk, VO_D_well_density_fluid,
        VO_D_well_last_true_density_fluid, VO_D_well_sat_liquid,
        VO_D_well_sat_vapor, VO_D_well_sat_halite, VO_D_well_mf_liquid,
        VO_D_well_mf_vapor, VO_D_well_mf_fluid, VO_D_well_mf_halite,
        VO_D_well_friction_factor, VO_D_well_cross_area,

        VO_D_well_viscosity_liquid, VO_D_well_viscosity_vapor,
        VO_D_well_viscosity_fluid,

        VO_D_well_Re,

        VO_D_well_salt_mass,

        VO_D_Psat, VO_D_Tsat,

        VO_D_well_air_saturation, VO_D_previous_well_air_saturation,
        VO_D_NR_previous_well_air_saturation, VO_D_time_to_edge,

        VO_D_time_out, VO_D_time_in,

        //// Previous version of above
        VO_D_previous_well_pressure, VO_D_previous_interp_well_pressure,
        VO_D_previous_well_velocity_liquid, VO_D_previous_well_velocity_vapor,
        VO_D_previous_well_velocity_fluid, VO_D_previous_well_s_enthalpy_liquid,
        VO_D_previous_well_s_enthalpy_vapor, VO_D_previous_well_s_enthalpy_halite,
        VO_D_previous_well_s_enthalpy_bulk, VO_D_previous_well_s_enthalpy_fluid,
        VO_D_previous_well_smf_liquid, VO_D_previous_well_smf_vapor,
        VO_D_previous_well_smf_fluid, VO_D_previous_well_smf_bulk,

        VO_D_previous_well_temperature, VO_D_previous_well_density_liquid,
        VO_D_previous_well_density_vapor, VO_D_previous_well_density_halite,
        VO_D_previous_well_density_bulk, VO_D_previous_well_density_fluid,
        VO_D_previous_well_sat_liquid, VO_D_previous_well_sat_vapor,
        VO_D_previous_well_sat_halite, VO_D_previous_well_mf_liquid,
        VO_D_previous_well_mf_vapor, VO_D_previous_well_mf_halite,
        VO_D_previous_well_friction_factor, VO_D_previous_well_cross_area,

        VO_D_previous_well_viscosity_liquid, VO_D_previous_well_viscosity_vapor,
        VO_D_previous_well_viscosity_fluid,

        //// NR Previous version of above

        VO_D_NR_previous_well_pressure, VO_D_NR_previous_interp_well_pressure,
        VO_D_NR_previous_well_s_enthalpy_liquid,
        VO_D_NR_previous_well_s_enthalpy_vapor,
        VO_D_NR_previous_well_s_enthalpy_halite,
        VO_D_NR_previous_well_s_enthalpy_bulk,
        VO_D_NR_previous_well_s_enthalpy_fluid, VO_D_NR_previous_well_smf_liquid,
        VO_D_NR_previous_well_smf_vapor, VO_D_NR_previous_well_smf_fluid,
        VO_D_NR_previous_well_smf_bulk,

        VO_D_NR_previous_well_density_liquid, VO_D_NR_previous_well_density_vapor,
        VO_D_NR_previous_well_density_halite, VO_D_NR_previous_well_density_bulk,
        VO_D_NR_previous_well_density_fluid, VO_D_NR_previous_well_sat_liquid,
        VO_D_NR_previous_well_sat_vapor, VO_D_NR_previous_well_sat_halite,
        VO_D_NR_previous_well_mf_liquid, VO_D_NR_previous_well_mf_vapor,
        VO_D_NR_previous_well_mf_halite,

        VO_D_NR_previous_well_viscosity_liquid,
        VO_D_NR_previous_well_viscosity_vapor,
        VO_D_NR_previous_well_viscosity_fluid,

        // radial heat
        V_temperature, VO_temperature, VO_D_temperature, VO_radial_heat,
        VO_E_radial_heat, VO_D_radial_heat,

        VO_D_pressure, // for visualization only

        // Derivatives for Jacobian
        VO_D_drhofdp, VO_D_drhofdhb, VO_D_drhofdsmfb,

        VO_D_drhobdp, VO_D_drhobdhb, VO_D_drhobdsmfb,

        // below, old:
        VO_D_drhodp, VO_D_drhodh, VO_D_drhodsmf,

        // below, unused:
        VO_D_dhdp, VO_D_dhdsmf,

        VO_D_dsmfdp, VO_D_dsmfdh,

        // Well geometry
        V_element_length, V_completion_length, V_equivalent_radius,
        V_element_inclination, VO_element_length, VO_completion_length,
        VO_equivalent_radius, VO_element_inclination, VO_D_segment_length,
        VO_D_segment_inclination, VO_D_segment_heel_depth, V_node_depth,
        V_element_depth,

        // Reservoir variables
        //// On mesh nodes
        V_pressure, V_s_enthalpy_liquid, V_s_enthalpy_vapor, V_density_liquid,
        V_density_vapor, V_density_bulk, V_viscosity_liquid, V_viscosity_vapor,
        V_sat_liquid, V_sat_vapor, V_sat_halite, V_smf_liquid, V_smf_vapor,
        V_fluid_state,

        VO_pressure, VO_pressure_at_rest, VO_well_node_pressure,
        VO_s_enthalpy_liquid, VO_s_enthalpy_vapor, VO_density_liquid,
        VO_density_vapor, VO_density_bulk, VO_viscosity_liquid,
        VO_viscosity_vapor, VO_sat_liquid, VO_sat_vapor, VO_sat_halite,
        VO_smf_liquid, VO_smf_vapor, VO_fluid_state,

        // Exchange
        V_well_index, V_relperm_visc_liquid, V_relperm_visc_vapor,
        V_formation_permeability, V_skin, V_bulk_volume, V_element_heel_depth,
        V_well_completion, VO_well_index, VO_well_injectivity,
        VO_relperm_visc_liquid, VO_relperm_visc_vapor, VO_formation_permeability,
        VO_skin, VO_bulk_volume, VO_element_heel_depth, VO_well_completion,
        VO_well_rate_factor, VO_previous_well_rate_factor, VO_mass_transfer_rate,
        VO_energy_transfer_rate_l, VO_energy_transfer_rate_v,
        VO_energy_transfer_rate, VO_salt_mass_transfer_rate,
        VO_D_mass_transfer_rate, VO_D_previous_mass_transfer_rate,
        VO_D_NR_previous_mass_transfer_rate, VO_D_energy_transfer_rate_l,
        VO_D_energy_transfer_rate_v, VO_D_energy_transfer_rate,
        VO_D_previous_energy_transfer_rate, VO_D_NR_previous_energy_transfer_rate,
        VO_D_salt_mass_transfer_rate, VO_D_previous_salt_mass_transfer_rate,
        VO_D_NR_previous_salt_mass_transfer_rate,

        VO_D_mass_Inj_source, VO_D_energy_Inj_source,

        // N-R
        Residuals, Sol_increment, Sol_increment_pressure,
        Sol_increment_velocity_fluid, Sol_increment_s_enthalpy_bulk,
        Sol_increment_smf_bulk, Sol_increment_s_enthalpy_fluid,
        Sol_increment_smf_fluid, Jacobian_segment_momentum,
        Jacobian_segment_continuity, Jacobian_segment_energy,
        Jacobian_segment_salt_mass;

    SparseMatrix Jacobian;

    bool top_pressure_BC, with_grav_pot_energy, with_kinetic_energy,
        with_radial_heat, with_inertial_terms_in_momentum, with_friction,
        pressure_at_heads, use_center_pressure,
        /** Homogeneous flow velocity model. The alternative — a Drift Flux Model —
            is NOT implemented: the else branch throws. Kept as the hook for that
            work. (`multiphase` and `with_salt` were removed: the first had no
            effect, the second left salt transported but not exchanged.) */
        homogeneous_flow,
        well_solution_converged, density_converged, exit_and_cut_dt,
        /// Kept in step with control_ == WellControl::TargetRate; the physics
        /// still reads it. Do not set it directly — use Set_target_rate().
        target_rate_active,
        injection_mode,
        water_table_depth_has_converged, with_water_table_displacement,
        thermal_eq_succeeded, thermal_eq_crashed_this_NR, spread_source,
        shut_down_clogged_segment, pause_at_end;
};
} // end namespace csmp

#endif