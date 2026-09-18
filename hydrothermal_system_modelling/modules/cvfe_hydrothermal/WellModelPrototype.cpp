// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "WellModelPrototype.h"
#include "CSMP_physical_constants.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include <sstream>

using namespace std;
namespace csmp {


// ── WellConfiguration ───────────────────────────────────────────────────────

void WellConfiguration::Validate() const {
    if (completions.empty())
        throw csmp::Exception(ERROR, "WellConfiguration::Validate",
                              "no completion intervals: the well would have no "
                              "connection to the reservoir. Give at least one.");

    for (size_t i = 0; i < completions.size(); ++i) {
        const CompletionInterval &ci = completions[i];
        if (ci.bottom_depth >= ci.top_depth) {
            std::ostringstream msg;
            msg << "completion " << i << ": bottom_depth (" << ci.bottom_depth
                << ") must be BELOW top_depth (" << ci.top_depth
                << "), i.e. more negative.";
            throw csmp::Exception(ERROR, "WellConfiguration::Validate", msg.str());
        }
    }

    if (radius_values.empty())
        throw csmp::Exception(ERROR, "WellConfiguration::Validate",
                              "radius_values is empty: give at least one radius.");

    if (radius_values.size() != radius_depths.size() + 1) {
        std::ostringstream msg;
        msg << "radius_values has " << radius_values.size() << " entries and "
            << "radius_depths has " << radius_depths.size()
            << "; expected radius_values.size() == radius_depths.size() + 1 "
               "(the last radius runs from the deepest break to the toe).";
        throw csmp::Exception(ERROR, "WellConfiguration::Validate", msg.str());
    }

    if (angle_fraction <= 0. || angle_fraction > 1.) {
        std::ostringstream msg;
        msg << "angle_fraction (" << angle_fraction << ") must be in (0, 1]. "
               "It is 1.0 for any well in the interior of the model; only a well "
               "lying on a symmetry plane needs less.";
        throw csmp::Exception(ERROR, "WellConfiguration::Validate", msg.str());
    }

    for (size_t i = 0; i < radius_values.size(); ++i)
        if (radius_values[i] <= 0.) {
            std::ostringstream msg;
            msg << "radius_values[" << i << "] = " << radius_values[i]
                << "; radii must be positive.";
            throw csmp::Exception(ERROR, "WellConfiguration::Validate", msg.str());
        }

    for (size_t i = 1; i < radius_depths.size(); ++i)
        if (radius_depths[i] >= radius_depths[i - 1]) {
            std::ostringstream msg;
            msg << "radius_depths must be strictly decreasing (negative downward), "
                << "but entry " << i << " (" << radius_depths[i] << ") is not below "
                << "entry " << i - 1 << " (" << radius_depths[i - 1] << ").";
            throw csmp::Exception(ERROR, "WellConfiguration::Validate", msg.str());
        }
}

const CompletionInterval *WellConfiguration::CompletionAtDepth(double depth) const {
    // Overlapping intervals are allowed; the first match wins.
    for (const CompletionInterval &ci : completions)
        if (depth <= ci.top_depth && depth >= ci.bottom_depth)
            return &ci;
    return nullptr;
}

double WellConfiguration::RadiusAtDepth(double depth) const {
    // radius_values[0] from the well top down to radius_depths[0],
    // radius_values[i] from radius_depths[i-1] down to radius_depths[i],
    // radius_values.back() from the last break to the toe.
    for (size_t i = 0; i < radius_depths.size(); ++i)
        if (depth > radius_depths[i])
            return radius_values[i];
    return radius_values.back();
}

template <uint32_t dim>
WellModelPrototype<dim>::WellModelPrototype(Model<dim> &model,
                                            const std::string &well_name,
                                            SolverKind well_solver_kind,
                                            const WellConfiguration &config)
    : // 1D well N-R solver
      // Only works properly within a 3D mesh

      well_name(well_name),

      fluid_equilibrator(m_fluid_, wt_, tp_, p_, H_current_),

      // Fluid object from which to fetch thermodynamic properties
      fluid(t_, p_, x_, h_fluid_, 0., 0., 1., false),
      // brine        ( t_, p_, smf_ ),

      well_bundle_(well_solver_kind), config_(config),

      model_ref_(model), prop_ref_(model.Database()),

      well_pressure_key_(prop_ref_.StorageKey("well fluid pressure")),
      wellhead_pressure_key_(prop_ref_.StorageKey("wellhead pressure")),
      well_velocity_liquid_key_(prop_ref_.StorageKey("well velocity liquid")),
      well_velocity_vapor_key_(prop_ref_.StorageKey("well velocity vapor")),
      well_velocity_fluid_key_(prop_ref_.StorageKey("well velocity fluid")),
      well_s_enthalpy_liquid_key_(prop_ref_.StorageKey("well enthalpy liquid")),
      well_s_enthalpy_vapor_key_(prop_ref_.StorageKey("well enthalpy vapor")),
      well_s_enthalpy_fluid_key_(prop_ref_.StorageKey("well enthalpy fluid")),
      well_smf_liquid_key_(prop_ref_.StorageKey("well salt fraction liquid")),
      well_smf_vapor_key_(prop_ref_.StorageKey("well salt fraction vapor")),
      well_smf_fluid_key_(prop_ref_.StorageKey("well salt fraction fluid")),

      pressure_key_(prop_ref_.StorageKey("fluid pressure")),
      s_enthalpy_liquid_key_(prop_ref_.StorageKey("enthalpy liquid")),
      s_enthalpy_vapor_key_(prop_ref_.StorageKey("enthalpy vapor")),
      smf_liquid_key_(prop_ref_.StorageKey("salt fraction liquid")),
      smf_vapor_key_(prop_ref_.StorageKey("salt fraction vapor")),

      viscosity_liquid_key_(prop_ref_.StorageKey("viscosity liquid")),
      viscosity_vapor_key_(prop_ref_.StorageKey("viscosity vapor")),
      sat_liquid_key_(prop_ref_.StorageKey("saturation liquid")),
      sat_vapor_key_(prop_ref_.StorageKey("saturation vapor")),
      sat_halite_key_(prop_ref_.StorageKey("saturation halite")),
      fluid_state_key_(prop_ref_.StorageKey("fluid state")),
      density_liquid_key_(prop_ref_.StorageKey("density liquid")),
      density_vapor_key_(prop_ref_.StorageKey("density vapor")),

      well_density_liquid_key_(prop_ref_.StorageKey("well density liquid")),
      well_density_vapor_key_(prop_ref_.StorageKey("well density vapor")),
      well_density_fluid_key_(prop_ref_.StorageKey("well density fluid")),
      well_viscosity_fluid_key_(prop_ref_.StorageKey("well viscosity fluid")),
      well_Re_key_(prop_ref_.StorageKey("well Re number")),
      well_sat_liquid_key_(prop_ref_.StorageKey("well saturation liquid")),
      well_sat_vapor_key_(prop_ref_.StorageKey("well saturation vapor")),
      well_sat_halite_key_(prop_ref_.StorageKey("well saturation halite")),
      well_friction_factor_key_(prop_ref_.StorageKey("well friction factor")),
      well_temperature_key_(prop_ref_.StorageKey("well temperature")),
      temperature_key_(prop_ref_.StorageKey("temperature")),

      node_depth_key_(prop_ref_.StorageKey("coordinate y")),
      element_depth_key_(prop_ref_.StorageKey("barycenter depth")),

      well_cross_area_key_(prop_ref_.StorageKey("well cross area")),
      well_completion_key_(prop_ref_.StorageKey("well completion")),
      well_completion_length_key_(
          prop_ref_.StorageKey("well completion length")),
      well_equivalent_radius_key_(
          prop_ref_.StorageKey("well effective radius")),

      well_injectivity_key_(prop_ref_.StorageKey("well injectivity")),
      mass_transfer_rate_key_(prop_ref_.StorageKey("well mass transfer rate")),
      energy_transfer_rate_l_key_(
          prop_ref_.StorageKey("well energy transfer rate liquid")),
      energy_transfer_rate_v_key_(
          prop_ref_.StorageKey("well energy transfer rate vapor")),
      salt_mass_transfer_rate_key_(
          prop_ref_.StorageKey("well salt mass transfer rate")),
      well_index_key_(prop_ref_.StorageKey("well index")),
      radial_heat_key_(prop_ref_.StorageKey("well reservoir radial heat")),
      rel_perm_visc_liquid_key_(prop_ref_.StorageKey(
                                    "relperm viscosity liquid")), // reservoir relative permeability
      // divided by viscosity (liquid)
      rel_perm_visc_vapor_key_(prop_ref_.StorageKey(
                                   "relperm viscosity vapor")), // reservoir relative permeability
      // divided by viscosity (vapor)
      formation_permeability_key_(prop_ref_.StorageKey(
                                      "nodal permeability")), // reservoir permeability extrapolated to
      // nodes SHOULD USE HORIZ PERM FOR VERT WELL
      permeability_key_(prop_ref_.StorageKey(
                            "permeability")), // reservoir scalar permeability
      bulk_volume_key_(prop_ref_.StorageKey("bulk volume")), // control volume
      well_skin_key_(prop_ref_.StorageKey("well skin")),

      worst_relative_residual(1.e8), previous_worst_relative_residual(1.e8),
      largest_increment(1.e8), DAMP(0.2), DAMP_JACOBIAN(1.), target_rate(0.),
      dt(0.), lowest_ever_time_step(1.e10),
      highest_ever_number_of_iterations_to_converge(0),

      //////////////////////////////////////
      pipe_roughness(config.pipe_roughness), // see WellConfiguration; historic note:
      // 0.02mm!!!!!!!!!!!!!!!!!
      water_table_depth(0.), // in m INITIAL WATER TABLE DEPTH
      previous_water_table_depth(water_table_depth),

      top_velocity(0.), // BC top velocity for flux target case, in that case it
      // replaces BC top/bottom pressure
      top_s_enthalpy(
          0.), // BC top enthalpy if top well segment behave as an injector,
      // might be better to have temperature instead
      top_smf(0.), // BC top salt mass fraction for injector
      top_temperature(config.top_temperature), top_pressure(config.top_pressure), mass_injected(0.),
      total_rate(0.), previous_total_rate(0.),

      grav_acc(
          9.81), // Constant, although it might change with depth and latitude!

      //////////////////////////////////////
      target_segment_length(
          config.target_segment_length), // subdivide the well elements into segments of the
      // target length, the average of the segments is then used to
      // store back on elements for visualization.
      // this could be variable depending on element's lengths or even be able
      // to change dynamically during the simulation

      n_of_eq(
          4), // number of equations, eventually, this should not be hardcoded
      converged_count(0), COUNTER(0), COUNTER_GOOD_CONVERGENCE(0),
      total_iterations(0),

      with_grav_pot_energy(config.with_grav_pot_energy), // Option to add grav. pot. energy to energy
      // conservation (generally important)
      with_kinetic_energy(config.with_kinetic_energy), // add kinetic energy to
      // energy conservation (generally negligible)
      with_radial_heat(config.with_radial_heat),
      with_inertial_terms_in_momentum(config.with_inertial_terms_in_momentum),
      // ^ generally very negligible
      with_friction(config.with_friction),

      use_center_pressure(true), // for thermal equilibration and derivatives (
      // not sure if the latter is very consistent..)

      homogeneous_flow(true), // Option to use homogeneous flow velocity model,
      // other option would be Drift Flux Model
      target_rate_active(false), injection_mode(false),
      spread_source(config.spread_source), shut_down_clogged_segment(true),
      // OFF until Set_water_table() switches it on, which is also where the
      // starting depth comes from. Not a configuration field: one source of
      // truth for a piece of state that changes during a run.
      with_water_table_displacement(false)

    // - Drift Flux Model implementation,

    // Other To do:

    // - improve the calculation of dRho_dP, dRho_dh, dh_dp

    // - add a function to set/modify BC externally.

    // - make various variable checks in constructor.

    // - add flux target mode.

    // - add bottom pressure BC.

    // - Optimize the N-R procedure

{
    Liquid.InitToZero(); // Useful?
    Vapor.InitToZero();  // Useful?
    Bulk.InitToZero();   // Useful?
    Salt.InitToZero();   // Useful?
    fluid_equilibrator.AttemptToSurviveFluidPropertiesError(true);

    // ── Linear solver for the Newton Jacobian ───────────────────────────────
    // Backend chosen by the caller (Eigen by default). Library defaults are
    // used deliberately, and are appropriate here: the per-well system is a few
    // unknowns per segment — small, dense-ish and non-symmetric — so a direct
    // factorisation is both fast and exact, and there is nothing for an
    // iterative method's tolerances to trade off.
    //
    //   Eigen  no settings exist; EigenSolver always does a direct SparseLU.
    //   PETSc  switched to a direct solve below ("preonly" + "lu") so it
    //          matches Eigen rather than iterating on a tiny system. Note
    //          PETSc rebuilds its Mat/Vec/KSP on every call, so per-iteration
    //          overhead may exceed the solve itself here.
    //   SAMG   defaults only. Before using it, check the no_unknowns/nsys
    //          handling: SolveWell passes 4 (pressure, velocity, enthalpy,
    //          salt per segment), and SAMG then builds its iu array with the
    //          unknown-based pattern [1,1,..,2,2,..], whereas this Jacobian is
    //          ordered point-wise [1,2,3,4,1,2,3,4,..]. That mismatch needs
    //          resolving (napproach, or an explicit iu) before trusting it.
    well_bundle_.UseDirectFactorisation();   // no-op unless the backend is PETSc

} // end WellModelPrototype

template <uint32_t dim>
WellModelPrototype<dim>::~WellModelPrototype() {} // end ~WellModelPrototype

// The main task:
template <uint32_t dim>
bool WellModelPrototype<dim>::Apply(int well_solves_count) {
    cerr << endl
         << "********************************" << endl
         << "Start of well solve, well name: " << well_name << endl
         << "********************************";

    // re-initialize some variables and booleans
    exit_and_cut_dt = false; // bool to stop computations, reset and cut the
    // timestep if various checks are not passed
    pause_at_end = false;

    if (!with_water_table_displacement)
        water_table_depth_has_converged = true;

    // water table location
    wt_index = -1;

    //////////////////////////////////////
    std::string filename("Well_segments_");
    char Incr[200];
    sprintf(Incr, "%ld", static_cast<long>(converged_count));
    filename += Incr;
    //////////////////////////////////////

    message = "";
    message_at_end = "";

    // cerr<<"Starting state: **************************"<<endl;

    if (well_solves_count == 0) {
        Compute_radial_heat();

        for (unsigned n = 1; n < number_of_well_nodes;
             n++) // we skip n = 0 because surely there will be no exchange term on
            // top virtual segment
        {
            // VO_well_rate_factor [ n ] = 3.5 - VO_well_sat_halite [ n - 1 ] * 5.; //
            // factor 0 to 1 when sat h goes from 0.5 to 0.7 VO_well_rate_factor [ n ]
            // = 3. - VO_well_sat_halite [ n - 1 ] * 5.; // factor 0 to 1 when sat h
            // goes from 0.4 to 0.6 VO_well_rate_factor [ n ] = 1.25 -
            // VO_well_sat_halite [ n - 1 ] * 2.5; // factor 0 to 1 when sat h goes
            // from 0.1 to 0.5 // A option
            VO_well_rate_factor[n] =
                    1.75 -
                    VO_well_sat_halite[n - 1] *
                    2.5; // factor 0 to 1 when sat h goes from 0.3 to 0.7 // B option
            // VO_well_rate_factor [ n ] = 1.5 - VO_well_sat_halite [ n - 1 ] * 5./3.;
            // // factor 0 to 1 when sat h goes from 0.3 to 0.9 // C option
            // VO_well_rate_factor [ n ] = 1.;
            if (VO_well_rate_factor[n] > 1.)
                VO_well_rate_factor[n] = 1.;
            // if( VO_well_rate_factor [ n ] < 0.01 ) VO_well_rate_factor [ n ] =
            // 0.01;

            // if(VO_previous_well_rate_factor [ n ] < 0.010000000001)
            // VO_well_rate_factor [ n ] = 0.01;//makes blockage irreversible...

            VO_well_rate_factor[n] = 1.; // disable well rate factor
        }
    }

    // FEB 2025 conmment all below
    Display_short();

    // Velocity_model();

    // //if( with_water_table_displacement ) Find_Water_Table(true); //IS THIS
    // NECESSARY? Find_Water_Table(with_water_table_displacement ); //IS THIS
    // NECESSARY?

    // Interpolate_pressure();

    // T_equilibration();

    // if( with_friction ) ComputeFrictionFactor();

    // ComputeDensityDerivatives();

    // Compute_Inflow_Outflow();

    // Weighted_Next_Guess(0.5, 0.5);//ADD FEB 2022

    NR_iteration(); // N-R loop

    if (exit_and_cut_dt) {
        cerr << endl
             << "WELL MODEL EXIT AND RESET!! The timestep was: " << dt << endl;
        DAMP = 0.1;
        // DAMP_JACOBIAN = 1.;
        if (dt < 1.) {
            PAUSE();
            pause_at_end = true;
        }
    }

    else {
        // Display_short();
        // Display();
        converged_count++;

        Check_mass_balance();

        if (number_of_iterations > highest_ever_number_of_iterations_to_converge)
            highest_ever_number_of_iterations_to_converge = number_of_iterations;

        // cerr<<endl<<"Highest_ever_number_of_iterations_to_converge :
        // "<<highest_ever_number_of_iterations_to_converge;
        // cerr<<endl<<"Lowest_ever_time_step : "<<lowest_ever_time_step;

        // DAMP_JACOBIAN *= 1.1;
        // DAMP_JACOBIAN = 1.;
        DAMP *= 1.005;
        // if( DAMP < 0.1 ) DAMP = 0.1; //TEST FEB 2022

        if (thermal_eq_crashed_this_NR) {
            cerr << "The well solution converged but T-eq did not converge in at "
                    "least one segment!!!"
                 << endl;
            // PAUSE();
        }
    }

    if (pause_at_end) {
        cerr << endl << endl << "Pause at end: " << message_at_end << endl;
        Display();
        PAUSE();
    }

    return (exit_and_cut_dt); // If returned true to the reservoir CVFEM, we skip
    // further CVFEM calculations and reset the
    // reservoir properties.
}

template <uint32_t dim>
void WellModelPrototype<dim>::Initialize_Well() {

    top_smf = config_.init_bulk_smf;
    // Well segments modelled above the domain top boundary, for a mesh that does
    // not reach the surface. UNTESTED: this has been 0 for a long time and
    // nothing exercises those branches. Before trusting it, check that the Ramey
    // radial-heat term means anything on a segment with no formation around it,
    // and that virtual segments are excluded from the loops that index reservoir
    // arrays. (The hardcoded per-case value that used to sit here is now
    // WellConfiguration::number_of_virtual_top_segments.)
    number_of_virtual_top_segments = config_.number_of_virtual_top_segments;

    // Determine the number of elements, nodes and segments of Well region
    number_of_well_elements = model_ref_.Region(well_name).Cells();
    number_of_well_nodes = model_ref_.Region(well_name).Nodes();
    //

    number_of_segment_divisions.resize(number_of_well_elements);
    index_of_last_segment.resize(number_of_well_elements);
    element_ordering_vector.resize(number_of_well_elements);
    element_deordering_vector.resize(number_of_well_elements);
    node_ordering_vector.resize(number_of_well_nodes);
    node_deordering_vector.resize(number_of_well_nodes);

    // Well geometry
    V_element_length.resize(number_of_well_elements);
    VO_element_length.resize(number_of_well_elements);
    V_element_heel_depth.resize(number_of_well_elements);
    V_element_inclination.resize(number_of_well_elements);
    VO_element_inclination.resize(number_of_well_elements);
    V_element_depth.resize(number_of_well_elements);
    V_node_depth.resize(number_of_well_nodes);

    top_of_well = -1.e8;
    double segment_end_depth1, segment_end_depth2;

    unsigned n = 0;
    for (typename std::vector<Element<dim> *>::const_iterator eit =
         model_ref_.Region(well_name).CellsBegin();
         eit != model_ref_.Region(well_name).CellsEnd(); ++eit) {
        segment_end_depth1 = (*eit)->N(0)->y();
        segment_end_depth2 = (*eit)->N(1)->y();

        if (segment_end_depth1 > top_of_well)
            top_of_well = segment_end_depth1;
        if (segment_end_depth2 > top_of_well)
            top_of_well = segment_end_depth2;
        n++;
    }

    cerr << "Top of well, i.e. top of top well element: " << top_of_well << endl;

    if (top_of_well >= 0. && number_of_virtual_top_segments > 0) {
        cerr << "The top of the well seems to already be at the surface and a "
                "positive number of top virtual well segments was set"
             << endl
             << "Setting it back to 0!! NOTE: this is assuming that the surface is "
                "at a depth (i.e. y coordinate) of 0 m!!!!!"
             << endl;

        number_of_virtual_top_segments = 0;
        PAUSE();
        PAUSE();
    }

    if (number_of_virtual_top_segments > 0.1)
        VO_element_heel_depth.resize(
                    number_of_well_nodes); // Add a "well element" that will host virtual
    // top segments !!!
    else
        VO_element_heel_depth.resize(number_of_well_elements);

    n = 0;
    for (typename std::vector<Element<dim> *>::const_iterator eit =
         model_ref_.Region(well_name).CellsBegin();
         eit != model_ref_.Region(well_name).CellsEnd(); ++eit) {
        segment_end_depth1 = (*eit)->N(0)->y();
        segment_end_depth2 = (*eit)->N(1)->y();

        V_element_length[n] =
                ((*eit)->N(0)->Coordinate().DistanceTo((*eit)->N(1)->Coordinate()));
        V_element_heel_depth[n] = std::min(segment_end_depth1, segment_end_depth2);
        V_element_inclination[n] = (acos(
                                        abs(segment_end_depth1 - segment_end_depth2) /
                                        ((*eit)->N(0)->Coordinate().DistanceTo((*eit)->N(1)->Coordinate()))));

        if (segment_end_depth1 > top_of_well)
            top_of_well = segment_end_depth1;
        if (segment_end_depth2 > top_of_well)
            top_of_well = segment_end_depth2;
        n++;
    }

    // Defining ordering vectors for later use. So that the well variables are
    // ordered from heel to toe ( which would not be the case if we simply loop
    // through region elements/nodes using for example: for( typename std::vector<
    // Element<dim>* >::const_iterator eit = ref.CellsBegin(); eit  !=
    // ref.CellsEnd(); ++eit )

    // Elements ordering vector
    n = 0;
    for (typename std::vector<Element<dim> *>::const_iterator eit =
         model_ref_.Region(well_name).CellsBegin();
         eit != model_ref_.Region(well_name).CellsEnd(); ++eit) {
        V_element_depth[n] = ((*eit)->Read(element_depth_key_));
        n++;
    }

    int x = 0;
    std::iota(element_ordering_vector.begin(), element_ordering_vector.end(),
              x++); // Initializing
    stable_sort(
                element_ordering_vector.begin(), element_ordering_vector.end(),
                [&](int i, int j) {
        return V_element_depth[i] > V_element_depth[j];
    });

    x = 0;
    std::iota(element_deordering_vector.begin(), element_deordering_vector.end(),
              x++); // Initializing
    stable_sort(element_deordering_vector.begin(),
                element_deordering_vector.end(), [&](int i, int j) {
        return element_ordering_vector[i] < element_ordering_vector[j];
    });

    n = 0;
    for (typename std::vector<Node<dim> *>::const_iterator nit =
         model_ref_.Region(well_name).NodesBegin();
         nit != model_ref_.Region(well_name).NodesEnd(); ++nit) {
        V_node_depth[n] = ((*nit)->Read(node_depth_key_));
        n++;
    }

    x = 0;
    std::iota(node_ordering_vector.begin(), node_ordering_vector.end(),
              x++); // Initializing
    stable_sort(node_ordering_vector.begin(), node_ordering_vector.end(),
                [&](int i, int j) {
        return V_node_depth[i] > V_node_depth[j];
    });

    x = 0;
    std::iota(node_deordering_vector.begin(), node_deordering_vector.end(),
              x++); // Initializing
    stable_sort(node_deordering_vector.begin(), node_deordering_vector.end(),
                [&](int i, int j) {
        return node_ordering_vector[i] < node_ordering_vector[j];
    });

    cerr << "\n ***************************" << endl;

    // We immediately order and write lengths and inclinations ( we assume they
    // will be constant... because no mechanical deformations )
    for (unsigned n = 0; n < number_of_well_elements; n++) {
        VO_element_length[n] = V_element_length[element_ordering_vector[n]];
        VO_element_inclination[n] =
                V_element_inclination[element_ordering_vector[n]];

        if (number_of_virtual_top_segments > 0.1)
            VO_element_heel_depth[n + 1] =
                    V_element_heel_depth[element_ordering_vector[n]];
        else
            VO_element_heel_depth[n] =
                    V_element_heel_depth[element_ordering_vector[n]];
    }

    if (number_of_virtual_top_segments > 0.1)
        VO_element_heel_depth[0] = top_of_well;

    cerr << "Ordered ele length and HEEL depth" << endl;
    for (unsigned n = 0; n < number_of_well_elements; n++) {
        if (number_of_virtual_top_segments > 0.1)
            cerr << "Element" << n << " has a length of " << VO_element_length[n]
                    << " and a heel depth of " << VO_element_heel_depth[n + 1] << ", "
                    << endl;
        else
            cerr << "Element" << n << " has a length of " << VO_element_length[n]
                    << " and a heel depth of " << VO_element_heel_depth[n] << ", "
                    << endl;
    }

    // Set completion nodes given the depth range provided
    // MAKE CHECKS!!! Elements associated with completions should not be
    // divided!!!

    // Can the last element be associated with a completion?
    // Yes but with a particular mesh: with an added line element at the end of
    // the well that does not belong to the WELL region! Otherwise the well
    // equivalent radius will likely be wrong!

    VO_well_completion.resize(number_of_well_nodes);
    VO_well_rate_factor.resize(number_of_well_nodes);
    VO_previous_well_rate_factor.resize(number_of_well_nodes);
    V_completion_length.resize(number_of_well_nodes);  // NEW, not used
    VO_completion_length.resize(number_of_well_nodes); // NEW

    VO_angle_fraction.assign(number_of_well_nodes, config_.angle_fraction);

    // Skin per completion node, from the interval it falls in. Kept separate from
    // VO_skin, which is not sized until much later and is filled every step from
    // the "well skin" property — see the seeding block below.
    VO_configured_skin.assign(number_of_well_nodes, config_.default_skin);

    std::fill(VO_well_rate_factor.begin(), VO_well_rate_factor.end(), 1.);
    VO_previous_well_rate_factor = VO_well_rate_factor;

    // ── Host reservoir element: lower-dimensional? and its thickness ────────
    // A well node's parents include its own line elements and the reservoir
    // elements it sits in. A surface parent means the node lies in a
    // lower-dimensional region — a fault mid-region — whose `thickness` is the
    // aperture. Wells are never modelled in 2D, so IsLine() identifies the well
    // itself and is skipped here.
    {
        const csmp::Index thickness_key(prop_ref_.StorageKey("thickness"));
        V_host_is_lower_dimensional.assign(number_of_well_nodes, 0.);
        V_host_thickness.assign(number_of_well_nodes, 0.);

        V_nodes.assign(number_of_well_nodes, nullptr);

        unsigned nh = 0;
        for (typename vector<Node<dim> *>::const_iterator it(
                 model_ref_.Region(well_name).NodesBegin());
             it != model_ref_.Region(well_name).NodesEnd(); ++it) {

            V_nodes[nh] = *it;

            const size_t nParents = (*it)->Parents();
            for (size_t i = 0; i < nParents; i++) {
                if ((*it)->Parent(i)->IsLine())     continue;   // the well itself
                if (!(*it)->Parent(i)->IsSurface()) continue;   // a volume element
                V_host_is_lower_dimensional[nh] = 1.;
                V_host_thickness[nh] = (*it)->Parent(i)->Read(thickness_key);
                break;                                          // first match wins
            }
            nh++;
        }

        VO_host_is_lower_dimensional.resize(number_of_well_nodes);
        VO_host_thickness.resize(number_of_well_nodes);
        for (unsigned n = 0; n < number_of_well_nodes; n++) {
            VO_host_is_lower_dimensional[n] =
                    V_host_is_lower_dimensional[node_ordering_vector[n]];
            VO_host_thickness[n] = V_host_thickness[node_ordering_vector[n]];
        }
    }

    // ── Completions ─────────────────────────────────────────────────────────
    // A node is a completion if its depth falls inside ANY configured interval.
    // Only completion nodes exchange with the reservoir (see the VO_well_completion
    // test in Reservoir_exchange), so this decides where the well is connected.
    //
    // Completion length scales the well index linearly. The default rule gives a
    // node half of each adjoining well element; FromHostThickness overrides that
    // where the node's host is a lower-dimensional element (a fault mid-region),
    // using the element's `thickness`, i.e. the aperture.
    if (config_.angle_fraction < 1.)
        cerr << endl << "Radial flow angle fraction for this well: "
             << config_.angle_fraction << " (well on a symmetry boundary)";

    cerr << endl << "Completion intervals (absolute depth, negative downward):";
    for (size_t i = 0; i < config_.completions.size(); ++i)
        cerr << endl << "  [" << i << "] " << config_.completions[i].top_depth
             << " m to " << config_.completions[i].bottom_depth
             << " m, skin "
             << config_.completions[i].skin.value_or(config_.default_skin)
             << (config_.completions[i].skin ? "" : " (well default)")
;

    std::vector<unsigned> nodes_per_interval(config_.completions.size(), 0U);
    unsigned n_completions = 0U;
    unsigned nn = 0;

    for (unsigned n = 1; n < number_of_well_nodes; n++) { // we skip node 0
        if (number_of_virtual_top_segments > 0.1)
            nn = n;
        else
            nn = n - 1;

        const double node_depth = VO_element_heel_depth[nn];

        cerr << "Node number " << n << " has a depth of: " << node_depth << endl;

        const CompletionInterval *ci = config_.CompletionAtDepth(node_depth);
        if (ci == nullptr)
            continue;

        VO_well_completion[n] = 1;
        VO_configured_skin[n] = ci->skin.value_or(config_.default_skin);
        VO_angle_fraction[n] = config_.angle_fraction;
        ++n_completions;
        for (size_t i = 0; i < config_.completions.size(); ++i)
            if (&config_.completions[i] == ci) ++nodes_per_interval[i];

        cerr << " therefore it is considered to be a completion" << endl;

        // default rule: half of each adjoining well element
        if (n == number_of_well_nodes - 1)
            VO_completion_length[n] = VO_element_length[n - 1] / 2.;
        else
            VO_completion_length[n] =
                    (VO_element_length[n - 1] + VO_element_length[n]) / 2.;

        // override where the host is lower-dimensional and the rule asks for it
        if (config_.completion_length_rule == CompletionLengthRule::FromHostThickness
                && VO_host_is_lower_dimensional[n] > 0.5) {
            VO_completion_length[n] = VO_host_thickness[n];
            cerr << "   host element is lower-dimensional: completion length set to"
                    " its thickness (" << VO_completion_length[n] << " m)" << endl;
        }
    }

    // ── Seed the "well skin" property from the configuration ────────────────
    // The property, not VO_skin, is the runtime source: Read_reservoir_variables
    // re-reads it every step, so anything that writes it during the simulation
    // takes effect on the next step. The configuration only sets the starting
    // value — default_skin on the whole well, overridden per completion interval.
    {
        const csmp::Index skin_key(prop_ref_.StorageKey("well skin"));
        ScalarVariable skin_value(ANY, config_.default_skin);
        for (unsigned n = 0; n < number_of_well_nodes; n++) {
            Node<dim> *node = V_nodes[node_ordering_vector[n]];
            if (node == nullptr) continue;
            skin_value() = VO_configured_skin[n];
            node->Store(skin_key, skin_value);
        }
    }

    cerr << endl << "Completion nodes found: " << n_completions;
    for (size_t i = 0; i < nodes_per_interval.size(); ++i) {
        cerr << endl << "  interval [" << i << "] captured "
             << nodes_per_interval[i] << " node(s)";
        if (nodes_per_interval[i] == 0U)
            cerr << "   <-- WARNING: no node falls in this interval; check the"
                    " depths against the well's own extent reported above";
    }
    if (n_completions == 0U)
        throw csmp::Exception(ERROR, "WellModelPrototype::Initialize_Well",
                              well_name,
                              "no node falls inside any completion interval: this "
                              "well would never exchange with the reservoir.");

    // Dividing ordered element lenghts into sub segments based on target segment
    // length

    cerr << "Defining " << number_of_virtual_top_segments
         << " virtual top segments" << endl;

    for (unsigned n = 0; n < number_of_well_elements; n++) {

        number_of_segment_divisions[n] =
                floor(VO_element_length[n] / target_segment_length + 0.5);

        if (number_of_segment_divisions[n] < 1)
            number_of_segment_divisions[n] = 1;

        cerr << "Element " << n << " is divided in "
             << number_of_segment_divisions[n] << "segments" << endl;
    }

    number_of_well_segments =
            std::accumulate(number_of_segment_divisions.begin(),
                            number_of_segment_divisions.end(), 0) +
            number_of_virtual_top_segments;

    cerr << "number of elements: " << number_of_well_elements
         << ", number of nodes: " << number_of_well_nodes
         << ", number of well segments: " << number_of_well_segments << endl;

    // ── Well cross-sectional area from the configured radius profile ────────
    // ONE definition of the radius now feeds both arrays:
    //   VO_well_cross_area    (per ELEMENT)  -> rw in the well-index formula
    //   VO_D_well_cross_area  (per SEGMENT)  -> the flow equations
    // These used to disagree: the element array came from a hardcoded 0.15 m
    // while the segment array was overwritten by a hardcoded 0.17/0.12/0.11 m
    // depth profile, so exchange and flow saw different wells.
    config_.Validate();

    cerr << endl
         << "Well radius profile (absolute depths, negative downward):";
    if (config_.radius_depths.empty()) {
        cerr << endl << "  uniform, r = " << config_.radius_values.front() << " m";
    } else {
        double from = top_of_well;
        for (size_t i = 0; i < config_.radius_depths.size(); ++i) {
            cerr << endl << "  " << from << " m to " << config_.radius_depths[i]
                 << " m : r = " << config_.radius_values[i] << " m";
            from = config_.radius_depths[i];
        }
        cerr << endl << "  " << from << " m to the toe : r = "
             << config_.radius_values.back() << " m";
    }
    cerr << endl << "  well top = " << top_of_well << " m";

    VO_well_cross_area.resize(number_of_well_elements);
    for (unsigned n = 0; n < number_of_well_elements; n++) {
        const double r = config_.RadiusAtDepth(VO_element_heel_depth[n]);
        VO_well_cross_area[n] = CSMP_PI * r * r;
    }

    // radius of the virtual segments above the mesh top, if any
    const double virtual_top_radius = (config_.virtual_top_radius > 0.)
            ? config_.virtual_top_radius
            : config_.radius_values.front();
    const double crossarea = CSMP_PI * virtual_top_radius * virtual_top_radius;

    VO_D_segment_heel_depth.resize(number_of_well_segments); // NEW
    VO_D_segment_length.resize(number_of_well_segments);
    VO_D_segment_inclination.resize(number_of_well_segments);
    VO_D_well_cross_area.resize(number_of_well_segments);

    VO_well_segment_activated.resize(number_of_well_segments);
    VO_previous_well_segment_activated.resize(number_of_well_segments);

    // Downscaling of lengths and inclinations on sub-segments + extra top virtual
    // segment

    if (number_of_virtual_top_segments > 0.1)
        VO_D_segment_heel_depth[number_of_virtual_top_segments - 1] =
                VO_element_heel_depth[0]; // NEW

    unsigned cumul_index =
            number_of_virtual_top_segments; // we start from segment 1 and will add
    // "virtual" segment 0

    for (unsigned k = 0; k < number_of_well_elements; k++) {
        cerr << endl;
        for (unsigned l = cumul_index;
             l < cumul_index + number_of_segment_divisions[k]; l++) {
            VO_D_segment_length[l] +=
                    VO_element_length[k] / number_of_segment_divisions[k];
            VO_D_segment_inclination[l] +=
                    VO_element_inclination[k]; // Inclination of segments is the same as
            // the element they reside in.
            VO_D_well_cross_area[l] += VO_well_cross_area[k];

            if (number_of_virtual_top_segments > 0.1)
                VO_D_segment_heel_depth[l] +=
                        VO_D_segment_heel_depth[l - 1] -
                        VO_D_segment_length[l] * cos(VO_D_segment_inclination[l]); // NEW

            else {
                // if( cumul_index == 0 )
                if (l == 0) {
                    VO_D_segment_heel_depth[l] +=
                            top_of_well -
                            VO_D_segment_length[l] * cos(VO_D_segment_inclination[l]); // NEW

                }
                else
                    VO_D_segment_heel_depth[l] +=
                            VO_D_segment_heel_depth[l - 1] -
                            VO_D_segment_length[l] * cos(VO_D_segment_inclination[l]); // NEW
            }

            cerr << "Segment index " << l << " has a length of "
                 << VO_D_segment_length[l] << "m, and a heel depth of"
                 << VO_D_segment_heel_depth[l] << "m" << endl;

            index_of_last_segment[k] = l;
        }

        cumul_index += number_of_segment_divisions[k];
    }

    cerr << "index of last segments: " << endl;
    for (unsigned k = 0; k < number_of_well_elements; k++) {
        cerr << index_of_last_segment[k] << ", ";
    }

    cerr << endl << "Virtual top segments: " << endl;
    for (unsigned k = 0; k < number_of_virtual_top_segments; k++) {

        VO_D_segment_length[k] = abs(top_of_well) / number_of_virtual_top_segments;
        VO_D_segment_inclination[k] = 0.;
        VO_D_well_cross_area[k] = crossarea;
        VO_D_segment_heel_depth[k] =
                -abs(top_of_well) / number_of_virtual_top_segments * (k + 1);

        cerr << "Segment index " << k << " has a length of "
             << VO_D_segment_length[k] << "m, and a heel depth of"
             << VO_D_segment_heel_depth[k] << "m" << endl;
    }

    if (number_of_virtual_top_segments < 0.1)
        cerr << " none." << endl;

    // PAUSE();

    // Segment cross-sectional area from the SAME radius profile, and the
    // activation flags. The per-case overrides that used to live here (IDDP2,
    // KD68B, KD45, and the 0.17/0.12/0.11 lithium profile) are gone: they belong
    // in a WellConfiguration, not in the source.
    for (unsigned n = 0; n < number_of_well_segments; n++) {

        VO_well_segment_activated[n] = 1;
        VO_previous_well_segment_activated[n] = 1;

        const double r = config_.RadiusAtDepth(VO_D_segment_heel_depth[n]);
        VO_D_well_cross_area[n] = CSMP_PI * r * r;
    }

    // Well variables
    //// On depth-ordered mesh elements
    VO_interp_well_pressure.resize(number_of_well_elements);
    VO_well_velocity_liquid.resize(number_of_well_elements);
    VO_well_velocity_vapor.resize(number_of_well_elements);
    VO_well_velocity_fluid.resize(number_of_well_elements);
    VO_well_s_enthalpy_liquid.resize(number_of_well_elements);
    VO_well_s_enthalpy_vapor.resize(number_of_well_elements);
    VO_well_s_enthalpy_fluid.resize(number_of_well_elements);
    VO_well_s_enthalpy_halite.resize(number_of_well_elements);
    VO_well_s_enthalpy_bulk.resize(number_of_well_elements);
    VO_well_smf_liquid.resize(number_of_well_elements);
    VO_well_smf_vapor.resize(number_of_well_elements);
    VO_well_fluid_state.resize(number_of_well_elements);
    VO_well_smf_fluid.resize(number_of_well_elements);
    VO_well_smf_bulk.resize(number_of_well_elements);

    VO_well_temperature.resize(number_of_well_elements);
    VO_well_density_liquid.resize(number_of_well_elements);
    VO_well_density_vapor.resize(number_of_well_elements);
    VO_well_density_fluid.resize(number_of_well_elements);
    VO_well_density_halite.resize(number_of_well_elements);
    VO_well_density_bulk.resize(number_of_well_elements);
    VO_well_sat_liquid.resize(number_of_well_elements);
    VO_well_sat_vapor.resize(number_of_well_elements);
    VO_well_sat_halite.resize(number_of_well_elements);
    VO_previous_well_sat_halite.resize(number_of_well_elements);
    VO_well_mf_liquid.resize(number_of_well_elements);
    VO_well_mf_vapor.resize(number_of_well_elements);
    VO_well_mf_halite.resize(number_of_well_elements);
    VO_well_friction_factor.resize(number_of_well_elements);

    VO_well_viscosity_liquid.resize(number_of_well_elements);
    VO_well_viscosity_vapor.resize(number_of_well_elements);
    VO_well_viscosity_fluid.resize(number_of_well_elements);

    VO_well_Re.resize(number_of_well_elements);

    //// Ordered, downscaled to segments, + virtual top segment
    VO_D_well_pressure.resize(number_of_well_segments);
    VO_D_interp_well_pressure.resize(number_of_well_segments);
    VO_D_well_velocity_liquid.resize(number_of_well_segments);
    VO_D_well_velocity_vapor.resize(number_of_well_segments);
    VO_D_well_velocity_fluid.resize(number_of_well_segments);
    VO_D_well_s_enthalpy_liquid.resize(number_of_well_segments);
    VO_D_well_s_enthalpy_vapor.resize(number_of_well_segments);
    VO_D_well_s_enthalpy_fluid.resize(number_of_well_segments);
    VO_D_well_s_enthalpy_halite.resize(number_of_well_segments);
    VO_D_well_s_enthalpy_bulk.resize(number_of_well_segments);
    VO_D_well_smf_liquid.resize(number_of_well_segments);
    VO_D_well_smf_vapor.resize(number_of_well_segments);
    VO_D_well_smf_fluid.resize(number_of_well_segments);
    VO_D_well_smf_bulk.resize(number_of_well_segments);

    VO_D_well_fluid_state.resize(number_of_well_segments);
    VO_D_previous_well_fluid_state.resize(number_of_well_segments);
    VO_D_NR_previous_well_fluid_state.resize(number_of_well_segments);

    VO_D_well_temperature.resize(number_of_well_segments);
    VO_D_well_density_liquid.resize(number_of_well_segments);
    VO_D_well_density_vapor.resize(number_of_well_segments);
    VO_D_well_density_fluid.resize(number_of_well_segments);
    VO_D_well_last_true_density_fluid.resize(number_of_well_segments);
    VO_D_well_density_halite.resize(number_of_well_segments);
    VO_D_well_density_bulk.resize(number_of_well_segments);
    VO_D_well_sat_liquid.resize(number_of_well_segments);
    VO_D_well_sat_vapor.resize(number_of_well_segments);

    VO_D_well_sat_halite.resize(number_of_well_segments);
    VO_D_well_mf_liquid.resize(number_of_well_segments);
    VO_D_well_mf_fluid.resize(number_of_well_segments);
    VO_D_well_mf_vapor.resize(number_of_well_segments);
    VO_D_well_mf_halite.resize(number_of_well_segments);
    VO_D_well_friction_factor.resize(number_of_well_segments);
    // VO_D_well_roughness.                    resize( number_of_well_segments );

    VO_D_well_viscosity_liquid.resize(number_of_well_segments);
    VO_D_well_viscosity_vapor.resize(number_of_well_segments);
    VO_D_well_viscosity_fluid.resize(number_of_well_segments);

    VO_D_well_Re.resize(number_of_well_segments);

    VO_D_well_salt_mass.resize(number_of_well_segments);

    VO_D_well_water_table.resize(number_of_well_segments);
    VO_D_previous_well_water_table.resize(number_of_well_segments);
    VO_D_NR_previous_well_water_table.resize(number_of_well_segments);
    VO_D_well_air_saturation.resize(number_of_well_segments);
    VO_D_previous_well_air_saturation.resize(number_of_well_segments);
    VO_D_NR_previous_well_air_saturation.resize(number_of_well_segments);

    //// Previous version of above
    VO_D_previous_well_pressure.resize(number_of_well_segments);
    VO_D_previous_interp_well_pressure.resize(number_of_well_segments);
    VO_D_previous_well_velocity_liquid.resize(number_of_well_segments);
    VO_D_previous_well_velocity_vapor.resize(number_of_well_segments);
    VO_D_previous_well_velocity_fluid.resize(number_of_well_segments);
    VO_D_previous_well_s_enthalpy_liquid.resize(number_of_well_segments);
    VO_D_previous_well_s_enthalpy_vapor.resize(number_of_well_segments);
    VO_D_previous_well_s_enthalpy_fluid.resize(number_of_well_segments);
    VO_D_previous_well_s_enthalpy_halite.resize(number_of_well_segments);
    VO_D_previous_well_s_enthalpy_bulk.resize(number_of_well_segments);

    VO_D_previous_well_smf_liquid.resize(number_of_well_segments);
    VO_D_previous_well_smf_vapor.resize(number_of_well_segments);
    VO_D_previous_well_smf_fluid.resize(number_of_well_segments);
    VO_D_previous_well_smf_bulk.resize(number_of_well_segments);

    VO_D_previous_well_temperature.resize(number_of_well_segments);
    VO_D_previous_well_density_liquid.resize(number_of_well_segments);
    VO_D_previous_well_density_vapor.resize(number_of_well_segments);
    VO_D_previous_well_density_fluid.resize(number_of_well_segments);
    VO_D_previous_well_density_halite.resize(number_of_well_segments);
    VO_D_previous_well_density_bulk.resize(number_of_well_segments);
    VO_D_previous_well_sat_liquid.resize(number_of_well_segments);
    VO_D_previous_well_sat_vapor.resize(number_of_well_segments);
    VO_D_previous_well_sat_halite.resize(number_of_well_segments);
    VO_D_previous_well_mf_liquid.resize(number_of_well_segments);
    VO_D_previous_well_mf_vapor.resize(number_of_well_segments);
    VO_D_previous_well_mf_halite.resize(number_of_well_segments);
    VO_D_previous_well_friction_factor.resize(number_of_well_segments);
    VO_D_previous_well_cross_area.resize(number_of_well_segments);

    VO_D_previous_well_viscosity_liquid.resize(number_of_well_segments);
    VO_D_previous_well_viscosity_vapor.resize(number_of_well_segments);
    VO_D_previous_well_viscosity_fluid.resize(number_of_well_segments);

    //// NR Previous version of above

    VO_D_NR_previous_well_pressure.resize(number_of_well_segments);
    VO_D_NR_previous_interp_well_pressure.resize(number_of_well_segments);

    VO_D_NR_previous_well_s_enthalpy_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_s_enthalpy_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_s_enthalpy_fluid.resize(number_of_well_segments);
    VO_D_NR_previous_well_s_enthalpy_halite.resize(number_of_well_segments);
    VO_D_NR_previous_well_s_enthalpy_bulk.resize(number_of_well_segments);

    VO_D_NR_previous_well_smf_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_smf_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_smf_fluid.resize(number_of_well_segments);
    VO_D_NR_previous_well_smf_bulk.resize(number_of_well_segments);

    VO_D_NR_previous_well_density_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_density_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_density_fluid.resize(number_of_well_segments);
    VO_D_NR_previous_well_density_halite.resize(number_of_well_segments);
    VO_D_NR_previous_well_density_bulk.resize(number_of_well_segments);
    VO_D_NR_previous_well_sat_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_sat_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_sat_halite.resize(number_of_well_segments);
    VO_D_NR_previous_well_mf_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_mf_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_mf_halite.resize(number_of_well_segments);

    VO_D_NR_previous_well_viscosity_liquid.resize(number_of_well_segments);
    VO_D_NR_previous_well_viscosity_vapor.resize(number_of_well_segments);
    VO_D_NR_previous_well_viscosity_fluid.resize(number_of_well_segments);

    // radial heat
    V_temperature.resize(number_of_well_nodes);
    VO_temperature.resize(number_of_well_nodes);
    VO_D_temperature.resize(number_of_well_segments);
    VO_D_pressure.resize(number_of_well_segments);

    VO_radial_heat.resize(number_of_well_nodes);
    VO_E_radial_heat.resize(number_of_well_elements);
    VO_D_radial_heat.resize(number_of_well_segments);

    // Derivatives for Jacobian
    VO_D_drhofdp.resize(number_of_well_segments);
    VO_D_drhofdhb.resize(number_of_well_segments);
    VO_D_drhofdsmfb.resize(number_of_well_segments);
    VO_D_drhobdp.resize(number_of_well_segments);
    VO_D_drhobdhb.resize(number_of_well_segments);
    VO_D_drhobdsmfb.resize(number_of_well_segments);

    VO_D_dhdp.resize(number_of_well_segments);
    VO_D_dhdsmf.resize(number_of_well_segments);
    VO_D_dsmfdh.resize(number_of_well_segments);
    VO_D_dsmfdp.resize(number_of_well_segments);

    // Well geometry
    V_equivalent_radius.resize(number_of_well_nodes);  // NEW
    VO_equivalent_radius.resize(number_of_well_nodes); // NEW

    // Reservoir variables
    //// On mesh nodes
    V_pressure.resize(number_of_well_nodes);
    V_s_enthalpy_liquid.resize(number_of_well_nodes);
    V_s_enthalpy_vapor.resize(number_of_well_nodes);
    V_density_liquid.resize(number_of_well_nodes);
    V_density_vapor.resize(number_of_well_nodes);
    V_viscosity_liquid.resize(number_of_well_nodes);
    V_viscosity_vapor.resize(number_of_well_nodes);
    V_sat_liquid.resize(number_of_well_nodes);
    V_sat_vapor.resize(number_of_well_nodes);
    V_sat_halite.resize(number_of_well_nodes);
    V_smf_liquid.resize(number_of_well_nodes);
    V_smf_vapor.resize(number_of_well_nodes);
    V_fluid_state.resize(number_of_well_nodes);
    VO_pressure.resize(number_of_well_nodes);
    VO_pressure_at_rest.resize(number_of_well_nodes);
    VO_well_node_pressure.resize(number_of_well_nodes);
    VO_s_enthalpy_liquid.resize(number_of_well_nodes);
    VO_s_enthalpy_vapor.resize(number_of_well_nodes);
    VO_density_liquid.resize(number_of_well_nodes);
    VO_density_vapor.resize(number_of_well_nodes);
    VO_viscosity_liquid.resize(number_of_well_nodes);
    VO_viscosity_vapor.resize(number_of_well_nodes);
    VO_sat_liquid.resize(number_of_well_nodes);
    VO_sat_vapor.resize(number_of_well_nodes);
    VO_sat_halite.resize(number_of_well_nodes);
    VO_smf_liquid.resize(number_of_well_nodes);
    VO_smf_vapor.resize(number_of_well_nodes);
    VO_fluid_state.resize(number_of_well_nodes);

    // Exchange, evaluated at completions / elements's heels, then interpolated to
    // segments if any (element heel correspond to reservoir node)
    V_well_index.resize(number_of_well_nodes);
    V_relperm_visc_liquid.resize(number_of_well_nodes);
    V_relperm_visc_vapor.resize(number_of_well_nodes);
    V_formation_permeability.resize(number_of_well_nodes);
    V_skin.resize(number_of_well_nodes);
    V_bulk_volume.resize(number_of_well_nodes);
    V_well_completion.resize(number_of_well_nodes);
    VO_well_index.resize(number_of_well_nodes);
    VO_well_injectivity.resize(number_of_well_nodes);
    VO_relperm_visc_liquid.resize(number_of_well_nodes);
    VO_relperm_visc_vapor.resize(number_of_well_nodes);
    VO_formation_permeability.resize(number_of_well_nodes);
    VO_drainage_radius.resize(number_of_well_nodes);
    VO_D_drainage_radius.resize(number_of_well_segments);
    VO_skin.resize(number_of_well_nodes);
    VO_bulk_volume.resize(number_of_well_nodes);
    VO_mass_transfer_rate.resize(number_of_well_nodes);
    VO_energy_transfer_rate_l.resize(number_of_well_nodes);
    VO_energy_transfer_rate_v.resize(number_of_well_nodes);
    VO_energy_transfer_rate.resize(number_of_well_nodes);
    VO_salt_mass_transfer_rate.resize(number_of_well_nodes);
    VO_D_mass_transfer_rate.resize(number_of_well_segments);
    VO_D_energy_transfer_rate_l.resize(number_of_well_segments);
    VO_D_energy_transfer_rate_v.resize(number_of_well_segments);
    VO_D_energy_transfer_rate.resize(number_of_well_segments);
    VO_D_salt_mass_transfer_rate.resize(number_of_well_segments);
    VO_D_NR_previous_mass_transfer_rate.resize(number_of_well_segments);
    VO_D_NR_previous_energy_transfer_rate.resize(number_of_well_segments);
    VO_D_NR_previous_salt_mass_transfer_rate.resize(number_of_well_segments);
    VO_D_previous_mass_transfer_rate.resize(number_of_well_segments);
    VO_D_previous_energy_transfer_rate.resize(number_of_well_segments);
    VO_D_previous_salt_mass_transfer_rate.resize(number_of_well_segments);

    VO_D_mass_Inj_source.resize(number_of_well_segments);
    VO_D_energy_Inj_source.resize(number_of_well_segments);

    // N-R
    Residuals.resize(number_of_well_segments * n_of_eq);
    // Sol_increment.                          resize( number_of_well_segments *
    // n_of_eq );
    Sol_increment_pressure.resize(number_of_well_segments * n_of_eq);
    Sol_increment_velocity_fluid.resize(number_of_well_segments * n_of_eq);
    Sol_increment_s_enthalpy_bulk.resize(number_of_well_segments * n_of_eq);
    Sol_increment_smf_bulk.resize(number_of_well_segments * n_of_eq);
    Sol_increment_s_enthalpy_fluid.resize(number_of_well_segments * n_of_eq);
    Sol_increment_smf_fluid.resize(number_of_well_segments * n_of_eq);
    Jacobian_segment_momentum.resize(n_of_eq * n_of_eq);
    Jacobian_segment_continuity.resize(n_of_eq * n_of_eq);
    Jacobian_segment_energy.resize(n_of_eq * n_of_eq);
    Jacobian_segment_salt_mass.resize(n_of_eq * n_of_eq);
    Jacobian.Resize(number_of_well_segments * n_of_eq);

    for (typename std::vector<Node<dim> *>::const_iterator nit =
         model_ref_.Region(well_name).PerimeterNodesBegin();
         nit != model_ref_.Region(well_name).NodesEnd(); ++nit) {
        (*nit)->Status(temperature_key_,
                       ANY); // Set the reservoir node at the top of the well to
        // ANY temperature (used to be DIRICH)
    }

    /////////////////////////
    // Define and store some initial Well variables

    double init_well_velocity_fluid(1.e-8);

    std::fill(VO_D_well_smf_bulk.begin(), VO_D_well_smf_bulk.end(),
              config_.init_bulk_smf);
    std::fill(VO_D_well_velocity_fluid.begin(), VO_D_well_velocity_fluid.end(),
              init_well_velocity_fluid);

    // (VO_skin is filled from the "well skin" node property in
    // Read_reservoir_variables, which runs every step; Initialize_Well seeds that
    // property from the configuration. Do not set VO_skin directly here — it
    // would be overwritten on the next read.)
    std::fill(V_skin.begin(), V_skin.end(), 0.);   // not final

    /////////////////////////
    /*if(with_water_table_displacement) */
    Initialize_Water_Table(); // NEW,
    // DEV+TESTING
    // NEEDED
    Read_reservoir_variables(true);

    Initial_T_equilibration(config_.initial_temperature);
    Velocity_model();

    // if( with_water_table_displacement ) Find_Water_Table();//ADD FEB 2022
    // Compute_Inflow_Outflow2();
    // Compute_Inflow_Outflow();
    // previous_total_rate = total_rate;

    // (The permeability>1e-13 override that used to sit here is gone: completion
    // length now comes from CompletionLengthRule, decided per node by whether the
    // host reservoir element is lower-dimensional. See Initialize_Well.)

    Display();
    Write_results();

    /////////////////////////
}

template <uint32_t dim>
void WellModelPrototype<dim>::Read_reservoir_variables(bool at_initialization) {
    unsigned n(0);
    for (typename std::vector<Node<dim> *>::const_iterator nit =
         model_ref_.Region(well_name).NodesBegin();
         nit != model_ref_.Region(well_name).NodesEnd(); ++nit) {
        V_pressure[n] = ((*nit)->Read(pressure_key_));
        V_s_enthalpy_liquid[n] = ((*nit)->Read(s_enthalpy_liquid_key_));
        V_s_enthalpy_vapor[n] = ((*nit)->Read(s_enthalpy_vapor_key_));
        V_smf_liquid[n] = ((*nit)->Read(smf_liquid_key_));
        V_smf_vapor[n] = ((*nit)->Read(smf_vapor_key_));
        V_fluid_state[n] = ((*nit)->Read(fluid_state_key_));
        V_density_liquid[n] = ((*nit)->Read(density_liquid_key_));
        V_density_vapor[n] = ((*nit)->Read(density_vapor_key_));
        V_viscosity_liquid[n] = ((*nit)->Read(viscosity_liquid_key_));
        V_viscosity_vapor[n] = ((*nit)->Read(viscosity_vapor_key_));
        V_sat_liquid[n] = ((*nit)->Read(sat_liquid_key_));
        V_sat_vapor[n] = ((*nit)->Read(sat_vapor_key_));
        V_sat_halite[n] = ((*nit)->Read(sat_halite_key_));
        V_temperature[n] = ((*nit)->Read(temperature_key_));

        V_relperm_visc_liquid[n] = ((*nit)->Read(rel_perm_visc_liquid_key_));
        V_relperm_visc_vapor[n] = ((*nit)->Read(rel_perm_visc_vapor_key_));
        V_formation_permeability[n] = ((*nit)->Read(formation_permeability_key_));
        V_skin[n] = ((*nit)->Read(well_skin_key_));
        V_bulk_volume[n] = ((*nit)->Read(bulk_volume_key_));

        n++;
    }

    for (unsigned n = 0; n < number_of_well_nodes; n++) {
        VO_pressure[n] = V_pressure[node_ordering_vector[n]];
        VO_s_enthalpy_liquid[n] = V_s_enthalpy_liquid[node_ordering_vector[n]];
        VO_s_enthalpy_vapor[n] = V_s_enthalpy_vapor[node_ordering_vector[n]];
        VO_smf_liquid[n] = V_smf_liquid[node_ordering_vector[n]];
        VO_smf_vapor[n] = V_smf_vapor[node_ordering_vector[n]];
        VO_fluid_state[n] = V_fluid_state[node_ordering_vector[n]];
        VO_density_liquid[n] = V_density_liquid[node_ordering_vector[n]];
        VO_density_vapor[n] = V_density_vapor[node_ordering_vector[n]];
        VO_viscosity_liquid[n] = V_viscosity_liquid[node_ordering_vector[n]];
        VO_viscosity_vapor[n] = V_viscosity_vapor[node_ordering_vector[n]];
        VO_sat_liquid[n] = V_sat_liquid[node_ordering_vector[n]];
        VO_sat_vapor[n] = V_sat_vapor[node_ordering_vector[n]];
        VO_sat_halite[n] = V_sat_halite[node_ordering_vector[n]];
        VO_temperature[n] = V_temperature[node_ordering_vector[n]];

        VO_relperm_visc_liquid[n] = V_relperm_visc_liquid[node_ordering_vector[n]];
        VO_relperm_visc_vapor[n] = V_relperm_visc_vapor[node_ordering_vector[n]];
        VO_formation_permeability[n] =
                V_formation_permeability[node_ordering_vector[n]];

        VO_skin[n] = V_skin[node_ordering_vector[n]];
        VO_bulk_volume[n] = V_bulk_volume[node_ordering_vector[n]];

        // Drainage radius of the reservoir control volume — the same r_e the well
        // index uses, so mass and heat cross the same annulus. MUST come after
        // VO_bulk_volume is assigned, just above: computing it earlier used the
        // previous step's volume, and uninitialised memory on the first pass.
        {
            const double len = (n == 0) ? VO_element_length[0]
                                        : ((n == number_of_well_nodes - 1)
                                           ? VO_element_length[n - 1]
                                           : 0.5 * (VO_element_length[n - 1] + VO_element_length[n]));
            VO_drainage_radius[n] = (len > 0. && VO_bulk_volume[n] > 0.)
                    ? sqrtf(VO_bulk_volume[n] / len / CSMP_PI)
                    : 0.;
        }

        if (at_initialization)
            VO_pressure_at_rest[n] =
                    VO_pressure[n]; // Used for calculating injectivity

        //        // OLD bug: if only vapor in formation, relperm_visc_liquid = 0
        //        and liquid cannot be injected!
    }

    /////////////////////////////////
    // we interpolate reservoir temperature to the segment centers for the
    // calculation of radial heat

    unsigned cumul_index =
            number_of_virtual_top_segments; // we start from segment 1 and will add
    // "virtual" segment 0

    for (unsigned k = 0; k < number_of_well_elements; k++) {
        // cerr<<endl<<"Interpolating temperature from nodes to segment centers";
        // cerr<<"Temperature at index "<<k<<" and "<<k+1<<": "<<VO_temperature [ k
        // ]<<" and "<<VO_temperature [ k + 1 ]<<" interpolated at segments centers:
        // ";
        unsigned m = 1;
        for (unsigned l = cumul_index;
             l < cumul_index + number_of_segment_divisions[k]; l++) {

            VO_D_temperature[l] =
                    VO_temperature[k] + (VO_temperature[k + 1] - VO_temperature[k]) /
                    (2 * number_of_segment_divisions[k]) * m;
            VO_D_pressure[l] =
                    VO_pressure[k] + (VO_pressure[k + 1] - VO_pressure[k]) /
                    (2 * number_of_segment_divisions[k]) * m;

            // Interpolation of the drainage radius to segment centres.
            VO_D_drainage_radius[l] =
                    VO_drainage_radius[k] +
                    (VO_drainage_radius[k + 1] - VO_drainage_radius[k]) /
                    (2 * number_of_segment_divisions[k]) * m;
            // cerr<<endl<<"Segment index "<<l<<" : "<<VO_D_temperature [ l ]<<" ,";
            m += 2;
        }

        cumul_index += number_of_segment_divisions[k];
    }

    // VO_D_temperature [ 0 ] = ( VO_temperature [ 0 ] + top_temperature ) / 2;//
    // top temperature would be the ground surface. !! the length of segment 0
    // must then be such that its head is at the surface !!

    // cerr<<endl<<"Virtual top segments, temperature between
    // "<<top_temperature<<" and "<<VO_temperature [ 0 ];

    //    unsigned m = 1;
    //    for( unsigned k = 0; k < number_of_virtual_top_segments; k++ )
    //    {
    //        VO_D_temperature [ k ] = top_temperature + ( VO_temperature [ 0 ] -
    //        top_temperature ) / ( 2 * number_of_virtual_top_segments ) * m;
    //        cerr<<endl<<"Segment index "<<k<<" : "<<VO_D_temperature [ k ]<<"
    //        ,"; m+=2;
    //    }

    // NEW VERSION:
    for (unsigned k = 0; k < number_of_virtual_top_segments; k++) {
        VO_D_temperature[k] =
                top_temperature +
                (VO_temperature[0] - top_temperature) /
                abs(VO_D_segment_heel_depth[number_of_virtual_top_segments - 1]) *
                (abs(VO_D_segment_heel_depth[k]) - 0.5 * VO_D_segment_length[k]);

        VO_D_pressure[k] =
                top_pressure +
                (VO_pressure[0] - top_pressure) /
                abs(VO_D_segment_heel_depth[number_of_virtual_top_segments - 1]) *
                (abs(VO_D_segment_heel_depth[k]) - 0.5 * VO_D_segment_length[k]);
        // cerr<<endl<<"NEW VERSION Segment index "<<k<<" : "<<VO_D_temperature [ k
        // ]<<" ,";
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Equivalent_Radius() // better way to
// compute well
// equivalent radius??
{
    //! calculates well pressure from the reservoir pressure solved numerically,
    //! using the formula: p_well = p_res +
    // ════════════════════════════════════════════════════════════════════════
    // TRANSMISSIBILITY-BASED EQUIVALENT RADIUS — NOT USED BY THE WELL INDEX
    //
    // Read this before trusting anything below. (Benoit DD/MM/YYYY)
    //
    // 1. The result is output only. This block fills `well equivalent radius`,
    //    V_equivalent_radius and VO_equivalent_radius, and NOTHING reads them.
    //    The live well index (see Reservoir_exchange) computes its own `re` from
    //    the control volume:  re = sqrt(bulk_volume / element_length / PI).
    //    Both formulations are Peaceman's; they differ only in how re is found.
    //
    // 2. Why this one could be better: it weights by finite-element
    //    transmissibilities, so it carries permeability anisotropy and the actual
    //    mesh connectivity, which the purely geometric re does not. On a regular
    //    isotropic mesh the two agree closely; they diverge exactly where the
    //    fault gives strong local anisotropy or the mesh is stretched.
    //
    // 3. Why it is NOT yet valid where you would most want it: `Volume()` returns
    //    an element's own measure — a volume for a tetrahedron but an AREA for a
    //    triangle — and the thickness attribute is not applied (see the comment
    //    in the loop). With lower-dimensional parents present, m^3 and m^2 are
    //    summed into the same accumulator, so a fault element's transmissibility
    //    is under-weighted by its aperture. Multiplying Volume by the element's
    //    `thickness` for non-volume parents would fix it.
    //
    // 4. `k_avg` harmonically averages over ALL parents, including the well's own
    //    line elements, which carry the well-region permeability rather than a
    //    formation value. Harmless today because k_avg is unused; it would matter
    //    if this path were switched on.
    //
    // 5. `angle_mult` is pinned at 1.0 here. The equivalent concept now lives on
    //    WellConfiguration::angle_fraction, where it actually reaches the well
    //    index.
    // ════════════════════════════════════════════════════════════════════════

    //! q*mu/(angle_mult*2*CSMP_PI*k*mt*h)*ln(r_eff/r_well) where q is the mass
    //! rate [kg/s], mu is the viscosity, k is the element permeability, mt is the
    //! total fluid mass in the control volume, h is the length of a well segment,
    //! r_well is the well radius, r_eff is the effective well radius from the
    //! Peaceman's well model extension for FE in 3D r_eff = exp(( sum_i (T_0i ln
    //! r_0i) - 2 CSMP_PI)/sum_i T_0i), sum over all neighboring nodes of the well
    //! node, where T_0i is the transmissivity, and r_0i is the edge length
    //! between the well node and the neighboring node

    ScalarVariable q, mu, mt, p_res, p_well, re;

    // TODO: how to calculate angle_mult, if model is not box shaped?
    // if well node is on the model edge angle_mult=0.25, if well node is on the
    // side boundary angle_mult=0.5, if well node is at the model corner
    // angle_mult = 0.125

    // the well node is surrounded by mesh elements from all sides
    double angle_mult = 1.;

    double log_r_eff, nom, denom, r0j, T0j, k_avg, h;
    DenseMatrix<DM_MIN> B;
    Point<dim> p0, pj;
    size_t nParents;

    unsigned n(0);

    for (typename vector<Node<dim> *>::const_iterator it(
             model_ref_.Region(well_name).NodesBegin());
         it != model_ref_.Region(well_name).NodesEnd(); ++it) {
        nom = 0., denom = 0.;
        k_avg = 0.;
        h = 0.;

        // coordinates of the current well node
        p0 = (*it)->Coordinate();

        // number of parent elements of the current well node
        nParents = (*it)->Parents();

        // for all parent elements
        for (size_t i = 0; i < nParents; i++) {
            // Volume of the current parent element
            double Volume = (*it)->Parent(i)->FE()->Volume();

            // harmonic average of parent element permeabilities is used to calculate
            // permeability on the node
            double k = (*it)->Parent(i)->Read(permeability_key_);
            k_avg += 1. / k;

            // the layer thickness is calculated from half-lengths of well line
            // segments
            if ((*it)->Parent(i)->IsLine()) {
                // for a line element Volume = Length
                h += Volume / 2;
            }
            else
                // currently thickness attribute is not used, and flow in 1D elements is
                // not considered that's why transmissibility is only calculated for
                // volumetric elements
            {
                // matrix of basis/shape functions derivatives wrt to x,y,z; size (dim,
                // Nodes of the element)
                (*it)->Parent(i)->dN(B);

                size_t j0 = (*it)->ParentNodeNumber(i);
                // for all the nodes of each parent element
                for (size_t j = 0U; j < (*it)->Parent(i)->Nodes(); j++) {
                    // if this node is not the one from which we started
                    if (j != j0) {
                        // current parent node coordinate
                        pj = (*it)->Parent(i)->N(j)->Coordinate();
                        // parent element edge length
                        r0j = p0.DistanceTo(pj);

                        // transbissibility calculation for the current parent element edge
                        T0j = 0.;
                        // scalar product of shape function gradients for the current well
                        // node and current parent node
                        for (size_t d = 0; d < dim; d++)
                            T0j += B(d, j0) * B(d, j);

                        // multiplied by -1/2 of the parent element volume

                        T0j *= -Volume / 2;

                        // assembling nominator and denominator of the effective radius
                        // logarithm
                        nom += T0j * log(r0j);
                        denom += T0j;
                    }
                }
            }
        }

        log_r_eff = (nom - angle_mult * 2 * CSMP_PI * h) / denom;
        re() = exp(log_r_eff);

        k_avg = nParents / k_avg;

        // cerr<<"\n new re: "<<re();
        // cerr<<"\n new log_re: "<<log_r_eff;

        (*it)->Store(well_equivalent_radius_key_, re);

        V_equivalent_radius[n] = re();
        n++;
        // cerr<<"\n new k_avg: "<<k_avg;

        // p_well() = p_res() +
        // q()*mu()/(angle_mult*2*CSMP_PI*k_avg*mt()*h)*(log_r_eff - log(r_well));
    }

    for (unsigned n = 0; n < number_of_well_nodes; n++) {
        VO_equivalent_radius[n] = V_equivalent_radius[node_ordering_vector[n]];
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Inflow_Outflow() // inflow performance
// relationship
{
    // VERSION WITH TOTAL RELATIVE MOBILITY OF THE RESERVOIR

    // Compute_Equivalent_Radius();//TEST computing equivalent radius with Alina's
    // code

    bool verbose(true);

    total_rate = 0.;
    true_total_rate = 0.;

    double total_injectivity = 0.;

    Read_reservoir_variables(false);
    Element_average(); // is it still useful here?

    std::fill(VO_mass_transfer_rate.begin(), VO_mass_transfer_rate.end(),
              0.); // [kg/s]
    std::fill(VO_energy_transfer_rate_l.begin(), VO_energy_transfer_rate_l.end(),
              0.); // [J/s]
    std::fill(VO_energy_transfer_rate_v.begin(), VO_energy_transfer_rate_v.end(),
              0.); // [J/s]
    std::fill(VO_energy_transfer_rate.begin(), VO_energy_transfer_rate.end(),
              0.); // [J/s]
    std::fill(VO_salt_mass_transfer_rate.begin(),
              VO_salt_mass_transfer_rate.end(), 0.); // [kg/s]

    std::fill(VO_D_mass_transfer_rate.begin(), VO_D_mass_transfer_rate.end(),
              0.); // [kg/s]
    std::fill(VO_D_energy_transfer_rate_l.begin(),
              VO_D_energy_transfer_rate_l.end(), 0.); // [J/s]
    std::fill(VO_D_energy_transfer_rate_v.begin(),
              VO_D_energy_transfer_rate_v.end(), 0.); // [J/s]
    std::fill(VO_D_energy_transfer_rate.begin(), VO_D_energy_transfer_rate.end(),
              0.); // [J/s]
    std::fill(VO_D_salt_mass_transfer_rate.begin(),
              VO_D_salt_mass_transfer_rate.end(), 0.); // [kg/s]

    for (unsigned n = 1; n < number_of_well_nodes;
         n++) // we skip n = 0 because surely there will be no exchange term on
        // top virtual segment
    {
        if (shut_down_clogged_segment) {
            //            if( VO_well_completion [ n ] > 0 &&
            //            VO_well_segment_activated [ n ] > 0 && VO_well_sat_halite [
            //            n - 1 ] > 0.5 )
            //                VO_well_segment_activated [ n ] = 0;

            //            if( VO_well_completion [ n ] > 0 &&
            //            VO_well_segment_activated [ n ] == 0 && VO_well_sat_halite [
            //            n - 1 ] < 0.45 )
            //                VO_well_segment_activated [ n ] = 1;
        }

        double node_depth = 1.e8;
        if (number_of_virtual_top_segments > 0.1)
            node_depth = VO_element_heel_depth[n];
        else
            node_depth = VO_element_heel_depth[n - 1];

        if (VO_well_completion[n] > 0 && (water_table_depth > node_depth)) {

            unsigned K =
                    index_of_last_segment[n - 1]; // The index of the segment associated
            // with this completion;

            double mass_transfer_rate_l(0.), mass_transfer_rate_v(0.);
            double relative_perm_vapor(0.), relative_perm_liquid(0.),
                    effective_sat_vapor(0.), effective_sat_liquid(0.);
            double relative_mobility_liquid(0.), relative_mobility_vapor(0.);
            double pressure_difference(VO_pressure[n] - VO_well_node_pressure[n]);

            double pressure_difference_with_reservoir_at_rest(
                        VO_pressure_at_rest[n] - VO_well_node_pressure[n]);

            double exchange_rate_total(0.), exchange_rate_liquid(0.),
                    exchange_rate_vapor(0.); // m3/s

            // Calculation of the well index:
            double cross_area =
                    VO_bulk_volume[n] /
                    ((VO_element_length[n - 1] + VO_element_length[n]) / 2.);

            if (n == number_of_well_nodes - 1) {
                cross_area = VO_bulk_volume[n] / ((VO_element_length[n - 1]));
            }
            // Careful!!, there is no VO_element_length [ n ] for n = last node
            // there cannot be a completion on the last node!! Unless we use a special
            // mesh...

            double re =
                    sqrtf(cross_area /
                          CSMP_PI); // this equation seems to provide the right result..!

            double rw =
                    sqrtf(VO_well_cross_area[n - 1] * (1. - VO_well_sat_halite[n - 1]) /
                          CSMP_PI); // DOES ADDING HALITE IS RIGHT?

            //                        cerr<<"Node "<<n<<" is a completion node. The
            //                        halite saturation of the element attached to
            //                        this node ( element: "<<n-1<<" ) is "<<
            //                              VO_well_sat_halite [ n - 1 ]<<". The rate
            //                              factor is "<<VO_well_rate_factor [ n
            //                              ]<<endl;
            // the well radius at the element's heel is calculated based on the
            // element crossarea

            // Peaceman well index. `angle_fraction` (per well, from the
            // configuration) is the share of the full 2*PI radial flow angle
            // available: 1 for an interior well — the normal case — and less only
            // for a well lying on a symmetry plane.
            // It is NOT the same as the control-volume truncation, which `re`
            // already accounts for through bulk_volume — the two are separate
            // effects and correcting one does not correct the other.
            VO_well_index[n] = VO_angle_fraction[n] * 2 * CSMP_PI *
                    VO_completion_length[n] *
                    VO_formation_permeability[n] /
                    (logf(re / rw) + VO_skin[n]);

            // The "formation permeability" is an extrapolation at the nodes of the
            // elements permeability (weighted by volume so as to not account for the
            // well elements low permeability), we could maybe use something else...
            // 2 * CSMP_PI is the angle of the well connected to the reservoir (it can
            // be less)

            // Calculation of exchange with the Peaceman formula (phase by phase),
            // assumes radial steady flow around well
            // - There is no gravity terms since they are not necessary with the CVFEM
            // formulation
            // - This formulation do not account for permeability anisotropy

            // index_of_last_segment [ n - 1 ] : we use the well segment values,
            // segment which heel coincide with the mesh node Should we use well
            // values interpolated to nodes??? Then the completion could be any length
            // centered on the node... Although I suppose that the completion length
            // cannot be larger than the sum of half of surrounding segments/elements
            // lengths Also, keep in mind: the source term will only be applied to ONE
            // segment, the one which heel coincide with the mesh node. Also: There is
            // no capillary pressure

            effective_sat_liquid =
                    VO_sat_liquid[n] /
                    (1 - VO_sat_halite[n]); // WILL CRASH IF sat_halite = 1
            effective_sat_vapor = VO_sat_vapor[n] / (1 - VO_sat_halite[n]);

            // Hard coded residual saturation liquid of 30%!!!!
            relative_perm_liquid = 1. / 0.7 * effective_sat_liquid - 0.3 / 0.7;
            relative_perm_vapor = 1. / 0.7 * effective_sat_vapor;

            if (effective_sat_liquid < 0.3) {
                relative_perm_liquid = 0.;
                relative_perm_vapor = 1.;
            }

            if (effective_sat_liquid > 0.) {
                relative_mobility_liquid =
                        relative_perm_liquid / VO_viscosity_liquid[n];
            }

            if (effective_sat_vapor > 0.) {
                relative_mobility_vapor = relative_perm_vapor / VO_viscosity_vapor[n];
            }

            // Here we assumed Kr = S linear.
            // Here we DO NOT upwind mobility, we always use the reservoir's.

            // Since we do not upwind mobility we use a total relative permeability
            // formulation for injection. Otherwise we might end up having no-flow
            // problems ( i.e. when there is only liquid in the reservoir and gas in
            // the well or the opposite ) For example if the reservoir is full of gas,
            // mobility_l = 0 therefore a liquid filled well would not inject anything

            exchange_rate_liquid = VO_well_index[n] * relative_mobility_liquid *
                    pressure_difference *
                    VO_well_rate_factor[n]; // Up to here we have m3/s
            exchange_rate_vapor = VO_well_index[n] * relative_mobility_vapor *
                    pressure_difference *
                    VO_well_rate_factor[n]; // Up to here we have m3/s
            exchange_rate_total = exchange_rate_liquid + exchange_rate_vapor;

            // Upwinding density, energy and salt (Reminder: mobility is not upwinded)
            if (pressure_difference > 0.) {
                // PRODUCTION (INFLOW), done phase by phase

                if (effective_sat_liquid < 0.3) {
                    mass_transfer_rate_l = 0.;
                    mass_transfer_rate_v = exchange_rate_vapor *
                            VO_density_vapor[n]; // Now we have kg/s vapor
                }

                else {
                    mass_transfer_rate_l =
                            exchange_rate_liquid *
                            VO_density_liquid[n]; // Now we have kg/s liquid
                    mass_transfer_rate_v = exchange_rate_vapor *
                            VO_density_vapor[n]; // Now we have kg/s vapor
                }

                VO_mass_transfer_rate[n] = mass_transfer_rate_l + mass_transfer_rate_v;

                VO_well_injectivity[n] =
                        abs(VO_mass_transfer_rate[n] /
                            pressure_difference_with_reservoir_at_rest);
                total_injectivity += VO_well_injectivity[n];

                VO_energy_transfer_rate_l[n] =
                        mass_transfer_rate_l * VO_s_enthalpy_liquid[n];
                VO_energy_transfer_rate_v[n] =
                        mass_transfer_rate_v * VO_s_enthalpy_vapor[n];

                VO_energy_transfer_rate[n] =
                        VO_energy_transfer_rate_l[n] + VO_energy_transfer_rate_v[n];

                // if(with_salt)
                {
                    VO_salt_mass_transfer_rate[n] =
                            mass_transfer_rate_l * VO_smf_liquid[n];
                    VO_salt_mass_transfer_rate[n] +=
                            mass_transfer_rate_v * VO_smf_vapor[n];
                }
            }

            else {
                // INJECTION (OUTFLOW), only caring about the total relative mobility of
                // the reservoir

                mass_transfer_rate_l =
                        exchange_rate_total * VO_D_well_density_liquid[K] *
                        VO_D_well_sat_liquid[K] /
                        (VO_D_well_sat_liquid[K] +
                         VO_D_well_sat_vapor[K]); // Now we have kg/s liquid
                mass_transfer_rate_v =
                        exchange_rate_total * VO_D_well_density_vapor[K] *
                        VO_D_well_sat_vapor[K] /
                        (VO_D_well_sat_liquid[K] +
                         VO_D_well_sat_vapor[K]); // Now we have kg/s vapor
                // IMPORTANT, IN CASE THERE IS "IMMOBILE" HALITE, WE NEED TO NORMALIZE
                // THE SATURATIONS!

                VO_mass_transfer_rate[n] = mass_transfer_rate_l + mass_transfer_rate_v;

                VO_well_injectivity[n] =
                        abs(VO_mass_transfer_rate[n] /
                            pressure_difference_with_reservoir_at_rest);
                total_injectivity += VO_well_injectivity[n];

                VO_energy_transfer_rate_l[n] =
                        mass_transfer_rate_l * VO_D_well_s_enthalpy_liquid[K];
                VO_energy_transfer_rate_v[n] =
                        mass_transfer_rate_v * VO_D_well_s_enthalpy_vapor[K];

                VO_energy_transfer_rate[n] =
                        VO_energy_transfer_rate_l[n] + VO_energy_transfer_rate_v[n];

                // if(with_salt)
                {
                    VO_salt_mass_transfer_rate[n] =
                            mass_transfer_rate_l * VO_D_well_smf_liquid[K];
                    VO_salt_mass_transfer_rate[n] +=
                            mass_transfer_rate_v * VO_D_well_smf_vapor[K];
                }
            }

            total_rate += mass_transfer_rate_l + mass_transfer_rate_v;

            if (total_injectivity > 25) {
                cerr << endl << "**********************************";
                cerr << endl << "At node: " << n;
                cerr << endl
                     << "P res: " << VO_pressure[n]
                        << ", P well: " << VO_well_node_pressure[n]
                           << ", P diff: " << pressure_difference;
                cerr << endl
                     << "P res at rest: " << VO_pressure_at_rest[n]
                        << ", P diff at rest: "
                        << pressure_difference_with_reservoir_at_rest;
                cerr << endl
                     << "Cross area: " << cross_area << ", re: " << re << ", rw: " << rw
                     << ", well index: " << VO_well_index[n];
                cerr << endl
                     << "Effective sat vapor: " << effective_sat_vapor
                     << ", liquid: " << effective_sat_liquid;
                cerr << endl
                     << "Relative perm vapor: " << relative_perm_vapor
                     << ", liquid: " << relative_perm_liquid;
                cerr << endl << "Well rate factor: " << VO_well_rate_factor[n];
                cerr << endl
                     << "Exchange rate vapor: " << exchange_rate_vapor
                     << ", liquid: " << exchange_rate_liquid
                     << ", total: " << exchange_rate_total;
                cerr << endl
                     << "Mass transfer rate vapor: " << mass_transfer_rate_v
                     << ", liquid: " << mass_transfer_rate_l
                     << ", total: " << VO_mass_transfer_rate[n];
                cerr << endl
                     << "Energy transfer rate vapor: " << VO_energy_transfer_rate_v[n]
                        << ", liquid: " << VO_energy_transfer_rate_l[n]
                           << ", total: " << VO_energy_transfer_rate[n];
                cerr << endl << "Well injectivity: " << VO_well_injectivity[n] * 100000;
                cerr << endl << "Total injectivity: " << total_injectivity * 100000;
                cerr << endl << "**********************************";
                // PAUSE(); PAUSE();
            }
        } // end completions
    } // end well loop

    for (unsigned n = 1; n < number_of_well_nodes;
         n++) // we skip n = 0 because surely there will be no exchange term on
        // top virtual segment
    {
        VO_well_injectivity[n] = total_injectivity * 100000; // convert to "kg/s/bar
    }

    true_total_rate = total_rate;

    //        cerr << endl << "\033[31m Mass transfer rate [kg/s]              :
    //        "; for (unsigned n = 0; n < number_of_well_nodes; n++)
    //            cerr << VO_mass_transfer_rate[ n ] << ", ";
    //        cerr<<endl<<" Well Index : ";
    //        for (unsigned n = 0; n < number_of_well_nodes; n++)
    //            cerr << VO_well_index [ n ] << ", ";
    //            cerr<<endl<< " Energy transfer rate [J/s]              : ";
    //            for (unsigned n = 0; n < number_of_well_nodes; n++)
    //                cerr << VO_energy_transfer_rate[ n ] << ", ";
    //            cerr<<endl<< " Salt transfer rate [Kg/s]              : ";
    //            for (unsigned n = 0; n < number_of_well_nodes; n++)
    //                cerr << VO_salt_mass_transfer_rate[ n ] << ", ";

    //    cerr<<"\033[0m";

    if (spread_source) {
        // Distributing evenly to segments composing the element

        unsigned cumul_index = number_of_virtual_top_segments;

        for (unsigned k = 0; k < number_of_well_elements; k++) {
            // cerr<<endl;
            double element_fluid_mass_rate, element_e_rate, element_liquid_e_rate,
                    element_vapor_e_rate, element_salt_mass_rate;

            if (k + 1 == number_of_well_nodes - 1) {
                element_fluid_mass_rate =
                        VO_mass_transfer_rate[k] / 2 + VO_mass_transfer_rate[k + 1];
                element_e_rate =
                        VO_energy_transfer_rate[k] / 2 + VO_energy_transfer_rate[k + 1];
                element_liquid_e_rate =
                        VO_energy_transfer_rate_l[k] / 2 + VO_energy_transfer_rate_l[k + 1];
                element_vapor_e_rate =
                        VO_energy_transfer_rate_v[k] / 2 + VO_energy_transfer_rate_v[k + 1];
                element_salt_mass_rate = VO_salt_mass_transfer_rate[k] / 2 +
                        VO_salt_mass_transfer_rate[k + 1];

            }
            else {
                element_fluid_mass_rate =
                        (VO_mass_transfer_rate[k] + VO_mass_transfer_rate[k + 1]) / 2.;
                element_e_rate =
                        (VO_energy_transfer_rate[k] + VO_energy_transfer_rate[k + 1]) / 2.;
                element_liquid_e_rate =
                        (VO_energy_transfer_rate_l[k] + VO_energy_transfer_rate_l[k + 1]) /
                        2.;
                element_vapor_e_rate =
                        (VO_energy_transfer_rate_v[k] + VO_energy_transfer_rate_v[k + 1]) /
                        2.;
                element_salt_mass_rate = (VO_salt_mass_transfer_rate[k] +
                                          VO_salt_mass_transfer_rate[k + 1]) /
                        2.;
            }
            //                cerr<<endl<<"Element index "<<k<<" is bounded by nodes
            //                "<<k<<" and "<<k+1; cerr<<endl<<"VO_mass_transfer_rate
            //                [ k ] = "<<VO_mass_transfer_rate      [ k ]<<",
            //                VO_mass_transfer_rate      [ k + 1 ] =
            //                "<<VO_mass_transfer_rate      [ k + 1 ]; cerr<<endl<<"
            //                The element has a mass_transfer_rate of
            //                "<<element_fluid_mass_rate;

            for (unsigned l = cumul_index;
                 l < cumul_index + number_of_segment_divisions[k]; l++) {
                VO_D_mass_transfer_rate[l] =
                        element_fluid_mass_rate / number_of_segment_divisions[k];
                VO_D_energy_transfer_rate[l] =
                        element_e_rate / number_of_segment_divisions[k];
                VO_D_energy_transfer_rate_l[l] =
                        element_liquid_e_rate / number_of_segment_divisions[k];
                VO_D_energy_transfer_rate_v[l] =
                        element_vapor_e_rate / number_of_segment_divisions[k];
                VO_D_salt_mass_transfer_rate[l] =
                        element_salt_mass_rate / number_of_segment_divisions[k];

                // cerr<<endl<<"Segment index "<<l<<" has a mass_transfer_rate of
                // "<<VO_D_mass_transfer_rate      [ l ];

                index_of_last_segment[k] = l;
            }

            cumul_index += number_of_segment_divisions[k];
        }
    }

    else {
        // Attributing to last segment of well element:
        for (unsigned n = 1; n < number_of_well_nodes; n++) {
            VO_D_mass_transfer_rate[index_of_last_segment[n - 1]] =
                    VO_mass_transfer_rate[n];
            VO_D_energy_transfer_rate_l[index_of_last_segment[n - 1]] =
                    VO_energy_transfer_rate_l[n];
            VO_D_energy_transfer_rate_v[index_of_last_segment[n - 1]] =
                    VO_energy_transfer_rate_v[n];
            VO_D_energy_transfer_rate[index_of_last_segment[n - 1]] =
                    VO_energy_transfer_rate[n];
            VO_D_salt_mass_transfer_rate[index_of_last_segment[n - 1]] =
                    VO_salt_mass_transfer_rate[n];
        }
    }

    ////////////////////////////////

    unsigned n = 0; // Necessary to store now as it is used in the external
    // Well-Reservoir convergence loop
    for (typename std::vector<Node<dim> *>::const_iterator nit =
         model_ref_.Region(well_name).NodesBegin();
         nit != model_ref_.Region(well_name).NodesEnd(); ++nit) {
        well_injectivity() =
                VO_well_injectivity[node_deordering_vector[n]]; // convert to "kg/s/bar
        // for output
        mass_transfer_rate() = VO_mass_transfer_rate[node_deordering_vector[n]];
        energy_transfer_rate_l() =
                VO_energy_transfer_rate_l[node_deordering_vector[n]];
        energy_transfer_rate_v() =
                VO_energy_transfer_rate_v[node_deordering_vector[n]];
        energy_transfer_rate() = VO_energy_transfer_rate[node_deordering_vector[n]];
        salt_mass_transfer_rate() =
                VO_salt_mass_transfer_rate[node_deordering_vector[n]];
        well_completion() = VO_well_completion[node_deordering_vector[n]];
        well_completion_length() = VO_completion_length[node_deordering_vector[n]];
        well_index() = VO_well_index[node_deordering_vector[n]];

        (*nit)->Store(well_injectivity_key_, well_injectivity);
        (*nit)->Store(mass_transfer_rate_key_, mass_transfer_rate);
        (*nit)->Store(energy_transfer_rate_l_key_, energy_transfer_rate_l);
        (*nit)->Store(energy_transfer_rate_v_key_, energy_transfer_rate_v);
        // (*nit)->Store( energy_transfer_rate_key_ ,    energy_transfer_rate );
        (*nit)->Store(salt_mass_transfer_rate_key_, salt_mass_transfer_rate);
        (*nit)->Store(well_completion_key_,
                      well_completion); // never updated after initialization,
        // should store in initialization?
        (*nit)->Store(well_completion_length_key_, well_completion_length);
        (*nit)->Store(well_index_key_, well_index);
        n++;
    }
}

/** Conductive heat exchange between the well and the formation, Q_loss [W] per
    segment.

    ── Known limitations, for whoever improves this next (Benoit DD/MM/YYYY) ──

    1. kappa is a single configured number, not the local formation conductivity.
       To make it follow the model: create a "nodal thermal conductivity" property
       and fill it with InterpolateCellToNodeProperty from the element values — the
       same trick "nodal permeability" already uses — then read it here per node
       and interpolate to segments alongside the drainage radius. Worth doing if a
       well crosses units with contrasting conductivity; unnecessary while the
       transient run uses a uniform value.

    2. ShapeFactor assumes the reservoir node temperature represents the formation
       at the drainage radius, and that conduction near the well is radial and
       steady within a timestep. Both are the same assumptions the Peaceman well
       index already makes for mass, which is the point — but they are assumptions,
       and they weaken on a coarse mesh or over a very short timestep.

    3. Neither model represents the casing, the cement, or the temperature drop
       across them; r_w is the flowing radius, whereas the physical interface is
       the cement-formation contact. The original code carried a comment noting
       this and it remains true.

    4. Ramey is kept for meshes that do not resolve the near-well formation, but
       it is fed `dt` rather than the time since the well started flowing. That is
       wrong in principle — see RadialHeatModel — and is why it returns zero below
       its validity limit instead of a negative heat loss. If Ramey is ever needed
       in earnest, thread an elapsed-time variable through and use it.
 Two formulations, selected by WellConfiguration::radial_heat_model —
    see RadialHeatModel for why ShapeFactor is the default.

    Both are of the form  Q = 2 pi kappa (T_well - T_formation) / D  [W/m], then
    multiplied by the segment length. Only the denominator D differs:

      ShapeFactor : D = ln(r_e / r_w)                  steady, mesh-based
      Ramey       : D = ln(2 sqrt(alpha t) / r_w) - 0.29   transient line source

    kappa is the FORMATION THERMAL CONDUCTIVITY [W/(m.K)] — Eq. (6) of
    Lamy-Chappuis et al. (2022). It is read from the reservoir by default; the
    variable here was previously named `overall_heat_transfer_coeff`, carried a
    comment guessing at W/(m2.K), and was hardcoded to 8 against the 2 W/(m.K)
    the model actually uses. (Benoit DD/MM/YYYY) */
template <uint32_t dim> void WellModelPrototype<dim>::Compute_radial_heat() {
    double full_well_radius(0.);
    double denominator(0.);

    // Compute radial heat terms on segments, using interpolated reservoir
    // temperature

    // cerr<<endl<<"Computing radial heat at node: ";

    std::fill(VO_D_radial_heat.begin(), VO_D_radial_heat.end(), 0.);

    // (A `dt > 1.` guard used to sit here, to dodge the Ramey function going
    // negative at small timesteps. The denominator check below covers that case
    // properly, and ShapeFactor has no time dependence at all.)
    if (with_radial_heat)
    {

        for (unsigned n = 0; n < number_of_well_segments; n++) {
            if (VO_D_well_water_table[n] == 1 && VO_D_well_temperature[n] < 1.e-3) {
                cerr << endl
                     << "In Compute_radial_heat, segment " << n
                     << " is filling up, its temperature is "
                     << VO_D_well_temperature[n];

                int k = n;
                while (VO_D_well_temperature[k] < 1.e-3) {
                    k = k + 1;
                    cerr << endl
                         << "We set it to the temperature: " << VO_D_well_temperature[k]
                            << ", on segment " << k
                            << ", for the computation of radial heat";

                    VO_D_well_temperature[n] = VO_D_well_temperature[k];
                }
                // PAUSE();
            }

            full_well_radius = sqrtf(VO_D_well_cross_area[n] / CSMP_PI);

            // Formation thermal conductivity [W/(m.K)], from the configuration.
            //
            // NOT read from the reservoir: "thermal conductivity" is an ELEMENT
            // property, and this model only has nodal values available (which is
            // why it reads "nodal permeability" rather than "permeability"). There
            // is no nodal thermal conductivity in the property database, so adding
            // one — interpolated from the element values the way nodal permeability
            // is — would be the way to make this follow the model automatically.
            // (Benoit DD/MM/YYYY)
            const double kappa = config_.formation_thermal_conductivity;

            if (config_.radial_heat_model == RadialHeatModel::ShapeFactor) {
                // Steady shape factor to the first reservoir control volume, the
                // thermal analogue of the Peaceman well index. r_e > r_w always, so
                // the denominator is positive and needs no guard.
                const double r_e = VO_D_drainage_radius[n];
                denominator = (r_e > full_well_radius)
                        ? log(r_e / full_well_radius)
                        : 0.;
            }
            else {
                // Ramey (1962), Eq. (6) of Lamy-Chappuis et al. (2022). t is the
                // time the formation has been exchanging heat with the well; dt is
                // used here, which is only meaningful when a step is long compared
                // with alpha t / r_w^2 ~ 1. Below that the function is negative and
                // the exchange is suppressed rather than reversed.
                denominator =
                        log((2. * sqrt(config_.formation_thermal_diffusivity * dt)) /
                            full_well_radius) - 0.29;
                if (denominator <= 0.) denominator = 0.;
            }

            if (denominator <= 0.) {           // outside the model's validity
                VO_D_radial_heat[n] = 0.;
                continue;
            }

            VO_D_radial_heat[n] = 2 * CSMP_PI * kappa *
                    (VO_D_well_temperature[n] - VO_D_temperature[n]); // W/m

            VO_D_radial_heat[n] /= denominator;

            VO_D_radial_heat[n] *= VO_D_segment_length[n];            // W

            // No conduction through an air-filled segment.
            VO_D_radial_heat[n] *= (1 - VO_D_well_air_saturation[n]);

            // cerr<<endl<<"Radial heat at index "<<n<<" is "<<VO_D_radial_heat [ n ];
        }

        // VALUE ON SEG 0 IS NOT ACTUALLY USED IN CURRENT VERSION SINCE FOR WT_INDEX
        // = 1: H=H_NEXT_SEG
    }

    // Adding the segment heat terms to compute the "interpolated" nodal ones
    //
    std::fill(VO_radial_heat.begin(), VO_radial_heat.end(), 0.);
    std::fill(VO_E_radial_heat.begin(), VO_E_radial_heat.end(), 0.);

    unsigned cumul_index =
            number_of_virtual_top_segments; // we start from segment:
    // number_of_virtual_top_segments
    // cerr<<endl<<"Computing radial heat on elements";
    for (unsigned k = 0; k < number_of_well_elements; k++) {
        // cerr<<endl<<"Summing for index "<<k;
        for (unsigned l = cumul_index;
             l < cumul_index + number_of_segment_divisions[k]; l++) {
            // cerr<<endl<<"index "<<l<<" value "<<VO_D_radial_heat [ l ];
            VO_E_radial_heat[k] += VO_D_radial_heat[l];
        }
        // cerr<<endl<<" Sum: "<<VO_E_radial_heat [ k ];
        cumul_index += number_of_segment_divisions[k];
    }

    // cerr<<endl<<"Contribution from top virtual segments to node 0 :";
    double top_virtual_segments_radial_heat(0.);
    for (unsigned k = 0; k < number_of_virtual_top_segments; k++) {
        // cerr<<endl<<"index "<<k<<" value "<<VO_D_radial_heat [ k ];
        top_virtual_segments_radial_heat += VO_D_radial_heat[k];
    }
    // cerr<<endl<<" Sum: "<<top_virtual_segments_radial_heat;

    // cerr<<endl<<"Interpolate on nodes";
    for (unsigned k = 0; k < number_of_well_nodes; k++) {
        if (k < number_of_well_nodes - 1) {
            VO_radial_heat[k] += VO_E_radial_heat[k] / 2;
        }

        if (k > 0) {
            VO_radial_heat[k] += VO_E_radial_heat[k - 1] / 2;
        }

        if (k == 0) {
            VO_radial_heat[k] += top_virtual_segments_radial_heat / 2;
        }

        // cerr<<endl<<"Index "<<k<<" value "<<VO_radial_heat [ k ];
    }

    VO_radial_heat[0] = 0.; // If the top reservoir node has Dirichlet
    // temperature, it is pointless to have a heat term !!
    // attempt to remove DIRICH from this node !!

    // Order and store
    unsigned n = 0;
    for (typename std::vector<Node<dim> *>::const_iterator nit =
         model_ref_.Region(well_name).NodesBegin();
         nit != model_ref_.Region(well_name).NodesEnd(); ++nit) {
        radial_heat() = VO_radial_heat[node_deordering_vector[n]];
        (*nit)->Store(radial_heat_key_, radial_heat);
        n++;
    }
    // We might need to change the CVFEM procedure if we want to iteratively
    // converge the Well-Res temperatures using the heat terms We might also
    // simply do explicit coupling...
}

template <uint32_t dim> void WellModelPrototype<dim>::T_equilibration() {
    // cerr<<endl<<"Thermal equilibration";

    thermal_eq_crashed_this_NR = false;

    for (unsigned k = 0; k < 1;
         k++) // DO TWICE SINCE AT FIRST ITERATION IT WILL USE NR OLD DENSITIES
        // and mass fractions!!! DISABLED
    {
        // for( unsigned n = 0; n < number_of_well_segments; n++ )
        for (int n = number_of_well_segments - 1; n >= 0; --n) {

            bool verbose(false);
            // if( n<3) verbose = true;
            // if( (n>=wt_index-1 && n<wt_index + 2 ) /*or n<4*/) verbose = true;

            if (VO_D_well_water_table[n] == 0) { // segment filled with air
                // cerr<<endl<<endl<<"Segment: "<<n<<" filled with air";

                //! Might want to allow halite to remain!!!!!

                VO_D_well_temperature[n] = 0.; // We could use air temp: reservoir temp?

                VO_D_well_sat_liquid[n] = 0.;
                VO_D_well_sat_vapor[n] = 0.;
                VO_D_well_sat_halite[n] = 0.;

                VO_D_well_mf_liquid[n] = 0.;
                VO_D_well_mf_fluid[n] = 0.;
                VO_D_well_mf_vapor[n] = 0.;
                VO_D_well_mf_halite[n] = 0.; //!

                VO_D_well_density_liquid[n] = 0.;
                VO_D_well_density_vapor[n] = 0.;
                VO_D_well_density_halite[n] = 0.; //!
                VO_D_well_density_fluid[n] = 0.;
                VO_D_well_density_bulk[n] = 0.; //!

                VO_D_well_s_enthalpy_liquid[n] = 0.;
                VO_D_well_s_enthalpy_vapor[n] = 0.;
                VO_D_well_s_enthalpy_halite[n] = 0.; //!
                VO_D_well_s_enthalpy_fluid[n] = 0.;
                VO_D_well_s_enthalpy_bulk[n] = 0.; //!

                VO_D_well_viscosity_liquid[n] = 0.;
                VO_D_well_viscosity_vapor[n] = 0.;
                VO_D_well_viscosity_fluid[n] = 0.;

                VO_D_well_smf_liquid[n] = 0.;
                VO_D_well_smf_vapor[n] = 0.;
                VO_D_well_smf_fluid[n] = 0.;
                VO_D_well_smf_bulk[n] = 0.; //!

                VO_D_well_salt_mass[n] = 0.; //!

                VO_D_well_fluid_state[n] = 0;
            }

            else { // segment is at least partially filled with water

                if (!(VO_D_well_s_enthalpy_fluid[n] > 0.)) {
                    cerr << endl
                         << "Segment: " << n
                         << " filled up or is the wt segment, using properties from "
                            "below the water table segment in the previous NR iteration"
                         << endl;

                    VO_D_well_smf_bulk[n] = smf_ = smf_below_wt;
                    wt_ = smf_ * 100.;
                    x_ = Weight2XNaCl(wt_);
                    t_ = tp_ = std::max(10., temp_below_wt);

                    p_ = VO_D_well_pressure[n];
                    if (use_center_pressure or n == 0)
                        p_ = VO_D_interp_well_pressure[n]; // FORCED FOR  n = 0!

                    VO_D_well_s_enthalpy_bulk[n] = h_fluid_ = enthalpy_below_wt;
                }

                else {
                    smf_ = VO_D_well_smf_bulk[n];
                    wt_ = smf_ * 100.;
                    x_ = Weight2XNaCl(wt_);
                    tp_ = std::max(VO_D_previous_well_temperature[n], 10.);

                    p_ = VO_D_well_pressure[n];
                    if (use_center_pressure or n == 0)
                        p_ = VO_D_interp_well_pressure[n]; // FORCED FOR  n = 0!

                    h_fluid_ = VO_D_well_s_enthalpy_bulk[n];
                }

                m_fluid_ = 1000.;

                H_current_ = m_fluid_ * h_fluid_;

                thermal_eq_succeeded =
                        fluid_equilibrator.ThreePhaseProperties(Bulk, Liquid, Vapor, Salt);

                // Here we update the thermodynamic properties
                VO_D_well_temperature[n] = Bulk.t;

                if (verbose) {
                    cerr << endl
                         << "Segment " << n << " T-eq............................"
                         << endl;
                    //                    cerr<<endl<<"Updating Segment: "<<n<<"
                    //                          <<", press_: "<<VO_D_well_pressure [ n
                    //                          ]<<", previous_press_:
                    //                          "<<VO_D_previous_well_pressure [ n ]
                    //                            <<", enthalpy_fluid_: "<<h_fluid_<<",
                    //                            previous_enthalpy_fluid_:
                    //                            "<<VO_D_previous_well_s_enthalpy_fluid [
                    //                            n ];

                    cerr << endl
                         << "VO_D_well_temperature [ n ]: " << VO_D_well_temperature[n];
                    cerr << endl
                         << "VO_D_well_s_enthalpy_fluid [ n ]: "
                         << VO_D_well_s_enthalpy_fluid[n]
                            << ", VO_D_well_smf_fluid [ n ]: " << VO_D_well_smf_fluid[n]
                               << ", VO_D_well_mf_halite [ n ]: " << VO_D_well_mf_halite[n]
                                  << endl;
                    cerr << endl
                         << "VO_D_well_s_enthalpy_bulk [ n ]: "
                         << VO_D_well_s_enthalpy_bulk[n]
                            << ", VO_D_well_smf_bulk [ n ]: " << VO_D_well_smf_bulk[n]
                               << endl;
                    cerr << endl
                         << "m_fluid_: " << m_fluid_ << ", h_fluid_: " << h_fluid_
                         << ", p_: " << p_ << ", smf_: " << smf_ << endl
                         << endl;
                }

                VO_D_well_sat_liquid[n] = Liquid.s;
                VO_D_well_sat_vapor[n] = Vapor.s;
                VO_D_well_sat_halite[n] = Salt.s;

                if (verbose) {
                    cerr << "sat liquid: " << VO_D_well_sat_liquid[n]
                            << ", sat vapor: " << VO_D_well_sat_vapor[n]
                               << ", sat halite: " << VO_D_well_sat_halite[n] << endl;
                }

                VO_D_well_mf_liquid[n] = Liquid.mf;
                VO_D_well_mf_vapor[n] = Vapor.mf;
                VO_D_well_mf_fluid[n] = Liquid.mf + Vapor.mf;
                VO_D_well_mf_halite[n] = Salt.mf;

                if (verbose) {
                    cerr << "mf liquid: " << VO_D_well_mf_liquid[n]
                            << ", mf vapor: " << VO_D_well_mf_vapor[n]
                               << ", mf fluid: " << VO_D_well_mf_fluid[n]
                                  << ", mf halite: " << VO_D_well_mf_halite[n] << endl;
                }

                if (VO_D_well_mf_fluid[n] + VO_D_well_mf_halite[n] > 1.) {
                    cerr << endl << "Sum of mass fractions above 1!!!";
                    cerr << endl
                         << endl
                         << "Segment " << n << " T-eq............................"
                         << endl;
                    cerr << endl
                         << "VO_D_well_s_enthalpy_fluid [ n ]: "
                         << VO_D_well_s_enthalpy_fluid[n]
                            << ", VO_D_well_smf_fluid [ n ]: " << VO_D_well_smf_fluid[n]
                               << ", VO_D_well_mf_halite [ n ]: " << VO_D_well_mf_halite[n]
                                  << endl;
                    cerr << endl
                         << "m_fluid_: " << m_fluid_ << ", h_fluid_: " << h_fluid_
                         << ", p_: " << p_ << ", smf_: " << smf_ << endl;
                    cerr << "mf liquid: " << VO_D_well_mf_liquid[n]
                            << ", mf vapor: " << VO_D_well_mf_vapor[n]
                               << ", mf fluid: " << VO_D_well_mf_fluid[n]
                                  << ", mf halite: " << VO_D_well_mf_halite[n] << endl;
                    // PAUSE();
                }

                VO_D_well_density_liquid[n] = Liquid.rho;
                VO_D_well_density_vapor[n] = Vapor.rho;
                VO_D_well_density_halite[n] = Salt.rho;

                if (VO_D_well_density_liquid[n] > 0. && VO_D_well_density_vapor[n] > 0.)
                    VO_D_well_density_fluid[n] =
                            (VO_D_well_mf_liquid[n] + VO_D_well_mf_vapor[n]) /
                            (VO_D_well_mf_liquid[n] / VO_D_well_density_liquid[n] +
                             VO_D_well_mf_vapor[n] / VO_D_well_density_vapor[n]);
                else if (VO_D_well_density_liquid[n] > 0.)
                    VO_D_well_density_fluid[n] = VO_D_well_density_liquid[n];
                else if (VO_D_well_density_vapor[n] > 0.)
                    VO_D_well_density_fluid[n] = VO_D_well_density_vapor[n];
                else
                    VO_D_well_density_fluid[n] = 0.;

                VO_D_well_density_bulk[n] = Bulk.rho;

                if (verbose) {
                    cerr << "previous_density liquid: "
                         << VO_D_previous_well_density_liquid[n]
                            << ", previous_density vapor: "
                            << VO_D_previous_well_density_vapor[n]
                               << ", previous_density fluid; "
                               << VO_D_previous_well_density_fluid[n]
                                  << ", previous_density bulk; "
                                  << VO_D_previous_well_density_bulk[n]
                                     << ", previous_density halite: "
                                     << VO_D_previous_well_density_halite[n] << endl;

                    cerr << "density liquid: " << VO_D_well_density_liquid[n]
                            << ", density vapor: " << VO_D_well_density_vapor[n]
                               << ", density fluid; " << VO_D_well_density_fluid[n]
                                  << ", density bulk; " << VO_D_well_density_bulk[n]
                                     << ", density halite; " << VO_D_well_density_halite[n] << endl;
                }

                VO_D_well_s_enthalpy_liquid[n] = Liquid.h;
                VO_D_well_s_enthalpy_vapor[n] = Vapor.h;
                VO_D_well_s_enthalpy_halite[n] = Salt.h;

                //////////////////
                VO_D_well_s_enthalpy_fluid[n] =
                        VO_D_well_s_enthalpy_liquid[n] * VO_D_well_mf_liquid[n] /
                        VO_D_well_mf_fluid[n] +
                        VO_D_well_s_enthalpy_vapor[n] * VO_D_well_mf_vapor[n] /
                        VO_D_well_mf_fluid[n];
                //////////////////

                // VO_D_well_s_enthalpy_bulk  [ n ] = Bulk.h; NO UPDATE

                if (verbose) {
                    cerr << "s_enthalpy liquid: " << VO_D_well_s_enthalpy_liquid[n]
                            << ", s_enthalpy vapor: " << VO_D_well_s_enthalpy_vapor[n]
                               << ", s_enthalpy fluid; " << VO_D_well_s_enthalpy_fluid[n]
                                  << ", s_enthalpy bulk; " << VO_D_well_s_enthalpy_bulk[n]
                                     << ", s_enthalpy halite; " << VO_D_well_s_enthalpy_halite[n]
                                        << endl;
                }

                VO_D_well_viscosity_liquid[n] = Liquid.mu;
                VO_D_well_viscosity_vapor[n] = Vapor.mu;
                VO_D_well_viscosity_fluid[n] =
                        pow(VO_D_well_viscosity_liquid[n],
                            VO_D_well_sat_liquid[n] / (1. - VO_D_well_sat_halite[n])) *
                        pow(VO_D_well_viscosity_vapor[n],
                            VO_D_well_sat_vapor[n] / (1. - VO_D_well_sat_halite[n]));

                VO_D_well_smf_liquid[n] = Liquid.smf;
                VO_D_well_smf_vapor[n] = Vapor.smf;

                //////////////////
                VO_D_well_smf_fluid[n] =
                        VO_D_well_smf_liquid[n] * VO_D_well_mf_liquid[n] /
                        VO_D_well_mf_fluid[n] +
                        VO_D_well_smf_vapor[n] * VO_D_well_mf_vapor[n] /
                        VO_D_well_mf_fluid[n];
                //////////////////

                // VO_D_well_smf_bulk    [ n ] = Bulk.smf; NO UPDATE

                if (verbose) {
                    cerr << "smf liquid: " << VO_D_well_smf_liquid[n]
                            << ", smf vapor: " << VO_D_well_smf_vapor[n] << ", smf fluid; "
                            << VO_D_well_smf_fluid[n] << ", smf bulk; "
                            << VO_D_well_smf_bulk[n] << endl;
                }

                // VO_D_well_salt_mass [ n ] = VO_D_well_density_bulk [ n ] *
                // VO_D_well_smf_bulk [ n ] * VO_D_segment_length [ n ] *
                // VO_D_well_cross_area [ n ] * ( 1. - VO_D_well_air_saturation [ n ] );

                VO_D_well_salt_mass[n] =
                        VO_D_segment_length[n] * VO_D_well_cross_area[n] *
                        (1. - VO_D_well_air_saturation[n]) *
                        ((1. - VO_D_well_sat_halite[n]) * VO_D_well_density_fluid[n] *
                         VO_D_well_smf_fluid[n] +
                         VO_D_well_sat_halite[n] * VO_D_well_density_halite[n]);

                if (verbose) {
                    cerr << "mass salt: " << VO_D_well_salt_mass[n]
                            << ", previous salt mass; "
                            << VO_D_segment_length[n] * VO_D_previous_well_cross_area[n] *
                               (1. - VO_D_previous_well_air_saturation[n]) *
                               ((1. - VO_D_previous_well_sat_halite[n]) *
                                VO_D_previous_well_density_fluid[n] *
                                VO_D_previous_well_smf_fluid[n] +
                                VO_D_previous_well_sat_halite[n] *
                                VO_D_previous_well_density_halite[n])
                            << endl;
                }

                VO_D_well_fluid_state[n] = Bulk.state;

                if (!thermal_eq_succeeded) {
                    thermal_eq_crashed_this_NR = true;

                    cerr << endl
                         << endl
                         << "Segment: " << n << " , T-eq did not converge";
                    cerr << endl << "h_fluid_: " << h_fluid_ << ", p_: " << p_;

                    if (VO_D_NR_previous_well_density_bulk[n] >
                            0.) // checking that if in the previous NR, the segment was not
                        // empty, else exit and cut dt!!!!
                    {
                        cerr << endl
                             << "previous_density liquid: "
                             << VO_D_previous_well_density_liquid[n]
                                << ", previous_density vapor: "
                                << VO_D_previous_well_density_vapor[n]
                                   << ", previous_density fluid; "
                                   << VO_D_previous_well_density_fluid[n]
                                      << ", previous_density bulk; "
                                      << VO_D_previous_well_density_bulk[n];
                        cerr << endl
                             << "NRprevious_density liquid: "
                             << VO_D_NR_previous_well_density_liquid[n]
                                << ", NRprevious_density vapor: "
                                << VO_D_NR_previous_well_density_vapor[n]
                                   << ", NRprevious_density fluid; "
                                   << VO_D_NR_previous_well_density_fluid[n]
                                      << ", NRprevious_density bulk; "
                                      << VO_D_NR_previous_well_density_bulk[n];
                        cerr << endl
                             << "density liquid: " << VO_D_well_density_liquid[n]
                                << ", density vapor: " << VO_D_well_density_vapor[n]
                                   << ", density fluid; " << VO_D_well_density_fluid[n]
                                      << ", density bulk; " << VO_D_well_density_bulk[n];

                        cerr << endl
                             << "previous_s_enthalpy liquid: "
                             << VO_D_previous_well_s_enthalpy_liquid[n]
                                << ", previous_s_enthalpy vapor: "
                                << VO_D_previous_well_s_enthalpy_vapor[n]
                                   << ", previous_s_enthalpy fluid; "
                                   << VO_D_previous_well_s_enthalpy_fluid[n]
                                      << ", previous_s_enthalpy bulk; "
                                      << VO_D_previous_well_s_enthalpy_bulk[n];
                        cerr << endl
                             << "NRprevious_s_enthalpy liquid: "
                             << VO_D_NR_previous_well_s_enthalpy_liquid[n]
                                << ", NRprevious_s_enthalpy vapor: "
                                << VO_D_NR_previous_well_s_enthalpy_vapor[n]
                                   << ", NRprevious_s_enthalpy fluid; "
                                   << VO_D_NR_previous_well_s_enthalpy_fluid[n]
                                      << ", NRprevious_s_enthalpy bulk; "
                                      << VO_D_NR_previous_well_s_enthalpy_bulk[n];
                        cerr << endl
                             << "s_enthalpy liquid: " << VO_D_well_s_enthalpy_liquid[n]
                                << ", s_enthalpy vapor: " << VO_D_well_s_enthalpy_vapor[n]
                                   << ", s_enthalpy fluid; " << VO_D_well_s_enthalpy_fluid[n]
                                      << ", s_enthalpy bulk; " << VO_D_well_s_enthalpy_bulk[n];
                        PAUSE();

                        // VO_D_well_temperature      [ n ] =
                        // VO_D_NR_previous_well_temperature      [ n ];

                        VO_D_well_sat_liquid[n] = VO_D_NR_previous_well_sat_liquid[n];
                        VO_D_well_sat_vapor[n] = VO_D_NR_previous_well_sat_vapor[n];
                        VO_D_well_sat_halite[n] = VO_D_NR_previous_well_sat_halite[n];

                        VO_D_well_mf_liquid[n] = VO_D_NR_previous_well_mf_liquid[n];
                        VO_D_well_mf_vapor[n] = VO_D_NR_previous_well_mf_vapor[n];
                        VO_D_well_mf_fluid[n] = VO_D_NR_previous_well_mf_vapor[n] +
                                VO_D_NR_previous_well_mf_liquid[n];
                        VO_D_well_mf_halite[n] = VO_D_NR_previous_well_mf_halite[n];

                        VO_D_well_density_liquid[n] =
                                VO_D_NR_previous_well_density_liquid[n];
                        VO_D_well_density_vapor[n] = VO_D_NR_previous_well_density_vapor[n];
                        VO_D_well_density_halite[n] =
                                VO_D_NR_previous_well_density_halite[n];
                        VO_D_well_density_fluid[n] = VO_D_NR_previous_well_density_fluid[n];
                        VO_D_well_density_bulk[n] = VO_D_NR_previous_well_density_bulk[n];

                        VO_D_well_s_enthalpy_liquid[n] =
                                VO_D_NR_previous_well_s_enthalpy_liquid[n];
                        VO_D_well_s_enthalpy_vapor[n] =
                                VO_D_NR_previous_well_s_enthalpy_vapor[n];
                        VO_D_well_s_enthalpy_halite[n] =
                                VO_D_NR_previous_well_s_enthalpy_halite[n];

                        //////////
                        VO_D_well_s_enthalpy_fluid[n] =
                                VO_D_NR_previous_well_s_enthalpy_fluid[n]; // not sure if thats
                        // gonna converge
                        //////////

                        VO_D_well_s_enthalpy_bulk[n] =
                                VO_D_NR_previous_well_s_enthalpy_bulk[n];

                        VO_D_well_viscosity_liquid[n] =
                                VO_D_NR_previous_well_viscosity_liquid[n];
                        VO_D_well_viscosity_vapor[n] =
                                VO_D_NR_previous_well_viscosity_vapor[n];
                        VO_D_well_viscosity_fluid[n] =
                                VO_D_NR_previous_well_viscosity_fluid[n];

                        VO_D_well_smf_liquid[n] = VO_D_NR_previous_well_smf_liquid[n];
                        VO_D_well_smf_vapor[n] = VO_D_NR_previous_well_smf_vapor[n];

                        //////////
                        VO_D_well_smf_fluid[n] = VO_D_NR_previous_well_smf_fluid[n];
                        //////////

                        VO_D_well_smf_bulk[n] = VO_D_NR_previous_well_smf_bulk[n];

                        //////////
                        VO_D_well_salt_mass[n] =
                                VO_D_segment_length[n] * VO_D_well_cross_area[n] *
                                ((1. - VO_D_well_air_saturation[n]) *
                                 (1. - VO_D_well_sat_halite[n]) *
                                 VO_D_well_density_fluid[n] * VO_D_well_smf_fluid[n] +
                                 VO_D_well_sat_halite[n] * VO_D_well_density_halite[n]);
                        //////////

                        VO_D_well_fluid_state[n] = VO_D_NR_previous_well_fluid_state[n];
                    }

                    else
                        exit_and_cut_dt = true;
                }

                // VO_D_well_last_true_density_fluid    [ n ] = VO_D_well_density_fluid
                // [ n ];
            }

            VO_D_well_last_true_density_fluid[n] = VO_D_well_density_fluid[n];

            if (VO_D_well_sat_halite[n] > 0.95 && VO_well_segment_activated[n] > 0) {
                // VO_well_segment_activated[ n ] = 0;
            }
        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Initial_T_equilibration(const InitialWellTemperature &t0) {
    // We make a first guess for the density in the well:
    std::fill(VO_D_well_density_fluid.begin(), VO_D_well_density_fluid.end(),
              1000.);
    std::fill(VO_D_well_pressure.begin(), VO_D_well_pressure.end(), top_pressure);

    if (t0.mode == InitialWellTemperatureMode::ReservoirWithMinimum) {
        top_temperature = t0.value;

        // double maxResT = 0.;
        // for( unsigned n = 0; n < number_of_well_segments; n++ )
        // {
        //     if ( VO_D_temperature [ n ] > maxResT )  maxResT = VO_D_temperature [
        //     n ];
        // }

        // double T_over_Res = delta_T - VO_D_temperature [ 0 ];
        // for( unsigned n = 0; n < number_of_well_segments; n++ )
        // {
        //     VO_D_well_temperature [ n ] = VO_D_temperature [ n ] + T_over_Res;
        //     if( VO_D_well_temperature [ n ] > maxResT )  VO_D_well_temperature [
        //     n ] = maxResT;
        // }

        for (unsigned n = 0; n < number_of_well_segments; n++) {
            VO_D_well_temperature[n] = VO_D_temperature[n];
            if (VO_D_well_temperature[n] < t0.value)
                VO_D_well_temperature[n] = t0.value;
        }
    }
    else if (t0.mode == InitialWellTemperatureMode::Uniform) {
        top_temperature = t0.value;
        for (unsigned n = 0; n < number_of_well_segments; n++)
            VO_D_well_temperature[n] = t0.value;
    }
    else { // FollowReservoir
        top_temperature = VO_D_temperature[0];
        VO_D_well_temperature = VO_D_temperature; // Twell = Tres
    }

    // The three branches are now mutually exclusive. Previously they were three
    // sequential `if`s, so follow_reservoir_T silently overrode init_T, which
    // silently overrode delta_T — and init_T was only applied above a hardcoded
    // threshold of 5 degrees. (Benoit DD/MM/YYYY)

    // std::fill( VO_D_well_temperature.begin(),   VO_D_well_temperature.end(),
    // 250. ); // option to put homogeneous T everywhere in the well

    // std::fill( VO_D_well_smf_bulk.begin(),   VO_D_well_smf_bulk.end(), 0.01 );

    //    //    //COLD INJECTION INITIAL: 20C down to 3400m,
    //            for( unsigned n = 0; n < number_of_well_segments; n++ )
    //            {
    //                if(VO_D_segment_heel_depth [ n ] > -3400. )
    //                {
    //                    VO_D_well_temperature [ n ] = 20. ;
    //                }

    //                if(VO_D_segment_heel_depth [ n ] < -3400. &&
    //                VO_D_segment_heel_depth [ n ] > -4300. )
    //                {
    //                    VO_D_well_temperature [ n ] = 20. +
    //                    abs(VO_D_segment_heel_depth [ n ] + 3400.) * 0.15;
    //                }
    //            }

    //    for( unsigned n = 0; n < number_of_well_segments; n++ )
    //    {
    //        if(VO_D_segment_heel_depth [ n ] > -1500. )
    //        {
    //            VO_D_well_temperature [ n ] = VO_D_temperature [n] * 0.5;
    //        }
    //        else
    //        {
    //            VO_D_well_temperature [ n ] = VO_D_temperature [n] * 1.5;
    //        }

    //        if(VO_D_well_temperature [ n ]>600.) VO_D_well_temperature [ n
    //        ]=600.; if(VO_D_well_temperature [ n ]<6.) VO_D_well_temperature [ n
    //        ]=6.;
    //    }

    // We perform an initialization loop (basically converging pressure-density
    // solution)
    for (unsigned k = 0; k < 5; k++) {
        // We compute the pressure in the well:
        Initial_P();

        // We compute thermodynamic properties
        for (unsigned n = 0; n < number_of_well_segments; n++) {

            if (VO_D_well_water_table[n] == 0) { // segment filled with air
                VO_D_well_temperature[n] = 0.; // We could use air temp: reservoir temp?

                VO_D_well_sat_liquid[n] = 0.;
                VO_D_well_sat_vapor[n] = 0.;
                VO_D_well_sat_halite[n] =
                        0.; //! I suppose there could still be halite..

                VO_D_well_mf_liquid[n] = 0.;
                VO_D_well_mf_vapor[n] = 0.;
                VO_D_well_mf_fluid[n] = 0.;
                VO_D_well_mf_halite[n] = 0.; //!

                VO_D_well_density_liquid[n] = 0.;
                VO_D_well_density_vapor[n] = 0.;
                VO_D_well_density_halite[n] = 0.; //!
                VO_D_well_density_fluid[n] = 0.;
                VO_D_well_density_bulk[n] = 0.; //!

                VO_D_well_s_enthalpy_liquid[n] = 0.;
                VO_D_well_s_enthalpy_vapor[n] = 0.;
                VO_D_well_s_enthalpy_halite[n] = 0.; //!
                VO_D_well_s_enthalpy_fluid[n] = 0.;
                VO_D_well_s_enthalpy_bulk[n] = 0.; //!

                VO_D_well_viscosity_liquid[n] = 0.;
                VO_D_well_viscosity_vapor[n] = 0.;
                VO_D_well_viscosity_fluid[n] = 0.;

                VO_D_well_smf_liquid[n] = 0.;
                VO_D_well_smf_vapor[n] = 0.;
                VO_D_well_smf_fluid[n] = 0.;
                VO_D_well_smf_bulk[n] = 0.; //!

                VO_D_well_salt_mass[n] = 0.; //!

                VO_D_well_fluid_state[n] = 0;
            }

            else { // segment is at least partially filled with water
                smf_ = VO_D_well_smf_bulk[n];
                wt_ = smf_ * 100.;
                x_ = Weight2XNaCl(wt_);
                p_ = VO_D_well_pressure[n];
                if (use_center_pressure or n == 0)
                    p_ = VO_D_interp_well_pressure[n]; // FORCED FOR  n = 0!
                t_ = tp_ = VO_D_well_temperature[n];

                // h_fluid_ = VO_D_well_s_enthalpy_bulk [ n ] = brine.Enthalpy();
                h_fluid_ = VO_D_well_s_enthalpy_bulk[n] = fluid.BulkProperties().h;

                // double rho = brine.Density();
                double rho = fluid.BulkProperties().rho;

                if (rho > 0)
                    VO_D_well_density_bulk[n] = rho;

                m_fluid_ = 1000.;

                H_current_ = m_fluid_ * h_fluid_;

                Bulk = fluid.BulkProperties();
                Liquid = fluid.ReportLiquidProperties();
                Vapor = fluid.ReportVaporProperties();
                Salt = fluid.ReportSaltProperties();

                // cerr<<endl<<"The expected T is: "<<t_<<", the bulk T is:
                // "<<Bulk.t<<endl; // Careful, Bulk.t needs to be written here
                // cerr<<"Segment: "<<n<<", h_fluid_: "<<h_fluid_<<", m_fluid_:
                // "<<m_fluid_<<", H_current: "<<H_current_<<", p_: "<<p_<<", smf:
                // "<<smf_<<endl;

                // Here we update thermodynamic properties
                VO_D_well_sat_liquid[n] = Liquid.s;
                VO_D_well_sat_vapor[n] = Vapor.s;
                VO_D_well_sat_halite[n] = Salt.s;

                // cerr<<"sat liquid: "<<VO_D_well_sat_liquid [ n ]<<", sat vapor:
                // "<<VO_D_well_sat_vapor [ n ]<<", sat halite: "<<VO_D_well_sat_halite
                // [ n ]<<endl;

                VO_D_well_mf_liquid[n] = Liquid.mf;
                VO_D_well_mf_vapor[n] = Vapor.mf;
                VO_D_well_mf_fluid[n] = Liquid.mf + Vapor.mf;
                VO_D_well_mf_halite[n] = Salt.mf;
                // cerr<<"mf liquid: "<<VO_D_well_mf_liquid [ n ]<<", mf vapor:
                // "<<VO_D_well_mf_vapor [ n ]<<", mf fluid: "<<VO_D_well_mf_vapor [ n
                // ]<<", mf halite: "<<VO_D_well_mf_halite [ n ]<<endl;

                VO_D_well_density_liquid[n] = Liquid.rho;
                VO_D_well_density_vapor[n] = Vapor.rho;
                VO_D_well_density_halite[n] = Salt.rho;

                if (VO_D_well_density_liquid[n] > 0. && VO_D_well_density_vapor[n] > 0.)
                    VO_D_well_density_fluid[n] =
                            (VO_D_well_mf_liquid[n] + VO_D_well_mf_vapor[n]) /
                            (VO_D_well_mf_liquid[n] / VO_D_well_density_liquid[n] +
                             VO_D_well_mf_vapor[n] / VO_D_well_density_vapor[n]);
                else if (VO_D_well_density_liquid[n] > 0.)
                    VO_D_well_density_fluid[n] = VO_D_well_density_liquid[n];
                else if (VO_D_well_density_vapor[n] > 0.)
                    VO_D_well_density_fluid[n] = VO_D_well_density_vapor[n];
                else
                    VO_D_well_density_fluid[n] = 0.;

                VO_D_well_density_bulk[n] = Bulk.rho;
                // cerr<<"density liquid: "<<VO_D_well_density_liquid [ n ]<<", density
                // vapor: "<<VO_D_well_density_vapor [ n ]<<", density fluid;
                // "<<VO_D_well_density_fluid [ n ]<<", density bulk; "<<
                // VO_D_well_density_bulk [ n ]<<endl;

                VO_D_well_s_enthalpy_liquid[n] = Liquid.h;
                VO_D_well_s_enthalpy_vapor[n] = Vapor.h;
                VO_D_well_s_enthalpy_halite[n] = Salt.h;
                VO_D_well_s_enthalpy_fluid[n] =
                        VO_D_well_s_enthalpy_liquid[n] * VO_D_well_mf_liquid[n] /
                        VO_D_well_mf_fluid[n] +
                        VO_D_well_s_enthalpy_vapor[n] * VO_D_well_mf_vapor[n] /
                        VO_D_well_mf_fluid[n];
                VO_D_well_s_enthalpy_bulk[n] = Bulk.h;
                // cerr<<"s_enthalpy liquid: "<<VO_D_well_s_enthalpy_liquid [ n ]<<",
                // s_enthalpy vapor: "<<VO_D_well_s_enthalpy_vapor [ n ]<<", s_enthalpy
                // fluid; "<<VO_D_well_s_enthalpy_fluid [ n ]<<", s_enthalpy bulk; "<<
                // VO_D_well_s_enthalpy_bulk [ n ]<<endl;

                VO_D_well_viscosity_liquid[n] = Liquid.mu;
                VO_D_well_viscosity_vapor[n] = Vapor.mu;
                VO_D_well_viscosity_fluid[n] =
                        pow(VO_D_well_viscosity_liquid[n],
                            VO_D_well_sat_liquid[n] / (1. - VO_D_well_sat_halite[n])) *
                        pow(VO_D_well_viscosity_vapor[n],
                            VO_D_well_sat_vapor[n] / (1. - VO_D_well_sat_halite[n]));

                VO_D_well_smf_liquid[n] = Liquid.smf;
                VO_D_well_smf_vapor[n] = Vapor.smf;
                VO_D_well_smf_fluid[n] =
                        VO_D_well_smf_liquid[n] * VO_D_well_mf_liquid[n] /
                        VO_D_well_mf_fluid[n] +
                        VO_D_well_smf_vapor[n] * VO_D_well_mf_vapor[n] /
                        VO_D_well_mf_fluid[n];
                VO_D_well_smf_bulk[n] = Bulk.smf;
                // cerr<<"smf liquid: "<<VO_D_well_smf_liquid [ n ]<<", smf vapor:
                // "<<VO_D_well_smf_vapor [ n ]<<", smf fluid: "<<VO_D_well_smf_fluid [
                // n ]<<", smf bulk: "<< VO_D_well_smf_bulk [ n ]<<endl;

                VO_D_well_salt_mass[n] = VO_D_segment_length[n] *
                        VO_D_well_cross_area[n] *
                        (1 - VO_D_well_air_saturation[n]) *
                        VO_D_well_smf_bulk[n];

                VO_D_well_fluid_state[n] = Bulk.state;
            }
        }

        if (with_friction)
            ComputeFrictionFactor(); // Should be zero at this point since we
        // initialized velocities to zero
        // Display();
    }

    t_ = top_temperature;
    p_ = top_pressure;
    smf_ = top_smf;
    wt_ = smf_ * 100.;
    x_ = Weight2XNaCl(wt_);

    // top_s_enthalpy = brine.Enthalpy();
    top_s_enthalpy = fluid.BulkProperties().h;

    cerr << endl
         << "Initialized top of well with: top_temperature: " << top_temperature
         << ", top_pressure: " << top_pressure << ", top_smf: " << top_smf
         << ", top_s_enthalpy: " << top_s_enthalpy;
    // (will only apply when there is downflow from the top)
}

template <uint32_t dim> void WellModelPrototype<dim>::Initial_P() {
    // Initial hydrostatic pressure
    double pressure(0.);

    // hydrostatic pressure at segments heads
    // VO_D_well_pressure [ 0 ] = VO_D_well_pressure [ 1 ]  = top_pressure;
    VO_D_well_pressure[0] = top_pressure;

    for (unsigned n = 1; n < number_of_well_segments; n++) {
        pressure = VO_D_well_pressure[n - 1] +
                VO_D_segment_length[n - 1] *
                cos(VO_D_segment_inclination[n - 1]) * grav_acc *
                VO_D_well_density_fluid[n - 1] *
                (1. - VO_D_well_air_saturation[n - 1]);
        VO_D_well_pressure[n] = pressure;
    }

    Interpolate_pressure();
}

template <uint32_t dim> void WellModelPrototype<dim>::Interpolate_pressure() {
    double interp_pressure(0.);

    // cerr<<endl<<"Interpolate pressure";

    // cerr<<endl<<"node index 0, value at segment head index
    // "<<number_of_virtual_top_segments<<" is "<<VO_D_well_pressure [
    // number_of_virtual_top_segments ];
    VO_well_node_pressure[0] = VO_D_well_pressure[number_of_virtual_top_segments];

    // pressure at heels (on nodes, in order to compare with reservoir pressure)
    for (unsigned n = 1; n < number_of_well_nodes - 1; n++) {
        VO_well_node_pressure[n] =
                VO_D_well_pressure[index_of_last_segment[n - 1] + 1];
        // cerr<<endl<<"node index "<<n<<", value at segment head
        // index"<<index_of_last_segment [ n - 1 ] + 1<<" is:
        // "<<VO_well_node_pressure [ n ]<<", ";
    }

    // interpolate pressure at segments centers
    // cerr<<endl<<endl<<"Interpolate pressure at segment centers:";
    for (unsigned n = 0; n < number_of_well_segments - 1; n++) {
        interp_pressure = 0.5 * (VO_D_well_pressure[n] + VO_D_well_pressure[n + 1]);
        VO_D_interp_well_pressure[n] = interp_pressure;
        // cerr<<endl<<"segment "<<n<<" center pressure: "<<interp_pressure<<", ";
    }

    unsigned n(number_of_well_segments - 1);
    // APPROXIMATION of bottom pressure ( no inertia ), using pressure at heels
    // formulation is better in this regard
    bottom_pressure = VO_D_well_pressure[n] +
            VO_D_segment_length[n] * cos(VO_D_segment_inclination[n]) *
            grav_acc * VO_D_well_density_fluid[n];

    if (VO_D_well_velocity_fluid[n] != 0) {
        bottom_pressure -=
                VO_D_well_friction_factor[n] * VO_D_well_density_fluid[n] *
                VO_D_well_velocity_fluid[n] * VO_D_well_velocity_fluid[n] *
                VO_D_well_velocity_fluid[n] * VO_D_segment_length[n] /
                abs(VO_D_well_velocity_fluid[n]);
    }

    VO_D_interp_well_pressure[n] =
            0.5 * (VO_D_well_pressure[n] + bottom_pressure);

    VO_well_node_pressure[number_of_well_nodes - 1] = bottom_pressure;

    // if the water table is in the segment then center pressure is at the center
    // of the fluid column!
}

template <uint32_t dim> void WellModelPrototype<dim>::Velocity_model() {
    if (homogeneous_flow) {
        // easy, all phases move at the same velocity (fluid velocity)

        for (unsigned n = 0; n < number_of_well_segments; n++) { // NEW n = 0
            VO_D_well_velocity_liquid[n] = VO_D_well_velocity_fluid[n];
            VO_D_well_velocity_vapor[n] = VO_D_well_velocity_fluid[n];
        }
    }

    else {
        // The Drift Flux Model has never been implemented. This branch used to be
        // empty, so homogeneous_flow=false silently left the phase velocities at
        // whatever they held. Fail loudly instead. (Benoit DD/MM/YYYY)
        throw csmp::Exception(ERROR, "WellModelPrototype::Continuity_residual",
                              "homogeneous_flow == false selects the Drift Flux "
                              "Model, which is not implemented.");
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::ComputeDensityDerivatives() { // Currently
    // disabled

    // double half_delta = 10000, rhof_below, rhof_above, rhob_below, rhob_above,
    // t_below, t_above, h_below, h_above;

    std::fill(VO_D_drhofdp.begin(), VO_D_drhofdp.end(), 0.);
    std::fill(VO_D_drhofdhb.begin(), VO_D_drhofdhb.end(), 0.);
    std::fill(VO_D_drhofdsmfb.begin(), VO_D_drhofdsmfb.end(), 0.);

    std::fill(VO_D_drhobdp.begin(), VO_D_drhobdp.end(), 0.);
    std::fill(VO_D_drhobdhb.begin(), VO_D_drhobdhb.end(), 0.);
    std::fill(VO_D_drhobdsmfb.begin(), VO_D_drhobdsmfb.end(), 0.);

    std::fill(VO_D_dhdp.begin(), VO_D_dhdp.end(), 0.);

    //Benoit: code below leads to problems so it is disabled, new ComputeDensityDerivatives() function, further down, partly works but it does not seem to lead to great gains:
    // this may be due to the otherwise very conservative solutio incrementation in the NR loop, with proper density derivatives we could be more agressive,
    // we may also be able to reach tighter convergence criteria

    //    for( unsigned n = 0; n < number_of_well_segments; n++ )
    //    {
    //        if( VO_D_well_water_table [ n ] == 2 )
    //        {
    //            // define the point p, t, h, smf around which we want to compute
    //            derivatives: p_  = VO_D_well_pressure [ n ]; if (
    //            use_center_pressure ) p_ = VO_D_interp_well_pressure [ n ]; t_ =
    //            tp_ = VO_D_well_temperature [ n ]; smf_ = VO_D_well_smf_bulk [ n
    //            ]; wt_ = smf_ * 100.; x_  = Weight2XNaCl(wt_); //added

    //            m_fluid_    = 1000.;
    //            h_fluid_    = VO_D_well_s_enthalpy_bulk [ n ];
    //            H_current_  = m_fluid_ * h_fluid_;

    //            //            cerr<<endl<<"--------------------------";
    //            //            cerr<<endl<<"Segment: "<<n<<" *** Starting
    //            conditions for derivatives calculations ***";
    //            //            cerr<<endl<<"T: "<<t_<<", P: "<<p_<<", Smf:
    //            "<<smf_<<", h: "<<h_fluid_;

    //            ////////////////////////////////////////////////////
    //            //first: drhofdp (isenthalpic), ////we could expand to drhodT
    //            (isenthalpic) * dTdp (isenthalpic)?
    //            //cerr<<endl<<"drhofdp (isenthalpic): ";
    //            p_ -= half_delta;

    //            fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //            t_below = Bulk.t;

    //            rhof_below = (Liquid.rho * Liquid.s + Vapor.rho *
    //            Vapor.s)/(Liquid.s+Vapor.s); rhob_below = Bulk.rho;

    //            p_ += 2*half_delta;

    //            fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //            t_above = Bulk.t;

    //            rhof_above = (Liquid.rho * Liquid.s + Vapor.rho *
    //            Vapor.s)/(Liquid.s+Vapor.s); rhob_above = Bulk.rho;

    //            VO_D_drhofdp [ n ] = ( rhof_above - rhof_below ) / ( 2 *
    //            half_delta); VO_D_drhobdp [ n ] = ( rhob_above - rhob_below ) /
    //            ( 2 * half_delta);

    //            p_ -= half_delta;// re-initialize pressure

    //            //////////////////////////////////////////////////////
    //            //second: drhofdhb (isobaric), we expand to drhodT (isobaric) *
    //            dTdh (isobaric) h_fluid_ -= half_delta; H_current_  = m_fluid_ *
    //            h_fluid_;

    //            fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //            t_below = Bulk.t;

    //            rhof_below = (Liquid.rho * Liquid.s + Vapor.rho *
    //            Vapor.s)/(Liquid.s+Vapor.s); rhob_below = Bulk.rho;

    //            h_fluid_ += 2*half_delta;
    //            H_current_  = m_fluid_ * h_fluid_;

    //            fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //            t_above = Bulk.t;

    //            rhof_above = (Liquid.rho * Liquid.s + Vapor.rho *
    //            Vapor.s)/(Liquid.s+Vapor.s); rhob_above = Bulk.rho;

    //            VO_D_drhofdhb [ n ] = ( rhof_above - rhof_below ) / ( 2 *
    //            half_delta ); VO_D_drhobdhb [ n ] = ( rhob_above - rhob_below )
    //            / ( 2 * half_delta );

    //            h_fluid_ -= half_delta;// re-initialize h
    //            H_current_  = m_fluid_ * h_fluid_;

    //            //////////////////////////////////////////////////
    //            ////third: drhofdsmfb (cst p and h?)
    //            if ( smf_ > 0.0002 )
    //            {
    //                //cerr<<endl<<"drhofdsmfb: ";
    //                double delta = 0.0001;
    //                smf_ += delta;
    //                wt_ = smf_ * 100.;
    //                x_  = Weight2XNaCl(wt_);
    //                fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //                t_above = Bulk.t;

    //                rhof_above = (Liquid.rho * Liquid.s + Vapor.rho *
    //                Vapor.s)/(Liquid.s+Vapor.s); rhob_above = Bulk.rho;

    //                smf_ -= 2*delta;
    //                wt_ = smf_ * 100.;
    //                x_  = Weight2XNaCl(wt_);
    //                fluid_equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
    //                t_below = Bulk.t;

    //                rhof_below = (Liquid.rho * Liquid.s + Vapor.rho *
    //                Vapor.s)/(Liquid.s+Vapor.s); rhob_below = Bulk.rho;

    //                VO_D_drhofdsmfb [ n ] = ( rhof_above - rhof_below ) / ( 2 *
    //                delta ); VO_D_drhobdsmfb [ n ] = ( rhob_above - rhob_below )
    //                / ( 2 * delta );

    //                smf_ += delta; //re-initialize smf
    //            }

    //            // ////////////////////////////////////////////////////
    //            //            //fourth: dhdp (isothermal, not isochoric, not
    //            sure...)

    //            wt_ = smf_ * 100.;
    //            x_  = Weight2XNaCl(wt_);

    //            half_delta = 1000.;
    //            p_ -= half_delta;
    //            h_below = brine.Enthalpy();
    //            p_ += 2*half_delta;
    //            h_above = brine.Enthalpy();

    //            //cerr<<endl<<"hfluid : "<<h_fluid_<<", h_below : "<<h_below<<",
    //            h_above : "<<h_above<<", dhdp_below : "<<dhdp_below<<",
    //            dhdp_above : "<<dhdp_above;

    //            VO_D_dhdp [ n ] = ( h_above - h_below ) / ( 2 * half_delta );
    //            //isothermal, not isochoric
    //}
    //}

    //    cerr << endl << endl << "DRHObDP                                 : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_drhobdp[ n ] << ", ";
    //    cerr << endl << endl << "DRHObDHb                                 : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_drhobdhb[ n ] << ", ";
    //        cerr << endl << endl << "DRHObDSMFb               : ";
    //        for (unsigned n = 0; n < number_of_well_segments; n++)
    //            cerr << VO_D_drhobdsmfb[ n ] << ", ";
    //    cerr << endl << endl << "DHDP               : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_dhdp[ n ] << ", ";
}

// template <uint32_t dim>
// void WellModelPrototype<dim>::ComputeDensityDerivatives() {//NEW, UNTESTED! partly tested in Aug 2026

//     // -------------------------------------------------------------------------
//     // PURPOSE
//     // -------------------------------------------------------------------------
//     // Fill VO_D_drhofdp, VO_D_drhofdhb, VO_D_drhofdsmfb  (fluid density)
//     //      VO_D_drhobdp, VO_D_drhobdhb, VO_D_drhobdsmfb  (bulk density)
//     // for every water-filled segment (water_table == 2).
//     //
//     // These are used in the Jacobian in Fill_*_Jacobian_part_salt().
//     // All six are derivatives of the density definitions used in the residuals:
//     //
//     //   rhof = (mf_l + mf_v) / (mf_l/rho_l + mf_v/rho_v)   [mass-fraction harmonic mean]
//     //   rhob = Bulk.rho                                       [from ThreePhaseProperties]
//     //
//     // IMPORTANT: rhof here MUST be computed the same way as in T_equilibration,
//     // i.e. using the mass-fraction harmonic mean, NOT the saturation arithmetic
//     // mean (Liquid.rho*Liquid.s + Vapor.rho*Vapor.s)/(Liquid.s+Vapor.s).
//     // The two are NOT equivalent in two-phase regions, and using the wrong one
//     // would make the Jacobian inconsistent with the residuals.
//     //
//     // METHOD: central finite differences using fluid_equilibrator.ThreePhaseProperties().
//     // We cannot use the `fluid` object here because `fluid` takes (T, p, x) as
//     // inputs, and we do not know T — we only know (p, h, smf). Only the
//     // equilibrator maps (p, h, smf) -> T and all phase properties. The `fluid`
//     // object is appropriate only when T is already known (e.g. for the isothermal
//     // dhdp derivative).
//     //
//     // PHASE-BOUNDARY SAFETY
//     // -------------------------------------------------------------------------
//     // The main reason this function was originally disabled is that a fixed
//     // perturbation delta can cross a phase boundary (liquid<->two-phase<->vapor,
//     // or fluid<->halite), making the finite difference meaningless or explosive.
//     // We guard against this by:
//     //   1. Choosing adaptive deltas that are small relative to the current state.
//     //   2. After each ThreePhaseProperties call, checking whether the phase state
//     //      (Bulk.state) changed relative to the base state. If it did, we fall
//     //      back to a one-sided (forward) difference with a smaller step, and if
//     //      that still crosses, we return zero for that derivative at that segment
//     //      (the solver will fall back to approximate Jacobian behaviour, which is
//     //      safe given the existing damping).
//     //   3. Restoring ALL equilibrator input variables unconditionally after each
//     //      segment, so that a failure in segment n cannot corrupt segment n+1.
//     //      The inputs shared with the equilibrator via references are:
//     //      p_, wt_, tp_, m_fluid_, H_current_  (and derived: smf_, x_, h_fluid_).
//     // -------------------------------------------------------------------------

//     std::fill(VO_D_drhofdp.begin(),    VO_D_drhofdp.end(),    0.);
//     std::fill(VO_D_drhofdhb.begin(),   VO_D_drhofdhb.end(),   0.);
//     std::fill(VO_D_drhofdsmfb.begin(), VO_D_drhofdsmfb.end(), 0.);
//     std::fill(VO_D_drhobdp.begin(),    VO_D_drhobdp.end(),    0.);
//     std::fill(VO_D_drhobdhb.begin(),   VO_D_drhobdhb.end(),   0.);
//     std::fill(VO_D_drhobdsmfb.begin(), VO_D_drhobdsmfb.end(), 0.);
//     std::fill(VO_D_dhdp.begin(),       VO_D_dhdp.end(),       0.);

//     //Benoit comment for test
//     // Helper lambda: compute rhof the same way T_equilibration does,
//     // i.e. mass-fraction harmonic mean of liquid and vapor densities.
//     // This MUST stay consistent with the formula in T_equilibration.
//     auto compute_rhof = [](const Fluidproperties& L, const Fluidproperties& V) -> double {
//         const double mf_fluid = L.mf + V.mf;
//         if (mf_fluid <= 0.) return 0.;
//         const double inv_rho = (L.rho > 0. ? L.mf / L.rho : 0.)
//                 + (V.rho > 0. ? V.mf / V.rho : 0.);
//         return (inv_rho > 0.) ? mf_fluid / inv_rho : 0.;
//     };

//      for (unsigned n = 0; n < number_of_well_segments; n++) {

//         // Only compute derivatives for water-filled segments.
//         // Empty (0) and water-table (1) segments use simplified BCs; their
//         // Jacobian rows don't need density derivatives.
//         if (VO_D_well_water_table[n] != 2)
//             continue;

//         // ---------------------------------------------------------------
//         // Save the full equilibrator state so we can restore it no matter
//         // what happens below.  The equilibrator reads these class members
//         // by reference: p_, wt_, tp_, m_fluid_, H_current_.
//         // We also save the derived quantities so the loop's own state is clean.
//         // ---------------------------------------------------------------
//         const double saved_p         = p_;
//         const double saved_wt        = wt_;
//         const double saved_smf       = smf_;
//         const double saved_x         = x_;
//         const double saved_tp        = tp_;
//         const double saved_m_fluid   = m_fluid_;
//         const double saved_h_fluid   = h_fluid_;
//         const double saved_H_current = H_current_;

//         // Set base-state inputs for this segment.
//         smf_    = VO_D_well_smf_bulk[n];
//         wt_     = smf_ * 100.;
//         x_      = Weight2XNaCl(wt_);
//         tp_     = std::max(VO_D_previous_well_temperature[n], 10.);
//         m_fluid_ = 1000.;
//         h_fluid_ = VO_D_well_s_enthalpy_bulk[n];
//         H_current_ = m_fluid_ * h_fluid_;

//         p_ = use_center_pressure ? VO_D_interp_well_pressure[n]
//                                    : VO_D_well_pressure[n];
//         if (n == 0)
//             p_ = VO_D_interp_well_pressure[n]; // always use centre pressure at top

//         // Record the base-state phase flag (Bulk.state integer from the
//         // equilibrator) so we can detect phase-boundary crossings.
//         const int base_state = static_cast<int>(VO_D_well_fluid_state[n]);

//         // Base-state densities (already computed by T_equilibration earlier
//         // in this NR iteration — we read them directly rather than re-calling
//         // the equilibrator, which saves cost and avoids any state mutation).
//         const double rhof_base = VO_D_well_density_fluid[n];
//         const double rhob_base = VO_D_well_density_bulk[n];

//         // ---------------------------------------------------------------
//         // Adaptive delta sizes.
//         // For pressure: 0.1% of current pressure, but at least 500 Pa and
//         //   at most 50000 Pa.  The original hard-coded 10000 Pa was fine at
//         //   ~10 MPa but too large near boiling at low pressure, and too
//         //   small at supercritical pressures.
//         // For enthalpy: 0.1% of current specific enthalpy.
//         // For smf: absolute 0.0001 (0.01 wt%), only if smf > 0.0002.
//         // ---------------------------------------------------------------
//         const double dp   = std::min(std::max(0.001 * p_,       500.), 50000.);
//         const double dh   = std::max(0.001 * std::abs(h_fluid_), 100.);
//         const double dsmf = 0.0001;

//         // Local copies of Bulk/Liquid/Vapor/Salt so we don't permanently
//         // overwrite the class members during the FD calls.
//         Fluidproperties Blo, Llo, Vlo, Slo;  // below/minus perturbation
//         Fluidproperties Bhi, Lhi, Vhi, Shi;  // above/plus perturbation

//         // ---------------------------------------------------------------
//         // 1.  drhofdp  and  drhobdp   (isenthalpic: h fixed, p varied)
//         // ---------------------------------------------------------------
//         {
//             bool use_central = true;

//             // --- lower perturbation ---
//             p_ -= dp;
//             bool ok_lo = fluid_equilibrator.ThreePhaseProperties(Blo, Llo, Vlo, Slo);
//             if (!ok_lo || static_cast<int>(Blo.state) != base_state) {
//                 use_central = false;  // phase boundary crossed going down
//             }
//             p_ += dp;  // restore

//             // --- upper perturbation ---
//             p_ += dp;
//             bool ok_hi = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//             if (!ok_hi || static_cast<int>(Bhi.state) != base_state) {
//                 use_central = false;
//             }
//             p_ -= dp;  // restore

//             if (use_central) {
//                 const double rhof_lo = compute_rhof(Llo, Vlo);
//                 const double rhof_hi = compute_rhof(Lhi, Vhi);
//                 VO_D_drhofdp[n] = (rhof_hi - rhof_lo) / (2. * dp);
//                 VO_D_drhobdp[n] = (Bhi.rho - Blo.rho) / (2. * dp);
//             } else {
//                 // Try a forward-only difference with half the step.
//                 // If that also crosses a phase boundary, leave at zero.
//                 const double dp_small = 0.5 * dp;
//                 p_ += dp_small;
//                 bool ok_fwd = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//                 p_ -= dp_small;
//                 if (ok_fwd && static_cast<int>(Bhi.state) == base_state) {
//                     VO_D_drhofdp[n] = (compute_rhof(Lhi, Vhi) - rhof_base) / dp_small;
//                     VO_D_drhobdp[n] = (Bhi.rho - rhob_base)                / dp_small;
//                 }
//                 // else: leave at zero — approximate Jacobian, solver handles it
//             }
//         }

//         //Benoit: this part is commented because it causes problems, need to debug
//         // // ---------------------------------------------------------------
//         // // 2.  drhofdhb  and  drhobdhb   (isobaric: p fixed, h varied)
//         // // ---------------------------------------------------------------
//         // {
//         //     bool use_central = true;

//         //     h_fluid_ -= dh;
//         //     H_current_ = m_fluid_ * h_fluid_;
//         //     bool ok_lo = fluid_equilibrator.ThreePhaseProperties(Blo, Llo, Vlo, Slo);
//         //     if (!ok_lo || static_cast<int>(Blo.state) != base_state)
//         //         use_central = false;
//         //     h_fluid_ += dh;
//         //     H_current_ = m_fluid_ * h_fluid_;

//         //     h_fluid_ += dh;
//         //     H_current_ = m_fluid_ * h_fluid_;
//         //     bool ok_hi = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//         //     if (!ok_hi || static_cast<int>(Bhi.state) != base_state)
//         //         use_central = false;
//         //     h_fluid_ -= dh;
//         //     H_current_ = m_fluid_ * h_fluid_;

//         //     if (use_central) {
//         //         const double rhof_lo = compute_rhof(Llo, Vlo);
//         //         const double rhof_hi = compute_rhof(Lhi, Vhi);
//         //         VO_D_drhofdhb[n] = (rhof_hi - rhof_lo) / (2. * dh);
//         //         VO_D_drhobdhb[n] = (Bhi.rho - Blo.rho) / (2. * dh);
//         //     } else {
//         //         const double dh_small = 0.5 * dh;
//         //         h_fluid_ += dh_small;
//         //         H_current_ = m_fluid_ * h_fluid_;
//         //         bool ok_fwd = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//         //         h_fluid_ -= dh_small;
//         //         H_current_ = m_fluid_ * h_fluid_;
//         //         if (ok_fwd && static_cast<int>(Bhi.state) == base_state) {
//         //             VO_D_drhofdhb[n] = (compute_rhof(Lhi, Vhi) - rhof_base) / dh_small;
//         //             VO_D_drhobdhb[n] = (Bhi.rho - rhob_base)                / dh_small;
//         //         }
//         //     }
//         // }

//         //Benoit: this part is not tested
//     //     // ---------------------------------------------------------------
//     //     // 3.  drhofdsmfb  and  drhobdsmfb   (p and h fixed, smf varied)
//     //     //     Only meaningful when there is actually some salt.
//     //     // ---------------------------------------------------------------
//     //     if (smf_ > 0.0002) {
//     //         bool use_central = true;

//     //         auto set_smf = [&](double new_smf) {
//     //             smf_ = new_smf;
//     //             wt_  = smf_ * 100.;
//     //             x_   = Weight2XNaCl(wt_);
//     //         };

//     //         set_smf(smf_ - dsmf);
//     //         bool ok_lo = fluid_equilibrator.ThreePhaseProperties(Blo, Llo, Vlo, Slo);
//     //         if (!ok_lo || static_cast<int>(Blo.state) != base_state)
//     //             use_central = false;
//     //         set_smf(smf_ + dsmf);  // restore to base

//     //         set_smf(smf_ + dsmf);
//     //         bool ok_hi = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//     //         if (!ok_hi || static_cast<int>(Bhi.state) != base_state)
//     //             use_central = false;
//     //         set_smf(smf_ - dsmf);  // restore to base

//     //         if (use_central) {
//     //             VO_D_drhofdsmfb[n] = (compute_rhof(Lhi, Vhi) - compute_rhof(Llo, Vlo)) / (2. * dsmf);
//     //             VO_D_drhobdsmfb[n] = (Bhi.rho - Blo.rho) / (2. * dsmf);
//     //         } else {
//     //             set_smf(smf_ + dsmf);
//     //             bool ok_fwd = fluid_equilibrator.ThreePhaseProperties(Bhi, Lhi, Vhi, Shi);
//     //             set_smf(smf_ - dsmf);
//     //             if (ok_fwd && static_cast<int>(Bhi.state) == base_state) {
//     //                 VO_D_drhofdsmfb[n] = (compute_rhof(Lhi, Vhi) - rhof_base) / dsmf;
//     //                 VO_D_drhobdsmfb[n] = (Bhi.rho - rhob_base)                / dsmf;
//     //             }
//     //         }
//     //     }
//     //     // else: smf derivatives stay at zero — correct, salt effect is negligible

//         // ---------------------------------------------------------------
//         // 4.  dhdp  (isothermal: T fixed at current value, p varied)
//         //     Here we CAN use the `fluid` object because T is already known.
//         //     The `fluid` object reads t_, p_, x_ by reference.
//         //     This derivative is used in the energy Jacobian for the
//         //     enthalpy-pressure coupling term.
//         // ---------------------------------------------------------------
//         {
//             // t_ is already set to the current temperature via T_equilibration
//             // (VO_D_well_temperature[n]). We just need to perturb p_.
//             // Note: we only do a forward difference here because the isothermal
//             // enthalpy is smooth and we want to avoid an extra equilibrator call.
//             const double p_base = p_;
//             t_ = VO_D_well_temperature[n];  // ensure t_ is correct

//             p_ = p_base + dp;
//             const double h_hi = fluid.BulkProperties().h;

//             p_ = p_base - dp;
//             const double h_lo = fluid.BulkProperties().h;

//             p_ = p_base;  // restore
//             VO_D_dhdp[n] = (h_hi - h_lo) / (2. * dp);
//         }

//         // ---------------------------------------------------------------
//         // Unconditional state restore.
//         // This guarantees that even if any of the branches above returned
//         // early (via the ok_* checks), the equilibrator inputs are left in
//         // a clean state for the next segment.
//         // ---------------------------------------------------------------
//         p_         = saved_p;
//         wt_        = saved_wt;
//         smf_       = saved_smf;
//         x_         = saved_x;
//         tp_        = saved_tp;
//         m_fluid_   = saved_m_fluid;
//         h_fluid_   = saved_h_fluid;
//         H_current_ = saved_H_current;

//     } // end segment loop
// }

template <uint32_t dim> void WellModelPrototype<dim>::NR_iteration() {
    cerr << endl;
    inflexion_pressure = 1.e8;
    well_solution_converged = false;
    density_converged = false;
    exit_and_cut_dt = false;

    unsigned n(1);
    number_of_iterations = 1;

    double previous_cycle_total_rate(0.);
    double previous_cycle_top_pressure(0.);
    double current_cycle_last_converged_total_rate(0.);
    double current_cycle_last_converged_top_pressure(0.);
    double total_injectivity(0.);

    for (unsigned n = 1; n < number_of_well_nodes;
         n++) // we skip n = 0 because surely there will be no exchange term on
        // top virtual segment
    {
        if (VO_well_completion[n] > 0) {
            total_injectivity += VO_well_injectivity[n];
        }
    }

    previous_cycle_total_rate = total_rate;
    previous_cycle_top_pressure = top_pressure;

    if (target_rate_active) {
        abs_p_min = 1.01325e5;
        // abs_p_min = 2.e5;
        abs_p_max = 1.e8;

        //        double pressure_change_rate = 0.01; //bars per min
        //        double max_pressure_correction = pressure_change_rate * dt / 60.
        //        * 1.e5;
        double proposed_pressure_correction = 0.;
        double proposed_pressure_correction2 = 0.;
        double proposed_pressure = 0.;
        double pressure_correction = 0.;

        // abs_p_max = previous_top_pressure + max_pressure_correction;
        // abs_p_min = previous_top_pressure - max_pressure_correction;

        double max_pressure_increment = 50000. / 60 * dt;
        // if (max_pressure_increment > 300000)
        //     max_pressure_increment = 300000;

        abs_p_max = previous_top_pressure + max_pressure_increment;
        abs_p_min = previous_top_pressure - (1.1 * max_pressure_increment);

        if (abs_p_min < 1.01325e5)
            abs_p_min = 1.01325e5;
        // if( abs_p_min <  2.e5 ) abs_p_min = 2.e5;

        number_of_target_prod_rate_iterations = 1;

        // OVERALL LOOP for target rate
        while (!exit_and_cut_dt &&
               !(well_solution_converged && water_table_depth_has_converged)) {
            well_solution_converged = false;

            n = 1; // reset number of iterations

            ////////////////////////////

            cerr << endl << "----------------------------------";
            cerr << endl << "----------------------------------";
            cerr << endl
                 << "In Target prod rate, PRE ITERATION: "
                 << number_of_target_prod_rate_iterations;
            cerr << endl << "TR MODE N-R iteration number " << n;
            cerr << endl << "Timestep: " << dt;
            cerr << endl
                 << "Previous cycle               : Top pressure: "
                 << previous_cycle_top_pressure
                 << ", Rate: " << previous_cycle_total_rate;
            cerr << endl
                 << "Current cycle, last converged: Top pressure: "
                 << current_cycle_last_converged_top_pressure
                 << ", Rate: " << current_cycle_last_converged_total_rate;
            cerr << endl
                 << "Current cycle                : Top pressure: " << top_pressure
                 << ", Rate: " << total_rate;
            cerr << endl
                 << "Inflexion_pressure: " << inflexion_pressure
                 << ", pmax: " << abs_p_max << ", pmin: " << abs_p_min;
            cerr << endl << "Target rate: " << target_rate << endl;

            ////////////////////////////
            // make the top pressure correction
            if (number_of_target_prod_rate_iterations == 1) {
                proposed_pressure_correction = pressure_correction = 0.;
            }

            top_pressure += pressure_correction;

            cerr << "NEW Top pressure: " << top_pressure << endl;

            if (top_pressure < abs_p_min) {
                top_pressure = abs_p_min;
                cerr << "Total pressure decrease exeeds max: Top pressure = p_min = "
                     << top_pressure << endl;
            }

            if (top_pressure > abs_p_max) {
                top_pressure = abs_p_max;
                cerr << "Total pressure increase exeeds max: Top pressure = p_max = "
                     << top_pressure << endl;
            }

            cerr << "NEW Top pressure: " << top_pressure << endl;

            ///// calculate the actual pressure correction applied
            if (number_of_target_prod_rate_iterations > 1) {
                pressure_correction =
                        top_pressure - current_cycle_last_converged_top_pressure;
            }

            ////////////////////////////
            ////////////////////////////
            // Now go into classic NR loop
            while (!exit_and_cut_dt &&
                   !(well_solution_converged && water_table_depth_has_converged)) {
                well_solution_converged = false;

                number_of_iterations = n;
                // cerr<<n;

                //                cerr<<endl<<"\033[1m\033[36mN-R iteration number
                //                "<<number_of_iterations<<", going to Solve and update.
                //                Damp = "<<Damp<<", counter = "<<COUNTER<<
                //                      ", counter good conv:
                //                      "<<COUNTER_GOOD_CONVERGENCE<<"\033[0m"<<endl;

                if (!exit_and_cut_dt) {
                    cerr << endl << "----------------------------------";
                    cerr << endl
                         << "In Target prod rate, ITERATION: "
                         << number_of_target_prod_rate_iterations;
                    cerr << ". TR MODE N-R iteration number " << number_of_iterations;
                    cerr << endl << "Timestep: " << dt;
                    if (VO_D_well_sat_vapor[0] > 0.)
                        cerr << endl << "Vapor present at top of well!!";
                    cerr << endl
                         << "Previous cycle               : Top pressure: "
                         << previous_cycle_top_pressure
                         << ", Rate: " << previous_cycle_total_rate;
                    cerr << endl
                         << "Current cycle, last converged: Top pressure: "
                         << current_cycle_last_converged_top_pressure
                         << ", Rate: " << current_cycle_last_converged_total_rate;
                    cerr << endl
                         << "Current cycle                : Top pressure: "
                         << top_pressure << ", Rate: " << total_rate;
                    cerr << endl
                         << "Inflexion_pressure: " << inflexion_pressure
                         << ", pmax: " << abs_p_max << ", pmin: " << abs_p_min;
                    cerr << endl << "Target rate: " << target_rate;
                    cerr << endl
                         << "++++++++++++++++"
                         << "Pressure correction: " << pressure_correction
                         << "++++++++++++++++";
                    cerr << endl << "NEW Top pressure: " << top_pressure;
                    cerr << endl << "water table depth: " << water_table_depth;
                    cerr << endl
                         << "previous water table depth: " << previous_water_table_depth
                         << endl;

                    if (!exit_and_cut_dt) {
                        SetupResidualVectorAndJacobian();

                        if (!density_converged) {
                            well_solution_converged = false;
                        }

                        if (!(well_solution_converged && water_table_depth_has_converged) &&
                                !exit_and_cut_dt) {
                            NR_Advance_well_variables(); // ADDED JUNE 2021

                            SolveWell();

                            NR_update_solution();

                            if (!exit_and_cut_dt) {
                                Velocity_model();

                                Interpolate_pressure();

                                Compute_Inflow_Outflow();

                                // if( with_water_table_displacement ) Find_Water_Table(true);
                                Find_Water_Table(with_water_table_displacement);

                                T_equilibration();

                                // Compute_radial_heat();

                                if (with_friction)
                                    ComputeFrictionFactor();

                                ComputeDensityDerivatives();

                                Weighted_Next_Guess(0.5, 0.5);

                                density_converged = Check_density_and_rate_convergence();
                            }
                        }
                    }
                }

                if (number_of_iterations > 10000) {
                    if (!well_solution_converged)
                        exit_and_cut_dt = true;
                }

                n++;
                total_iterations++;

                cerr << endl << message;
                cerr << endl << "well_solution_converged: " << well_solution_converged;
                cerr << endl
                     << "water_table_depth_has_converged: "
                     << water_table_depth_has_converged;
                cerr << endl
                     << "density_and_rate_has_converged: " << density_converged << endl;

            } // exit or converged (well)
            ////////////////////////////
            ////////////////////////////

            cerr << endl << " !!!!! Stopped for top pressure: " << top_pressure;
            COUNTER = 0;
            COUNTER_GOOD_CONVERGENCE = 0;

            // if( (target_rate > 0. && total_rate < 0.) or (target_rate < 0. &&
            // total_rate > 0.) )
            // {
            //     cerr<<endl<<"Exiting due to reversed rate"<<endl;
            //     exit_and_cut_dt = true;
            //     well_solution_converged = false;
            // }

            if (exit_and_cut_dt) {
                cerr << endl << "CRASH!!, Reseting";
                cerr << endl
                     << "In Target prod rate, ITERATION: "
                     << number_of_target_prod_rate_iterations;
                cerr << endl << "TR MODE N-R iteration number " << number_of_iterations;

                cerr << endl << "Timestep: " << dt;
                cerr << endl
                     << "Previous cycle               : Top pressure: "
                     << previous_cycle_top_pressure
                     << ", Rate: " << previous_cycle_total_rate;
                cerr << endl
                     << "Current cycle, last converged: Top pressure: "
                     << current_cycle_last_converged_top_pressure
                     << ", Rate: " << current_cycle_last_converged_total_rate;
                cerr << endl
                     << "Current cycle, not converged: Top pressure: " << top_pressure
                     << ", Rate: " << total_rate;
                cerr << endl
                     << "Inflexion_pressure: " << inflexion_pressure
                     << ", pmax: " << abs_p_max << ", pmin: " << abs_p_min;
                cerr << endl << "Target rate: " << target_rate;
                cerr << endl
                     << "++++++++++++++++"
                     << "Pressure correction: " << pressure_correction
                     << "++++++++++++++++";
                cerr << endl << "water table depth: " << water_table_depth << endl;
            }

            if (well_solution_converged) { //
                cerr << " !!!!! Well solution converged for top pressure: "
                     << top_pressure << endl;

                // DAMP_JACOBIAN *= 1.001;
                DAMP_JACOBIAN = 1.;
                // DAMP *= 1.0001;
                // if( DAMP < 0.05 ) DAMP = 0.05; //TEST FEB 2022

                well_solution_converged =
                        false; // we reset to false until we make further checks on the
                // overlying convergence loop

                /////////////////////////////////////////////////////////
                // Updating bounds
                if (number_of_target_prod_rate_iterations > 2) {
                    if ((total_rate > target_rate) &&
                            (current_cycle_last_converged_total_rate < target_rate))
                        inflexion_pressure = current_cycle_last_converged_top_pressure;

                    if ((total_rate < target_rate) &&
                            (current_cycle_last_converged_total_rate > target_rate))
                        inflexion_pressure = current_cycle_last_converged_top_pressure;
                }

                ///////////////////////////////////////////////////////////
                //  convergence loop termination conditions

                if (/*abs( ( total_rate - target_rate ) / target_rate ) < 0.1 or*/ abs(
                            total_rate - target_rate) < 0.3) {
                    cerr << "Total rate is close enough to target. EXIT" << endl;
                    well_solution_converged = true;
                }

                if (number_of_target_prod_rate_iterations < 4 &&
                        abs(total_rate - target_rate) > 0.05) {
                    cerr << "Less than 4 target rate iterations done and we could be "
                            "closer to target, we keep on going..."
                         << endl;
                    well_solution_converged = false;
                }
                if (number_of_target_prod_rate_iterations < 2) {
                    cerr << "Less than 2 target rate iterations done, we keep on going..."
                         << endl;
                    well_solution_converged = false;
                }

                //                if( total_rate < target_rate && ( abs( top_pressure
                //                - 1.01325e5 ) < 200. ) )
                //                {
                //                    cerr<<"Pressure cannot go below 1 bar.
                //                    EXIT"<<endl;  well_solution_converged = true;
                //                }
                /////////////////////////////

                cerr << " !!!!!Target rate converged for top pressure: " << top_pressure
                     << endl;

                // cerr<<endl<<"In Target prod rate, ITERATION:
                // "<<number_of_target_prod_rate_iterations; cerr<<endl<<"TR MODE N-R
                // iteration number "<<number_of_iterations;

                // cerr<<endl<<"Timestep: "<<dt;
                // cerr<<endl<<"Previous cycle               : Top pressure:
                // "<<previous_cycle_top_pressure<<", Rate:
                // "<<previous_cycle_total_rate; cerr<<endl<<"Current cycle, last
                // converged: Top pressure:
                // "<<current_cycle_last_converged_top_pressure<<", Rate:
                // "<<current_cycle_last_converged_total_rate; cerr<<endl<<"Current
                // cycle, just converged: Top pressure: "<<top_pressure<<", Rate:
                // "<<total_rate; cerr<<endl<<"Inflexion_pressure:
                // "<<inflexion_pressure<<", pmax: "<<abs_p_max<<", pmin: "<<abs_p_min;
                // cerr<<endl<<"Target rate: "<<target_rate;
                // cerr<<endl<<"++++++++++++++++"<<"Pressure correction:
                // "<<pressure_correction<<"++++++++++++++++"; cerr<<endl<<"water table
                // depth: "<<water_table_depth<<endl;

                // determine a pressure correction to reach target rate
                proposed_pressure_correction2 =
                        0.7 * (total_rate - target_rate) / total_injectivity * 100000;
                proposed_pressure_correction = proposed_pressure_correction2;
                cerr << "Proposed pressure correction based on index: "
                     << proposed_pressure_correction2 << endl;

                if ((abs(top_pressure - abs_p_min) < 50. &&
                     proposed_pressure_correction2 < 0.) or
                        (abs(top_pressure - abs_p_max) < 50. &&
                         proposed_pressure_correction2 > 0.)) {
                    cerr << "Close to min or max allowed pressure. EXIT" << endl;
                    well_solution_converged = true;
                }

                if (number_of_target_prod_rate_iterations > 250) {
                    cerr << "More than 250 target rate iterations done, we stop. EXIT"
                         << endl;
                    well_solution_converged = true;
                }

                // if( number_of_target_prod_rate_iterations == 1 )
                // {
                //     if( ( ( total_rate > previous_cycle_total_rate ) && (
                //     top_pressure < previous_cycle_top_pressure ) ) or
                //             ( ( total_rate < previous_cycle_total_rate ) && (
                //             top_pressure > previous_cycle_top_pressure ) ) )
                //     {// we are approaching the target
                //         proposed_pressure = 0.1 * ( previous_cycle_top_pressure -
                //         top_pressure ) / ( previous_cycle_total_rate - total_rate ) *
                //         ( target_rate - total_rate ) + top_pressure;
                //         proposed_pressure_correction = proposed_pressure -
                //         top_pressure;

                //         cerr<<"Pressure correction based on extrapolation:
                //         "<<proposed_pressure_correction<<endl; if(
                //         abs(proposed_pressure_correction) > abs(
                //         proposed_pressure_correction2 ) )
                //             proposed_pressure_correction =
                //             proposed_pressure_correction2;
                //     }
                // }

                // else
                // {
                //     if( ( ( total_rate > current_cycle_last_converged_total_rate ) &&
                //     ( top_pressure < current_cycle_last_converged_top_pressure ) ) or
                //             ( ( total_rate < current_cycle_last_converged_total_rate
                //             ) && ( top_pressure >
                //             current_cycle_last_converged_top_pressure ) ) )
                //     {// we are approaching the target
                //         proposed_pressure = 0.1 * (
                //         current_cycle_last_converged_top_pressure - top_pressure ) /
                //         ( current_cycle_last_converged_total_rate - total_rate ) * (
                //         target_rate - total_rate ) + top_pressure;
                //         proposed_pressure_correction = proposed_pressure -
                //         top_pressure;

                //         cerr<<"Pressure correction based on extrapolation:
                //         "<<proposed_pressure_correction<<endl; if(
                //         abs(proposed_pressure_correction) > abs(
                //         proposed_pressure_correction2 ) )
                //             proposed_pressure_correction =
                //             proposed_pressure_correction2;
                //     }
                // }

                pressure_correction = proposed_pressure_correction;

                // if ( total_rate > current_cycle_last_converged_total_rate &&
                // pressure_correction > 50000 )  pressure_correction = 50000; else if (
                // pressure_correction > 1000 )  pressure_correction = 1000;

                // if ( pressure_correction < -10000)  pressure_correction = -10000;

                //                if ( pressure_correction > 500/60*dt )
                //                {
                //                    pressure_correction = 500/60*dt;
                //                    cerr<<"Pressure correction limited by
                //                    timestep"<<endl;
                //                }
                //                if ( pressure_correction < -1500/60*dt)
                //                {
                //                    pressure_correction = -1500/60*dt;
                //                    cerr<<"Pressure correction limited by
                //                    timestep"<<endl;
                //                }
                ////////////////////////////////////////////////////////////////

                // if( top_pressure > inflexion_pressure && top_pressure - 10000. <
                // inflexion_pressure)
                // {
                //     if ( pressure_correction < -250)  pressure_correction = -250;
                // }

                // if( top_pressure < inflexion_pressure && top_pressure + 10000. >
                // inflexion_pressure)
                // {
                //     if ( pressure_correction > 250 )  pressure_correction = 250;
                // }

                if (VO_D_well_sat_vapor[0] > 0.) {
                    if (pressure_correction < -500)
                        pressure_correction = -500;
                    if (pressure_correction > 500)
                        pressure_correction = 500;
                }
                else {
                    if (pressure_correction < -5000)
                        pressure_correction = -5000;
                    if (pressure_correction > 5000)
                        pressure_correction = 5000;
                }
                if (number_of_target_prod_rate_iterations > 10) {
                    if (pressure_correction < -2500)
                        pressure_correction = -2500;
                    if (pressure_correction > 2500)
                        pressure_correction = 2500;
                }

                if(number_of_target_prod_rate_iterations > 2  && (abs(pressure_correction) < 100.))
                {
                    cerr << "Pressure correction < 100. EXIT" << endl;
                    well_solution_converged = true;
                }

                cerr << "++++++++++++++++"
                     << "NEW Pressure correction: " << pressure_correction
                     << "++++++++++++++++" << endl;

                current_cycle_last_converged_total_rate = total_rate;
                current_cycle_last_converged_top_pressure = top_pressure;

                if (number_of_iterations >
                        highest_ever_number_of_iterations_to_converge)
                    highest_ever_number_of_iterations_to_converge = number_of_iterations;

            } // end converged target

            if (well_solution_converged && water_table_depth_has_converged &&
                    exit_and_cut_dt) {
                cerr << endl << "converged && exit_and_cut_dt" << endl;
                // PAUSE();
                exit_and_cut_dt = true;
            }

            number_of_target_prod_rate_iterations++;
        } // exit or converged target

        cerr << endl
             << "----------------------------------" << endl
             << "Target prod rate FINISHED at iteration  "
             << number_of_target_prod_rate_iterations - 1;
        cerr << endl
             << "Current rate: " << total_rate << ", target rate: " << target_rate;
        cerr << endl
             << "Previous cycle               : Top pressure: "
             << previous_cycle_top_pressure
             << ", Rate: " << previous_cycle_total_rate;
        cerr << endl
             << "Current cycle, last converged: Top pressure: "
             << current_cycle_last_converged_top_pressure
             << ", Rate: " << current_cycle_last_converged_total_rate;
        cerr << endl
             << "Current cycle, FINAL         : Top pressure: " << top_pressure
             << ", Rate: " << total_rate;
        cerr << endl << "Actual top pressure: " << VO_D_well_pressure[0];
        cerr << endl << "Target rate: " << target_rate << endl;
    }

    else {
        number_of_target_prod_rate_iterations = 99;

        while (!exit_and_cut_dt &&
               !(well_solution_converged && water_table_depth_has_converged)) {
            number_of_iterations = n;
            // cerr<<n<<" ";

            //            cerr<<"\033[1m\033[36mN-R iteration number
            //            "<<number_of_iterations<<", going to Solve and update. Damp
            //            = "<<Damp<<", counter = "<<COUNTER<<
            //                  ", counter good conv:
            //                  "<<COUNTER_GOOD_CONVERGENCE<<"\033[0m"<<endl;

            if (!exit_and_cut_dt) {
                SetupResidualVectorAndJacobian();

                if (!density_converged) {
                    well_solution_converged = false;
                }

                if (!(well_solution_converged && water_table_depth_has_converged) &&
                        !exit_and_cut_dt) {
                    NR_Advance_well_variables(); // ADDED JUNE 2021

                    SolveWell();

                    NR_update_solution();

                    if (!exit_and_cut_dt) {
                        Velocity_model();

                        Interpolate_pressure();

                        Compute_Inflow_Outflow();

                        // if( with_water_table_displacement) Find_Water_Table(true);
                        Find_Water_Table(with_water_table_displacement);

                        T_equilibration();

                        // Compute_radial_heat();

                        if (with_friction)
                            ComputeFrictionFactor();

                        ComputeDensityDerivatives();

                        Weighted_Next_Guess(0.5, 0.5);

                        density_converged = Check_density_and_rate_convergence();
                    }
                }
            }

            if (number_of_iterations > 10000) {
                if (!well_solution_converged)
                    exit_and_cut_dt = true;
            }

            n++;
            total_iterations++;
            // cerr<<endl<<message;
            //             cerr<<endl<<"well_solution_converged:
            //             "<<well_solution_converged;
            //             cerr<<endl<<"water_table_depth_has_converged:
            //             "<<water_table_depth_has_converged;
            //             cerr<<endl<<"density_and_rate_has_converged:
            //             "<<density_converged<<endl;
        }

        COUNTER = 0;
        COUNTER_GOOD_CONVERGENCE = 0;
    }

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    //if (exit_and_cut_dt) {
    cerr << endl
         << "N-R well solve EXIT after " << number_of_iterations
         << " iterations..." << endl;
    cerr << endl << "Damp : " << DAMP;
    cerr << endl
         << "Worst momentum,    abs : " << worst_value_p
         << " // rel : " << worst_value_p_;
    cerr << endl
         << "Worst continuity,  abs : " << worst_value_v
         << " // rel : " << worst_value_v_;
    cerr << endl
         << "Worst energy,      abs : " << worst_value_h
         << " // rel : " << worst_value_h_;
    cerr << endl
         << "Worst salt,        abs : " << worst_value_smf
         << " // rel : " << worst_value_smf_;
    //}
}

/** Solves the linearised well system J * dx = -R for one Newton iteration.

    The backend is selected at construction (Eigen SparseLU by default) — see
    CVFEM_SolverChoice.h. The old inline block of commented-out SAMG tuning
    calls was removed: SAMG parameters are documented and set for real in
    CVFEM_PHX_Scheme's constructor, and duplicating them
    here as dead comments only invited the two copies to drift apart.

    Benoit DD/MM/YYYY */
template <uint32_t dim> void WellModelPrototype<dim>::SolveWell() {

    // Four unknowns per well segment: pressure, fluid velocity, enthalpy, salt.
    const int nunknowns(4);

    well_bundle_.Get().Solve(Jacobian, Residuals, Sol_increment, nunknowns);
}

template <uint32_t dim> void WellModelPrototype<dim>::NR_update_solution() {
    std::fill(Sol_increment_pressure.begin(), Sol_increment_pressure.end(), 0.);
    std::fill(Sol_increment_velocity_fluid.begin(),
              Sol_increment_velocity_fluid.end(), 0.);
    std::fill(Sol_increment_s_enthalpy_bulk.begin(),
              Sol_increment_s_enthalpy_bulk.end(), 0.);
    std::fill(Sol_increment_smf_bulk.begin(), Sol_increment_smf_bulk.end(), 0.);
    //    std::fill( Sol_increment_s_enthalpy_fluid.begin(),
    //    Sol_increment_s_enthalpy_fluid.end(), 0. ); std::fill(
    //    Sol_increment_smf_fluid.begin(), Sol_increment_smf_fluid.end(), 0. );

    if (previous_worst_relative_residual > worst_relative_residual)
        COUNTER_GOOD_CONVERGENCE++;
    else
        COUNTER_GOOD_CONVERGENCE = 0;

    if (COUNTER_GOOD_CONVERGENCE % 30 == 0 && COUNTER_GOOD_CONVERGENCE > 1) {
        DAMP *= 1.001;
    }

    COUNTER++;
    if (COUNTER == 100) {
        if (COUNTER_GOOD_CONVERGENCE < 5) {
            DAMP *= 0.95;
        } // Oct 2024: Was 0.9

        COUNTER = 0;
    }

    if (DAMP > 0.2)
        DAMP = 0.2;
    if (DAMP < 0.005)
        DAMP = 0.005; // Oct 2024 change from 0.001

    //        if( DAMP > 0.25 ) DAMP = 0.25; //TEST NOV 2021
    //        if( number_of_iterations > 500 && DAMP > 0.2 ) DAMP = 0.2;
    //        if( number_of_iterations > 1000 && DAMP > 0.15 ) DAMP = 0.15;
    //        if( number_of_iterations > 2000 && DAMP > 0.1 ) DAMP = 0.1;
    //        if( number_of_iterations > 4000 && DAMP > 0.05 ) DAMP = 0.05;
    //        if( number_of_iterations > 6000 && DAMP > 0.01 ) DAMP = 0.01;
    //        if( number_of_iterations > 10000 && DAMP > 0.005 ) DAMP = 0.005;

    //    if( number_of_iterations < 100 && DAMP < 0.05 )
    //        Damp = 0.05;
    //    else
    Damp = DAMP;
    largest_increment = 0.;

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        unsigned k(n_of_eq * n);
        Sol_increment_pressure[n] = Sol_increment[k];
        Sol_increment_velocity_fluid[n] = Sol_increment[k + 1];
        Sol_increment_s_enthalpy_bulk[n] = Sol_increment[k + 2];
        Sol_increment_smf_bulk[n] = Sol_increment[k + 3];

        if (abs(Sol_increment_pressure[n]) > largest_increment)
            largest_increment = abs(Sol_increment_pressure[n]);
        // if(abs(Sol_increment_velocity_fluid [ n ]) > largest_increment)
        // largest_increment = abs(Sol_increment_velocity_fluid [ n ]);
        if (abs(Sol_increment_s_enthalpy_bulk[n]) > largest_increment)
            largest_increment = abs(Sol_increment_s_enthalpy_bulk[n]);
        // if(abs(Sol_increment_smf_bulk [ n ]) > largest_increment)
        // largest_increment = abs(Sol_increment_smf_bulk [ n ]);

        // TEST MARCH 2022

        if (Sol_increment_s_enthalpy_bulk[n] > 0.05 * VO_D_well_s_enthalpy_bulk[n])
            Sol_increment_s_enthalpy_bulk[n] = 0.05 * VO_D_well_s_enthalpy_bulk[n];
        if (Sol_increment_s_enthalpy_bulk[n] < -0.05 * VO_D_well_s_enthalpy_bulk[n])
            Sol_increment_s_enthalpy_bulk[n] = -0.05 * VO_D_well_s_enthalpy_bulk[n];

        if (Sol_increment_pressure[n] > 0.05 * VO_D_well_pressure[n])
            Sol_increment_pressure[n] = 0.05 * VO_D_well_pressure[n];
        if (Sol_increment_pressure[n] < -0.05 * VO_D_well_pressure[n])
            Sol_increment_pressure[n] = -0.05 * VO_D_well_pressure[n];

        if (Sol_increment_smf_bulk[n] > 0.005)
            Sol_increment_smf_bulk[n] = 0.005;
        if (Sol_increment_smf_bulk[n] < -0.005)
            Sol_increment_smf_bulk[n] = -0.005;
    }

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        VO_D_well_velocity_fluid[n] += Sol_increment_velocity_fluid[n];
    }

    for (int n = number_of_well_segments - 1; n >= 0; --n) {
        bool verbose(false);

        // if(n<3) verbose = true;

        if (verbose) {
            cerr << endl << "--------------------------";
            cerr << endl << "Previous to correction, Segment: " << n;
            cerr << endl
                 << "p = "
                 << VO_D_well_pressure[n] /*<<", v = "<<VO_D_well_velocity_fluid   [ n
                                                       ]*/
                 << ", h = " << VO_D_well_s_enthalpy_bulk[n]
                    << ", smf = " << VO_D_well_smf_bulk[n];
        }

        // PRESSURE
        VO_D_well_pressure[n] += Sol_increment_pressure[n] * Damp;

        // ENTHALPY
        VO_D_well_s_enthalpy_bulk[n] += Sol_increment_s_enthalpy_bulk[n] * Damp;

        // SALT
        VO_D_well_smf_bulk[n] += Sol_increment_smf_bulk[n] * 0.1;

        /////////////////////
        if (VO_D_well_smf_bulk[n] < 0.)
            VO_D_well_smf_bulk[n] = 0.;

        if (VO_D_well_smf_bulk[n] > 0.9999999) {
            cerr << endl
                 << " The bulk salt mass fraction is above  0.9999999 in segment "
                 << n << ". VO_D_well_smf_bulk [ n ]: " << VO_D_well_smf_bulk[n]
                    << endl;
            VO_D_well_smf_bulk[n] = 0.9999999;
            // pause_at_end = true;
            message_at_end = "The salt mass fraction has been set to 0.9999999";
        }

        if (VO_D_well_sat_halite[n] > 0.) {
            if (VO_D_well_s_enthalpy_bulk[n] < VO_D_well_s_enthalpy_halite[n]) {
                // cerr<<endl<<" The bulk enthalpy is below halite enthalpy in segment
                // "<<n<<". VO_D_well_s_enthalpy_bulk [ n ]:
                // "<<VO_D_well_s_enthalpy_bulk [ n ]; pause_at_end = true;
                VO_D_well_s_enthalpy_bulk[n] = VO_D_well_s_enthalpy_halite[n];
                message_at_end = "The bulk enthalpy has been below halite enthalpy";
            }
        }

        // BLOCKED SEGMENT
        if (VO_well_segment_activated[n] < 1) {
            VO_D_well_velocity_fluid[n] = 0.;
        }

        // WT SEGMENT
        if (VO_D_well_water_table[n] == 1) {
            VO_D_well_pressure[n] = top_pressure;
            VO_D_well_velocity_fluid[n] = velocity_wt;
            // VO_D_well_velocity_fluid   [ n ]  = VO_D_well_velocity_fluid   [ n + 1
            // ]; VO_D_well_s_enthalpy_bulk  [ n ] = VO_D_well_s_enthalpy_bulk  [ n +
            // 1 ];
            VO_D_well_smf_bulk[n] = VO_D_well_smf_bulk[n + 1];
        }

        if (VO_D_well_water_table[n] == 1 or
                VO_D_previous_well_water_table[n] < 2) { // TEST AUGUST 2024
            if (injection_mode or VO_D_well_velocity_fluid[n] > 0.) {
                VO_D_well_s_enthalpy_bulk[n] = top_s_enthalpy;
                VO_D_well_smf_bulk[n] = top_smf;
            }
        }

        // EMPTY SEGMENTS
        if (VO_D_well_water_table[n] == 0) {
            VO_D_well_velocity_fluid[n] = velocity_wt;
            // VO_D_well_velocity_fluid   [ n ]  = VO_D_well_velocity_fluid   [ n + 1
            // ]; VO_D_well_s_enthalpy_fluid  [ n ] = 0.; VO_D_well_smf_fluid [ n ] =
            // 0.;
            VO_D_well_s_enthalpy_bulk[n] = 0.; // not sure...
            VO_D_well_smf_bulk[n] = 0.;        // not sure...
            VO_D_well_pressure[n] = top_pressure;
        }

        if (verbose) {
            cerr << endl << endl << "Segment: " << n;
            cerr << endl
                 << "Density fluid before: " << VO_D_previous_well_density_fluid[n]
                    << ", density fluid: " << VO_D_well_density_fluid[n];
            cerr << endl
                 << "p = " << VO_D_well_pressure[n]
                    << ", v = " << VO_D_well_velocity_fluid[n]
                       << ", h = " << VO_D_well_s_enthalpy_bulk[n]
                          << ", smf = " << VO_D_well_smf_bulk[n];
            cerr << endl
                 << "Incr_p = " << Sol_increment_pressure[n]
                    << ", Incr_v = " << Sol_increment_velocity_fluid[n]
                       << ", Incr_h = " << Sol_increment_s_enthalpy_bulk[n]
                          << ", Incr smf = " << Sol_increment_smf_bulk[n];
        }

        if (VO_D_well_pressure[n] < 1.e5) {
            cerr << endl
                 << "Segment : " << n
                 << " Pressure inferior to 1 bar in the N-R loop.";
            VO_D_well_pressure[n] = 1.01325e5;
            // exit_and_cut_dt = true;
            // pause_at_end = true;
        }

        if (VO_D_well_smf_fluid[n] > 1.) {
            cerr << endl
                 << "Segment : " << n
                 << " Salt mass fraction above 100% in the N-R loop, resetting..";
            exit_and_cut_dt = true;
        }

        smf_ = VO_D_well_smf_bulk[n];
        wt_ = smf_ * 100.;
        x_ = Weight2XNaCl(wt_);

        p_ = VO_D_well_pressure[n];

        h_fluid_ = VO_D_well_s_enthalpy_bulk[n];

        t_ = 5.; // HARD LOW TEMP

        // double hmin=brine.Enthalpy();
        double hmin = fluid.BulkProperties().h;

        if (h_fluid_ < hmin && VO_D_well_water_table[n] > 0) {
            cerr << endl
                 << "Segment : " << n << "Looks like temp below min temp (" << t_
                 << ") in the N-R loop, resetting.." << endl;
            exit_and_cut_dt = true;
        }

        t_ = 1000.; // HARD MAX TEMP

        // double hmax=brine.Enthalpy();
        double hmax = fluid.BulkProperties().h;
        if (h_fluid_ > hmax) {
            cerr << endl
                 << "Segment : " << n << "Looks like temp above max temp (" << t_
                 << ") in the N-R loop, resetting.." << endl;
            VO_D_well_s_enthalpy_bulk[n] = hmax;

            // if(number_of_iterations >500)
            exit_and_cut_dt = true;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// RESIDUALS /////////////////////////////////////////////////////////////////
template <uint32_t dim>
void WellModelPrototype<dim>::SetupResidualVectorAndJacobian() {
    //    if (DAMP_JACOBIAN > 0.8) DAMP_JACOBIAN = 0.8;
    //    if( number_of_iterations > 500 && DAMP_JACOBIAN > 0.7 ) DAMP_JACOBIAN =
    //    0.7; if( number_of_iterations > 1000 && DAMP_JACOBIAN > 0.6 )
    //    DAMP_JACOBIAN = 0.6; if( number_of_iterations > 2000 && DAMP_JACOBIAN >
    //    0.5 ) DAMP_JACOBIAN = 0.5;
    DAMP_JACOBIAN = 1.;
    //    if( number_of_iterations > 10000 ) DAMP_JACOBIAN = 0.5;
    //    else DAMP_JACOBIAN = 1.;

    std::fill(Residuals.begin(), Residuals.end(), 0.);

    Jacobian.Zero();
    Jacobian.Resize(number_of_well_segments * n_of_eq);

    // init
    drhofdp_ = drhofdhb_ = drhofdsmfb_ = dhdp_ = drhofdp_prev_seg_ =
            drhofdhb_prev_seg_ = drhofdsmfb_prev_seg_ = dhdp_prev_seg_ =
            drhofdp_next_seg_ = drhofdhb_next_seg_ = drhofdsmfb_next_seg_ =
            dhdp_next_seg_ = 0.;

    drhobdp_ = drhobdhb_ = drhobdsmfb_ = drhobdp_prev_seg_ = drhobdhb_prev_seg_ =
            drhobdsmfb_prev_seg_ = drhobdp_next_seg_ = drhobdhb_next_seg_ =
            drhobdsmfb_next_seg_ = 0.;

    // unused:
    dhdsmf = dhdsmf_prev_seg = dhdsmf_next_seg = dsmfdp = dsmfdp_prev_seg =
            dsmfdp_next_seg = dsmfdh = dsmfdh_prev_seg = dsmfdh_next_seg = 0.;
    ///////////

    double test_value(0.), eps_(1.e-4), worst_absolute_residual(1.e8);
    double max_error_P(1000.), max_error_V(0.005), max_prct_error_H(0.005),
            max_error_S(0.05);

    unsigned worst_node_P_rel(0), worst_node_V_rel(0), worst_node_H_rel(0),
            worst_node_S_rel(0);
    unsigned worst_seg_P_abs(0), worst_seg_V_abs(0), worst_seg_H_abs(0),
            worst_seg_S_abs(0);
    worst_value_p = worst_value_v = worst_value_h = worst_value_smf =
            worst_value_p_ = worst_value_v_ = worst_value_h_ = worst_value_smf_ = 0.;

    if (number_of_iterations > 1000) {
        eps_ *= 5.;         // used to be *5
        max_error_P *= 10.; // used to be *10
        max_error_V *= 2.;
        max_prct_error_H *= 2.;
        max_error_S *= 2.;
    }

    for (int n = number_of_well_segments - 1; n >= 0; --n) {
        test_p = 100.;
        test_v = 100.;
        test_h = 100.;
        test_smf = 100.;

        Momentum_residual = 0.;
        Continuity_residual = 0.;
        Energy_residual = 0.;
        Salt_mass_residual = 0.;

        seg_air_saturation_scaling = 1. - VO_D_well_air_saturation[n];

        p = VO_D_well_pressure[n];
        vm = VO_D_well_velocity_fluid[n];

        hb = VO_D_well_s_enthalpy_bulk[n];
        hf = VO_D_well_s_enthalpy_fluid[n];
        hh = VO_D_well_s_enthalpy_halite[n];

        smf_f = VO_D_well_smf_fluid[n];
        smf_b = VO_D_well_smf_bulk[n];

        sath = VO_D_well_sat_halite[n];
        sath_prev_time = VO_D_previous_well_sat_halite[n];

        rhob = VO_D_well_density_bulk[n];
        rhof = VO_D_well_density_fluid[n];
        rhoh = VO_D_well_density_halite[n];

        seg_length = VO_D_segment_length[n];
        seg_crossarea =
                VO_D_well_cross_area[n]; // should use effective area in some places!
        seg_inclination = VO_D_segment_inclination[n];

        radius = sqrtf(VO_D_well_cross_area[n] *
                       (1 - VO_D_previous_well_sat_halite[n]) / CSMP_PI);
        friction_factor = VO_D_well_friction_factor[n] / (4 * radius);

        drhofdp_ = VO_D_drhofdp[n];
        drhofdhb_ = VO_D_drhofdhb[n];
        drhofdsmfb_ = VO_D_drhofdsmfb[n];
        drhobdp_ = VO_D_drhobdp[n];
        drhobdhb_ = VO_D_drhobdhb[n];
        drhobdsmfb_ = VO_D_drhobdsmfb[n];
        dhdp_ = VO_D_dhdp[n];

        p_prev_time = VO_D_previous_well_pressure[n];
        vm_prev_time = VO_D_previous_well_velocity_fluid[n];

        hb_prev_time = VO_D_previous_well_s_enthalpy_bulk[n];
        hf_prev_time = VO_D_previous_well_s_enthalpy_fluid[n];
        hh_prev_time = VO_D_previous_well_s_enthalpy_halite[n];

        smf_f_prev_time = VO_D_previous_well_smf_fluid[n];
        smf_b_prev_time = VO_D_previous_well_smf_bulk[n];

        rhob_prev_time = VO_D_previous_well_density_bulk[n];
        rhof_prev_time = VO_D_previous_well_density_fluid[n];
        rhoh_prev_time = VO_D_previous_well_density_halite[n];

        prev_time_seg_crossarea = VO_D_previous_well_cross_area[n];

        if (n > 0) {
            p_prev_seg = VO_D_well_pressure[n - 1];
            vm_prev_seg = VO_D_well_velocity_fluid[n - 1];

            prev_seg_air_saturation_scaling = 1 - VO_D_well_air_saturation[n - 1];

            vm_prev_time_prev_seg = VO_D_previous_well_velocity_fluid[n - 1];

            hb_prev_seg = VO_D_well_s_enthalpy_bulk[n - 1];
            hf_prev_seg = VO_D_well_s_enthalpy_fluid[n - 1];

            smf_f_prev_seg = VO_D_well_smf_fluid[n - 1];

            sath_prev_seg = VO_D_well_sat_halite[n - 1];
            sath_prev_time_prev_seg = VO_D_previous_well_sat_halite[n - 1];

            rhof_prev_seg = VO_D_well_density_fluid[n - 1];
            rhof_prev_time_prev_seg = VO_D_previous_well_density_fluid[n - 1];

            radius = sqrtf(VO_D_well_cross_area[n - 1] *
                           (1 - VO_D_previous_well_sat_halite[n - 1]) / CSMP_PI);
            friction_factor_prev_seg =
                    VO_D_well_friction_factor[n - 1] / (4 * radius);

            prev_seg_length = VO_D_segment_length[n - 1];
            prev_seg_crossarea = VO_D_well_cross_area[n - 1];
            prev_seg_inclination = VO_D_segment_inclination[n - 1];

            drhofdp_prev_seg_ = VO_D_drhofdp[n - 1];
            drhofdhb_prev_seg_ = VO_D_drhofdhb[n - 1];
            drhofdsmfb_prev_seg_ = VO_D_drhofdsmfb[n - 1];
            drhobdp_prev_seg_ = VO_D_drhobdp[n - 1];
            drhobdhb_prev_seg_ = VO_D_drhobdhb[n - 1];
            drhobdsmfb_prev_seg_ = VO_D_drhobdsmfb[n - 1];
            dhdp_prev_seg_ = VO_D_dhdp[n - 1];
        }

        if (n < number_of_well_segments - 1) {
            p_next_seg = VO_D_well_pressure[n + 1];
            vm_next_seg = VO_D_well_velocity_fluid[n + 1];

            hf_next_seg = VO_D_well_s_enthalpy_fluid[n + 1];
            hb_next_seg = VO_D_well_s_enthalpy_bulk[n + 1];

            smf_f_next_seg = VO_D_well_smf_fluid[n + 1];

            sath_next_seg = VO_D_well_sat_halite[n + 1];
            sath_prev_time_next_seg = VO_D_previous_well_sat_halite[n + 1];

            rhof_next_seg = VO_D_well_density_fluid[n + 1];

            next_seg_length = VO_D_segment_length[n + 1];
            next_seg_crossarea = VO_D_well_cross_area[n + 1];
            next_seg_inclination = VO_D_segment_inclination[n + 1];

            drhofdp_next_seg_ = VO_D_drhofdp[n + 1];
            drhofdhb_next_seg_ = VO_D_drhofdhb[n + 1];
            drhofdsmfb_next_seg_ = VO_D_drhofdsmfb[n + 1];
            drhobdp_next_seg_ = VO_D_drhobdp[n + 1];
            drhobdhb_next_seg_ = VO_D_drhobdhb[n + 1];
            drhobdsmfb_next_seg_ = VO_D_drhobdsmfb[n + 1];
            dhdp_next_seg_ = VO_D_dhdp[n + 1];
        }

        // velocity signs
        velocity_sign = velocity_sign_prev_seg = 0.;
        if (vm != 0.)
            velocity_sign = vm / (abs(vm));
        if (vm_prev_seg != 0.)
            velocity_sign_prev_seg = vm_prev_seg / (abs(vm_prev_seg));

        bool verbose(false);
        // if(n>45) verbose = true;

        // if(n<3) verbose = true;

        if (verbose) {
            cerr << endl << "--------------------------" << endl;
            cerr << "Segment: " << n << endl;
            //            cerr<<endl<<"WATER TABLE INDEX: "<<VO_D_well_water_table [ n
            //            ]; cerr<<endl<<"PREVIOUS WATER TABLE INDEX:
            //            "<<VO_D_previous_well_water_table [ n ]; cerr<<endl<<"NR
            //            PREVIOUS WATER TABLE INDEX:
            //            "<<VO_D_NR_previous_well_water_table [ n ];
        }
        Compute_Momentum_Residuals(n, false);
        Compute_Continuity_Residuals(n, false);
        Compute_Energy_Residuals(n, false);
        Compute_Salt_Mass_Residuals(n, verbose);

        m_fluid_ = seg_length * seg_crossarea * seg_air_saturation_scaling *
                ((1. - sath) * rhof + sath * rhoh);
        H_current_ = m_fluid_ * hb;

        double internal_energy(0.), water_mass(0.);
        internal_energy =
                seg_length * seg_crossarea * seg_air_saturation_scaling *
                ((1. - sath) * (rhof * (hf + vm * vm * 0.5) - p) + sath * rhoh * hh);
        water_mass = seg_length * seg_crossarea * seg_air_saturation_scaling *
                (1. - sath) * rhof;

        // Check convergence
        ///////////////////////////////
        ///// Normalized residuals criterion:
        //

        if (abs(Momentum_residual / VO_D_well_pressure[n]) < eps_)
            test_p = 0.;

        if (m_fluid_ > 0.) {
            if (abs(Continuity_residual / water_mass) < eps_)
                test_v = 0.;
            if (abs(Energy_residual / internal_energy) < eps_)
                test_h = 0.;
        }

        else {
            test_v = 0.;
            test_h = 0.;
        }

        if (VO_D_well_salt_mass[n] > 0.) {
            if (abs(Salt_mass_residual / VO_D_well_salt_mass[n]) < 0.005)
                test_smf = 0.;
        }

        else
            test_smf = 0.;

        ///////////////////////////////
        //// Absolute residuals criterion:
        if (abs(Momentum_residual) < max_error_P)
            test_p = 0.;
        if (abs(Salt_mass_residual) < max_error_S)
            test_smf = 0.;

        if (m_fluid_ > 0.) {
            if (abs(Continuity_residual) < max_error_V)
                test_v = 0.;
            if (abs(Energy_residual) < max_prct_error_H * hb)
                test_h = 0.;
        }

        ///////////////////////////////
        //  solution at water table always "correct":

        if (VO_D_well_water_table[n] == 1) {
            test_p = 0.;
            test_v = 0.;
            test_smf = 0.;
            test_h = 0.;
        }

        //////////////////////////////

        if (test_p > test_value)
            test_value = test_p;
        if (test_v > test_value)
            test_value = test_v;
        if (test_h > test_value)
            test_value = test_h;
        if (test_smf > test_value)
            test_value = test_smf;

        // Check worst segments
        // abs
        if (abs(Momentum_residual) > abs(worst_value_p) && n > 0) {
            worst_value_p = Momentum_residual;
            worst_seg_P_abs = n;
        }
        if (m_fluid_ > 0. && abs(Continuity_residual) > abs(worst_value_v) &&
                n > 0) {
            worst_value_v = Continuity_residual;
            worst_seg_V_abs = n;
        }
        if (m_fluid_ > 0. && abs(Energy_residual) > abs(worst_value_h)) {
            worst_value_h = Energy_residual;
            worst_seg_H_abs = n;
        }
        if (VO_D_well_salt_mass[n] > 0. &&
                abs(Salt_mass_residual) > abs(worst_value_smf) && n > 0) {
            worst_value_smf = Salt_mass_residual;
            worst_seg_S_abs = n;
        }

        // rel
        if (abs(Momentum_residual / VO_D_well_pressure[n]) > abs(worst_value_p_) &&
                n > 0) {
            worst_value_p_ = (Momentum_residual / VO_D_well_pressure[n]);
            worst_node_P_rel = n;
        }
        if (m_fluid_ > 0. &&
                abs(Continuity_residual / water_mass) > abs(worst_value_v_) && n > 0) {
            worst_value_v_ = (Continuity_residual / water_mass);
            worst_node_V_rel = n;
        }
        if (m_fluid_ > 0. &&
                abs(Energy_residual / internal_energy) > abs(worst_value_h_)) {
            worst_value_h_ = (Energy_residual / internal_energy);
            worst_node_H_rel = n;
        }
        if (VO_D_well_salt_mass[n] > 0. &&
                abs(Salt_mass_residual / VO_D_well_salt_mass[n]) >
                abs(worst_value_smf_) &&
                n > 0) {
            worst_value_smf_ = (Salt_mass_residual / VO_D_well_salt_mass[n]);
            worst_node_S_rel = n;
        }

        //        if ( verbose )
        //        {
        //            cerr<<endl<<endl<<"res_m = "<<res_p<<", res_c = "<<res_v<<",
        //            res_e = "<<res_h<<", res_smf = "<<res_smf; cerr<<endl<<"p =
        //            "<<VO_D_well_pressure[ n ]<<", mass fluid = "<<m_fluid_ <<",
        //            energy fluid = "<<internal_energy<<", mass salt =
        //            "<<VO_D_well_salt_mass [ n ]; cerr<<endl<<"RES M =
        //            "<<Momentum_residual<<", RES C = "<<Continuity_residual<<",
        //            RES E = "<<Energy_residual<<", RES SMF =
        //            "<<Salt_mass_residual; cerr<<endl<<"R_M   =
        //            "<<Momentum_residual / VO_D_well_pressure [ n ]<<
        //                  ", R_C = "<<Continuity_residual / water_mass<<
        //                  ", R_E = "<<Energy_residual / internal_energy<<
        //                  ", R_SMF = "<<Salt_mass_residual / VO_D_well_salt_mass [
        //                  n ];
        //        }

        ////////////////////////////////////////////

        unsigned k(n_of_eq * n);
        Residuals[k] = -Momentum_residual;
        Residuals[k + 1] = -Continuity_residual;
        Residuals[k + 2] = -Energy_residual;
        Residuals[k + 3] = -Salt_mass_residual;

        Fill_momentum_Jacobian_part_salt(n);
        Fill_continuity_Jacobian_part_salt(n);
        Fill_energy_Jacobian_part_salt(n);
        Fill_salt_mass_Jacobian_part(n);

        verbose = false;
        if (verbose) {
            cerr << endl << "--------------------------" << endl;
            cerr << "Segment: " << n << endl;

            //            cerr<<endl<<"Jacobian segment momentum:     "<<endl;
            //            for ( unsigned  n = 0; n < Jacobian_segment_momentum.size();
            //            n++ )
            //            {
            //                cerr<<", "<<Jacobian_segment_momentum [ n ];
            //            }
            cerr << endl << "Jacobian segment continuity:   " << endl;
            for (unsigned n = 0; n < Jacobian_segment_continuity.size(); n++) {
                cerr << ", " << Jacobian_segment_continuity[n];
            }
            //            cerr<<endl<<"Jacobian segment energy:       "<<endl;
            //            for ( unsigned  n = 0; n < Jacobian_segment_energy.size();
            //            n++ )
            //            {
            //                cerr<<", "<<Jacobian_segment_energy [ n ];
            //            }
            cerr << endl << "Jacobian segment salt mass:    " << endl;
            for (unsigned n = 0; n < Jacobian_segment_salt_mass.size(); n++) {
                cerr << ", " << Jacobian_segment_salt_mass[n];
            }
        }
    }

    //    cerr<<"Worst abs momentum segment: "<<worst_seg_P_abs<<", residual
    //    value: "<<worst_value_p<<" // Worst rel momentum segment:
    //    "<<worst_node_P_rel<<", relative residual value:
    //    "<<worst_value_p_<<endl; cerr<<"Worst abs continuity segment:
    //    "<<worst_seg_V_abs<<", residual value: "<<worst_value_v<<" // Worst rel
    //    continuity node: "<<worst_node_V_rel<<", relative residual value:
    //    "<<worst_value_v_<<", density:
    //    "<<VO_D_well_density_fluid[worst_node_V_rel]<<endl; cerr<<"Worst abs
    //    energy segment: "<<worst_seg_H_abs<<", residual value:
    //    "<<worst_value_h<<" // Worst rel energy segment: "<<worst_node_H_rel<<",
    //    relative residual value: "<<worst_value_h_<<", density:
    //    "<<VO_D_well_density_fluid[worst_node_H_rel]<<endl; cerr<<"Worst abs
    //    salt segment: "<<worst_seg_S_abs<<", residual value:
    //    "<<worst_value_smf<<" // Worst rel salt segment: "<<worst_node_S_rel<<",
    //    relative residual value: "<<worst_value_smf_<<endl;

    previous_worst_relative_residual = worst_relative_residual;
    worst_absolute_residual =
            std::max(std::max(abs(worst_value_p), abs(worst_value_v)),
                     std::max(abs(worst_value_h), abs(worst_value_smf)));
    worst_relative_residual =
            std::max(std::max(abs(worst_value_p_), abs(worst_value_v_)),
                     std::max(abs(worst_value_h_), abs(worst_value_smf_)));

    // cerr<<"In well: "<<well_name<<". Worst relative residual:
    // "<<worst_relative_residual; cerr<<", Worst absolute residual:
    // "<<worst_absolute_residual<<endl;

    // if( worst_value_p > 1.e12 ) exit_and_cut_dt = true;

    if (largest_increment < 100.)
        cerr << "Solution of pressure and energy not changing but potentially bad "
                "convergence"
             << endl;

    if (test_value < 1.e-10 or largest_increment < 100.) {
        cerr << endl
             << "Stopping criteria met at iteration: " << number_of_iterations
             << " // Current rate: " << total_rate
             << " // Total injectivity: " << VO_well_injectivity[1];
        if (number_of_iterations > 1) {
            well_solution_converged = true;
            // cerr<<endl<<"Damp : "<<DAMP;
            // cerr<<endl<<"Worst momentum,    abs : "<<worst_value_p<<" // rel :
            // "<<worst_value_p_; cerr<<endl<<"Worst continuity,  abs :
            // "<<worst_value_v<<" // rel : "<<worst_value_v_; cerr<<endl<<"Worst
            // energy,      abs : "<<worst_value_h<<" // rel : "<<worst_value_h_;
            // cerr<<endl<<"Worst salt,        abs : "<<worst_value_smf<<" // rel :
            // "<<worst_value_smf_;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// RESIDUALS //////////////////////////////////////////////////////////////////

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Momentum_Residuals(unsigned n,
                                                         bool verbose) {

    if (VO_D_well_water_table[n] < 2) { // BC
        Momentum_residual = p - top_pressure;
    }

    else { // gravity
        Momentum_residual += (p - p_prev_seg) -
                rhof_prev_seg * grav_acc * cos(prev_seg_inclination) *
                prev_seg_length * prev_seg_air_saturation_scaling;

        if (abs(vm_prev_seg) > 1.e-6 && n > 0) {
            if (with_friction) {
                // friction
                Momentum_residual += velocity_sign_prev_seg * friction_factor_prev_seg *
                        rhof_prev_seg * vm_prev_seg * vm_prev_seg *
                        prev_seg_air_saturation_scaling * prev_seg_length;
            }

            if (with_inertial_terms_in_momentum && rhof > 0.) {
                // inertia
                if (rhof_prev_time_prev_seg > 0.) {
                    Momentum_residual +=
                            (rhof_prev_seg * vm_prev_seg -
                             rhof_prev_time_prev_seg * vm_prev_time_prev_seg) /
                            dt * prev_seg_length;
                }

                if (rhof_prev_seg > 0.) {
                    Momentum_residual +=
                            (rhof * vm * vm - rhof_prev_seg * vm_prev_seg * vm_prev_seg);
                }
            }
        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Continuity_Residuals(
        unsigned n, bool verbose) { // Need multi phase version with DFM

    // Multiphase Drift Flux Model requires the continuity equation to be
    // expressed in multiphase form Homogeneous flow model can still be expressed
    // in term of fluid
    if (verbose)
        cerr << endl << endl << "--------------------------";
    if (verbose)
        cerr << endl << "compute continuity residual";

    // (Was `if ((multiphase && homogeneous_flow) || !multiphase)`. That is true
    // unless multiphase && !homogeneous_flow, and homogeneous_flow is always
    // true, so the test never excluded anything. `multiphase` removed.
    // Benoit DD/MM/YYYY)
    if (homogeneous_flow) {

        if (VO_well_segment_activated[n] < 1) {
            Continuity_residual = vm;
        }

        else {
            if (VO_D_well_water_table[n] < 2 or
                    (VO_D_NR_previous_well_water_table[n] == 1 &&
                     VO_D_well_water_table[n] == 2)) {
                Continuity_residual = vm - velocity_wt; // override the above
                // Continuity_residual = vm - vm_next_seg;
            }

            else {
                // Accumulation
                if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022
                    Continuity_residual += seg_length *
                            (seg_crossarea * (1. - sath) * rhof -
                             prev_time_seg_crossarea *
                             (1. - sath_prev_time) * rhof_prev_time) /
                            dt;
                }

                // Source
                Continuity_residual -= VO_D_mass_transfer_rate[n];

                if (n == number_of_well_segments - 1) {
                    if (vm >= 0.) {
                        Continuity_residual -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg * vm; // in
                    }

                    if (vm <= 0.) {
                        Continuity_residual -=
                                seg_crossarea * (1. - sath_prev_time) * rhof * vm; // out
                    }
                }

                else {
                    if (vm >= 0. && vm_next_seg >= 0.) {
                        if (n == 0)
                            Continuity_residual -= 0. * seg_crossarea *
                                    (1. - sath_prev_time) * 1000. *
                                    vm; // in !!

                        else
                            Continuity_residual -=
                                    prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                    rhof_prev_seg * vm; // in !!would give a zero in Jacobian
                        // diagonal for water table segment!!
                        Continuity_residual += seg_crossarea * (1. - sath_prev_time) *
                                rhof * vm_next_seg; // out
                    }

                    if (vm <= 0. && vm_next_seg <= 0.) {
                        Continuity_residual -= seg_crossarea * (1. - sath_prev_time) *
                                rhof *
                                vm; // out !!could give a zero in Jacobian
                        // diagonal  for water table segment!!
                        Continuity_residual += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                rhof_next_seg * vm_next_seg; // in
                    }

                    if (vm >= 0. && vm_next_seg <= 0.) {
                        if (n == 0)
                            Continuity_residual -= 0. * seg_crossarea *
                                    (1. - sath_prev_time) * 1000. *
                                    vm; // in !!

                        else
                            Continuity_residual -=
                                    prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                    rhof_prev_seg * vm; // in !!would give a zero in Jacobian
                        // diagonal for water table segment!!
                        Continuity_residual += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                rhof_next_seg * vm_next_seg; // in
                    }

                    if (vm <= 0. && vm_next_seg >= 0.) {
                        Continuity_residual -= seg_crossarea * (1. - sath_prev_time) *
                                rhof *
                                vm; // out !!could give a zero in Jacobian
                        // diagonal for water table segment!!

                        Continuity_residual += seg_crossarea * (1. - sath_prev_time) *
                                rhof * vm_next_seg; // out
                    }
                }
            }
        }
    }

    else // Multiphase Drift Flux Model, we probably want to skip DFM if not
        // multiphase, of course this needs Drift Flux Model implemented
        // elsewhere otherwise it will end up being the same as above
    {
        //
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Energy_Residuals(
        unsigned n, bool verbose) { // Need multi phase version with DFM
    if (VO_D_well_water_table[n] == 0) {
        // zero energy if empty segment
        // Energy_residual = hf;
        Energy_residual = hb;

        // What if halite is present?
    }

    else if (VO_D_well_water_table[n] == 1 or
             VO_D_previous_well_water_table[n] < 2) { // TEST MARCH 2022
        if (injection_mode or vm > 0.) {
            // Energy_residual = hf - top_s_enthalpy;
            Energy_residual = hb - top_s_enthalpy;
            // What if halite was present???
        }

        else {
            // Energy_residual = hf - hf_next_seg;
            Energy_residual = hb - hb_next_seg; // SWITCHED FEB 2022
            // What if halite is already present???

            // Approx
            if (with_grav_pot_energy) {
                Energy_residual +=
                        grav_acc * (next_seg_length * cos(next_seg_inclination) / 2. +
                                    seg_length * cos(seg_inclination) *
                                    seg_air_saturation_scaling / 2.);
            }
        }
    }

    else {
        // Accumulation
        if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022
            //            // fluid phase
            //            Energy_residual += seg_length * ( seg_crossarea * ( 1. -
            //            sath ) * ( rhof * hf - p ) -
            //                                              prev_time_seg_crossarea *
            //                                              ( 1. - sath_prev_time ) *
            //                                              ( rhof_prev_time *
            //                                              hf_prev_time - p_prev_time
            //                                              ) ) / dt;

            //            // halite phase:
            //            Energy_residual += seg_length * ( seg_crossarea * sath *
            //            rhoh * hh -
            //                                              prev_time_seg_crossarea *
            //                                              sath_prev_time *
            //                                              rhoh_prev_time *
            //                                              hh_prev_time ) / dt;

            // Above equivalent to:
            Energy_residual +=
                    seg_length *
                    (seg_crossarea * (rhob * hb - (1. - sath) * p) -
                     prev_time_seg_crossarea * (rhob_prev_time * hb_prev_time -
                                                (1. - sath_prev_time) * p_prev_time)) /
                    dt;

            // POTENTIAL PROBLEM: halite will be removed if air is pushed back into
            // the well

            if (with_kinetic_energy) {
                Energy_residual +=
                        seg_length *
                        (seg_crossarea * (1. - sath) * (rhof * 0.5 * vm * vm) -
                         prev_time_seg_crossarea * (1. - sath_prev_time) *
                         (rhof_prev_time * 0.5 * vm_prev_time * vm_prev_time)) /
                        dt;
            }
        }

        // source
        Energy_residual -= VO_D_energy_transfer_rate[n];

        if (n == number_of_well_segments - 1) {
            if (vm >= 0.) {
                Energy_residual -= prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                        rhof_prev_seg * hf_prev_seg * vm; // in
            }

            if (vm <= 0.) {
                Energy_residual -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf * vm; // out
            }
        }

        else {
            if (vm >= 0. && vm_next_seg >= 0.) {
                if (n == wt_index) {
                    Energy_residual -= 0. * seg_crossarea * (1. - sath_prev_time) * 1000 *
                            vm * top_s_enthalpy; // in !!
                    cerr << "SHOULD NOT REACH HERE!!!!";
                }

                else
                    Energy_residual -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                            hf_prev_seg * vm; // in
                Energy_residual += seg_crossarea * (1. - sath_prev_time) * rhof * hf *
                        vm_next_seg; // out
            }

            if (vm <= 0. && vm_next_seg <= 0.) {
                Energy_residual -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf * vm; // out
                Energy_residual += next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * hf_next_seg * vm_next_seg; // in
            }

            if (vm >= 0. && vm_next_seg <= 0.) {
                if (n == wt_index) {
                    Energy_residual -= 0. * seg_crossarea * (1. - sath_prev_time) * 1000 *
                            vm * top_s_enthalpy; // in !!
                    cerr << "SHOULD NOT REACH HERE!!!!";
                }

                else
                    Energy_residual -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                            hf_prev_seg * vm; // in
                Energy_residual += next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * hf_next_seg * vm_next_seg; // in
            }

            if (vm <= 0. && vm_next_seg >= 0.) {
                Energy_residual -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf * vm; // out
                Energy_residual += seg_crossarea * (1. - sath_prev_time) * rhof * hf *
                        vm_next_seg; // out
            }
        }

        if (with_grav_pot_energy) {
            if (vm >= 0.) {
                Energy_residual -= seg_crossarea * (1. - sath_prev_time) * seg_length *
                        rhof * vm * grav_acc * cos(seg_inclination);
            }
            if (vm <= 0.) {
                Energy_residual -= seg_crossarea * (1. - sath_prev_time) * seg_length *
                        rhof * vm * grav_acc * cos(seg_inclination);
            }
        }

        if (with_kinetic_energy) {
            if (n == number_of_well_segments - 1) {
                if (vm >= 0.) {
                    Energy_residual -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                            0.5 * vm_prev_seg * vm_prev_seg * vm; // in
                }

                if (vm <= 0.) {
                    Energy_residual -= seg_crossarea * (1. - sath_prev_time) * rhof *
                            0.5 * vm * vm * vm; // out
                }
            }

            else {
                if (vm >= 0. && vm_next_seg >= 0.) {
                    if (n == wt_index)
                        Energy_residual -= 0. * seg_crossarea * (1. - sath_prev_time) *
                                1000 * vm * 0.5 * vm * vm; // in !!

                    else
                        Energy_residual -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                                0.5 * vm_prev_seg * vm_prev_seg * vm; // in
                    Energy_residual += seg_crossarea * (1. - sath_prev_time) * rhof *
                            0.5 * vm * vm * vm_next_seg; // out
                }

                if (vm <= 0. && vm_next_seg <= 0.) {
                    Energy_residual -= seg_crossarea * (1. - sath_prev_time) * rhof *
                            0.5 * vm * vm * vm; // out
                    Energy_residual += next_seg_crossarea *
                            (1. - sath_prev_time_next_seg) * rhof_next_seg *
                            0.5 * vm_next_seg * vm_next_seg *
                            vm_next_seg; // in
                }

                if (vm >= 0. && vm_next_seg <= 0.) {
                    if (n == wt_index)
                        Energy_residual -= 0. * seg_crossarea * (1. - sath_prev_time) *
                                1000 * vm * 0.5 * vm * vm; // in !!

                    else
                        Energy_residual -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                                0.5 * vm_prev_seg * vm_prev_seg * vm; // in
                    Energy_residual += next_seg_crossarea *
                            (1. - sath_prev_time_next_seg) * rhof_next_seg *
                            0.5 * vm_next_seg * vm_next_seg *
                            vm_next_seg; // in
                }

                if (vm <= 0. && vm_next_seg >= 0.) {
                    Energy_residual -= seg_crossarea * (1. - sath_prev_time) * rhof *
                            0.5 * vm * vm * vm; // out
                    Energy_residual += seg_crossarea * (1. - sath_prev_time) * rhof *
                            0.5 * vm * vm * vm_next_seg; // out
                }
            }
        }

        // We finally add heat loss ( it does not depend on the flow direction)
        if (with_radial_heat &&
                VO_D_well_water_table[n] > 0) { // > 0 : not empty segment
            Energy_residual += VO_D_radial_heat[n];
            if (verbose)
                cerr << endl
                     << "After conduction, ENERGY RESIDUAL = " << Energy_residual;
        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Compute_Salt_Mass_Residuals(
        unsigned n, bool verbose) { // Need multi phase version with DFM

    if (VO_D_well_water_table[n] == 0) { // empty segment
        if (verbose)
            cerr << endl << "empty segment";
        // Salt_mass_residual = smf_f;
        Salt_mass_residual = smf_b; // Switched FEB 2022
        // What if halite is present???
    }

    else if (VO_D_well_water_table[n] == 1 or
             VO_D_previous_well_water_table[n] < 2) { // TEST MARCH 2022
        if (verbose)
            cerr << endl << "wt segment or previously empty";

        if (injection_mode or vm > 0.) {
            Salt_mass_residual = smf_b - top_smf; // Oct 2023
        }

        else {
            // Salt_mass_residual = smf_f - smf_f_next_seg;//SWITCHED FEB 2022
            Salt_mass_residual = smf_b - smf_b_next_seg; // SWITCHED FEB 2022
            // What if halite is already present???
        }
    }

    else {
        // Accumulation
        if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022

            //            //fluid phase:
            //            Salt_mass_residual +=  seg_length * ( seg_crossarea * ( 1. -
            //            sath ) * rhof * smf_f -
            //                                                  prev_time_seg_crossarea
            //                                                  * ( 1. -
            //                                                  sath_prev_time ) *
            //                                                  rhof_prev_time *
            //                                                  smf_f_prev_time ) /
            //                                                  dt;
            //            //halite phase:
            //            Salt_mass_residual +=  seg_length * ( seg_crossarea * sath *
            //            rhoh  -
            //                                                  prev_time_seg_crossarea
            //                                                  * sath_prev_time *
            //                                                  rhoh_prev_time ) / dt;
            // Above is equivalent to:
            Salt_mass_residual +=
                    seg_length *
                    (seg_crossarea * rhob * smf_b -
                     prev_time_seg_crossarea * rhob_prev_time * smf_b_prev_time) /
                    dt;

            if (verbose)
                cerr << endl << "Accumulation term: " << Salt_mass_residual;
        }
        // PROBLEM: halite will be removed if air is pushed back into the well

        // Source. (This was gated by a `with_salt` flag. n_of_eq is fixed at 4, so
        // the salt equation is always assembled and solved; the flag only removed
        // the EXCHANGE term, leaving salt transported but not exchanged — worse
        // than either setting. Historically with_salt=false dropped the equation
        // entirely, with n_of_eq=3; that path is long gone. Benoit DD/MM/YYYY)
        Salt_mass_residual -= VO_D_salt_mass_transfer_rate[n];

        if (n == number_of_well_segments - 1) {
            if (vm >= 0.) {
                Salt_mass_residual -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                        smf_f_prev_seg * vm; // in
            }

            if (vm <= 0.) {
                Salt_mass_residual -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f * vm; // out
            }
        }

        else {
            if (vm >= 0. && vm_next_seg >= 0.) {
                if (verbose)
                    cerr << endl << "flow down segment";

                if (n == 0)
                    Salt_mass_residual -= 0. * seg_crossarea * (1. - sath_prev_time) *
                            1000. * top_smf * vm; // in !! //top_smf=0..

                else
                    Salt_mass_residual -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                            smf_f_prev_seg *
                            vm; // in !!would give a zero in Jacobian
                // diagonal for water table segment!!
                Salt_mass_residual += seg_crossarea * (1. - sath_prev_time) * rhof *
                        smf_f * vm_next_seg; // out

                if (verbose)
                    cerr << endl
                         << "In term: "
                         << -prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg *smf_f_prev_seg *vm;
                if (verbose)
                    cerr << endl
                         << "Out term: "
                         << seg_crossarea * (1. - sath_prev_time) * rhof *smf_f *
                            vm_next_seg;

                if (verbose)
                    cerr << endl << "Final: " << Salt_mass_residual;
            }

            if (vm <= 0. && vm_next_seg <= 0.) {
                if (verbose)
                    cerr << endl << "flow up  segment";

                Salt_mass_residual -= seg_crossarea * (1. - sath_prev_time) * rhof *
                        smf_f * vm; // out !!could give a zero in Jacobian
                // diagonal  for water table segment!!
                Salt_mass_residual += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) * rhof_next_seg *
                        smf_f_next_seg * vm_next_seg; // in
            }

            if (vm >= 0. && vm_next_seg <= 0.) {
                if (verbose)
                    cerr << endl << "flow converge segment";
                if (n == 0)
                    Salt_mass_residual -= 0. * seg_crossarea * (1. - sath_prev_time) *
                            1000. * top_smf * vm; // in !! //top_smf=0..

                else
                    Salt_mass_residual -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * rhof_prev_seg *
                            smf_f_prev_seg *
                            vm; // in !!would give a zero in Jacobian
                // diagonal for water table segment!!
                Salt_mass_residual += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) * rhof_next_seg *
                        smf_f_next_seg * vm_next_seg; // in
            }

            if (vm <= 0. && vm_next_seg >= 0.) {
                if (verbose)
                    cerr << endl << "flow diverge segment";
                Salt_mass_residual -= seg_crossarea * (1. - sath_prev_time) * rhof *
                        smf_f * vm; // out !!could give a zero in Jacobian
                // diagonal for water table segment!!
                Salt_mass_residual += seg_crossarea * (1. - sath_prev_time) * rhof *
                        smf_f * vm_next_seg; // out
            }
        }
    }

    if (verbose)
        cerr << endl << "Salt mass resid: " << Salt_mass_residual;
}

//////////////////////////////////////////////////////////////////////////////
// JACOBIAN //////////////////////////////////////////////////////////////////

template <uint32_t dim>
void WellModelPrototype<dim>::Fill_momentum_Jacobian_part_salt(
        unsigned
        n) { // always expressed in term of fluid, ONLY FLUID DERIVATIVES NEEDED

    // Jacobian_segment_momentum =
    // [ dRp_n/dp_(n-1),    0
    //   dRp_n/dv_(n-1),    1
    //   dRp_n/dh_(n-1),    2
    //   dRp_n/dsmf_(n-1),  3
    //   dRp_n/dp_n,        4
    //   dRp_n/dv_n,        5
    //   dRp_n/dh_n,        6
    //   dRP_n/dsmf_n,      7
    //   dRP_n/dp_(n+1),    8 ]

    // TESTING DERIVATIVES
    drhofdp = drhofdhb = drhofdsmfb = dhdp =
            drhofdp_prev_seg = drhofdhb_prev_seg = drhofdsmfb_prev_seg = dhdp_prev_seg =
            drhofdp_next_seg = drhofdhb_next_seg = drhofdsmfb_next_seg = dhdp_next_seg = 0.;

    drhobdp = drhobdhb = drhobdsmfb =
            drhobdp_prev_seg = drhobdhb_prev_seg = drhobdsmfb_prev_seg =
            drhobdp_next_seg = drhobdhb_next_seg = drhobdsmfb_next_seg = 0.;


    drhofdp = drhofdp_;
    drhofdp_next_seg = drhofdp_next_seg_;
    drhofdp_prev_seg = drhofdp_prev_seg_;

    drhofdhb = drhofdhb_;
    drhofdhb_next_seg = drhofdhb_next_seg_;
    drhofdhb_prev_seg = drhofdhb_prev_seg_;

    drhofdsmfb = drhofdsmfb_;
    drhofdsmfb_next_seg = drhofdsmfb_next_seg_;
    drhofdsmfb_prev_seg = drhofdsmfb_prev_seg_;

    drhobdp = drhobdp_;
    drhobdp_next_seg = drhobdp_next_seg_;
    drhobdp_prev_seg = drhobdp_prev_seg_;

    drhobdhb = drhobdhb_;
    drhobdhb_next_seg = drhobdhb_next_seg_;
    drhobdhb_prev_seg = drhobdhb_prev_seg_;

    drhobdsmfb = drhobdsmfb_;
    drhobdsmfb_next_seg = drhobdsmfb_next_seg_;
    drhobdsmfb_prev_seg = drhobdsmfb_prev_seg_;

    dhdp = dhdp_;
    dhdp_next_seg = dhdp_next_seg_;
    dhdp_prev_seg = dhdp_prev_seg_;


    std::fill(Jacobian_segment_momentum.begin(), Jacobian_segment_momentum.end(),
              0.);

    if (VO_D_well_water_table[n] < 2) {
        // Momentum_residual = p - top_pressure;
        Jacobian_segment_momentum[4] = 1.;
    }

    else {
        // Momentum_residual += ( p - p_prev_seg ) - rhof_prev_seg * grav_acc * cos(
        // prev_seg_inclination ) * prev_seg_length *
        // prev_seg_air_saturation_scaling;

        Jacobian_segment_momentum[0] -= 1.;
        //        Jacobian_segment_momentum       [ 0 ] -= grav_acc * cos(
        //        prev_seg_inclination ) * prev_seg_length *
        //        prev_seg_air_saturation_scaling * drhofdp_prev_seg;
        //        Jacobian_segment_momentum       [ 2 ] -= grav_acc * cos(
        //        prev_seg_inclination ) * prev_seg_length *
        //        prev_seg_air_saturation_scaling * drhofdhb_prev_seg;
        //        Jacobian_segment_momentum       [ 3 ] -= grav_acc * cos(
        //        prev_seg_inclination ) * prev_seg_length *
        //        prev_seg_air_saturation_scaling * drhofdsmfb_prev_seg;
        Jacobian_segment_momentum[4] += 1.;

        if (abs(vm_prev_seg) > 1.e-6 && n > 0) {
            if (with_friction) {
                // Momentum_residual += velocity_sign_prev_seg *
                // friction_factor_prev_seg * rhof_prev_seg * vm_prev_seg * vm_prev_seg
                // * prev_seg_air_saturation_scaling * prev_seg_length

                Jacobian_segment_momentum[0] +=
                        velocity_sign_prev_seg * friction_factor_prev_seg * vm_prev_seg *
                        vm_prev_seg * drhofdp_prev_seg * prev_seg_air_saturation_scaling *
                        prev_seg_length;
                Jacobian_segment_momentum[1] +=
                        velocity_sign_prev_seg * friction_factor_prev_seg * rhof_prev_seg *
                        2 * abs(vm_prev_seg) * prev_seg_air_saturation_scaling *
                        prev_seg_length;
                // abs(vm_prev_seg) coming from friction term

                Jacobian_segment_momentum[2] +=
                        velocity_sign_prev_seg * friction_factor_prev_seg * vm_prev_seg *
                        vm_prev_seg * drhofdhb_prev_seg * prev_seg_air_saturation_scaling *
                        prev_seg_length;
                Jacobian_segment_momentum[3] +=
                        velocity_sign_prev_seg * friction_factor_prev_seg * vm_prev_seg *
                        vm_prev_seg * drhofdsmfb_prev_seg *
                        prev_seg_air_saturation_scaling * prev_seg_length;
            }

            if (with_inertial_terms_in_momentum && rhof > 0.) {
                // Momentum_residual += ( rhof_prev_seg * vm_prev_seg -
                // rhof_prev_time_prev_seg * vm_prev_time_prev_seg ) / dt *
                // prev_seg_length* dt; Momentum_residual += ( rhof * vm * vm -
                // rhof_prev_seg * vm_prev_seg * vm_prev_seg )* dt;

                if (rhof_prev_time_prev_seg > 0.) {
                    Jacobian_segment_momentum[0] +=
                            drhofdp_prev_seg * vm_prev_seg / dt * prev_seg_length;
                    Jacobian_segment_momentum[1] += rhof_prev_seg / dt * prev_seg_length;
                    Jacobian_segment_momentum[2] +=
                            drhofdhb_prev_seg * vm_prev_seg / dt * prev_seg_length;
                    Jacobian_segment_momentum[3] +=
                            drhofdsmfb_prev_seg * vm_prev_seg / dt * prev_seg_length;
                }

                if (rhof_prev_seg > 0.) {
                    Jacobian_segment_momentum[0] -=
                            drhofdp_prev_seg * vm_prev_seg * vm_prev_seg;
                    Jacobian_segment_momentum[1] -= rhof_prev_seg * 2 * abs(vm_prev_seg);
                    Jacobian_segment_momentum[2] -=
                            drhofdhb_prev_seg * vm_prev_seg * vm_prev_seg;
                    Jacobian_segment_momentum[3] -=
                            drhofdsmfb_prev_seg * vm_prev_seg * vm_prev_seg;
                    Jacobian_segment_momentum[4] += drhofdp * vm * vm;
                    Jacobian_segment_momentum[5] += rhof * 2 * abs(vm);
                    Jacobian_segment_momentum[6] += drhofdhb * vm * vm;
                    Jacobian_segment_momentum[7] += drhofdsmfb * vm * vm;
                }
            }
        }

        //        //        //DAMPING
        //        Jacobian_segment_momentum   [ 0 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 1 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 2 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 3 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 4 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 5 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 6 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 7 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_momentum   [ 8 ] *= DAMP_JACOBIAN;
    }

    // Assign
    unsigned m = n_of_eq * n;
    unsigned k = n_of_eq * (n - 1); // n - 1: we start from prev n segment

    if (n == 0) {
        // no prev segment
        for (unsigned l = n_of_eq; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k + n_of_eq, Jacobian_segment_momentum[l]);
            k++;
        }
    }

    else {
        for (unsigned l = 0; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k, Jacobian_segment_momentum[l]);
            k++;
        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Fill_continuity_Jacobian_part_salt(unsigned n) {
    // MANY MODIFS NECESSARY IF WE USE THE FORMULATION WITH SATURATIONS
    // NEW JACOBIAN TERMS NEED TO BE CALCULATED, instead of drho/dx we need
    // d(sat*rho)/dx

    // Jacobian_segment_continuity =
    // [ dRc_n/dp_(n-1),        0
    //   dRc_n/dv_(n-1),        1
    //   dRc_n/dh_(n-1),        2
    //   dRc_n/dsmf_(n-1),      3
    //   dRc_n/dp_n,        0->4
    //   dRc_n/dv_n,        1->5
    //   dRc_n/dh_n,        2->6
    //   dRc_n/dsmf_n,      3->7
    //   dRc_n/dp_(n+1),    4->8
    //   dRc_n/dv_(n+1),    5->9
    //   dRc_n/dh_(n+1),    6->10
    //   dRc_n/dsmf_(n+1)   7->11 ]

    // TESTING DERIVATIVES
    drhofdp = drhofdhb = drhofdsmfb = dhdp =
            drhofdp_prev_seg = drhofdhb_prev_seg = drhofdsmfb_prev_seg = dhdp_prev_seg =
            drhofdp_next_seg = drhofdhb_next_seg = drhofdsmfb_next_seg = dhdp_next_seg = 0.;

    drhobdp = drhobdhb = drhobdsmfb =
            drhobdp_prev_seg = drhobdhb_prev_seg = drhobdsmfb_prev_seg =
            drhobdp_next_seg = drhobdhb_next_seg = drhobdsmfb_next_seg = 0.;

    drhofdp = drhofdp_;
    drhofdp_next_seg = drhofdp_next_seg_;
    drhofdp_prev_seg = drhofdp_prev_seg_;

    drhofdhb = drhofdhb_;
    drhofdhb_next_seg = drhofdhb_next_seg_;
    drhofdhb_prev_seg = drhofdhb_prev_seg_;

    drhofdsmfb = drhofdsmfb_;
    drhofdsmfb_next_seg = drhofdsmfb_next_seg_;
    drhofdsmfb_prev_seg = drhofdsmfb_prev_seg_;

    drhobdp = drhobdp_;
    drhobdp_next_seg = drhobdp_next_seg_;
    drhobdp_prev_seg = drhobdp_prev_seg_;

    drhobdhb = drhobdhb_;
    drhobdhb_next_seg = drhobdhb_next_seg_;
    drhobdhb_prev_seg = drhobdhb_prev_seg_;

    drhobdsmfb = drhobdsmfb_;
    drhobdsmfb_next_seg = drhobdsmfb_next_seg_;
    drhobdsmfb_prev_seg = drhobdsmfb_prev_seg_;

    dhdp = dhdp_;
    dhdp_next_seg = dhdp_next_seg_;
    dhdp_prev_seg = dhdp_prev_seg_;


    std::fill(Jacobian_segment_continuity.begin(),
              Jacobian_segment_continuity.end(), 0.);

    // (Was `if ((multiphase && homogeneous_flow) || !multiphase)`. That is true
    // unless multiphase && !homogeneous_flow, and homogeneous_flow is always
    // true, so the test never excluded anything. `multiphase` removed.
    // Benoit DD/MM/YYYY)
    if (homogeneous_flow) {
        if (VO_well_segment_activated[n] < 1) {
            // Continuity_residual = vm;
            Jacobian_segment_continuity[5] = 1.;
        }

        else {
            if (VO_D_well_water_table[n] < 2 or
                    (VO_D_NR_previous_well_water_table[n] == 1 &&
                     VO_D_well_water_table[n] == 2)) {
                // Continuity_residual = vm - WT_VELOCITY;//override the above
                // Continuity_residual = vm - vm_next_seg;
                Jacobian_segment_continuity[5] = 1.;
                // Jacobian_segment_continuity     [ 9 ] = -1.;
            }

            else {
                // Accumulation
                if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022
                    //                    Continuity_residual += seg_length * (
                    //                    seg_crossarea * ( 1. - sath ) * rhof -
                    //                                                          prev_time_seg_crossarea
                    //                                                          * ( 1. -
                    //                                                          sath_prev_time
                    //                                                          ) *
                    //                                                          rhof_prev_time
                    //                                                          ) / dt;

                    Jacobian_segment_continuity[4] +=
                            seg_length * seg_crossarea * (1. - sath) * drhofdp / dt;
                    Jacobian_segment_continuity[6] +=
                            seg_length * seg_crossarea * (1. - sath) * drhofdhb / dt;
                    Jacobian_segment_continuity[7] +=
                            seg_length * seg_crossarea * (1. - sath) * drhofdsmfb / dt;
                }

                if (n == number_of_well_segments - 1) {
                    if (vm >= 0.) {
                        // Continuity_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg
                        // ) * rhof_prev_seg * vm;//in
                        Jacobian_segment_continuity[0] -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) *
                                vm * drhofdp_prev_seg;
                        Jacobian_segment_continuity[5] -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg;
                        Jacobian_segment_continuity[2] -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) *
                                vm * drhofdhb_prev_seg;
                        Jacobian_segment_continuity[3] -= prev_seg_crossarea *
                                (1. - sath_prev_time_prev_seg) *
                                vm * drhofdsmfb_prev_seg;
                    }

                    if (vm <= 0.) {
                        // Continuity_residual -= seg_crossarea * ( 1. - sath ) * rhof *
                        // vm;//out
                        Jacobian_segment_continuity[4] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdp;
                        Jacobian_segment_continuity[5] -=
                                seg_crossarea * (1. - sath_prev_time) * rhof;
                        Jacobian_segment_continuity[6] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdhb;
                        Jacobian_segment_continuity[7] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdsmfb;
                    }
                }

                else {
                    if (vm >= 0. && vm_next_seg >= 0.) {
                        if (n == 0) {
                            // Continuity_residual -= 0. * seg_crossarea * ( 1. - sath ) *
                            // 1000. * vm;//in !!
                            Jacobian_segment_continuity[5] -=
                                    0. * seg_crossarea * (1. - sath_prev_time) * 1000.;
                        }
                        else {
                            // Continuity_residual -= prev_seg_crossarea * ( 1. -
                            // sath_prev_seg ) * rhof_prev_seg * vm;//in
                            Jacobian_segment_continuity[0] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdp_prev_seg;
                            Jacobian_segment_continuity[5] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    rhof_prev_seg;
                            Jacobian_segment_continuity[2] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdhb_prev_seg;
                            Jacobian_segment_continuity[3] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdsmfb_prev_seg;
                        }

                        // Continuity_residual += seg_crossarea * ( 1. - sath ) * rhof *
                        // vm_next_seg;//out
                        Jacobian_segment_continuity[4] +=
                                seg_crossarea * (1. - sath_prev_time) * vm_next_seg * drhofdp;
                        Jacobian_segment_continuity[9] +=
                                seg_crossarea * (1. - sath_prev_time) * rhof;
                        Jacobian_segment_continuity[6] +=
                                seg_crossarea * (1. - sath_prev_time) * vm_next_seg * drhofdhb;
                        Jacobian_segment_continuity[7] += seg_crossarea *
                                (1. - sath_prev_time) *
                                vm_next_seg * drhofdsmfb;
                    }

                    else if (vm <= 0. && vm_next_seg <= 0.) {
                        // Continuity_residual -= seg_crossarea * ( 1. - sath ) * rhof *
                        // vm;//out
                        Jacobian_segment_continuity[4] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdp;
                        Jacobian_segment_continuity[5] -=
                                seg_crossarea * (1. - sath_prev_time) * rhof;
                        Jacobian_segment_continuity[6] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdhb;
                        Jacobian_segment_continuity[7] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdsmfb;

                        // Continuity_residual += next_seg_crossarea * ( 1. - sath_next_seg
                        // ) * rhof_next_seg * vm_next_seg;//in
                        Jacobian_segment_continuity[8] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdp_next_seg;
                        Jacobian_segment_continuity[9] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                rhof_next_seg;
                        Jacobian_segment_continuity[10] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdhb_next_seg;
                        Jacobian_segment_continuity[11] +=
                                next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdsmfb_next_seg;
                    }

                    else if (vm >= 0. && vm_next_seg <= 0.) {
                        if (n == 0) {
                            // Continuity_residual -= 0. * seg_crossarea * ( 1. - sath ) *
                            // 1000. * vm;//in !!
                            Jacobian_segment_continuity[5] -=
                                    0. * seg_crossarea * (1. - sath_prev_time) * 1000.;
                        }
                        else {
                            // Continuity_residual -= prev_seg_crossarea * ( 1. -
                            // sath_prev_seg ) * rhof_prev_seg * vm;//in
                            Jacobian_segment_continuity[0] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdp_prev_seg;
                            Jacobian_segment_continuity[5] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    rhof_prev_seg;
                            Jacobian_segment_continuity[2] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdhb_prev_seg;
                            Jacobian_segment_continuity[3] -= prev_seg_crossarea *
                                    (1. - sath_prev_time_prev_seg) *
                                    vm * drhofdsmfb_prev_seg;
                        }

                        // Continuity_residual += next_seg_crossarea * ( 1. - sath_next_seg
                        // ) * rhof_next_seg * vm_next_seg;//in
                        Jacobian_segment_continuity[8] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdp_next_seg;
                        Jacobian_segment_continuity[9] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                rhof_next_seg;
                        Jacobian_segment_continuity[10] += next_seg_crossarea *
                                (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdhb_next_seg;
                        Jacobian_segment_continuity[11] +=
                                next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                                vm_next_seg * drhofdsmfb_next_seg;
                    }

                    else if (vm <= 0. && vm_next_seg >= 0.) {
                        // Continuity_residual -= seg_crossarea * ( 1. - sath ) * rhof *
                        // vm;//out
                        Jacobian_segment_continuity[4] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdp;
                        Jacobian_segment_continuity[5] -=
                                seg_crossarea * (1. - sath_prev_time) * rhof;
                        Jacobian_segment_continuity[6] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdhb;
                        Jacobian_segment_continuity[7] -=
                                seg_crossarea * (1. - sath_prev_time) * vm * drhofdsmfb;

                        // Continuity_residual += seg_crossarea * ( 1. - sath ) * rhof *
                        // vm_next_seg;//out
                        Jacobian_segment_continuity[4] +=
                                seg_crossarea * (1. - sath_prev_time) * vm_next_seg * drhofdp;
                        Jacobian_segment_continuity[9] +=
                                seg_crossarea * (1. - sath_prev_time) * rhof;
                        Jacobian_segment_continuity[6] +=
                                seg_crossarea * (1. - sath_prev_time) * vm_next_seg * drhofdhb;
                        Jacobian_segment_continuity[7] += seg_crossarea *
                                (1. - sath_prev_time) *
                                vm_next_seg * drhofdsmfb;
                    }
                }

                // DAMPING
                //             Jacobian_segment_continuity   [ 0 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 1 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 2 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 3 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 4 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 5 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 6 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 7 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 8 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 9 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 10 ] *= DAMP_JACOBIAN;
                //             Jacobian_segment_continuity   [ 11 ] *= DAMP_JACOBIAN;
            }
        }

        if (Jacobian_segment_continuity[5] == 0.) {
            cerr << endl << " Jacobian Diag term in continuity= 0 " << endl;
            PAUSE();
        }

        // Assign
        unsigned m = n_of_eq * n + 1;
        unsigned k = n_of_eq * (n - 1); // n - 1: we start from prev n segment

        if (n == 0) {
            // no prev segment
            for (unsigned l = n_of_eq; l < n_of_eq * n_of_eq; l++) {
                Jacobian.Assign(m, k + n_of_eq, Jacobian_segment_continuity[l]);
                k++;
            }
        }

        else if (n == number_of_well_segments - 1) {
            // no next segment
            for (unsigned l = 0; l < n_of_eq * (n_of_eq - 1); l++) {
                Jacobian.Assign(m, k, Jacobian_segment_continuity[l]);
                k++;
            }
        }

        else {
            for (unsigned l = 0; l < n_of_eq * n_of_eq; l++) {
                Jacobian.Assign(m, k, Jacobian_segment_continuity[l]);
                k++;
            }
        }
    }

    else // this needs Drift Flux Model implemented elsewhere otherwise it will
        // end up being the same as above
    {
        // need to write it...!!!!!!!!!!!!!!!!!!!!!!!!
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Fill_energy_Jacobian_part_salt(unsigned n) {
    // NEW JACOBIAN TERMS NEED TO BE CALCULATED, instead of drho/dx we need
    // d(sat*rho)/dx

    // Jacobian_segment_energy =
    // [ dRe_n/dp_(n-1),    0
    //   dRe_n/dv_(n-1),    1
    //   dRe_n/dh_(n-1),    2
    //   dRe_n/dsmf(n-1),   3
    //   dRe_n/dp_n,        4
    //   dRe_n/dv_n,        5
    //   dRe_n/dh_n,        6
    //   dRe_n/dsmf_n,      7
    //   dRe_n/dp_(n+1),    8
    //   dRe_n/dv_(n+1),    9
    //   dRe_n/dh_(n+1),    10
    //   dRe_n/dsmf_(n+1)   11 ]

    // TESTING DERIVATIVES
    drhofdp = drhofdhb = drhofdsmfb = dhdp =
            drhofdp_prev_seg = drhofdhb_prev_seg = drhofdsmfb_prev_seg = dhdp_prev_seg =
            drhofdp_next_seg = drhofdhb_next_seg = drhofdsmfb_next_seg = dhdp_next_seg = 0.;

    drhobdp = drhobdhb = drhobdsmfb =
            drhobdp_prev_seg = drhobdhb_prev_seg = drhobdsmfb_prev_seg =
            drhobdp_next_seg = drhobdhb_next_seg = drhobdsmfb_next_seg = 0.;

    drhofdp = drhofdp_;
    drhofdp_next_seg = drhofdp_next_seg_;
    drhofdp_prev_seg = drhofdp_prev_seg_;

    drhofdhb = drhofdhb_;
    drhofdhb_next_seg = drhofdhb_next_seg_;
    drhofdhb_prev_seg = drhofdhb_prev_seg_;

    drhofdsmfb = drhofdsmfb_;
    drhofdsmfb_next_seg = drhofdsmfb_next_seg_;
    drhofdsmfb_prev_seg = drhofdsmfb_prev_seg_;

    drhobdp = drhobdp_;
    drhobdp_next_seg = drhobdp_next_seg_;
    drhobdp_prev_seg = drhobdp_prev_seg_;

    drhobdhb = drhobdhb_;
    drhobdhb_next_seg = drhobdhb_next_seg_;
    drhobdhb_prev_seg = drhobdhb_prev_seg_;

    drhobdsmfb = drhobdsmfb_;
    drhobdsmfb_next_seg = drhobdsmfb_next_seg_;
    drhobdsmfb_prev_seg = drhobdsmfb_prev_seg_;

    dhdp = dhdp_;
    dhdp_next_seg = dhdp_next_seg_;
    dhdp_prev_seg = dhdp_prev_seg_;


    std::fill(Jacobian_segment_energy.begin(), Jacobian_segment_energy.end(), 0.);

    if (VO_D_well_water_table[n] == 0) {
        // zero energy if empty segment
        // Energy_residual = hf;
        // Energy_residual = hb;
        Jacobian_segment_energy[6] += 1.; // diag
    }

    else if (VO_D_well_water_table[n] == 1 or
             VO_D_previous_well_water_table[n] < 2) { // TEST MARCH 2022
        if (injection_mode or vm > 0.) {
            // Energy_residual = hf - top_s_enthalpy;
            // Energy_residual = hb - top_s_enthalpy;
            Jacobian_segment_energy[6] += 1.; // diag
        }
        else {
            // Energy_residual = hf - hf_next_seg; + pot energy
            // Energy_residual = hb - hb_next_seg; + pot energy
            Jacobian_segment_energy[6] += 1.;  // diag
            Jacobian_segment_energy[10] -= 1.; // diag
        }
    }

    else {
        // Accumulation
        if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022

            //        // fluid phase
            //        Energy_residual += seg_length * ( seg_crossarea * ( 1. - sath )
            //        * ( rhof * hf - p ) -
            //                                          prev_time_seg_crossarea * ( 1.
            //                                          - sath_prev_time ) * (
            //                                          rhof_prev_time * hf_prev_time
            //                                          - p_prev_time ) ) / dt;

            //        // halite phase:
            //        Energy_residual += seg_length * ( seg_crossarea * sath * rhoh *
            //        hh -
            //                                          prev_time_seg_crossarea *
            //                                          sath_prev_time *
            //                                          rhoh_prev_time * hh_prev_time
            //                                          ) / dt;

            //            Jacobian_segment_energy         [ 4 ] += seg_length * ( 1. -
            //            sath ) * seg_crossarea * ( hf * drhofdp + rhof * dhdp - 1 )
            //            / dt; Jacobian_segment_energy         [ 6 ] += seg_length *
            //            ( 1. - sath ) * seg_crossarea * ( hf * drhofdhb + rhof )/
            //            dt;//diag Jacobian_segment_energy         [ 7 ] +=
            //            seg_length * ( 1. - sath ) * seg_crossarea * ( hf *
            //            drhofdsmfb + rhof * dhdsmf )/ dt;
            // THIS NEGLECTS HALITE RELATED DERIVATIVES

            // Above equivalent to:
            // Energy_residual += seg_length * ( seg_crossarea * ( rhob * hb - ( 1. -
            // sath ) * p ) - prev_time_seg_crossarea * ( rhob_prev_time *
            // hb_prev_time - ( 1. - sath_prev_time ) * p_prev_time ) ) / dt;

            Jacobian_segment_energy[4] += seg_length * seg_crossarea *
                    (hb * drhobdp + rhob * dhdp - (1. - sath)) /
                    dt;
            Jacobian_segment_energy[6] +=
                    seg_length * seg_crossarea * (hb * drhobdhb + rhob) / dt; // diag
            Jacobian_segment_energy[7] +=
                    seg_length * seg_crossarea * (hb * drhobdsmfb + rhob * dhdsmf) / dt;
            // BULK DENSITY DERIVATIVES NEEDED ABOVE!

            if (with_kinetic_energy) {
                //                Energy_residual += seg_length * ( seg_crossarea * ( 1.
                //                - sath ) * ( rhof * 0.5 * vm * vm ) -
                //                                                  prev_time_seg_crossarea
                //                                                  * ( 1. -
                //                                                  sath_prev_time ) * (
                //                                                  rhof_prev_time * 0.5
                //                                                  * vm_prev_time *
                //                                                  vm_prev_time ) ) /
                //                                                  dt;

                Jacobian_segment_energy[4] += seg_length * seg_crossarea * (1. - sath) *
                        0.5 * vm * vm * drhofdp / dt;
                Jacobian_segment_energy[5] += seg_length * seg_crossarea * (1. - sath) *
                        rhof * abs(vm) /
                        dt; // abs(vm) coming from kinetic energy
                Jacobian_segment_energy[6] += seg_length * seg_crossarea * (1. - sath) *
                        0.5 * vm * vm * drhofdhb / dt; // diag
                Jacobian_segment_energy[7] += seg_length * seg_crossarea * (1. - sath) *
                        0.5 * vm * vm * drhofdsmfb / dt;
            }
        }

        if (n == number_of_well_segments - 1) {
            if (vm >= 0.) {
                // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                // rhof_prev_seg * hf_prev_seg * vm;//in
                Jacobian_segment_energy[0] -=
                        prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                        (hf_prev_seg * drhofdp_prev_seg + rhof_prev_seg * dhdp_prev_seg);
                Jacobian_segment_energy[5] -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) *
                        rhof_prev_seg * hf_prev_seg;
                Jacobian_segment_energy[2] -=
                        prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                        (hf_prev_seg * drhofdhb_prev_seg + rhof_prev_seg);
                Jacobian_segment_energy[3] -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) * vm *
                        (hf_prev_seg * drhofdsmfb_prev_seg +
                         rhof_prev_seg * dhdsmf_prev_seg);
            }

            if (vm <= 0.) {
                // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * hf *
                // vm;//out
                Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdp + rhof * dhdp);
                Jacobian_segment_energy[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf;
                Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdhb + rhof); // diag
                Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdsmfb + rhof * dhdsmf);
            }
        }

        else {
            if (vm >= 0. && vm_next_seg >= 0.) {
                if (n == wt_index) {
                    // Energy_residual -= 0. * seg_crossarea * ( 1. - sath ) * 1000 * vm *
                    // top_s_enthalpy;//in !!
                    Jacobian_segment_energy[5] -= 0. * seg_crossarea *
                            (1. - sath_prev_time) * 1000 *
                            top_s_enthalpy;
                }

                else {
                    // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                    // rhof_prev_seg * hf_prev_seg * vm;//in
                    Jacobian_segment_energy[0] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdp_prev_seg + rhof_prev_seg * dhdp_prev_seg);
                    Jacobian_segment_energy[5] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * hf_prev_seg;
                    Jacobian_segment_energy[2] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdhb_prev_seg + rhof_prev_seg);
                    Jacobian_segment_energy[3] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdsmfb_prev_seg +
                             rhof_prev_seg * dhdsmf_prev_seg);

                    // Energy_residual += seg_crossarea * ( 1. - sath ) * rhof * hf *
                    // vm_next_seg;//out
                    Jacobian_segment_energy[4] += seg_crossarea * (1. - sath_prev_time) *
                            vm_next_seg *
                            (hf * drhofdp + rhof * dhdp);
                    Jacobian_segment_energy[9] +=
                            seg_crossarea * (1. - sath_prev_time) * rhof * hf;
                    Jacobian_segment_energy[6] += seg_crossarea * (1. - sath_prev_time) *
                            vm_next_seg *
                            (hf * drhofdhb + rhof); // diag
                    Jacobian_segment_energy[7] += seg_crossarea * (1. - sath_prev_time) *
                            vm_next_seg *
                            (hf * drhofdsmfb + rhof * dhdsmf);
                }
            }

            if (vm <= 0. && vm_next_seg <= 0.) {
                // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * hf *
                // vm;//out
                Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdp + rhof * dhdp);
                Jacobian_segment_energy[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf;
                Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdhb + rhof); // diag
                Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdsmfb + rhof * dhdsmf);

                // Energy_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                // rhof_next_seg * hf_next_seg * vm_next_seg;//in
                Jacobian_segment_energy[8] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdp_next_seg + rhof_next_seg * dhdp_next_seg);
                Jacobian_segment_energy[9] += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * hf_next_seg;
                Jacobian_segment_energy[10] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdhb_next_seg + rhof_next_seg);
                Jacobian_segment_energy[11] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdsmfb_next_seg +
                         rhof_next_seg * dhdsmf_next_seg);
            }

            if (vm >= 0. && vm_next_seg <= 0.) {
                if (n == wt_index) {
                    // Energy_residual -= 0. * seg_crossarea * ( 1. - sath ) * 1000 * vm *
                    // top_s_enthalpy;//in !!
                    Jacobian_segment_energy[5] -= 0. * seg_crossarea *
                            (1. - sath_prev_time) * 1000 *
                            top_s_enthalpy;
                }

                else {
                    // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                    // rhof_prev_seg * hf_prev_seg * vm;//in
                    Jacobian_segment_energy[0] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdp_prev_seg + rhof_prev_seg * dhdp_prev_seg);
                    Jacobian_segment_energy[5] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * hf_prev_seg;
                    Jacobian_segment_energy[2] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdhb_prev_seg + rhof_prev_seg);
                    Jacobian_segment_energy[3] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (hf_prev_seg * drhofdsmfb_prev_seg +
                             rhof_prev_seg * dhdsmf_prev_seg);
                }

                // Energy_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                // rhof_next_seg * hf_next_seg * vm_next_seg;//in
                Jacobian_segment_energy[8] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdp_next_seg + rhof_next_seg * dhdp_next_seg);
                Jacobian_segment_energy[9] += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * hf_next_seg;
                Jacobian_segment_energy[10] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdhb_next_seg + rhof_next_seg);
                Jacobian_segment_energy[11] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (hf_next_seg * drhofdsmfb_next_seg +
                         rhof_next_seg * dhdsmf_next_seg);
            }

            if (vm <= 0. && vm_next_seg >= 0.) {
                // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * hf *
                // vm;//out
                Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdp + rhof * dhdp);
                Jacobian_segment_energy[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf;
                Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdhb + rhof); // diag
                Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (hf * drhofdsmfb + rhof * dhdsmf);

                // Energy_residual += seg_crossarea * ( 1. - sath ) * rhof * hf *
                // vm_next_seg;//out
                Jacobian_segment_energy[4] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (hf * drhofdp + rhof * dhdp);
                Jacobian_segment_energy[9] +=
                        seg_crossarea * (1. - sath_prev_time) * rhof * hf;
                Jacobian_segment_energy[6] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (hf * drhofdhb + rhof); // diag
                Jacobian_segment_energy[7] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (hf * drhofdsmfb + rhof * dhdsmf);
            }
        }

        if (with_grav_pot_energy) {
            if (vm >= 0.) {
                // Energy_residual -= seg_crossarea * ( 1. - sath ) * seg_length * rhof
                // * vm * grav_acc * cos( seg_inclination );
                Jacobian_segment_energy[4] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdp;
                Jacobian_segment_energy[5] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * rhof * grav_acc *
                        cos(seg_inclination);
                Jacobian_segment_energy[6] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdhb; // diag
                Jacobian_segment_energy[7] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdsmfb;
            }
            if (vm <= 0.) {
                // Energy_residual -= seg_crossarea * ( 1. - sath ) * seg_length * rhof
                // * vm * grav_acc * cos( seg_inclination );
                Jacobian_segment_energy[4] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdp;
                Jacobian_segment_energy[5] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * rhof * grav_acc *
                        cos(seg_inclination);
                Jacobian_segment_energy[6] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdhb; // diag
                Jacobian_segment_energy[7] -= seg_length * seg_crossarea *
                        (1. - sath_prev_time) * vm * grav_acc *
                        cos(seg_inclination) * drhofdsmfb;
            }
        }

        if (with_kinetic_energy) {
            if (n == number_of_well_segments - 1) {
                if (vm >= 0.) {
                    // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                    // rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg * vm;//in
                    Jacobian_segment_energy[0] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                            vm_prev_seg * vm_prev_seg * vm * drhofdp_prev_seg;
                    Jacobian_segment_energy[1] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * abs(vm_prev_seg) *
                            vm; // abs(vm_prev_seg) from kinetic energy
                    Jacobian_segment_energy[2] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                            vm_prev_seg * vm_prev_seg * vm * drhofdhb_prev_seg;
                    Jacobian_segment_energy[3] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                            vm_prev_seg * vm_prev_seg * vm * drhofdsmfb_prev_seg;
                    Jacobian_segment_energy[5] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg;
                }

                if (vm <= 0.) {
                    // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * 0.5 * vm
                    // * vm * vm;//out
                    Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdp;
                    Jacobian_segment_energy[5] -=
                            seg_crossarea * (1. - sath_prev_time) * 1.5 * rhof * vm * vm;
                    Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdhb; // diag
                    Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdsmfb;
                }
            }

            else {
                if (vm >= 0. && vm_next_seg >= 0.) {
                    if (n == wt_index) {
                        // Energy_residual -= 0. * seg_crossarea * ( 1. - sath ) * 1000. *
                        // vm * 0.5 * vm * vm;//in !!
                        Jacobian_segment_energy[5] -= 0. * seg_crossarea *
                                (1. - sath_prev_time) *
                                (1.5 * 1000. * vm * vm);
                    }

                    else {
                        // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                        // rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg * vm;//in
                        Jacobian_segment_energy[0] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdp_prev_seg;
                        Jacobian_segment_energy[1] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg * abs(vm_prev_seg) *
                                vm; // abs(vm_prev_seg) from kinetic energy
                        Jacobian_segment_energy[2] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdhb_prev_seg;
                        Jacobian_segment_energy[3] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdsmfb_prev_seg;
                        Jacobian_segment_energy[5] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg;
                    }

                    // Energy_residual += seg_crossarea * ( 1. - sath ) * rhof * 0.5 * vm
                    // * vm * vm_next_seg;//out
                    Jacobian_segment_energy[4] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg * drhofdp;
                    Jacobian_segment_energy[5] +=
                            seg_crossarea * (1. - sath_prev_time) * rhof * abs(vm) *
                            vm_next_seg; // abs(vm) from kinetic energy
                    Jacobian_segment_energy[6] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg *
                            drhofdhb; // diag
                    Jacobian_segment_energy[7] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg *
                            drhofdsmfb;
                    Jacobian_segment_energy[9] +=
                            seg_crossarea * (1. - sath_prev_time) * rhof * 0.5 * vm * vm;
                }

                if (vm <= 0. && vm_next_seg <= 0.) {
                    // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * 0.5 * vm
                    // * vm * vm;//out
                    Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdp;
                    Jacobian_segment_energy[5] -=
                            seg_crossarea * (1. - sath_prev_time) * 1.5 * rhof * vm * vm;
                    Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdhb; // diag
                    Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdsmfb;

                    // Energy_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                    // rhof_next_seg * 0.5 * vm_next_seg * vm_next_seg * vm_next_seg;//in
                    Jacobian_segment_energy[8] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdp_next_seg;
                    Jacobian_segment_energy[9] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 1.5 *
                            rhof_next_seg * vm_next_seg * vm_next_seg;
                    Jacobian_segment_energy[10] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdhb_next_seg;
                    Jacobian_segment_energy[11] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdsmfb_next_seg;
                }

                if (vm >= 0. && vm_next_seg <= 0.) {
                    if (n == wt_index) {
                        // Energy_residual -= 0. * seg_crossarea * ( 1. - sath ) * 1000. *
                        // vm * 0.5 * vm * vm;//in !!
                        Jacobian_segment_energy[5] -= 0. * seg_crossarea *
                                (1. - sath_prev_time) *
                                (1.5 * 1000. * vm * vm);
                    }

                    else {
                        // Energy_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                        // rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg * vm;//in
                        Jacobian_segment_energy[0] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdp_prev_seg;
                        Jacobian_segment_energy[1] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg * abs(vm_prev_seg) *
                                vm; // abs(vm_prev_seg) from kinetic energy
                        Jacobian_segment_energy[2] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdhb_prev_seg;
                        Jacobian_segment_energy[3] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) * 0.5 *
                                vm_prev_seg * vm_prev_seg * vm * drhofdsmfb_prev_seg;
                        Jacobian_segment_energy[5] -=
                                prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                                rhof_prev_seg * 0.5 * vm_prev_seg * vm_prev_seg;
                    }

                    // Energy_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                    // rhof_next_seg * 0.5 * vm_next_seg * vm_next_seg * vm_next_seg;//in
                    Jacobian_segment_energy[8] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdp_next_seg;
                    Jacobian_segment_energy[9] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 1.5 *
                            rhof_next_seg * vm_next_seg * vm_next_seg;
                    Jacobian_segment_energy[10] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdhb_next_seg;
                    Jacobian_segment_energy[11] +=
                            next_seg_crossarea * (1. - sath_prev_time_next_seg) * 0.5 *
                            vm_next_seg * vm_next_seg * vm_next_seg * drhofdsmfb_next_seg;
                }

                if (vm <= 0. && vm_next_seg >= 0.) {
                    // Energy_residual -= seg_crossarea * ( 1. - sath ) * rhof * 0.5 * vm
                    // * vm * vm;//out
                    Jacobian_segment_energy[4] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdp;
                    Jacobian_segment_energy[5] -=
                            seg_crossarea * (1. - sath_prev_time) * 1.5 * rhof * vm * vm;
                    Jacobian_segment_energy[6] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdhb; // diag
                    Jacobian_segment_energy[7] -= seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm * drhofdsmfb;

                    // Energy_residual += seg_crossarea * ( 1. - sath ) * rhof * 0.5 * vm
                    // * vm * vm_next_seg;//out
                    Jacobian_segment_energy[4] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg * drhofdp;
                    Jacobian_segment_energy[5] +=
                            seg_crossarea * (1. - sath_prev_time) * rhof * abs(vm) *
                            vm_next_seg; // abs(vm) from kinetic energy
                    Jacobian_segment_energy[6] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg *
                            drhofdhb; // diag
                    Jacobian_segment_energy[7] += seg_crossarea * (1. - sath_prev_time) *
                            0.5 * vm * vm * vm_next_seg *
                            drhofdsmfb;
                    Jacobian_segment_energy[9] +=
                            seg_crossarea * (1. - sath_prev_time) * rhof * 0.5 * vm * vm;
                }
            }
        }

        // DAMPING

        Jacobian_segment_energy[0] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[1] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[2] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[3] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[4] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[5] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[6] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[7] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[8] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[9] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[10] *= DAMP_JACOBIAN;
        Jacobian_segment_energy[11] *= DAMP_JACOBIAN;
    }

    if (Jacobian_segment_energy[6] == 0.) {
        cerr << endl
             << "segment " << n << ": Jacobian Diag term in energy= 0 " << endl;
        PAUSE();
    }

    // Assign
    unsigned m = n_of_eq * n + 2;
    unsigned k = n_of_eq * (n - 1); // n - 1: we start from prev n segment

    if (n == 0) {
        // no prev segment
        for (unsigned l = n_of_eq; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k + n_of_eq, Jacobian_segment_energy[l]);
            k++;
        }
    }

    else if (n == number_of_well_segments - 1) {
        // no next segment
        for (unsigned l = 0; l < n_of_eq * (n_of_eq - 1); l++) {
            Jacobian.Assign(m, k, Jacobian_segment_energy[l]);
            k++;
        }
    }

    else {
        for (unsigned l = 0; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k, Jacobian_segment_energy[l]);
            k++;
        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::Fill_salt_mass_Jacobian_part(unsigned n) {
    // MANY MODIFS NECESSARY IF WE USE THE FORMULATION WITH SATURATIONS
    // NEW JACOBIAN TERMS NEED TO BE CALCULATED, instead of drho/dx we need
    // d(sat*rho)/dx also we need to make the distinction between rho fluid and
    // rho bulk..

    // Jacobian_segment_salt_mass =
    // [ dRS_n/dp_(n-1),    0
    //   dRS_n/dv_(n-1),    1
    //   dRS_n/dh_(n-1),    2
    //   dRS_n/dsmf(n-1),   3
    //   dRS_n/dp_n,        4
    //   dRS_n/dv_n,        5
    //   dRS_n/dh_n,        6
    //   dRS_n/dsmf_n,      7
    //   dRS_n/dp_(n+1),    8
    //   dRS_n/dv_(n+1),    9
    //   dRS_n/dh_(n+1),    10
    //   dRS_n/dsmf_(n+1)   11 ]

    // TESTING DERIVATIVES
    drhofdp = drhofdhb = drhofdsmfb = dhdp = drhofdp_prev_seg =
            drhofdhb_prev_seg = drhofdsmfb_prev_seg = dhdp_prev_seg =
            drhofdp_next_seg = drhofdhb_next_seg = drhofdsmfb_next_seg =
            dhdp_next_seg = 0.;

    drhobdp = drhobdhb = drhobdsmfb = drhobdp_prev_seg = drhobdhb_prev_seg =
            drhobdsmfb_prev_seg = drhobdp_next_seg = drhobdhb_next_seg =
            drhobdsmfb_next_seg = 0.;

    drhofdp = drhofdp_;
    drhofdp_next_seg = drhofdp_next_seg_;
    drhofdp_prev_seg = drhofdp_prev_seg_;

    drhofdhb = drhofdhb_;
    drhofdhb_next_seg = drhofdhb_next_seg_;
    drhofdhb_prev_seg = drhofdhb_prev_seg_;

    drhofdsmfb = drhofdsmfb_;
    drhofdsmfb_next_seg = drhofdsmfb_next_seg_;
    drhofdsmfb_prev_seg = drhofdsmfb_prev_seg_;

    drhobdp = drhobdp_;
    drhobdp_next_seg = drhobdp_next_seg_;
    drhobdp_prev_seg = drhobdp_prev_seg_;

    drhobdhb = drhobdhb_;
    drhobdhb_next_seg = drhobdhb_next_seg_;
    drhobdhb_prev_seg = drhobdhb_prev_seg_;

    drhobdsmfb = drhobdsmfb_;
    drhobdsmfb_next_seg = drhobdsmfb_next_seg_;
    drhobdsmfb_prev_seg = drhobdsmfb_prev_seg_;

    dhdp = dhdp_;
    dhdp_next_seg = dhdp_next_seg_;
    dhdp_prev_seg = dhdp_prev_seg_;


    std::fill(Jacobian_segment_salt_mass.begin(),
              Jacobian_segment_salt_mass.end(), 0.);

    if (VO_D_well_water_table[n] == 0) { // empty segment
        // Salt_mass_residual = smf_f;
        // Salt_mass_residual = smf_b;
        Jacobian_segment_salt_mass[7] = 1.;
    }

    else if (VO_D_well_water_table[n] == 1 or
             VO_D_previous_well_water_table[n] < 2) { // TEST MARCH 2022
        if (injection_mode or vm > 0.) {
            // Salt_mass_residual = smf_b - top_smf;//Oct 2023
            Jacobian_segment_salt_mass[7] += 1.; // diag
        }

        else {
            // Salt_mass_residual = smf_f - smf_f_next_seg;
            // Salt_mass_residual = smf_b - smf_b_next_seg;
            Jacobian_segment_salt_mass[7] = 1.;  // diag
            Jacobian_segment_salt_mass[11] = -1.; // diag
        }
    }

    else {
        // Accumulation
        if (VO_D_previous_well_water_table[n] > 1) { // TEST MARCH 2022

            // Salt_mass_residual += seg_length * ( seg_crossarea * rhob * smf_b -
            // prev_time_seg_crossarea * rhob_prev_time * smf_b_prev_time ) / dt;

            Jacobian_segment_salt_mass[4] +=
                    seg_length * seg_crossarea * (drhobdp * smf_b + rhob * dsmfdp) / dt;
            Jacobian_segment_salt_mass[6] +=
                    seg_length * seg_crossarea * (drhobdhb * smf_b + rhob * dsmfdh) / dt;
            Jacobian_segment_salt_mass[7] +=
                    seg_length * seg_crossarea * (drhobdsmfb * smf_b + rhob) / dt;
        }
        // BULK DENSITY DERIVATIVES NEEDED ABOVE!

        if (n == number_of_well_segments - 1) {
            // assuming the wt segment will never reach here

            if (vm >= 0.) {
                // Salt_mass_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                // rhof_prev_seg * smf_f_prev_seg * vm;//in
                Jacobian_segment_salt_mass[0] -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) * vm *
                        (smf_f_prev_seg * drhofdp_prev_seg +
                         rhof_prev_seg * dsmfdp_prev_seg);
                Jacobian_segment_salt_mass[2] -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) * vm *
                        (smf_f_prev_seg * drhofdhb_prev_seg +
                         rhof_prev_seg * dsmfdh_prev_seg);
                Jacobian_segment_salt_mass[3] -=
                        prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                        (drhofdsmfb_prev_seg * smf_f_prev_seg + rhof_prev_seg) * vm;
                Jacobian_segment_salt_mass[5] -= prev_seg_crossarea *
                        (1. - sath_prev_time_prev_seg) *
                        rhof_prev_seg * smf_f_prev_seg;
            }

            if (vm <= 0.) {
                // Salt_mass_residual -= seg_crossarea * ( 1. - sath ) * rhof * smf_f *
                // vm;//out
                Jacobian_segment_salt_mass[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (drhofdp * smf_f + rhof * dsmfdp);
                Jacobian_segment_salt_mass[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f;
                Jacobian_segment_salt_mass[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm *
                        (drhofdhb * smf_f + rhof * dsmfdh);
                Jacobian_segment_salt_mass[7] -= seg_crossarea * (1. - sath_prev_time) *
                        (drhofdsmfb * smf_f + rhof) * vm;
            }
        }

        else {
            if (vm >= 0. && vm_next_seg >= 0.) {
                if (n == 0) {
                    // Salt_mass_residual -= 0. * seg_crossarea * ( 1. - sath ) * 1000. *
                    // top_smf * vm;//in !! //top_smf=0..
                    Jacobian_segment_salt_mass[5] -=
                            0. * seg_crossarea * (1. - sath_prev_time) * 1000. * top_smf;
                }

                else {
                    // Salt_mass_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                    // rhof_prev_seg * smf_f_prev_seg * vm;//in
                    Jacobian_segment_salt_mass[0] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (smf_f_prev_seg * drhofdp_prev_seg +
                             rhof_prev_seg * dsmfdp_prev_seg);
                    Jacobian_segment_salt_mass[2] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (smf_f_prev_seg * drhofdhb_prev_seg +
                             rhof_prev_seg * dsmfdh_prev_seg);
                    Jacobian_segment_salt_mass[3] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                            (drhofdsmfb_prev_seg * smf_f_prev_seg + rhof_prev_seg) * vm;
                    Jacobian_segment_salt_mass[5] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * smf_f_prev_seg;
                }

                // Salt_mass_residual += seg_crossarea * ( 1. - sath ) * rhof * smf_f *
                // vm_next_seg;//out
                Jacobian_segment_salt_mass[4] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (drhofdp * smf_f + rhof * dsmfdp);
                Jacobian_segment_salt_mass[6] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (drhofdhb * smf_f + rhof * dsmfdh);
                Jacobian_segment_salt_mass[7] += seg_crossarea * (1. - sath_prev_time) *
                        (drhofdsmfb * smf_f + rhof) *
                        vm_next_seg;
                Jacobian_segment_salt_mass[9] +=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f;
            }

            if (vm <= 0. && vm_next_seg <= 0.) {
                // Salt_mass_residual -= seg_crossarea * ( 1. - sath ) * rhof * smf_f *
                // vm;//out
                Jacobian_segment_salt_mass[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (drhofdp * smf_f + rhof * dsmfdp);
                Jacobian_segment_salt_mass[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f;
                Jacobian_segment_salt_mass[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm *
                        (drhofdhb * smf_f + rhof * dsmfdh);
                Jacobian_segment_salt_mass[7] -= seg_crossarea * (1. - sath_prev_time) *
                        (drhofdsmfb * smf_f + rhof) * vm;

                // Salt_mass_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                // rhof_next_seg * smf_f_next_seg * vm_next_seg;//in
                Jacobian_segment_salt_mass[8] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (smf_f_next_seg * drhofdp_next_seg +
                         rhof_next_seg * dsmfdp_next_seg);
                Jacobian_segment_salt_mass[9] += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * smf_f_next_seg;
                Jacobian_segment_salt_mass[10] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (smf_f_next_seg * drhofdhb_next_seg +
                         rhof_next_seg * dsmfdh_next_seg);
                Jacobian_segment_salt_mass[11] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                        (drhofdsmfb_next_seg * smf_f_next_seg + rhof_next_seg) *
                        vm_next_seg;
            }

            if (vm >= 0. && vm_next_seg <= 0.) {
                if (n == 0) {
                    // Salt_mass_residual -= seg_crossarea * ( 1. - sath ) * 1000. *
                    // top_smf * vm;//in !! //top_smf=0..
                    Jacobian_segment_salt_mass[5] -=
                            0. * seg_crossarea * (1. - sath_prev_time) * 1000. * top_smf;
                }

                else {
                    // Salt_mass_residual -= prev_seg_crossarea * ( 1. - sath_prev_seg ) *
                    // rhof_prev_seg * smf_f_prev_seg * vm;//in
                    Jacobian_segment_salt_mass[0] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (smf_f_prev_seg * drhofdp_prev_seg +
                             rhof_prev_seg * dsmfdp_prev_seg);
                    Jacobian_segment_salt_mass[5] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) *
                            rhof_prev_seg * smf_f_prev_seg;
                    Jacobian_segment_salt_mass[2] -= prev_seg_crossarea *
                            (1. - sath_prev_time_prev_seg) * vm *
                            (smf_f_prev_seg * drhofdhb_prev_seg +
                             rhof_prev_seg * dsmfdh_prev_seg);
                    Jacobian_segment_salt_mass[3] -=
                            prev_seg_crossarea * (1. - sath_prev_time_prev_seg) *
                            (drhofdsmfb_prev_seg * smf_f_prev_seg + rhof_prev_seg) * vm;
                }

                // Salt_mass_residual += next_seg_crossarea * ( 1. - sath_next_seg ) *
                // rhof_next_seg * smf_f_next_seg * vm_next_seg;//in
                Jacobian_segment_salt_mass[8] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (smf_f_next_seg * drhofdp_next_seg +
                         rhof_next_seg * dsmfdp_next_seg);
                Jacobian_segment_salt_mass[9] += next_seg_crossarea *
                        (1. - sath_prev_time_next_seg) *
                        rhof_next_seg * smf_f_next_seg;
                Jacobian_segment_salt_mass[10] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) * vm_next_seg *
                        (smf_f_next_seg * drhofdhb_next_seg +
                         rhof_next_seg * dsmfdh_next_seg);
                Jacobian_segment_salt_mass[11] +=
                        next_seg_crossarea * (1. - sath_prev_time_next_seg) *
                        (drhofdsmfb_next_seg * smf_f_next_seg + rhof_next_seg) *
                        vm_next_seg;
            }

            if (vm <= 0. && vm_next_seg >= 0.) {
                // Salt_mass_residual -= seg_crossarea * ( 1. - sath ) * rhof * smf_f *
                // vm;//out
                Jacobian_segment_salt_mass[4] -= seg_crossarea * (1. - sath_prev_time) *
                        vm * (drhofdp * smf_f + rhof * dsmfdp);
                Jacobian_segment_salt_mass[5] -=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f;
                Jacobian_segment_salt_mass[6] -= seg_crossarea * (1. - sath_prev_time) *
                        vm *
                        (drhofdhb * smf_f + rhof * dsmfdh);
                Jacobian_segment_salt_mass[7] -= seg_crossarea * (1. - sath_prev_time) *
                        (drhofdsmfb * smf_f + rhof) * vm;

                // Salt_mass_residual += seg_crossarea * ( 1. - sath ) * rhof * smf_f *
                // vm_next_seg;//out
                Jacobian_segment_salt_mass[4] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (drhofdp * smf_f + rhof * dsmfdp);
                Jacobian_segment_salt_mass[9] +=
                        seg_crossarea * (1. - sath_prev_time) * rhof * smf_f;
                Jacobian_segment_salt_mass[6] += seg_crossarea * (1. - sath_prev_time) *
                        vm_next_seg *
                        (drhofdhb * smf_f + rhof * dsmfdh);
                Jacobian_segment_salt_mass[7] += seg_crossarea * (1. - sath_prev_time) *
                        (drhofdsmfb * smf_f + rhof) *
                        vm_next_seg;
            }
        }

        //        //        //DAMPING
        //        Jacobian_segment_salt_mass   [ 0 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 1 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 2 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 3 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 4 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 5 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 6 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 7 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 8 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 9 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 10 ] *= DAMP_JACOBIAN;
        //        Jacobian_segment_salt_mass   [ 11 ] *= DAMP_JACOBIAN;
    }

    // Assign
    unsigned m = n_of_eq * n + 3;
    unsigned k = n_of_eq * (n - 1); // n - 1: we start from prev n segment

    if (n == 0) {
        // no prev segment
        for (unsigned l = n_of_eq; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k + n_of_eq, Jacobian_segment_salt_mass[l]);
            k++;
        }
    }

    else if (n == number_of_well_segments - 1) {
        // no next segment
        for (unsigned l = 0; l < n_of_eq * (n_of_eq - 1); l++) {
            Jacobian.Assign(m, k, Jacobian_segment_salt_mass[l]);
            k++;
        }
    }

    else {
        for (unsigned l = 0; l < n_of_eq * n_of_eq; l++) {
            Jacobian.Assign(m, k, Jacobian_segment_salt_mass[l]);
            k++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template <uint32_t dim> void WellModelPrototype<dim>::ComputeFrictionFactor() {

    double Re_number(0.);
    // double radius ( 0. );

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        if (VO_D_well_water_table[n] == 0) {
            VO_D_well_Re[n] = 0.;
            VO_D_well_friction_factor[n] = 0.;
        }

        else {
            if (abs(VO_D_well_velocity_fluid[n]) > 1.e-6) {
                radius = sqrtf(VO_D_well_cross_area[n] *
                               (1 - VO_D_previous_well_sat_halite[n]) / CSMP_PI);
                // !! EFFECTIVE CROSS AREA !! that account for the halite (assumed to
                // form an homogeneous crust on the surface of the well). Although it
                // does not account for a potential change in roughness!

                Re_number = VO_D_well_density_fluid[n] *
                        abs(VO_D_well_velocity_fluid[n]) * 2 * radius /
                        VO_D_well_viscosity_fluid[n];

                VO_D_well_Re[n] = Re_number;

                // cerr<<endl<<"segment: "<<n<<" , radius: "<<radius<<", density:
                // "<<VO_D_well_density_fluid[ n ]<<", velocity:
                // "<<abs(VO_D_well_velocity_fluid[ n ])<<", viscosity:
                // "<<VO_D_well_viscosity_fluid [ n ]<<", Re number: "<<Re_number;

                if (Re_number > 0.) {
                    if (Re_number < 4000.) {
                        VO_D_well_friction_factor[n] = 64. / Re_number;
                    }

                    else {
                        double X = log10(pipe_roughness / (3.7 * 2 * radius) -
                                         5.02 / Re_number *
                                         log10(pipe_roughness / (3.7 * 2 * radius) +
                                               13. / Re_number));
                        VO_D_well_friction_factor[n] = 1. / (4. * X * X);
                    }
                }
            }
        }
    }
}

template <uint32_t dim> void WellModelPrototype<dim>::Element_average() {
    // Averaging segments on elements

    std::fill(VO_interp_well_pressure.begin(), VO_interp_well_pressure.end(), 0);
    std::fill(VO_well_s_enthalpy_liquid.begin(), VO_well_s_enthalpy_liquid.end(),
              0);
    std::fill(VO_well_s_enthalpy_vapor.begin(), VO_well_s_enthalpy_vapor.end(),
              0);
    std::fill(VO_well_s_enthalpy_fluid.begin(), VO_well_s_enthalpy_fluid.end(),
              0);
    std::fill(VO_well_smf_liquid.begin(), VO_well_smf_liquid.end(), 0);
    std::fill(VO_well_smf_vapor.begin(), VO_well_smf_vapor.end(), 0);
    std::fill(VO_well_fluid_state.begin(), VO_well_fluid_state.end(), 0);
    std::fill(VO_well_smf_fluid.begin(), VO_well_smf_fluid.end(), 0);
    std::fill(VO_well_density_liquid.begin(), VO_well_density_liquid.end(), 0);
    std::fill(VO_well_density_vapor.begin(), VO_well_density_vapor.end(), 0);
    std::fill(VO_well_density_fluid.begin(), VO_well_density_fluid.end(), 0);
    std::fill(VO_well_sat_liquid.begin(), VO_well_sat_liquid.end(), 0);
    std::fill(VO_well_sat_vapor.begin(), VO_well_sat_vapor.end(), 0);
    std::fill(VO_well_sat_halite.begin(), VO_well_sat_halite.end(), 0);
    std::fill(VO_well_temperature.begin(), VO_well_temperature.end(), 0);
    std::fill(VO_well_velocity_liquid.begin(), VO_well_velocity_liquid.end(), 0);
    std::fill(VO_well_velocity_vapor.begin(), VO_well_velocity_vapor.end(), 0);
    std::fill(VO_well_velocity_fluid.begin(), VO_well_velocity_fluid.end(), 0);
    std::fill(VO_well_cross_area.begin(), VO_well_cross_area.end(), 0);
    std::fill(VO_well_friction_factor.begin(), VO_well_friction_factor.end(), 0);
    std::fill(VO_well_viscosity_fluid.begin(), VO_well_viscosity_fluid.end(), 0);
    std::fill(VO_well_Re.begin(), VO_well_Re.end(), 0);

    unsigned cumul_index = number_of_virtual_top_segments;
    for (unsigned k = 0; k < number_of_well_elements; k++) {
        for (unsigned l = cumul_index;
             l < cumul_index + number_of_segment_divisions[k]; l++) {
            VO_interp_well_pressure[k] +=
                    VO_D_interp_well_pressure[l] /
                    number_of_segment_divisions[k]; // INTERPOLATED at segment centers
            // then averaged at element center,
            // not great
            VO_well_s_enthalpy_liquid[k] +=
                    VO_D_well_s_enthalpy_liquid[l] / number_of_segment_divisions[k];
            VO_well_s_enthalpy_vapor[k] +=
                    VO_D_well_s_enthalpy_vapor[l] / number_of_segment_divisions[k];
            VO_well_s_enthalpy_fluid[k] +=
                    VO_D_well_s_enthalpy_fluid[l] / number_of_segment_divisions[k];
            VO_well_smf_liquid[k] +=
                    VO_D_well_smf_liquid[l] / number_of_segment_divisions[k];
            VO_well_smf_vapor[k] +=
                    VO_D_well_smf_vapor[l] / number_of_segment_divisions[k];
            VO_well_fluid_state[k] +=
                    VO_D_well_fluid_state[l] / number_of_segment_divisions[k];
            VO_well_smf_fluid[k] +=
                    VO_D_well_smf_fluid[l] / number_of_segment_divisions[k];
            VO_well_density_liquid[k] +=
                    VO_D_well_density_liquid[l] / number_of_segment_divisions[k];
            VO_well_density_vapor[k] +=
                    VO_D_well_density_vapor[l] / number_of_segment_divisions[k];
            VO_well_density_fluid[k] +=
                    VO_D_well_density_fluid[l] / number_of_segment_divisions[k];
            VO_well_sat_liquid[k] +=
                    VO_D_well_sat_liquid[l] / number_of_segment_divisions[k];
            VO_well_sat_vapor[k] +=
                    VO_D_well_sat_vapor[l] / number_of_segment_divisions[k];
            VO_well_sat_halite[k] +=
                    VO_D_well_sat_halite[l] / number_of_segment_divisions[k];
            VO_well_temperature[k] +=
                    VO_D_well_temperature[l] / number_of_segment_divisions[k];
            VO_well_velocity_liquid[k] +=
                    VO_D_well_velocity_liquid[l] / number_of_segment_divisions[k];
            VO_well_velocity_vapor[k] +=
                    VO_D_well_velocity_vapor[l] / number_of_segment_divisions[k];
            VO_well_velocity_fluid[k] +=
                    VO_D_well_velocity_fluid[l] / number_of_segment_divisions[k];
            VO_well_cross_area[k] +=
                    VO_D_well_cross_area[l] / number_of_segment_divisions[k];
            VO_well_friction_factor[k] +=
                    VO_D_well_friction_factor[l] / number_of_segment_divisions[k];
            VO_well_viscosity_fluid[k] +=
                    VO_D_well_viscosity_fluid[l] / number_of_segment_divisions[k];
            VO_well_Re[k] += VO_D_well_Re[l] / number_of_segment_divisions[k];
        }

        cumul_index += number_of_segment_divisions[k];
    }
}

template <uint32_t dim> void WellModelPrototype<dim>::Write_results() {
    Element_average();

    unsigned n = 0;
    for (typename std::vector<Element<dim> *>::const_iterator eit =
         model_ref_.Region(well_name).CellsBegin();
         eit != model_ref_.Region(well_name).CellsEnd(); ++eit) {
        well_pressure() = VO_interp_well_pressure[element_deordering_vector[n]];
        well_s_enthalpy_liquid() =
                VO_well_s_enthalpy_liquid[element_deordering_vector[n]];
        well_s_enthalpy_vapor() =
                VO_well_s_enthalpy_vapor[element_deordering_vector[n]];
        well_s_enthalpy_fluid() =
                VO_well_s_enthalpy_fluid[element_deordering_vector[n]];
        well_smf_liquid() = VO_well_smf_liquid[element_deordering_vector[n]];
        well_smf_vapor() = VO_well_smf_vapor[element_deordering_vector[n]];
        well_fluid_state() = VO_well_fluid_state[element_deordering_vector[n]];
        well_smf_fluid() = VO_well_smf_fluid[element_deordering_vector[n]];
        well_density_liquid() =
                VO_well_density_liquid[element_deordering_vector[n]];
        well_density_vapor() = VO_well_density_vapor[element_deordering_vector[n]];
        well_density_fluid() = VO_well_density_fluid[element_deordering_vector[n]];
        well_sat_liquid() = VO_well_sat_liquid[element_deordering_vector[n]];
        well_sat_vapor() = VO_well_sat_vapor[element_deordering_vector[n]];
        well_sat_halite() = VO_well_sat_halite[element_deordering_vector[n]];
        well_temperature() = VO_well_temperature[element_deordering_vector[n]];
        well_velocity_liquid() =
                VO_well_velocity_liquid[element_deordering_vector[n]];
        well_velocity_vapor() =
                VO_well_velocity_vapor[element_deordering_vector[n]];
        well_velocity_fluid() =
                VO_well_velocity_fluid[element_deordering_vector[n]];
        well_cross_area() = VO_well_cross_area[element_deordering_vector[n]];
        well_friction_factor() =
                VO_well_friction_factor[element_deordering_vector[n]];
        well_viscosity_fluid() =
                VO_well_viscosity_fluid[element_deordering_vector[n]];
        well_Re() = VO_well_Re[element_deordering_vector[n]];

        (*eit)->Store(well_pressure_key_, well_pressure);
        (*eit)->Store(well_s_enthalpy_liquid_key_, well_s_enthalpy_liquid);
        (*eit)->Store(well_s_enthalpy_vapor_key_, well_s_enthalpy_vapor);
        (*eit)->Store(well_s_enthalpy_fluid_key_, well_s_enthalpy_fluid);
        (*eit)->Store(well_smf_liquid_key_, well_smf_liquid);
        (*eit)->Store(well_smf_vapor_key_, well_smf_vapor);
        (*eit)->Store(well_fluid_state_key_, well_fluid_state);
        (*eit)->Store(well_smf_fluid_key_, well_smf_fluid);
        (*eit)->Store(well_density_liquid_key_, well_density_liquid);
        (*eit)->Store(well_density_vapor_key_, well_density_vapor);
        (*eit)->Store(well_density_fluid_key_, well_density_fluid);
        (*eit)->Store(well_sat_liquid_key_, well_sat_liquid);
        (*eit)->Store(well_sat_vapor_key_, well_sat_vapor);
        (*eit)->Store(well_sat_halite_key_, well_sat_halite);
        (*eit)->Store(well_temperature_key_, well_temperature);
        (*eit)->Store(well_velocity_liquid_key_, well_velocity_liquid);
        (*eit)->Store(well_velocity_vapor_key_, well_velocity_vapor);
        (*eit)->Store(well_velocity_fluid_key_, well_velocity_fluid);
        (*eit)->Store(well_cross_area_key_, well_cross_area);
        (*eit)->Store(well_friction_factor_key_, well_friction_factor);
        (*eit)->Store(well_viscosity_fluid_key_, well_viscosity_fluid);
        (*eit)->Store(well_Re_key_, well_Re);

        wellhead_pressure() = top_pressure;
        (*eit)->Store(wellhead_pressure_key_, wellhead_pressure);

        n++;
    }
}

template <uint32_t dim> void WellModelPrototype<dim>::Display() {
    cerr << defaultfloat;
    cerr.precision(6);

    cerr << endl << endl << "Timestep: " << dt << " seconds *****";

    //    cerr << endl << endl << "ordered Heel depth                     : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_element_heel_depth[ n ] << ", ";

    cerr << endl << endl << "Well seg mass rate  [Kg/s]      : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_mass_transfer_rate[n] << ", ";

    cerr << endl << endl << "Res nodes pressure [Pa]            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_pressure[n] << ", ";

    cerr << endl << endl << "Well-Reservoir nodes skin [-]            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_skin[n] << ", ";

    cerr << endl << endl << "Res nodes smf liquid            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_smf_liquid[n] << ", ";

    cerr << endl << endl << "Res nodes smf vapor            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_smf_vapor[n] << ", ";

    cerr << endl << endl << "Res nodes sat vapor            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_sat_vapor[n] << ", ";

    //    cerr << endl << endl << "well seg length                        : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_segment_length[ n ] << ", ";

    cerr << endl << endl << "Well seg s_enthalpy fluid  [J/kg]      : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_s_enthalpy_fluid[n] << ", ";

    cerr << endl << endl << "Well seg s_enthalpy bulk  [J/kg]       : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_s_enthalpy_bulk[n] << ", ";

    //    cerr << endl << endl << "Well seg s_enthalpy halite  [J/kg]     : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_s_enthalpy_halite[ n ] << ", ";

    //    cerr << endl << endl << "Well seg density liquid [kg/m3]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_density_liquid[ n ] << ", ";

    //    cerr << endl << endl << "Well seg density vapor [kg/m3]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_density_vapor[ n ] << ", ";

    //    cerr << endl << endl << "Well seg density bulk [kg/m3]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_density_bulk[ n ] << ", ";

    //    cerr << endl << endl << "Well seg saturation liquid [/]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_sat_liquid[ n ] << ", ";

    //    cerr << endl << endl << "Salt transfer rate [/s]               : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_salt_mass_transfer_rate[ n ] << ", ";

    //    cerr << endl << endl << "Well seg Re [/]                        : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_Re[ n ] << ", ";

    //    cerr << endl << endl << "Well seg viscosity fluid [Pa/s]       : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_viscosity_fluid[ n ] << ", ";

    //    cerr << endl << endl << "Well seg friction factor [?]           : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_friction_factor[ n ] << ", ";

    cerr << endl << endl << "\033[32m Well seg velocity fluid [m/s]        : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_velocity_fluid[n] << ", ";

    cerr << endl
         << endl
         << "\033[0m Well seg salt mass fraction (in fluid) [kg/kg]   : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_smf_fluid[n] << ", ";

    //    cerr << endl << endl << "Well mass salt [kg]                    : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_salt_mass [ n ] << ", ";

    //    cerr << endl << endl << "drhofdp                                 : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_drhofdp[ n ] << ", ";
    //    cerr << endl << endl << "drhofdhb                                 : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_drhofdhb[ n ] << ", ";

    //    cerr << endl << endl << "DHDP               : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_dhdp[ n ] << ", ";
    //    cerr << endl << endl << "drhofdsmfb               : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_drhofdsmfb[ n ] << ", ";

    //    cerr << endl << endl << "Well seg salt mass change term [kg/s]        :
    //    "; for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_segment_length[ n ] * VO_D_well_cross_area[ n ] / dt *
    //                (VO_D_well_density_bulk[ n ] * VO_D_well_smf_bulk[ n ] -
    //                 VO_D_previous_well_density_fluid[ n ] *
    //                 VO_D_previous_well_smf_bulk[ n ])
    //             << ", ";

    //    cerr << endl << endl << "Well seg salt mass Inflow [kg/s]             :
    //    "; for (unsigned n = 0; n < number_of_well_segments - 1; n++)
    //        cerr << VO_D_well_cross_area[n + 1] * (VO_D_well_density_fluid[n +
    //        1] *
    //                VO_D_well_velocity_fluid[n + 1] * VO_D_well_smf_fluid[n +
    //                1])
    //                << ", "; cerr<<"Source rate";

    //    cerr << endl << endl << "Well seg salt mass Outflow [kg/s]            :
    //    "; for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_cross_area[ n ] *
    //                (VO_D_well_density_fluid[ n ] * VO_D_well_velocity_fluid[ n
    //                ] * VO_D_well_smf_fluid[ n ])
    //             << ", ";

    //    cerr << endl << endl << "Nodal radial heat                      : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_radial_heat[ n ] << ", ";

    //    cerr << endl << endl << "well index                             : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_well_index[ n ] << ", ";

    //    cerr << endl << endl << "well index2                            : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_well_index2[ n ] << ", ";

    //    cerr << endl << endl << "relperm_visc_liquid                    : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_relperm_visc_liquid[ n ] << ", ";

    //    cerr << endl << endl << "relperm_visc_vapor                     : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_relperm_visc_vapor[ n ] << ", ";

    //    cerr << endl << endl << "formation_permeability                 : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_formation_permeability[ n ] << ", ";

    //    cerr << endl << endl << "bulk_volume                            : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_bulk_volume[ n ] << ", ";

    //    cerr << endl << endl << "Mass transfer rate [kg/s]              : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_mass_transfer_rate[ n ] << ", ";

    //    cerr << endl << endl << "Mass transfer rate [Kg/s]  (segments) : "; for
    //    (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_mass_transfer_rate[ n ]<< ", ";

    //    cerr << endl << endl << "\e[1m\033[31m Timestep: " << dt << " seconds
    //    *****";

    //    cerr << endl << endl << "Vol flow rate (Rm3/day)                : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_velocity_fluid[ n ]*VO_D_well_cross_area[ n
    //        ]*86400 << ", ";
    //    cerr << "\033[0m";

    //    cerr << endl << endl << "Nodal res T                            : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_temperature[ n ] << ", ";

    cerr << endl << endl << "Segments res T                         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_temperature[n] << ", ";

    cerr << endl << endl << "Well seg temperature [C]               : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_temperature[n] << ", ";

    cerr << endl << endl << "Segments radial heat                   : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_radial_heat[n] << ", ";

       cerr << endl << endl << "Segments Cross area                   : ";
       for (unsigned n = 0; n < number_of_well_segments; n++)
           cerr << VO_D_well_cross_area[ n ] << ", ";

    //    cerr << endl << endl << "Well seg heel depth [m]                : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_segment_heel_depth[ n ] << ", ";

    cerr << endl
         << endl
         << "previous water table depth: " << previous_water_table_depth
         << ", water table depth: " << water_table_depth;

    cerr << endl << endl << "Well seg water table index      : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_water_table[n] << ", ";

    cerr << endl << endl << "Well completions (nodes)  [/]          : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_well_completion[n] << ", ";

    cerr << endl << endl << "Well  completion length                 : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_completion_length[n] << ", ";

    //    cerr << endl << endl << "Well seg centers pressure [Pa]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_interp_well_pressure[ n ] << ", ";

    //    cerr << endl << endl << "Well nodes pressure [Pa]               : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_well_node_pressure[ n ] << ", ";

    cerr << endl << endl << "Well seg (head or heel) pressure [Pa]  : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_pressure[n] << ", ";

    //    cerr << endl << endl << "ele interp pressure                    : ";
    //    for (unsigned n = 0; n < number_of_well_elements; n++)
    //        cerr << VO_interp_well_pressure[ n ] << ", ";

    //    cerr << endl << endl << "time to edge                           : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_time_to_edge[ n ] << ", ";

    //    cerr << endl << endl << "time out                               : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_time_out[ n ] << ", ";
    //    cerr << endl << endl << "time in                                : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_time_in[ n ] << ", ";

    cerr << endl
         << endl
         << "Well seg Fluid state                                : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_fluid_state[n] << ", ";
    //    cerr << endl << endl << "Res node Fluid state : "; for (unsigned n = 0;
    //    n < number_of_well_nodes; n++)
    //        cerr << VO_fluid_state[ n ] << ", ";
    //    cerr << endl << endl << "well index3                            : ";
    //    for (unsigned n = 0; n < number_of_well_nodes; n++)
    //        cerr << VO_well_index3[ n ] << ", ";
    //    cerr << endl << endl << "Well air saturation                    : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_air_saturation[ n ] << ", ";

    cerr << endl << endl << "Well seg saturation halite  [/]     : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_sat_halite[n] << ", ";

    cerr << endl << endl << "Well seg bulk smf  [/]     : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_smf_bulk[n] << ", ";

    cerr << endl << endl << "Well seg fluid smf  [/]     : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_smf_fluid[n] << ", ";

    cerr << endl << endl << "Well seg saturation vapor [/]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_sat_vapor[n] << ", ";

    cerr << endl << endl << "Well seg density fluid [kg/m3]         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_fluid[n] << ", ";

    //    cerr << endl << endl << "Well seg activated         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_well_segment_activated[ n ] << ", ";

    //    cerr << endl << endl << "Well seg inclination        : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_segment_inclination[n] << ", ";

    cerr << endl;
}

template <uint32_t dim> void WellModelPrototype<dim>::Display_short() {
    cerr << defaultfloat;
    cerr.precision(6);

    cerr << endl << endl << "Well Elmt heels pressure [Pa]           : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_well_node_pressure[n] << ", ";

    cerr << endl << endl << "Well-Reservoir nodes skin [-]            : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_skin[n] << ", ";

    cerr << endl << endl << "Res Elmt heel pressure [Pa]             : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_pressure[n] << ", ";

    cerr << endl << endl << "Well seg pressure [Pa]               : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_pressure[n] << ", ";

    cerr << endl << endl << "Well seg temperature [C]               : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_temperature[n] << ", ";

    cerr << endl << endl << "Well seg density bulk [kg/m3]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_bulk[n] << ", ";

    cerr << endl << endl << "Well seg smf bulk [kg/kg]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_smf_bulk[n] << ", ";

    //    cerr << endl << endl << " Well seg velocity fluid [m/s]         : ";
    //    for (unsigned n = 0; n < number_of_well_segments; n++)
    //        cerr << VO_D_well_velocity_fluid[ n ] << ", ";

    cerr << endl << endl << "Mass transfer rate [kg/s]               : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_mass_transfer_rate[n] << ", ";

    cerr << endl << endl << "Completion length [m]               : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_completion_length[n] << ", ";

    //    cerr << endl << endl << "previous water table depth:
    //    "<<previous_water_table_depth<<", water table depth:
    //    "<<water_table_depth;

    cerr << endl;
}

template <uint32_t dim> void WellModelPrototype<dim>::Display_extended() {
    cerr << endl << "Converged result";
    cerr << endl << "Timestep: " << dt << " seconds *****";

    cerr << endl << endl << "Well seg velocity liquid [m/s]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_velocity_liquid[n] << ", ";

    cerr << endl << endl << "Well seg velocity vapor [m/s]           : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_velocity_vapor[n] << ", ";

    cerr << endl << endl << "Well seg velocity fluid change`[m/s]  : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_velocity_fluid[n] - VO_D_previous_well_velocity_fluid[n]
                << ", ";

    cerr << endl << endl << "Well seg temperature change [C]         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_temperature[n] - VO_D_previous_well_temperature[n]
                << ", ";

    cerr << endl << endl << "Well seg s_enthalpy change [J/kg]       : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_s_enthalpy_bulk[n] - VO_D_previous_well_s_enthalpy_bulk[n]
                << ", ";

    cerr << endl << endl << "Well seg density fluid change [kg/m3] : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_fluid[n] - VO_D_previous_well_density_fluid[n]
                << ", ";

    cerr << endl << endl << "Well seg mass change term [kg/s]        : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_segment_length[n] * VO_D_well_cross_area[n] / dt *
                (VO_D_well_density_fluid[n] -
                 VO_D_previous_well_density_fluid[n])
             << ", ";

    cerr << endl << endl << "Well seg mass Inflow [kg/s]             : ";
    for (unsigned n = 0; n < number_of_well_segments - 1; n++)
        cerr << VO_D_well_cross_area[n + 1] * (VO_D_well_density_fluid[n + 1] *
                                               VO_D_well_velocity_fluid[n + 1])
             << ", ";

    cerr << endl << endl << "Well seg mass Outflow [kg/s]            : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_cross_area[n] *
                (VO_D_well_density_fluid[n] * VO_D_well_velocity_fluid[n])
             << ", ";

    cerr << endl << endl << "Well seg drhofdp [kg/m3.Pa]              : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_drhofdp[n] << ", ";

    cerr << endl << endl << "Well seg drhofdhb [kg/m3.J]               : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_drhofdhb[n] << ", ";

    cerr << endl << endl << "Well seg dhdp [J/Pa]                    : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_dhdp[n] << ", ";

    cerr << endl << endl << "Well seg density liquid [kg/m3]         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_liquid[n] << ", ";

    cerr << endl << endl << "Well seg density vapor [kg/m3]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_vapor[n] << ", ";

    cerr << endl << endl << "Mass transfer rate [kg/s]               : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_mass_transfer_rate[n] << ", ";

    cerr << endl << endl << "Energy transfer rate [J/s] : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_energy_transfer_rate[n] << ", ";

    cerr << endl << "Converged result";
    cerr << endl << "Timestep: " << dt << " seconds *****";

    cerr << endl << endl << "Well seg velocity fluid [m/s]         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_velocity_fluid[n] << ", ";

    cerr << endl << endl << "Well seg heads pressure [Pa]            : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_pressure[n] << ", ";

    cerr << endl << endl << "Well seg centers presure [Pa]           : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_interp_well_pressure[n] << ", ";

    cerr << endl << endl << "Well Elmt heels pressure [Pa]           : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_well_node_pressure[n] << ", ";

    cerr << endl << endl << "Res Elmt heel pressure [Pa]             : ";
    for (unsigned n = 0; n < number_of_well_nodes; n++)
        cerr << VO_pressure[n] << ", ";

    cerr << endl << endl << "Well seg temperature [C]                : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_temperature[n] << ", ";

    cerr << endl << endl << "Well seg s_enthalpy fluid  [J/kg]     : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_s_enthalpy_fluid[n] << ", ";

    cerr << endl << endl << "Well seg density fluid [kg/m3]        : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_density_fluid[n] << ", ";

    cerr << endl << "Well seg saturation liquid [/]          : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_sat_liquid[n] << ", ";

    cerr << endl << endl << "Well seg saturation vapor [/]           : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_sat_vapor[n] << ", ";

    cerr << endl << endl << "Well seg Re [/]                         : ";
    for (unsigned n = 0; n < number_of_well_segments; n++)
        cerr << VO_D_well_Re[n] << ", ";

    cerr << endl;
}

template <uint32_t dim>
void WellModelPrototype<dim>::SetTimeStep(const double &time_step) {
    dt = time_step;

    if (dt < lowest_ever_time_step)
        lowest_ever_time_step = dt;

    if (dt < 1.e-4) {
        cerr << endl << "timestep too small";
        throw Exception(FATAL_ERROR, " WellModelPrototype ", "timestep too small");
    }
}

/** Switches the operational mode, keeping the legacy flags in step and reporting
    the change. The three modes drive the top of the well and are mutually
    exclusive: as one enum they cannot overlap, which is the point. */
template <uint32_t dim>
void WellModelPrototype<dim>::SwitchControl(WellControl to) {

    if (to == control_) return;

    cerr << endl << "Well " << well_name << ": control changes from "
         << WellControlName(control_) << " to " << WellControlName(to);

    if (to == WellControl::TargetRate && with_water_table_displacement)
        cerr << endl << "  WARNING: target-rate control together with water-table"
                        " displacement is not expected to behave sensibly — both"
                        " drive the top of the well. Water-table displacement is"
                        " intended for top injection.";

    control_ = to;

    // Keep the flags the physics still reads in step with the mode.
    target_rate_active = (control_ == WellControl::TargetRate);
    injection_mode                 = (control_ == WellControl::TopInjection);
}

template <uint32_t dim>
void WellModelPrototype<dim>::ReportControl() const {
    cerr << endl << "Well " << well_name << " is driven by " << WellControlName(control_);
    switch (control_) {
        case WellControl::WellheadPressure:
            cerr << ": " << top_pressure << " Pa at " << top_temperature << " oC";
            if (config_.wellhead_pressure_ramp > 0.)
                cerr << " (ramping at " << config_.wellhead_pressure_ramp << " Pa/s)";
            break;
        case WellControl::TargetRate:
            cerr << ": " << target_rate << " kg/s ("
                 << (target_rate >= 0. ? "production" : "injection")
                 << "), current rate " << total_rate << " kg/s";
            break;
        case WellControl::TopInjection:
            cerr << ": fluid at " << config_.injection_temperature << " oC, salt mass"
                    " fraction " << config_.injection_smf;
            break;
    }
    cerr << ", water-table displacement "
         << (with_water_table_displacement ? "ON" : "off");
}

template <uint32_t dim>
void WellModelPrototype<dim>::Set_well_top_pressure(const double &pressure) {

    SwitchControl(WellControl::WellheadPressure);

    if (config_.wellhead_pressure_ramp <= 0.) {
        top_pressure = pressure;   // immediate: the default
    }
    else {
        // Move toward the requested pressure at the configured rate [Pa/s],
        // stopping exactly on it. Replaces a hardcoded 0.001 bar/min that applied
        // whenever the old `forced` argument was false. (Benoit DD/MM/YYYY)
        const double step = config_.wellhead_pressure_ramp * dt;
        if (top_pressure > pressure)
            top_pressure = std::max(pressure, top_pressure - step);
        else if (top_pressure < pressure)
            top_pressure = std::min(pressure, top_pressure + step);
    }

    t_ = top_temperature;
    p_ = top_pressure;
    smf_ = top_smf;
    wt_ = smf_ * 100.;
    x_ = Weight2XNaCl(wt_);

    // top_s_enthalpy = brine.Enthalpy();
    top_s_enthalpy = fluid.BulkProperties().h;
    cerr << endl
         << "Reset top of well with: top_temperature: " << top_temperature
         << ", top_pressure: " << top_pressure << ", top_smf: " << top_smf
         << ", top_s_enthalpy: " << top_s_enthalpy;
}

template <uint32_t dim>
void WellModelPrototype<dim>::Set_well_top_temperature(
        const double &temperature) {

    top_temperature = temperature;

    t_ = top_temperature;
    p_ = top_pressure;
    smf_ = top_smf;
    wt_ = smf_ * 100.;
    x_ = Weight2XNaCl(wt_);

    // top_s_enthalpy = brine.Enthalpy();
    top_s_enthalpy = fluid.BulkProperties().h;
    cerr << endl
         << "Reset top of well with: top_temperature: " << top_temperature
         << ", top_pressure: " << top_pressure << ", top_smf: " << top_smf
         << ", top_s_enthalpy: " << top_s_enthalpy;
}

template <uint32_t dim>
void WellModelPrototype<dim>::Set_water_table(
        const bool &toggle_water_table_displacement,
        const double &initial_water_table_depth) {
    // The only place this toggle is set — deliberately not a configuration field,
    // so the state has a single source. (Benoit DD/MM/YYYY)
    with_water_table_displacement = toggle_water_table_displacement;
    water_table_depth = previous_water_table_depth = initial_water_table_depth;

    if (with_water_table_displacement && control_ == WellControl::TargetRate)
        cerr << endl << "  WARNING: water-table displacement switched on while this"
                        " well is under target-rate control; the two are not"
                        " expected to work together.";

    cerr << endl
         << "Setting the water table at " << water_table_depth
         << " and toggling its displacement to ";
    if (!with_water_table_displacement)
        cerr << "false.";
    else
        cerr << "true.";
    cerr << endl;
}
////////////////////
template <uint32_t dim> void WellModelPrototype<dim>::Advance_well_variables() {
    // Advance Well variables
    VO_D_previous_well_pressure = VO_D_well_pressure;
    VO_D_previous_interp_well_pressure = VO_D_interp_well_pressure;

    VO_D_previous_well_s_enthalpy_fluid = VO_D_well_s_enthalpy_fluid;
    VO_D_previous_well_s_enthalpy_liquid = VO_D_well_s_enthalpy_liquid;
    VO_D_previous_well_s_enthalpy_vapor = VO_D_well_s_enthalpy_vapor;
    VO_D_previous_well_s_enthalpy_halite = VO_D_well_s_enthalpy_halite;
    VO_D_previous_well_s_enthalpy_bulk = VO_D_well_s_enthalpy_bulk;

    VO_D_previous_well_smf_fluid = VO_D_well_smf_fluid;
    VO_D_previous_well_smf_liquid = VO_D_well_smf_liquid;
    VO_D_previous_well_smf_vapor = VO_D_well_smf_vapor;
    VO_D_previous_well_smf_bulk = VO_D_well_smf_bulk;

    VO_D_previous_well_velocity_fluid = VO_D_well_velocity_fluid;

    VO_D_previous_well_density_liquid = VO_D_well_density_liquid;
    VO_D_previous_well_density_vapor = VO_D_well_density_vapor;
    VO_D_previous_well_density_halite = VO_D_well_density_halite;
    VO_D_previous_well_density_fluid = VO_D_well_density_fluid;
    VO_D_previous_well_density_bulk = VO_D_well_density_bulk;

    VO_D_previous_well_viscosity_liquid = VO_D_well_viscosity_liquid;
    VO_D_previous_well_viscosity_vapor = VO_D_well_viscosity_vapor;
    VO_D_previous_well_viscosity_fluid = VO_D_well_viscosity_fluid;

    VO_D_previous_well_sat_liquid = VO_D_well_sat_liquid;
    VO_D_previous_well_sat_vapor = VO_D_well_sat_vapor;
    VO_D_previous_well_sat_halite = VO_D_well_sat_halite;

    VO_D_previous_well_mf_liquid = VO_D_well_mf_liquid;
    VO_D_previous_well_mf_vapor = VO_D_well_mf_vapor;
    VO_D_previous_well_mf_halite = VO_D_well_mf_halite;

    VO_D_previous_well_temperature = VO_D_well_temperature;
    VO_D_previous_well_cross_area = VO_D_well_cross_area;
    VO_D_previous_well_friction_factor = VO_D_well_friction_factor;

    VO_D_previous_well_air_saturation = VO_D_well_air_saturation;
    VO_D_previous_well_water_table = VO_D_well_water_table;
    previous_water_table_depth = water_table_depth;

    VO_D_previous_well_fluid_state = VO_D_well_fluid_state;

    VO_previous_well_segment_activated = VO_well_segment_activated;
    VO_previous_well_rate_factor = VO_well_rate_factor;
    VO_previous_well_sat_halite = VO_well_sat_halite;

    previous_top_pressure = top_pressure;
    previous_total_rate = total_rate;
    Previous_TOTAL_MASS = TOTAL_MASS;

    VO_D_previous_mass_transfer_rate = VO_D_mass_transfer_rate;
    VO_D_previous_energy_transfer_rate = VO_D_energy_transfer_rate;
    VO_D_previous_salt_mass_transfer_rate = VO_D_salt_mass_transfer_rate;
}

template <uint32_t dim>
void WellModelPrototype<dim>::NR_Advance_well_variables() {
    // Advance Well variables
    VO_D_NR_previous_well_pressure = VO_D_well_pressure;
    VO_D_NR_previous_interp_well_pressure = VO_D_interp_well_pressure;

    VO_D_NR_previous_well_s_enthalpy_fluid = VO_D_well_s_enthalpy_fluid;
    VO_D_NR_previous_well_s_enthalpy_liquid = VO_D_well_s_enthalpy_liquid;
    VO_D_NR_previous_well_s_enthalpy_vapor = VO_D_well_s_enthalpy_vapor;
    VO_D_NR_previous_well_s_enthalpy_halite = VO_D_well_s_enthalpy_halite;
    VO_D_NR_previous_well_s_enthalpy_bulk = VO_D_well_s_enthalpy_bulk;

    VO_D_NR_previous_well_smf_fluid = VO_D_well_smf_fluid;
    VO_D_NR_previous_well_smf_liquid = VO_D_well_smf_liquid;
    VO_D_NR_previous_well_smf_vapor = VO_D_well_smf_vapor;
    VO_D_NR_previous_well_smf_bulk = VO_D_well_smf_bulk;

    VO_D_NR_previous_well_density_liquid = VO_D_well_density_liquid;
    VO_D_NR_previous_well_density_vapor = VO_D_well_density_vapor;
    VO_D_NR_previous_well_density_halite = VO_D_well_density_halite;
    VO_D_NR_previous_well_density_fluid = VO_D_well_density_fluid;
    VO_D_NR_previous_well_density_bulk = VO_D_well_density_bulk;

    VO_D_NR_previous_well_viscosity_liquid = VO_D_well_viscosity_liquid;
    VO_D_NR_previous_well_viscosity_vapor = VO_D_well_viscosity_vapor;
    VO_D_NR_previous_well_viscosity_fluid = VO_D_well_viscosity_fluid;

    VO_D_NR_previous_well_sat_liquid = VO_D_well_sat_liquid;
    VO_D_NR_previous_well_sat_vapor = VO_D_well_sat_vapor;
    VO_D_NR_previous_well_sat_halite = VO_D_well_sat_halite;

    VO_D_NR_previous_well_mf_liquid = VO_D_well_mf_liquid;
    VO_D_NR_previous_well_mf_vapor = VO_D_well_mf_vapor;
    VO_D_NR_previous_well_mf_halite = VO_D_well_mf_halite;

    VO_D_NR_previous_well_air_saturation = VO_D_well_air_saturation;
    VO_D_NR_previous_well_water_table = VO_D_well_water_table;
    NR_previous_water_table_depth = water_table_depth;

    VO_D_NR_previous_well_fluid_state = VO_D_well_fluid_state;

    NR_previous_top_pressure = top_pressure;
    NR_previous_total_rate = total_rate;

    VO_D_NR_previous_mass_transfer_rate = VO_D_mass_transfer_rate;
    VO_D_NR_previous_energy_transfer_rate = VO_D_energy_transfer_rate;
    VO_D_NR_previous_salt_mass_transfer_rate = VO_D_salt_mass_transfer_rate;
}

template <uint32_t dim> void WellModelPrototype<dim>::Reset() {
    DAMP = 0.1;

    exit_and_cut_dt = false;

    VO_D_well_pressure = VO_D_NR_previous_well_pressure =
            VO_D_previous_well_pressure;
    VO_D_interp_well_pressure = VO_D_NR_previous_interp_well_pressure =
            VO_D_previous_interp_well_pressure;

    VO_D_well_s_enthalpy_liquid = VO_D_NR_previous_well_s_enthalpy_liquid =
            VO_D_previous_well_s_enthalpy_liquid;
    VO_D_well_s_enthalpy_vapor = VO_D_NR_previous_well_s_enthalpy_vapor =
            VO_D_previous_well_s_enthalpy_vapor;
    VO_D_well_s_enthalpy_fluid = VO_D_NR_previous_well_s_enthalpy_fluid =
            VO_D_previous_well_s_enthalpy_fluid;
    VO_D_well_s_enthalpy_halite = VO_D_NR_previous_well_s_enthalpy_halite =
            VO_D_previous_well_s_enthalpy_halite;
    VO_D_well_s_enthalpy_bulk = VO_D_NR_previous_well_s_enthalpy_bulk =
            VO_D_previous_well_s_enthalpy_bulk;

    VO_D_well_smf_fluid = VO_D_NR_previous_well_smf_fluid =
            VO_D_previous_well_smf_fluid;
    VO_D_well_smf_liquid = VO_D_NR_previous_well_smf_liquid =
            VO_D_previous_well_smf_liquid;
    VO_D_well_smf_vapor = VO_D_NR_previous_well_smf_vapor =
            VO_D_previous_well_smf_vapor;
    VO_D_well_smf_bulk = VO_D_NR_previous_well_smf_bulk =
            VO_D_previous_well_smf_bulk;

    VO_D_well_velocity_fluid = VO_D_previous_well_velocity_fluid;

    VO_D_well_density_liquid = VO_D_NR_previous_well_density_liquid =
            VO_D_previous_well_density_liquid;
    VO_D_well_density_vapor = VO_D_NR_previous_well_density_vapor =
            VO_D_previous_well_density_vapor;
    VO_D_well_density_halite = VO_D_NR_previous_well_density_halite =
            VO_D_previous_well_density_halite;
    VO_D_well_density_fluid = VO_D_NR_previous_well_density_fluid =
            VO_D_previous_well_density_fluid;
    VO_D_well_density_bulk = VO_D_NR_previous_well_density_bulk =
            VO_D_previous_well_density_bulk;

    VO_D_well_viscosity_liquid = VO_D_NR_previous_well_viscosity_liquid =
            VO_D_previous_well_viscosity_liquid;
    VO_D_well_viscosity_vapor = VO_D_NR_previous_well_viscosity_vapor =
            VO_D_previous_well_viscosity_vapor;
    VO_D_well_viscosity_fluid = VO_D_NR_previous_well_viscosity_fluid =
            VO_D_previous_well_viscosity_fluid;

    VO_D_well_sat_liquid = VO_D_NR_previous_well_sat_liquid =
            VO_D_previous_well_sat_liquid;
    VO_D_well_sat_vapor = VO_D_NR_previous_well_sat_vapor =
            VO_D_previous_well_sat_vapor;
    VO_D_well_sat_halite = VO_D_NR_previous_well_sat_halite =
            VO_D_previous_well_sat_halite;

    VO_D_well_mf_liquid = VO_D_NR_previous_well_mf_liquid =
            VO_D_previous_well_mf_liquid;
    VO_D_well_mf_vapor = VO_D_NR_previous_well_mf_vapor =
            VO_D_previous_well_mf_vapor;
    VO_D_well_mf_halite = VO_D_NR_previous_well_mf_halite =
            VO_D_previous_well_mf_halite;

    VO_D_well_temperature = VO_D_previous_well_temperature;
    VO_D_well_cross_area = VO_D_previous_well_cross_area;
    VO_D_well_friction_factor = VO_D_previous_well_friction_factor;

    VO_D_well_air_saturation = VO_D_NR_previous_well_air_saturation =
            VO_D_previous_well_air_saturation;
    VO_D_well_water_table = VO_D_NR_previous_well_water_table =
            VO_D_previous_well_water_table;
    water_table_depth = NR_previous_water_table_depth =
            previous_water_table_depth;

    VO_D_well_fluid_state = VO_D_previous_well_fluid_state;

    VO_well_segment_activated = VO_previous_well_segment_activated;
    VO_well_rate_factor = VO_previous_well_rate_factor;
    VO_well_sat_halite = VO_previous_well_sat_halite;

    top_pressure = NR_previous_top_pressure = previous_top_pressure;

    TOTAL_MASS = Previous_TOTAL_MASS;
    total_rate = previous_total_rate;

    VO_D_mass_transfer_rate = VO_D_NR_previous_mass_transfer_rate =
            VO_D_previous_mass_transfer_rate;
    VO_D_energy_transfer_rate = VO_D_NR_previous_energy_transfer_rate =
            VO_D_previous_energy_transfer_rate;
    VO_D_salt_mass_transfer_rate = VO_D_NR_previous_salt_mass_transfer_rate =
            VO_D_previous_salt_mass_transfer_rate;
}

///////////////////
template <uint32_t dim>
void WellModelPrototype<dim>::Weighted_Next_Guess(double w1, double w2) {
    // cerr<<endl<<"//////////////// MAKING NEXT GUESS ////////////////"<<endl;

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        if (VO_D_well_density_fluid[n] > 0 &&
                VO_D_NR_previous_well_density_fluid[n] > 0) {
            if (VO_D_well_density_fluid[n] / VO_D_NR_previous_well_density_fluid[n] >
                    1.05)
                VO_D_well_density_fluid[n] =
                        VO_D_NR_previous_well_density_fluid[n] * 1.05;

            if (VO_D_well_density_fluid[n] / VO_D_NR_previous_well_density_fluid[n] <
                    0.95)
                VO_D_well_density_fluid[n] =
                        VO_D_NR_previous_well_density_fluid[n] * 0.95;

            //            if( DAMP < 0.1 )
            //            {
            //                if ( VO_D_well_density_fluid [ n ] -
            //                VO_D_NR_previous_well_density_fluid [ n ] > 0.5 )
            //                    VO_D_well_density_fluid [ n ] =
            //                    VO_D_NR_previous_well_density_fluid [ n ] + 0.5;

            //                if ( VO_D_well_density_fluid [ n ] -
            //                VO_D_NR_previous_well_density_fluid [ n ] < -0.5 )
            //                    VO_D_well_density_fluid [ n ] =
            //                    VO_D_NR_previous_well_density_fluid [ n ] - 0.5;
            //            }
        }

        if (VO_D_well_sat_halite[n] > 0.1 &&
                VO_D_NR_previous_well_sat_halite[n] > 0) {
            if (VO_D_well_sat_halite[n] / VO_D_NR_previous_well_sat_halite[n] > 1.05)
                VO_D_well_sat_halite[n] = VO_D_NR_previous_well_sat_halite[n] * 1.05;

            if (VO_D_well_sat_halite[n] / VO_D_NR_previous_well_sat_halite[n] < 0.95)
                VO_D_well_sat_halite[n] = VO_D_NR_previous_well_sat_halite[n] * 0.95;
        }

        //        if( VO_D_well_water_table [ n ] > 0 &&
        //        VO_D_NR_previous_well_water_table [ n ] > 0 )
        //        {
        //            VO_D_mass_transfer_rate [ n ] = w1 * VO_D_mass_transfer_rate [
        //            n ] + w2 * VO_D_NR_previous_mass_transfer_rate [ n ];
        //            VO_D_salt_mass_transfer_rate [ n ] = w1 *
        //            VO_D_salt_mass_transfer_rate [ n ] + w2 *
        //            VO_D_NR_previous_salt_mass_transfer_rate [ n ];
        //            VO_D_energy_transfer_rate [ n ] = w1 *
        //            VO_D_energy_transfer_rate [ n ] + w2 *
        //            VO_D_NR_previous_energy_transfer_rate [ n ];
        //        }
    }
}

template <uint32_t dim>
void WellModelPrototype<dim>::WriteSolutionToFile(const std::string &filename) {
    std::ofstream SolutionFile(filename);
    cout << "Saving file: " << filename << endl;

    SolutionFile << "W_pressure" << "\t" << "W_temperature" << "\t"
                 << "R_pressure" << "\t" << "R_temperature" << "\t" << "velocity"
                 << "\t" << "h_vapor" << "\t" << "h_liquid" << "\t" << "h_fluid"
                 << "\t" << "smf_vapor" << "\t" << "smf_liquid" << "\t"
                 << "smf_fluid" << "\t" << "rho_vapor" << "\t" << "rho_liquid"
                 << "\t" << "rho_fluid" << "\t" << "sat_liquid" << "\t"
                 << "sat_vapor" << "\t" << "sat_air" << "\t" << "sat_halite" <<

                    "\t" << "friction_factor" << "\t" << "cross_area" << "\t" << "visc_liquid"
                 << "\t" << "visc_vapor" << "\t" << "visc_fluid" << "\t"
                 << "Re_number" <<

                    "\t" << "heel_depth" << "\t" << "water_table" << "\t" << "top_pressure"
                 << "\t"
                 << "total_rate"
                    "\n";

    for (size_t ix = 0; ix < number_of_well_segments; ix++) {
        SolutionFile << VO_D_well_pressure[ix] << "\t" << VO_D_well_temperature[ix]
                        << "\t" << VO_D_pressure[ix] << "\t" << VO_D_temperature[ix]
                           << "\t" << VO_D_well_velocity_fluid[ix] << "\t"
                           << VO_D_well_s_enthalpy_vapor[ix] << "\t"
                           << VO_D_well_s_enthalpy_liquid[ix] << "\t"
                           << VO_D_well_s_enthalpy_fluid[ix] << "\t"
                           << VO_D_well_smf_vapor[ix] << "\t" << VO_D_well_smf_liquid[ix]
                              << "\t" << VO_D_well_smf_fluid[ix] << "\t"
                              << VO_D_well_density_vapor[ix] << "\t"
                              << VO_D_well_density_liquid[ix] << "\t"
                              << VO_D_well_density_fluid[ix] << "\t"
                              << VO_D_well_sat_liquid[ix] << "\t" << VO_D_well_sat_vapor[ix]
                                 << "\t" << VO_D_well_air_saturation[ix] << "\t"
                                 << VO_D_well_sat_halite[ix] <<

                                    "\t" << VO_D_well_friction_factor[ix] << "\t"
                                 << VO_D_well_cross_area[ix] << "\t"
                                 << VO_D_well_viscosity_liquid[ix] << "\t"
                                 << VO_D_well_viscosity_vapor[ix] << "\t"
                                 << VO_D_well_viscosity_fluid[ix] << "\t" << VO_D_well_Re[ix] <<

                                    "\t" << VO_D_segment_heel_depth[ix] << "\t" << water_table_depth << "\t"
                                 << top_pressure << "\t" << total_rate;

        SolutionFile << std::endl;
    }
}

template <uint32_t dim> void WellModelPrototype<dim>::Initialize_Water_Table() {

    // water_table_depth = previous_water_table_depth;

    wt_index = 0;
    for (unsigned n = 0; n < number_of_well_segments; n++) {
        if ((water_table_depth == VO_D_segment_heel_depth[n]) or
                (water_table_depth > VO_D_segment_heel_depth[n])) {
            break;
        }
        wt_index++;
    }

    for (unsigned n = 0; n < wt_index + 1; n++) {
        // cerr<<endl<<"Loop below, index "<<n<<endl;
        if (water_table_depth < VO_D_segment_heel_depth[n]) {
            // cerr<<endl<<"water level is not inside"<<endl;
            VO_D_well_water_table[n] = 0;
            VO_D_well_air_saturation[n] = 1.;
        }

        else {
            if (water_table_depth == VO_D_segment_heel_depth[n]) {
                // cerr<<endl<<"Water level is at bottom"<<endl;
                VO_D_well_water_table[n] = 0;
                VO_D_well_air_saturation[n] = 1.;
            }
            else {
                // cerr<<endl<<"Water level is inside but not at bottom"<<endl;
                VO_D_well_water_table[n] = 1;
                VO_D_well_air_saturation[n] =
                        (1 -
                         (water_table_depth - VO_D_segment_heel_depth[n]) /
                         (VO_D_segment_length[n] * cos(VO_D_segment_inclination[n])));
            }
        }
    }

    for (unsigned n = wt_index + 1; n < number_of_well_segments; n++) {
        // cerr<<endl<<"Loop above, index "<<n<<endl;
        VO_D_well_water_table[n] = 2;
        VO_D_well_air_saturation[n] = 0.;
    }

    VO_D_previous_well_water_table = VO_D_well_water_table;
}

template <uint32_t dim>
void WellModelPrototype<dim>::Find_Water_Table(bool move) {
    double previous_total_mass = 0.;
    double mass_exchanged = 0.;
    double total_mass = 0.;

    double total_height(0.);
    double density_fluid(0.);

    double test_water_table_depth = water_table_depth;
    double target_water_table_depth;

    /////////////////
    // Find segment where the water table is currently located
    wt_index = 0;
    for (unsigned n = 0; n < number_of_well_segments; n++) {
        if ((water_table_depth == VO_D_segment_heel_depth[n]) or
                (water_table_depth > VO_D_segment_heel_depth[n])) {
            break;
        }
        wt_index++;
    }

    //////////////////////////////////////
    enthalpy_below_wt = VO_D_well_s_enthalpy_bulk[wt_index + 1];
    smf_below_wt = VO_D_well_smf_bulk[wt_index + 1];
    temp_below_wt = VO_D_well_temperature[wt_index + 1];

    //    cerr<<endl<<"Enthalpy below wt is: "<<enthalpy_below_wt;
    //    cerr<<endl<<"Smf below wt is: "<<smf_below_wt<<endl;
    /////////////////////////////////////

    /////////////////
    for (unsigned n = 0; n < number_of_well_segments; n++) {
        previous_total_mass += VO_D_previous_well_density_fluid[n] *
                (1. - VO_D_previous_well_air_saturation[n]) *
                VO_D_segment_length[n] *
                VO_D_previous_well_cross_area[n] *
                (1. - VO_D_previous_well_sat_halite[n]);

        mass_exchanged += VO_D_mass_transfer_rate[n] * dt;
    }

    if (mass_injected > 0.)
        mass_exchanged += mass_injected;

    //    cerr<<endl<<" The timestep is "<<dt;
    //    cerr<<endl<<" The mass exchanged is "<<mass_exchanged<<". In tons per
    //    hour: "<<mass_exchanged/dt*3.6<<"&&&&&&&&&&&&&&&&&&"; cerr<<endl<<" The
    //    variable total rate is: "<<total_rate<<". In tons per hour:
    //    "<<total_rate*3.6<<"&&&&&&&&&&&&&&&&&&"; cerr<<endl<<" Among which mass
    //    injected is "<<mass_injected;

    total_mass =
            previous_total_mass + mass_exchanged; // not caring for fluid moving out
    // of the well during upflow
    // cerr<<endl<<"The previous total mass is : "<<previous_total_mass<<", adding
    // mass exchanged: "<<total_mass;

    ////////////
    // Calculate fluid height from the bottom up
    int k(number_of_well_segments - 1);
    while (total_mass > 0. && k > 0) {
        density_fluid = 0.;
        int l = k;

        while (!(density_fluid > 0.)) {
            density_fluid = VO_D_well_density_fluid[l];
            l++;
        }

        if (VO_D_segment_length[k] <
                total_mass / (density_fluid * VO_D_well_cross_area[k] *
                              (1. - VO_D_previous_well_sat_halite[k]))) {
            total_height += VO_D_segment_length[k] * cos(VO_D_segment_inclination[k]);

            total_mass -= density_fluid * VO_D_well_cross_area[k] *
                    VO_D_segment_length[k] *
                    (1. - VO_D_previous_well_sat_halite[k]);

        }

        else {
            total_height += total_mass /
                    (density_fluid * VO_D_well_cross_area[k] *
                     (1. - VO_D_previous_well_sat_halite[k])) *
                    cos(VO_D_segment_inclination[k]);

            total_mass = 0.;
        }

        k--;
    }

    if (total_mass > 0.) { // there is mass left to fill segment 0
        density_fluid = 0.;
        int l = 0;
        while (!(density_fluid > 0.)) {
            density_fluid = VO_D_well_density_fluid[l];
            l++;
        }

        if (VO_D_segment_length[0] >
                total_mass / (density_fluid * VO_D_well_cross_area[0] *
                              (1. - VO_D_previous_well_sat_halite[0]))) {
            total_height += total_mass /
                    (density_fluid * VO_D_well_cross_area[0] *
                     (1. - VO_D_previous_well_sat_halite[0])) *
                    cos(VO_D_segment_inclination[0]);
        }

        else {
            total_height += total_mass /
                    (density_fluid * VO_D_well_cross_area[0] *
                     (1. - VO_D_previous_well_sat_halite[0])) *
                    cos(VO_D_segment_inclination[0]);
        }
    }

    water_table_depth =
            VO_D_segment_heel_depth[number_of_well_segments - 1] + total_height;
    target_water_table_depth = water_table_depth;

    double wt_disp = (water_table_depth - test_water_table_depth);

    if (!move)
        wt_disp = 0.;

    //    //        //TEST FEB 2022
    // if(wt_disp > 0.5) wt_disp = 0.5;
    // if(wt_disp < -0.5) wt_disp = -0.5;

    WT_VELOCITY = -(water_table_depth - previous_water_table_depth) / dt;
    // cerr<<endl<<"WT_VELOCITY (based on height): "<<WT_VELOCITY;
    velocity_wt = WT_VELOCITY;

    water_table_depth = test_water_table_depth + wt_disp;

    //    cerr<<endl<<"/////////////////////////"<<endl<<
    //          "target: "<<target_water_table_depth<<", test:
    //          "<<test_water_table_depth<<", actual: "<<water_table_depth<<
    //          endl<<"/////////////////////////"<<endl;

    if (!water_table_depth_has_converged) {
        if (abs(water_table_depth - target_water_table_depth) < 0.1)
            water_table_depth_has_converged = true;
        else if (abs(water_table_depth - target_water_table_depth) * density_fluid *
                 grav_acc <
                 1000.)
            water_table_depth_has_converged = true;
        else
            water_table_depth_has_converged = false;
    }

    //    // the water table cannot rise above the top of segment 0
    //    if ( water_table_depth > 0. or essentiallyEqual(water_table_depth,
    //    0., 1.e-6 )) //surface = 0 meters!!!!
    //    {
    //        water_table_depth = 0.;
    //    }

    // the water table cannot rise above the top of segment 0
    if (water_table_depth > top_of_well or
            essentiallyEqual(water_table_depth, top_of_well,
                             1.e-6)) { // surface = 0 meters!!!!
        water_table_depth = top_of_well;
    }

    ///////////////////
    // Update water table indexes and air saturation
    // Find segment where the water table is now located
    wt_index = 0;
    for (unsigned n = 0; n < number_of_well_segments; n++) {
        // cerr<<endl<<endl<<"Looping through segments, trying to find if the water
        // level is at the bottom or inside the segment"<<endl;
        if ((water_table_depth == VO_D_segment_heel_depth[n]) or
                (water_table_depth > VO_D_segment_heel_depth[n])) {
            break;
        }
        wt_index++;
    }

    for (unsigned n = 0; n < wt_index + 1; n++) {
        // cerr<<endl<<"Loop below, index "<<n<<endl;
        if (water_table_depth < VO_D_segment_heel_depth[n]) {
            // cerr<<endl<<"water level is not inside"<<endl;
            VO_D_well_water_table[n] = 0;
            VO_D_well_air_saturation[n] = 1.;
        }

        else {
            //            if (  water_table_depth == VO_D_segment_heel_depth  [ n ]  )
            //            {
            //                //cerr<<endl<<"Water level is at bottom"<<endl;
            //                VO_D_well_water_table [ n ] = 0;
            //                //VO_D_well_air_saturation [ n ] = 1.;
            //            }
            //            else
            //            {
            // cerr<<endl<<"Water level is inside but not at bottom"<<endl;
            VO_D_well_water_table[n] = 1;
            VO_D_well_air_saturation[n] =
                    (1 - (water_table_depth - VO_D_segment_heel_depth[n]) /
                     VO_D_segment_length[n]);
            //}
        }
    }

    for (unsigned n = wt_index + 1; n < number_of_well_segments; n++) {
        // cerr<<endl<<"Loop above, index "<<n<<endl;
        VO_D_well_water_table[n] = 2;
        VO_D_well_air_saturation[n] = 0.;
    }

    // cerr<<"wt_index: "<<wt_index<<", density there: "<<density_fluid<<endl;
    // if ( water_table_depth_has_converged ) cerr<<"The water table depth has
    // converged!"<<endl;
    //         cerr << "Well seg water table index      : ";
    //         for (unsigned n = 0; n < number_of_well_segments; n++)
    //             cerr << VO_D_well_water_table[ n ] << ", ";
}

template <uint32_t dim>
void WellModelPrototype<dim>::Top_Injection(double mass_rate_injected) {
    // OLD CODE , NEEDS  REVIEW, INCLINATION AND SAT HALITE MISSING IN PLACES
    // Injection is limited to filling the well!, no pressurization is done...

    if (mass_rate_injected <= 0.)
        throw csmp::Exception(ERROR, "WellModelPrototype::Top_Injection", well_name,
                              "the injection rate must be positive. To stop top "
                              "injection, put the well under another control mode "
                              "(Set_well_top_pressure or Set_target_rate) "
                              "— a non-positive rate is no longer a way of saying "
                              "'off'.");

    SwitchControl(WellControl::TopInjection);

    // Injected fluid properties, from the configuration (they used to be
    // hardcoded here as 10 oC and zero salinity). (Benoit DD/MM/YYYY)
    t_ = config_.injection_temperature;
    p_ = top_pressure;
    smf_ = config_.injection_smf;
    wt_ = smf_ * 100.;
    x_ = Weight2XNaCl(wt_);

    mass_injected = 0.;

    // top_s_enthalpy = brine.Enthalpy();
    // top_density = brine.Density();
    top_s_enthalpy = fluid.BulkProperties().h;
    top_density = fluid.BulkProperties().rho;

    injection_mode = true;
    mass_injected = mass_rate_injected * dt;

    std::fill(VO_D_mass_Inj_source.begin(), VO_D_mass_Inj_source.end(),
              0); // [kg]
    std::fill(VO_D_energy_Inj_source.begin(), VO_D_energy_Inj_source.end(),
              0); // [J]

    // DOES NOT TAKE INTO ACCOUNT INCLINATION!?

    volume_injected = mass_rate_injected * dt / top_density;
    enthalpy_injected = mass_rate_injected * dt * top_s_enthalpy;

    //////////////////////////////////////////////
    //   INJECTION

    double distance_injected_below = 0.;
    double distance_injected_at_wt = 0.;
    volume_injected_at_wt = 0.;
    volume_injected_below = 0.;

    cerr << endl
         << "Water table is at" << water_table_depth
         << " after normal displacement, volume injected is: "
         << volume_injected;

    // Find segment where the water table is located

    wt_index = 0;
    for (unsigned n = 0; n < number_of_well_segments; n++) {
        // cerr<<endl<<endl<<"Looping through segments, trying to find if the
        // water level is at the bottom or inside the segment"<<endl;
        if ((water_table_depth == VO_D_segment_heel_depth[n]) or
                (water_table_depth > VO_D_segment_heel_depth[n])) {
            break;
        }
        wt_index++;
    }

    cerr << endl
         << "The water table is in seg " << wt_index
         << ", now we add injected fluid on top";
    volume_injected_at_wt = volume_injected;
    distance_injected_at_wt = volume_injected / VO_D_well_cross_area[wt_index];

    double volume_left_to_inject = volume_injected;

    double next_seg_depth = 0.;

    for (int n = wt_index; n > -1; n--) {
        if (volume_left_to_inject > 0.) {
            if (n < 0)
                n = 0;

            volume_injected_below = 0.;

            cerr << endl << endl << "In seg " << n;
            double distance_to_inject =
                    volume_left_to_inject / VO_D_well_cross_area[n];

            if (n > 0)
                next_seg_depth = VO_D_segment_heel_depth[n - 1];
            else
                next_seg_depth = 0.;

            cerr << endl << "Distance to inject: " << distance_to_inject;

            if (water_table_depth + distance_to_inject > next_seg_depth) {
                cerr << endl << "Injection is filling up segment " << n;
                cerr << endl << "Above seg heel is at " << next_seg_depth;
                distance_injected_below =
                        next_seg_depth -
                        water_table_depth; // INCLINATION MISSING!!!!!!!!!!!!!!
                volume_injected_below =
                        distance_injected_below * VO_D_well_cross_area[n];

                water_table_depth += distance_injected_below;
                volume_left_to_inject -= volume_injected_below;
            }

            else {
                cerr << endl << "Injection is partially filling seg " << n;

                distance_injected_below =
                        volume_left_to_inject / VO_D_well_cross_area[n];
                volume_injected_below = volume_left_to_inject;

                water_table_depth += distance_injected_below;
                volume_left_to_inject = 0.;
            }

            if (volume_injected_below > 0.) {
                enthalpy_injected_below =
                        volume_injected_below * top_density * top_s_enthalpy;

                cerr << endl << "Volume injected: " << volume_injected_below;
                cerr << endl << "Water table depth: " << water_table_depth;
                cerr << endl << "Volume left to inject: " << volume_left_to_inject;

                double fluid_vol = VO_D_segment_length[n] * VO_D_well_cross_area[n] *
                        (1. - VO_D_well_air_saturation[n]) *
                        (1. - VO_D_well_sat_halite[n]);

                // smf_ =  VO_D_well_smf_fluid [ n ];
                // h_fluid_ = VO_D_well_s_enthalpy_fluid [ n ];

                smf_ = VO_D_well_smf_bulk[n];            // June 2022
                h_fluid_ = VO_D_well_s_enthalpy_bulk[n]; // June 2022
                m_fluid_ = VO_D_well_density_fluid[n] * fluid_vol;
                H_current_ = m_fluid_ * h_fluid_;
                double m_salt_f = m_fluid_ * smf_;

                H_current_ += enthalpy_injected_below;

                cerr << endl
                     << "On segment " << n << ", due to injection we add "
                     << enthalpy_injected_below << " Joules and "
                     << volume_injected_below *top_density << " kg.";

                cerr << endl << "m_fluid was " << m_fluid_;
                m_fluid_ += volume_injected_below * top_density;
                cerr << endl << "m_fluid is now " << m_fluid_;

                //                cerr<<endl<<"hf was "<<VO_D_well_s_enthalpy_fluid [
                //                n ]; VO_D_well_s_enthalpy_fluid [ n ] = H_current_ /
                //                m_fluid_; cerr<<endl<<"hf is now
                //                "<<VO_D_well_s_enthalpy_fluid [ n ];

                // June 2022
                cerr << endl << "hb was " << VO_D_well_s_enthalpy_bulk[n];
                VO_D_well_s_enthalpy_bulk[n] = H_current_ / m_fluid_;
                cerr << endl << "hb is now " << VO_D_well_s_enthalpy_bulk[n];

                cerr << endl << "densityf was " << VO_D_well_density_fluid[n];
                VO_D_well_density_fluid[n] =
                        m_fluid_ / (fluid_vol + volume_injected_below);
                cerr << endl << "densityf is now " << VO_D_well_density_fluid[n];

                cerr << endl << "smff was " << VO_D_well_smf_fluid[n];
                VO_D_well_smf_fluid[n] = m_salt_f / m_fluid_;
                cerr << endl << "smff is now " << VO_D_well_smf_fluid[n];
            }
        }
    }

    if (volume_left_to_inject > 1.e-3) {
        cerr << endl
             << "Volume left to inject is positive after filling up the well!, "
                "we increase top pressure";
        // top_pressure=2.e5 + (volume_left_to_inject/VO_D_well_cross_area [ 0
        // ])/10.*1.e5;
        top_pressure *= 1.05;
    }
    else
        top_pressure *= 0.99;

    if (top_pressure < 1.e5)
        top_pressure = 1.e5;

    Initialize_Water_Table();

    cerr << endl
         << "UPDATED WT_INDEX= " << wt_index << endl
         << "*****************";
    cerr << endl << "Water table is now at" << water_table_depth << endl;
}

template <uint32_t dim>
void WellModelPrototype<dim>::Set_target_rate(
        const bool &trigger, const double &ext_target_rate) {
    target_rate = ext_target_rate;
    SwitchControl(trigger ? WellControl::TargetRate
                          : WellControl::WellheadPressure);
    if (trigger)
        cerr << endl << "Well " << well_name << ": target rate set to " << target_rate
             << " kg/s (" << (target_rate >= 0. ? "production" : "injection") << ")";
}

template <uint32_t dim> void WellModelPrototype<dim>::PAUSE() {
    string dummy;
    cerr << " press ENTER to continue... " << endl;
    getline(cin, dummy);
}

template <uint32_t dim>
bool WellModelPrototype<dim>::Check_density_and_rate_convergence() {
    bool density_has_converged(false);

    double max_density_rel_change(0.);
    double density_rel_change(0.);
    double TOTAL_RATE(0.);
    double worst_density_diff(0.);

    unsigned worst_node(0);
    double worst_density_last(0.);
    double worst_density_current(0.);
    double worst_density_last_true(0.);

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        if (VO_D_well_density_fluid[n] > 0. &&
                VO_D_NR_previous_well_density_fluid[n] > 0.) {

            // TOTAL_RATE += VO_D_mass_transfer_rate[n];
            // cerr <<endl<< "Segment "<<n<<" density current: "<<
            // VO_D_well_density_fluid[ n ]<<", density previous NR:
            // "<<VO_D_NR_previous_well_density_fluid[ n ]<<", density truw:
            // "<<VO_D_well_last_true_density_fluid[ n ];

            density_rel_change = abs(VO_D_well_last_true_density_fluid[n] -
                                     VO_D_well_density_fluid[n]) /
                    VO_D_well_density_fluid[n];

            if (abs(VO_D_NR_previous_well_density_fluid[n] -
                    VO_D_well_density_fluid[n]) > abs(worst_density_diff)) {
                worst_density_diff =
                        VO_D_NR_previous_well_density_fluid[n] - VO_D_well_density_fluid[n];
                worst_node = n;
                worst_density_last = VO_D_NR_previous_well_density_fluid[n];
                worst_density_current = VO_D_well_density_fluid[n];
                worst_density_last_true = VO_D_well_last_true_density_fluid[n];
            }

            if (abs(VO_D_well_last_true_density_fluid[n] -
                    VO_D_well_density_fluid[n]) < 0.1) {
                density_rel_change = 0.;
            }

            if (density_rel_change > max_density_rel_change) {
                max_density_rel_change = density_rel_change;
            }
        }
    }

    // cerr<<endl<<"Worst density diff: "<<worst_density_diff<<", on node:
    // "<<worst_node<<". True density: "<<worst_density_last_true<<", current
    // density: "<<worst_density_current<<", last NR density:
    // "<<worst_density_last<<endl; cerr<<endl<<"Max density rel change =
    // "<<max_density_rel_change; cerr<<endl<<"Rate diff = "<<abs(
    // NR_previous_total_rate - total_rate );

    if (max_density_rel_change <
            0.05 /*&& abs( TOTAL_RATE - true_total_rate ) < 0.01*/) {
        density_has_converged = true;
    }

    return density_has_converged;
}

template <uint32_t dim> void WellModelPrototype<dim>::Check_mass_balance() {
    // NOT GOOD, USE FOR ROUGH ESTIMATE, cross area needs to be multiplied by
    // 1-Shalite in some places

    double total_mass_exchange(0.);
    double total_mass_out_at_head(0.);

    TOTAL_MASS = 0.;

    for (unsigned n = 0; n < number_of_well_segments; n++) {
        TOTAL_MASS += VO_D_well_density_fluid[n] *
                (1. - VO_D_well_air_saturation[n]) * VO_D_segment_length[n] *
                VO_D_well_cross_area[n];
        total_mass_exchange += VO_D_mass_transfer_rate[n] * dt;
    }

    if (VO_D_well_velocity_fluid[0] < 0. && !(water_table_depth < 0.)) {
        total_mass_out_at_head = VO_D_well_density_fluid[0] *
                VO_D_well_cross_area[0] *
                VO_D_well_velocity_fluid[0] * dt;
    }

    if (abs(TOTAL_MASS - Previous_TOTAL_MASS - total_mass_exchange -
            total_mass_out_at_head) > abs(WORST_MASS_MISMATCH) &&
            Previous_TOTAL_MASS > 0.)
        WORST_MASS_MISMATCH = TOTAL_MASS - Previous_TOTAL_MASS -
                total_mass_exchange - total_mass_out_at_head;

    //    cerr<<endl<<"////////////////////////////"<<endl<<"   TOTAL MASS BEFORE:
    //    "<<Previous_TOTAL_MASS<<", TOTAL MASS CURRENT: "<<TOTAL_MASS<<", MASS
    //    EXCHANGED: "<<total_mass_exchange<<", MASS EXCAPED AT HEAD:
    //    "<<total_mass_out_at_head<<endl; cerr<<"   MASS MISMATCH =
    //    "<<TOTAL_MASS-Previous_TOTAL_MASS-total_mass_exchange-total_mass_out_at_head<<",
    //    the timestep was: "<<dt<<", mismatch in Kg/s:
    //    "<<(TOTAL_MASS-Previous_TOTAL_MASS-total_mass_exchange-total_mass_out_at_head)/dt<<endl;
    //    cerr<<"   WORST MASS MISMATCH= "<<WORST_MASS_MISMATCH<<endl;
    //    cerr<<"////////////////////////////"<<endl;

    //    if(abs(TOTAL_MASS-Previous_TOTAL_MASS-total_mass_exchange-total_mass_out_at_head)
    //    > 200 && Previous_TOTAL_MASS> 0.)
    //        PAUSE();
}

template <uint32_t dim>
bool WellModelPrototype<dim>::Check_too_large_change_in_total_rate() {
    bool rate_changing_too_much = false;

    if (abs(total_rate - previous_total_rate) > 2.) {
        cerr << endl
             << "Total rate is changing a lot in one step!! Previous: "
             << previous_total_rate << ", current: " << total_rate;
        rate_changing_too_much = true;
    }
    else
        cerr << endl
             << "Total rate change ok. Previous: " << previous_total_rate
             << ", current: " << total_rate;

    return !rate_changing_too_much;
}

template class WellModelPrototype<1U>;
template class WellModelPrototype<2U>;
template class WellModelPrototype<3U>;

} // end namespace csmp