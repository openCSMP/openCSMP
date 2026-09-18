// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_PHX_VARIABLENAMES_H
#define CVFEM_PHX_VARIABLENAMES_H

#include "CSMP_definitions.h"

/*   Changelog
     February 2015, Philipp Weis:
     - initial port to CSMP++.
     2025, Benoit LC:
     - added air phase (lhs/rhs vectors, fv_transport, upwind, grav, mobility,
       transient/reset variable lists)
*/

using namespace std;

namespace csmp {

/**
     @class CVFEM_PHX_VariableNames CVFEM_PHX_VariableNames.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes

     @section motivation Motivation
      Class to bundle all variable names needed for CVFEM scheme.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
     Applies source term at the node.

     @endcode

     @section dependencies Dependencies
     Tailored for CVFEM_PHX_Scheme

     @section issues Known issues

     @section testing Testing
     testing was done in the period before publication in 2014.

  */
class CVFEM_PHX_VariableNames {
public:

    CVFEM_PHX_VariableNames(
        bool thermodymamic_density_in_gravity_term = true,
        bool with_zinc = false,
        bool with_lithium = false);

    std::vector<std::string>
        pres_grad_variables;

    std::vector<std::string>
        pore_visitor_variables;

    std::vector<std::string>
        lhs_liquid_vector,
        rhs_liquid_vector,
        lhs_vapor_vector,
        rhs_vapor_vector,
        lhs_air_vector,
        rhs_air_vector;

    std::vector<std::string>
        fv_transport_liquid_variables,
        fv_transport_vapor_variables,
        fv_transport_air_variables;


    std::vector<std::string>
        upwind_control_variables;

    std::vector<std::string>
        relperm_visc,
        densities,
        saturations,
        velocities,
        pore_velocities,
        time_control;

    std::vector<std::string>
        capacitance_lhs_variables,
        conductance_variables,
        capacitance_rhs_variables,
        heat_bottom_variables;

    std::vector<std::string>
        capacitance_lhs_p_variables,
        conductance_p_upwind_liquid_variables,
        conductance_p_upwind_vapor_variables,
        conductance_p_upwind_air_variables,
        capacitance_rhs_p_variables,
        grav_liq_variables,
        grav_vap_variables,
        grav_air_variables,
        source_p_variables,
        source_p2_variables;

    std::vector<std::string>
        transient_variables,
        full_reset_variables;

    // NOTE: AddAdvectionVariable adds a tracer species carried only by the H2O-NaCl
    // phases (liquid + vapor).  Air does not carry salt, magmatic mass, lithium, or
    // tracer species, so no air variant is needed here.
    //
    // (Benoit 18/08/2026) The new_rhs_liquid / new_rhs_vapor arguments name the
    // species mobility fields. They are registered and validated, but transport
    // DOES NOT READ THEM: every content property is advected by rebuilding its
    // facet flux from the carrier mass flux at the donor's current content/mass
    // ratio. The mobility fields are diagnostic/output only. Setting one to zero
    // will NOT stop the species advecting. See
    // ExplicitFiniteVolumeTransportPHX::DetermineFacetFluxRatesOnce and
    // StencilProcessorPHX::RebuildContentFluxFromMass.
    void AddAdvectionVariable( std::string balanced_variable,
                              std::string new_lhs_liquid, std::string new_rhs_liquid,
                              std::string new_lhs_vapor,  std::string new_rhs_vapor  );

    void ActivateAir(bool thermodymamic_density_in_gravity_term);

};

inline CVFEM_PHX_VariableNames::CVFEM_PHX_VariableNames(
    bool thermodymamic_density_in_gravity_term,
    bool with_zinc,
    bool with_lithium)
{

    // variable names for constructor arguments
    // please see individual constructors and
    // the CVFEM_PHX_Scheme-Class

    // pressure gradient
    pres_grad_variables.push_back("fluid pressure");
    pres_grad_variables.push_back("reference pressure");
    pres_grad_variables.push_back("permeability tensor");
    //pres_grad_variables.push_back("permeability");
    pres_grad_variables.push_back("KgradP");
    pres_grad_variables.push_back("gradP scaling");

    // pore volume
    pore_visitor_variables.push_back("nodal porosity");
    pore_visitor_variables.push_back("bulk volume");
    pore_visitor_variables.push_back("pore volume");

    // ── Mass and enthalpy advection variables per phase ───────────────────────
    // liquid
    lhs_liquid_vector.push_back("fluid mass liquid");
    rhs_liquid_vector.push_back("liquid mass mobility");
    lhs_liquid_vector.push_back("enthalpy content liquid");
    rhs_liquid_vector.push_back("liquid enthalpy mobility");

    // vapor
    lhs_vapor_vector.push_back("fluid mass vapor");
    rhs_vapor_vector.push_back("vapor mass mobility");
    lhs_vapor_vector.push_back("enthalpy content vapor");
    rhs_vapor_vector.push_back("vapor enthalpy mobility");

    // ── Air advection: mass and enthalpy only ─────────────────────────────────
    // NOT gated by with_air: fv_transport_air is a value member of Coupled_reservoir
    // and CVFEM_PHX_Scheme, so its constructor arguments must always be valid.
    // When with_air=false, fv_transport_air is constructed but never called
    // (ThreePhaseTransportPHX gates all calls on with_air_).
    // air carries no salt, magmatic mass, or tracers
    lhs_air_vector.push_back("fluid mass air");
    rhs_air_vector.push_back("air mass mobility");
    lhs_air_vector.push_back("enthalpy content air");
    rhs_air_vector.push_back("air enthalpy mobility");

    // ── Salt advection ────────────────────────────────────────────────────────
    lhs_liquid_vector.push_back("salt content liquid");
    rhs_liquid_vector.push_back("liquid salt mobility");
    lhs_vapor_vector.push_back("salt content vapor");
    rhs_vapor_vector.push_back("vapor salt mobility");
    // air: no salt

    // ── Magmatic mass advection ───────────────────────────────────────────────
    lhs_liquid_vector.push_back("magmatic fluid mass liquid");
    rhs_liquid_vector.push_back("magmatic liquid mass mobility");
    lhs_vapor_vector.push_back("magmatic fluid mass vapor");
    rhs_vapor_vector.push_back("magmatic vapor mass mobility");
    // air: no magmatic mass

    // ── Magmatic salt advection ───────────────────────────────────────────────
    lhs_liquid_vector.push_back("magmatic salt content liquid");
    rhs_liquid_vector.push_back("magmatic liquid salt mobility");
    lhs_vapor_vector.push_back("magmatic salt content vapor");
    rhs_vapor_vector.push_back("magmatic vapor salt mobility");
    // air: no magmatic salt

    // ── FV transport constructor arguments per phase ──────────────────────────
    // [0] = velocity field (KgradP shared by all phases — phase splitting is
    //       done in UpwindControlVisitor by multiplying with relperm_visc)
    // [1] = nodal source
    // [2] = phase density (sets which upwind matrix this instance uses)

    // liquid
    fv_transport_liquid_variables.push_back("KgradP");
    fv_transport_liquid_variables.push_back("nodal source liquid");
    fv_transport_liquid_variables.push_back(thermodymamic_density_in_gravity_term
                                                ? "density liquid"
                                                : "density liquid transport");

    // vapor
    fv_transport_vapor_variables.push_back("KgradP");
    fv_transport_vapor_variables.push_back("nodal source vapor");
    fv_transport_vapor_variables.push_back(thermodymamic_density_in_gravity_term
                                               ? "density vapor"
                                               : "density vapor transport");


    // air — NOT gated: see comment above on lhs_air_vector.
    // Same applies to conductance_p_upwind_air_variables and grav_air_variables below.
    fv_transport_air_variables.push_back("KgradP");
    fv_transport_air_variables.push_back("nodal source air");
    fv_transport_air_variables.push_back(thermodymamic_density_in_gravity_term
                                             ? "density air"
                                             : "density air transport");


    // ── UpwindControlVisitor phase-count vectors — GATED on with_air ─────────
    // These vectors control how many phases UpwindControlVisitor loops over:
    //   phases = densities.size()
    // If air entries are present, phases=3 and air is computed each timestep.
    // If absent, phases=2 and fv_transport_air is never touched by the upwind visitor.
    // This is the ONLY place where with_air actually gates runtime behavior.
    // NOTE: UpwindControlVisitor::uc_key is hardcoded for indices 0,1,2 —
    // it must check (phases==3) before accessing uc_key[2]. See UpwindControlVisitor.cpp.

    upwind_control_variables.push_back("vertical permeability");
    //upwind_control_variables.push_back("permeability");
    upwind_control_variables.push_back("porosity");

    densities.push_back(thermodymamic_density_in_gravity_term ? "density liquid"           : "density liquid transport");
    densities.push_back(thermodymamic_density_in_gravity_term ? "density vapor"            : "density vapor transport");

    relperm_visc.push_back("relperm viscosity liquid");
    relperm_visc.push_back("relperm viscosity vapor");

    saturations.push_back("saturation liquid");
    saturations.push_back("saturation vapor");

    time_control.push_back("courant liquid");
    time_control.push_back("courant vapor");

    velocities.push_back("velocity liquid");
    velocities.push_back("velocity vapor");

    pore_velocities.push_back("pore velocity liquid");
    pore_velocities.push_back("pore velocity vapor");

    // ── PDE operators — heat diffusion ────────────────────────────────────────
    capacitance_lhs_variables.push_back("nodal heat capacity");
    capacitance_lhs_variables.push_back("temperature");
    capacitance_lhs_variables.push_back("temperature");

    conductance_variables.push_back("thermal conductivity");
    conductance_variables.push_back("temperature");
    conductance_variables.push_back("temperature");

    capacitance_rhs_variables.push_back("nodal heat capacity");
    capacitance_rhs_variables.push_back("temperature");

    heat_bottom_variables.push_back("nodal heat flux bottom");
    heat_bottom_variables.push_back("temperature");

    // ── PDE operators — pressure equation ────────────────────────────────────
    capacitance_lhs_p_variables.push_back("nodal total compressibility");
    capacitance_lhs_p_variables.push_back("fluid pressure");
    capacitance_lhs_p_variables.push_back("fluid pressure");

    // conductance liquid
    conductance_p_upwind_liquid_variables.push_back("permeability tensor");
    //conductance_p_upwind_liquid_variables.push_back("permeability");
    conductance_p_upwind_liquid_variables.push_back("fluid pressure");
    conductance_p_upwind_liquid_variables.push_back("fluid pressure");
    conductance_p_upwind_liquid_variables.push_back("liquid mass mobility");
    conductance_p_upwind_liquid_variables.push_back(thermodymamic_density_in_gravity_term
                                                        ? "density liquid"
                                                        : "density liquid transport");

    // conductance vapor
    conductance_p_upwind_vapor_variables.push_back("permeability tensor");
    //conductance_p_upwind_vapor_variables.push_back("permeability");
    conductance_p_upwind_vapor_variables.push_back("fluid pressure");
    conductance_p_upwind_vapor_variables.push_back("fluid pressure");
    conductance_p_upwind_vapor_variables.push_back("vapor mass mobility");
    conductance_p_upwind_vapor_variables.push_back(thermodymamic_density_in_gravity_term
                                                       ? "density vapor"
                                                       : "density vapor transport");

    // conductance air
    conductance_p_upwind_air_variables.push_back("permeability tensor");
    //conductance_p_upwind_air_variables.push_back("permeability");
    conductance_p_upwind_air_variables.push_back("fluid pressure");
    conductance_p_upwind_air_variables.push_back("fluid pressure");
    conductance_p_upwind_air_variables.push_back("air mass mobility");
    conductance_p_upwind_air_variables.push_back(thermodymamic_density_in_gravity_term
                                                     ? "density air"
                                                     : "density air transport");

    capacitance_rhs_p_variables.push_back("nodal total compressibility");
    capacitance_rhs_p_variables.push_back("fluid pressure");

    // ── Gravity terms (each uses the VERTICAL permeability component) ─────────
    // liquid
    grav_liq_variables.push_back("vertical permeability");
    //grav_liq_variables.push_back("permeability");
    grav_liq_variables.push_back("fluid pressure");
    grav_liq_variables.push_back("liquid mass mobility density");
    grav_liq_variables.push_back(thermodymamic_density_in_gravity_term
                                     ? "density liquid"
                                     : "density liquid transport");

    // vapor
    grav_vap_variables.push_back("vertical permeability");
    //grav_vap_variables.push_back("permeability");
    grav_vap_variables.push_back("fluid pressure");
    grav_vap_variables.push_back("vapor mass mobility density");
    grav_vap_variables.push_back(thermodymamic_density_in_gravity_term
                                     ? "density vapor"
                                     : "density vapor transport");

    // air
    grav_air_variables.push_back("vertical permeability");
    //grav_air_variables.push_back("permeability");
    grav_air_variables.push_back("fluid pressure");
    grav_air_variables.push_back("air mass mobility density");
    grav_air_variables.push_back(thermodymamic_density_in_gravity_term
                                     ? "density air"
                                     : "density air transport");

    // source terms
    source_p_variables.push_back("nodal fluid volume source");
    source_p_variables.push_back("fluid pressure");

    source_p2_variables.push_back("fluid source rate");
    source_p2_variables.push_back("fluid pressure");

    // ── Transient variables ───────────────────────────────────────────────────
    transient_variables.push_back("temperature");
    transient_variables.push_back("fluid pressure");
    transient_variables.push_back("salinity");
    transient_variables.push_back("fluid mass liquid");
    transient_variables.push_back("fluid mass vapor");
    transient_variables.push_back("fluid density");
    transient_variables.push_back("enthalpy content liquid");
    transient_variables.push_back("enthalpy content vapor");
    transient_variables.push_back("salt content liquid");
    transient_variables.push_back("salt content vapor");
    transient_variables.push_back("salt content fluid");

    transient_variables.push_back("magmatic fluid mass");
    transient_variables.push_back("magmatic fluid mass liquid");
    transient_variables.push_back("magmatic fluid mass vapor");
    transient_variables.push_back("magmatic mass salt");
    transient_variables.push_back("magmatic salt content liquid");
    transient_variables.push_back("magmatic salt content vapor");
    transient_variables.push_back("magmatic salt content halite");
    transient_variables.push_back("salt flux integral");
    transient_variables.push_back("energy flux integral");
    transient_variables.push_back("fluid flux integral");
    transient_variables.push_back("magmatic fluid flux integral");
    transient_variables.push_back("shell magmatic fluid flux integral");
    transient_variables.push_back("nodal fluid volume source");
    transient_variables.push_back("nodal porosity");

    // transient_variables.push_back("liquidus temperature");
    // transient_variables.push_back("solidus temperature");

    transient_variables.push_back("nodal density rock");

    // air transient variables — NOT gated: these are zero-initialized at startup
    // (StoreInitialPropertiesAndFlags writes zeros unconditionally) so including
    // them in the reset list is harmless.
    transient_variables.push_back("fluid mass air");
    transient_variables.push_back("enthalpy content air");

    if(with_zinc)
    {
        transient_variables.push_back("zinc solid");
        transient_variables.push_back("zinc t1 solid");
        transient_variables.push_back("zinc t2 solid");
        transient_variables.push_back("zinc t3 solid");
        transient_variables.push_back("zinc t4 solid");
        transient_variables.push_back("total zinc mass out");
        transient_variables.push_back("total zinc mass out precip");
        transient_variables.push_back("total zinc t1 mass out");
        transient_variables.push_back("total zinc t1 mass out precip");
        transient_variables.push_back("total zinc t2 mass out");
        transient_variables.push_back("total zinc t2 mass out precip");
        transient_variables.push_back("total zinc t3 mass out");
        transient_variables.push_back("total zinc t3 mass out precip");
        transient_variables.push_back("total zinc t4 mass out");
        transient_variables.push_back("total zinc t4 mass out precip");
        transient_variables.push_back("zinc content liquid");//no need to add here
        transient_variables.push_back("zinc content vapor");//no need to add here
        transient_variables.push_back("zinc content fluid");//no need to add here
    }

    if(with_lithium)
    {
        transient_variables.push_back("lithium solid");
        transient_variables.push_back("lithium content liquid");//no need to add here
        transient_variables.push_back("lithium content vapor");//no need to add here
        transient_variables.push_back("lithium content fluid");//no need to add here
        transient_variables.push_back("tracer content liquid");//no need to add here
        transient_variables.push_back("tracer content vapor");//no need to add here
        transient_variables.push_back("tracer content fluid");//no need to add here
    }

    // ── Full reset variables ──────────────────────────────────────────────────
    // Benoit: these variables would probably need to be reset if reset comes after thermal equilibration
    // which is NOT the case at the moment
    // Other variables may need to be reset: previous total enthalpy, previous mass salt and previous fluid density
    //     ( see thermal equilibrator: PrepareVariablesForStorage())

    full_reset_variables.push_back("liquid mass mobility");
    full_reset_variables.push_back("vapor mass mobility");
    full_reset_variables.push_back("liquid enthalpy mobility");
    full_reset_variables.push_back("vapor enthalpy mobility");
    full_reset_variables.push_back("liquid salt mobility");
    full_reset_variables.push_back("vapor salt mobility");
    full_reset_variables.push_back("density liquid");
    full_reset_variables.push_back("density vapor");
    full_reset_variables.push_back("density liquid transport");
    full_reset_variables.push_back("density vapor transport");
    full_reset_variables.push_back("relperm viscosity liquid");
    full_reset_variables.push_back("relperm viscosity vapor");
    full_reset_variables.push_back("liquid mass mobility density");
    full_reset_variables.push_back("vapor mass mobility density");
    full_reset_variables.push_back("enthalpy liquid");
    full_reset_variables.push_back("enthalpy vapor");
    full_reset_variables.push_back("nodal heat capacity");
    full_reset_variables.push_back("nodal total compressibility");
    full_reset_variables.push_back("fluid source rate");
    full_reset_variables.push_back("fluid source h");
    full_reset_variables.push_back("fluid source wt");
    full_reset_variables.push_back("previous mass salt");
    full_reset_variables.push_back("previous fluid density");
    full_reset_variables.push_back("magmatic liquid mass mobility");
    full_reset_variables.push_back("magmatic vapor mass mobility");
    full_reset_variables.push_back("magmatic liquid salt mobility");
    full_reset_variables.push_back("magmatic vapor salt mobility");
    full_reset_variables.push_back("viscosity liquid");
    full_reset_variables.push_back("viscosity vapor");
    full_reset_variables.push_back("saturation liquid");
    full_reset_variables.push_back("saturation vapor");
    full_reset_variables.push_back("saturation halite");

    // air full-reset variables — NOT gated: same rationale as transient_variables above.
    // Resetting zero values to zero has low cost and no side effect.
    full_reset_variables.push_back("air mass mobility");
    full_reset_variables.push_back("air enthalpy mobility");
    full_reset_variables.push_back("density air");
    full_reset_variables.push_back("density air transport");
    full_reset_variables.push_back("relperm viscosity air");
    full_reset_variables.push_back("air mass mobility density");
    full_reset_variables.push_back("enthalpy air");
    full_reset_variables.push_back("viscosity air");
    full_reset_variables.push_back("saturation air");

    if(with_zinc)
    {
        full_reset_variables.push_back("liquid zinc mobility");
        full_reset_variables.push_back("vapor zinc mobility");
    }

    if(with_lithium)
    {
        full_reset_variables.push_back("liquid lithium mobility");
        full_reset_variables.push_back("vapor lithium mobility");
        full_reset_variables.push_back("liquid tracer mobility");
        full_reset_variables.push_back("vapor tracer mobility");
    }

    cerr << "\n=== CVFEM_PHX_VariableNames ===";
    cerr << "\n  with_lithium: " << (with_lithium ? "ON" : "OFF");
    cerr << "\n  with_zinc:    " << (with_zinc    ? "ON" : "OFF");
    cerr << "\n===============================\n";
};

inline void CVFEM_PHX_VariableNames::AddAdvectionVariable( std::string balanced_variable,
                                                          std::string new_lhs_liquid, std::string new_rhs_liquid,
                                                          std::string new_lhs_vapor,  std::string new_rhs_vapor  )
{
    // Adds a tracer species that travels with the H2O-NaCl phases only.
    // Air does not carry salt, magmatic mass, lithium, or tracer species,
    // so no air equivalent is registered here.
    transient_variables.push_back(balanced_variable);
    transient_variables.push_back(new_lhs_liquid);
    transient_variables.push_back(new_lhs_vapor);

    lhs_liquid_vector.push_back(new_lhs_liquid);
    rhs_liquid_vector.push_back(new_rhs_liquid);
    lhs_vapor_vector.push_back(new_lhs_vapor);
    rhs_vapor_vector.push_back(new_rhs_vapor);
};

inline void CVFEM_PHX_VariableNames::ActivateAir(bool thermodymamic_density_in_gravity_term)
{
    densities.push_back(thermodymamic_density_in_gravity_term
                            ? "density air" : "density air transport");
    relperm_visc.push_back("relperm viscosity air");
    saturations.push_back("saturation air");
    time_control.push_back("courant air");
    velocities.push_back("velocity air");
    pore_velocities.push_back("pore velocity air");
}

} // end namespace csmp

#endif