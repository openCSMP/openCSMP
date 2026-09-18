// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef   H2ONACLTHERMALEQUILIBRATOR_H
#define   H2ONACLTHERMALEQUILIBRATOR_H

#include "ErrorHandler.h"

#include "H2ONaClFluidProperties.h"
#include "Fluidproperties.h"
#include "H2OThermalEquilibrator.h"

#include "Rock.h"
#include "Air.h"


namespace csmp
{
class H2ONaClThermalEquilibrator
{

    //***TD*** cross-check the actual porosity meaning such that any implicit scaling with volume or so is clearly spelled out in documentation

public:
    H2ONaClThermalEquilibrator(const double& external_mass_rock,  // [kg]
                               const double& external_cp_rock,    // [J/kg/K]
                               const double& external_rho_rock,   // [kg/m^3]

                               const double& external_mass_air,  // [kg] NEW
                               const double& external_cp_air,    // [J/kg/K] NEW

                               const double& external_porosity,   // [ ]
                               const double& external_mass_fluid, // [kg]
                               const double& external_wt_current, // [wt% NaCl]
                               const double& external_t_previous, // [C]
                               const double& external_p_current,  // [bar]
                               const double& external_H_current,  // [J] NOT [J/kg]
                               const double& external_H_previous, // [J]
                               const bool&     fixed_external_t,
                               const double& t_fix,
                               const double& external_t_diff,
                               const bool& verbose);
    ~H2ONaClThermalEquilibrator();
    
    Fluidproperties            Liquid();
    Fluidproperties            Vapor();
    Fluidproperties            Bulk();
    Fluidproperties            Salt();

    bool                       ThreePhaseProperties(Fluidproperties& bulk_external,
                                                    Fluidproperties& liquid_external,
                                                    Fluidproperties& vapor_external,
                                                    Fluidproperties& salt_external);

    Fluidproperties            ReportLiquidProperties(const double& t,
                                                      const double& p,
                                                      const double& x,
                                                      const double& h);
    Fluidproperties            ReportVaporProperties( const double& t,
                                                      const double& p,
                                                      const double& x,
                                                      const double& h);
    Fluidproperties            ReportBulkProperties(  const double& t,
                                                      const double& p,
                                                      const double& x,
                                                      const double& h);
    Fluidproperties            ReportSaltProperties(  const double& t,
                                                      const double& p,
                                                      const double& x,
                                                      const double& h);

    double                   Resid();
    int                        n_iterations();
    bool                       Equilibrated();
    bool                       Fatal();

    // PW May 2016
    void                       SetConvergenceSpeedUpTo(bool boost);
    //BenoitLC add
    void                       WithRockLiquidusSolidus(bool with_rock_liquidus_solidus_);
    void                       WithAirPhase(bool enable_air_phase_);
    void                       AttemptToSurviveFluidPropertiesError(bool attempt_to_survive_fluid_properties_error_);
    void                       SetRockLiquidusSolidusTemperatures(double tl_, double ts_);
    void                       SetRockHeatCapacity( double mini_cp );
    void                       SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                                           double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_);

    double SaturationPressureFromT(double t_C);

protected:

    const double&            mass_rock;
    const double&            cp_rock;
    const double&            rho_rock;

    const double&            mass_air;
    const double&            cp_air;

    const double&            phi;
    const double&            mass_fluid;
    const double&            wt_current;
    const double&            t_previous;
    const double&            p_current;
    const double&            H_current;
    const double&            H_previous;
    const double&            t_fix;
    const double&            t_diff;
    const bool&              fixed_t;
    const bool&              verbose;
    //    const bool& lognode;
    bool                       equilibrated;
    bool                       fatal;

    // PW May 2016
    bool                       convergence_speed_up;
    //BenoitLC add
    bool with_rock_liquidus_solidus;
    double tl;
    double ts;
    bool enable_air_phase;

    double                 temp_correction;
    bool                   properties_check_fail;
    bool                   attempt_to_survive_fluid_properties_error;
    //---

    const double             convergence_criterion; // fraction of fluid enthalpy
    const double             minimum_dt;

    double                   t_eq;
    double                   x_current_eq;
    double                   h_fluid_eq;
    double                   h_fluid_test;

    double                   h_air_eq;
    double                   h_air_test;

    double                   tmin;
    double                   tmax;
    double                   resid;
    double                   Hmax;
    double                   Hmin;
    double                   mass_rock_eq;
    double                   cp_rock_eq;

    double                   mass_air_eq;
    double                   cp_air_eq;

    double                   mass_fluid_eq;
    double                   wt_current_eq;
    double                   t_previous_eq;
    double                   p_current_eq;
    double                   H_current_eq;
    double                   H_previous_eq;
    double                   H_test;
    double                   tdummy;
    double                   hdummy;
    double                   t_previous_ini;
    // new debug stuff
    double                   mass_rock_ini;
    double                   cp_rock_ini;
    double                   rho_rock_ini;

    double                   mass_air_ini;
    double                   cp_air_ini;

    double                   porosity_ini;
    double                   mass_fluid_ini;
    double                   wt_current_ini;
    double                   x_current_ini;
    double                   p_current_ini;
    double                   H_current_ini;
    double                   H_previous_ini;
    // end new debug stuff

    int                        icrit;
    int                        icrit_max;
    //id long myid;

    //    ofstream logfile,Hfile;

    Rock                       rock;
    Air                        air;
    H2ONaClFluidProperties     fluid;
    H2OThermalEquilibrator     water_equilibrator;

    double                   ComputeTotalEnthalpyAtTemperature( const double& t );
    void                       ErrorCheckHmin(const int& i);
    void                       ErrorCheckHmax(const int& i);
    void                       FluidPropertiesErrorCheck();
    void                       PrintStatusToCerr();
    void                       ErrorConditionsToScreen();
    void                       FindInitialValues();

    // PW May 2016
    void                       InitialValues();

    Fluidproperties            liquidprops;
    Fluidproperties            vaporprops;
    Fluidproperties            bulkprops;
    Fluidproperties            saltprops;

    void                       Equilibrate();
    void                       UpdateEquilibratedFluidproperties();
    void                       AssignAllPropertiesViaReport();
    void                       AssignAllPropertiesNormally();

    ErrorHandler&              csmp_error;
};


/**
     @class H2ONaClThermalEquilibrator.h "eos/h2o_nacl/H2ONaClThermalEquilibrator.h" (preliminary)

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch
     
     @section motivation Motivation
     H2ONaClThermalEquilibrator performs an isobaric thermal equilibration between a mass of rock and a mass of fluid (from the system H2O-NaCl) given a total enthalpy for this system.
     
     @section assumptions Assumptions
     Currently, we assume that we can consistently treat treat the rock's enthalpy as h_rock = cp_rock*t . As only differences are used, this principally wrong formula [rather, it'd be h0+integral(cp_rock(t) dt)], should work here.  It does, however, fail in the case of a t-dependent heat capacity.

     @section usage Usage
     Construct an instance of H2ONaClThermalEquilibrator with references to the following properties (units are given in square brackets) as constructor arguments:
     - const double& external_mass_rock  [kg]    : the mass of rock that shall be equilibrated with the fluid
     - const double& external_cp_rock    [J/kg/K]: the rock's isobaric heat capacity (assumed to be independent of temperature)
     - const double& external_rho_rock   [kg/m^3]: the rock's mass density
     - const double& external_porosity   [ ]     : (currently not used) the rock's porosity (this is allowed to be different from the actual fluid volume);
     - const double& external_mass_fluid [kg]    : the mass of "fluid" - this is kg H2O + kg NaCl (solid+dissolved); including the solid is essential here as salt may be dissolved during equilibration
     - const double& external_wt_current [ ]     : 100 * [kg NaCl (solid+dissolved)] / [kg NaCl (solid+dissolved) + kg H2O]
     - const double& external_t_previous [C]     : temperature at thermal equilibrium from the previous time step
     - const double& external_p_current  [bar]   : preesure ar which thermal equilibration will be performed
     - const double& external_H_current  [J]     : total enthalpy in system (total enthalpy of rock + total enthalpy of fluid); this is the result of a transport step, including advection and thermal diffusion
     - const double& external_H_previous [J]     : the above at  the previous time step
     - const bool&      fixed_external_t    [ ]     : a bool, mostly for use at fixed temperature boundary condition; if true, no thermal equilibration will be performed, the fluid properties and volume at p_current and t_fix (see next entry) will be calculated
     - const double& t_fix               [C]     : the fixed temperature to be used if fixed_external_t == true
     - const double& external_t_diff     [C]     : rock temperature after the thermal diffusion step and prior to equilibration. By default, if thermal equilibration fails, this will be assumed to be the best proxy to an equilibration result

     H2ONaClThermalEquilibrator is typically used by NaClH2OPropertiesVisitor.

     @code
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues

     @section todo To do list
     - Constitutive relation for rock enthalpies
     - make t_diff as default in non-equilibration case an option
     
     @section testing Testing
  */


}// csmp
#endif
