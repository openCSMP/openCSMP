// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef H2OTHERMALEQUILIBRATOR_H
#define H2OTHERMALEQUILIBRATOR_H

#include<cmath>
#include<iostream>
#include<string>
#include<sstream>

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

#include "H2OFluidProperties.h"
#include "Fluidproperties.h"
#include "States.h"

//new June 14, 2011
#include "Rock.h"

#include "Air.h"

using namespace std;


namespace csmp
{
  class H2OThermalEquilibrator
  {
  public:
    H2OThermalEquilibrator(const double& external_mass_rock,  // [kg]
                           const double& external_cp_rock,    // [J/kg/K]
                           const double& external_rho_rock,   // [kg/m^3]

                             const double& external_mass_air,  // [kg] NEW
                             const double& external_cp_air,    // [J/kg/K] NEW

                           const double& external_porosity,   // [dimensionless]
                           const double& external_mass_fluid, // [kg]
                           const double& external_t_previous, // [C]
                           const double& external_p_current,  // [Pa]
                           const double& external_H_current,  // [J] NOT [J/kg]
                           const double& external_H_previous, // [J]
                           const bool&      fixed_external_t,
                           const double& t_fix,
                           const double& external_t_diff);

    ~H2OThermalEquilibrator();
    
    Fluidproperties        Liquid();
    Fluidproperties        Vapor();
    Fluidproperties        Bulk();

    Fluidproperties        ReportLiquidProperties( const double& t, const double& p, const double& h );
    Fluidproperties        ReportVaporProperties(  const double& t, const double& p, const double& h );
    Fluidproperties        ReportBulkProperties(   const double& t, const double& p, const double& h );

    int                    n_iterations();
    int                    EqType();
    bool                   Equilibrated();
    bool                   Fatal();
    

	// PW May 2016 - two functions added for convergence speed-up
    bool                   TwoPhaseProperties(Fluidproperties& bulk_external,
		                                      Fluidproperties& liquid_external,
		                                      Fluidproperties& vapor_external);
	void                   SetConvergenceSpeedUpTo(bool boost);

    //BenoitLC add
    void                   WithRockLiquidusSolidus(bool with_rock_liquidus_solidus_);
    void                   WithAirPhase(bool enable_air_phase_);
    void                   AttemptToSurviveFluidPropertiesError(bool attempt_to_survive_fluid_properties_error_);
    void                   SetRockLiquidusSolidusTemperatures(double tl_, double ts_);
    void                   SetRockHeatCapacity( double mini_cp );
    void                   SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                                       double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_);
    //---

  protected:

    const double&        mass_rock;
    const double&        cp_rock;
    const double&        rho_rock;

    const double&            mass_air;
    const double&            cp_air;

    const double&        phi;
    const double&        mass_fluid;
    const double&        t_previous;
    const double&        p_current;
    const double&        H_current;
    const double&        H_previous;
    const double&        t_fix;
    const double&        t_diff;
    const bool&            fixed_t;
				   
    bool                   equilibrated;
    bool                   fatal;

	// PW May 2016 - boolian for convergence speed-up
	bool                  convergence_speed_up;

    //BenoitLC add
    bool with_rock_liquidus_solidus;
    double tl;
    double ts;
    bool enable_air_phase;


    bool                   attempt_to_survive_fluid_properties_error;
    bool                   properties_check_fail;
    double                 temp_correction;

    const double         convergence_criterion; // fraction of fluid enthalpy

    //---

    double               t_eq;
    double               h_fluid_eq;

    double               h_air_eq;

    double               tmin;
    double               tmax;
    double               resid;
    double               Hmax;
    double               Hmin;

    int                    icrit;
    int                    icrit_max;

    double               mass_rock_eq;
    double               cp_rock_eq;

    double                   mass_air_eq;
    double                   cp_air_eq;

    double               mass_fluid_eq;
    double               t_previous_eq;
    double               p_current_eq;
    double               H_current_eq;
    double               H_previous_eq;
    double               H_test;
    double               tdummy;
    double               hdummy;

    double               cpr_min;
    double               cpr_max;
    double               t_min;
    double               t_max;
    double               cpr_t_dep;

    bool                   t_dependent_cpr;

    Rock                   rock;
    Air                    air;
    H2OFluidProperties     fluid;

    Fluidproperties        liquidprops;
    Fluidproperties        vaporprops;
    Fluidproperties        bulkprops;

    double               ComputeTotalEnthalpyAtTemperature( const double& t );
    void                   ErrorCheckHmin(const int& i);
    void                   ErrorCheckHmax(const int& i);

    void                   Equilibrate();
    void                   FindInitialValues();
    void                   AssignAllPropertiesViaReport();
    void                   CheckForTminTmaxOutOfRange();
    void                   FluidPropertiesErrorCheck(); 

	// PW May 2016 - function for convergence speed-up
    void                   InitialValues();


    ErrorHandler&          csmp_error;
  };


  /**
     @class H2OThermalEquilibrator H2OThermalEquilibrator.h 

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch
     
     @section motivation Motivation
     H2OThermalEquilibrator finds the temperature at which rock and a pure H2O fluid are in thermal equilibrium at a given pressure. This is done by an isobaric bisection iteration through the water phase diagram, which compares the total enthalpy at a test temperature t_eq (H_test = mass_rock*enthalpy_rock(t_eq) + mass_fluid*enthalpy_fluid(t_eq) ) with the enthalpy value H_current_eq that was provided by the transport code (specifically: H2OPropertiesVisitor). 
 
     Currently (as of January 28, 2011) iteration starts from an upper temperature limt tmax at 1000C and a lower one at 5C and follows true bisection. The actual computation of fluid properties as well as determination of fluid phase state is done by an H2OFluidProperties object. 
     
     @section usage Usage
     Instantiate H2OThermalEquilibrator with references to various properties seen as constructor arguments. It is YOUR responsibility to make sure they are properly up to date when H2OThermalEquilibrator ist called.

     The main calls are Liquid(), Vapor(), and Bulk(). Calling one of them triggers a full equilibration (except if the values of those properties referenced in the constructor havn't changed, then the last equilibration results are assumed to be correct). It returns fluid property values as Fluidproperties objects (see documentation there for full property set). These fluidproperties refer to the phase indicate by the name of the function call.

     Another set of calls - ReportLiquidProperties(const double& t, const double& p, const double& h), ReportVaporProperties(const double& t, const double& p, const double& h), ReportBulkProperties(const double& t, const double& p, const double& h) - lets you query the properties for a desired combination of temperature, pressure and specific fluid enthalpy. Use with care because getting the logic for this right in a simulation context is non-trivial.
     @section dependencies Dependencies
     requires "H2OFluidProperties.h", "Fluidproperties.h", "States.h"
     
     @section issues Known issues 
     Pressure() still returns values in bars, a legacy from the SoWat development.
     Composition() is in mole fraction NaCl, to convert to other units, the functions in "ConvertConcentrationUnitsNaCl.h" may be useful
     
     @section testing Testing
     Testing is in progress. Besides general functionality, the following aspects need to be tested:

     - Convergence is assumed once the difference (H_test-H_current_eq) divided by the fluid mass is less than 10^-4 of the fluid's specific enthalpy. Wether a tighter criterion is necessary (and can successfully be fulfilled by this method) remains subject to further testing.

     - Care has been taken that the H2OFluidProperties object catches problems arising from numerical precision when being close to phase boundaries (e.g. when t_eq is within numerical precision identical to the boiling temperature and a vapor state is determined although the actual state is liquid). However, this has not been tested in great detail. Further details will be added to this section with time. 

     - The current bisection scheme is likely to be small. If a Newton-like component is added, the treatment of rock heat capacity / enthalpy and the prediction of tmin / tmax from t_eq may become non-trivial
  */



}// csmp
#endif
