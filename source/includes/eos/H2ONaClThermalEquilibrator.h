#ifndef   H2ONACLTHERMALEQUILIBRATOR_H
#define   H2ONACLTHERMALEQUILIBRATOR_H

#include "ErrorHandler.h"

#include "H2ONaClFluidProperties.h"
#include "Fluidproperties.h"
#include "H2OThermalEquilibrator.h"

//new June 14, 2011
// needs improvement/generalization
#include "Rock.h"


namespace csmp
{
  class H2ONaClThermalEquilibrator
  {

    //***TD*** cross-check the actual porosity meaning such that any implicit scaling with volume or so is clearly spelled out in documentation

  public: 
    H2ONaClThermalEquilibrator(const double64& external_mass_rock,  // [kg]
                               const double64& external_cp_rock,    // [J/kg/K]
                               const double64& external_rho_rock,   // [kg/m^3]
                               const double64& external_porosity,   // [ ]
                               const double64& external_mass_fluid, // [kg]
                               const double64& external_wt_current, // [wt% NaCl]
                               const double64& external_t_previous, // [C]
                               const double64& external_p_current,  // [bar]
                               const double64& external_H_current,  // [J] NOT [J/kg]
                               const double64& external_H_previous, // [J]
                               const bool&     fixed_external_t,
                               const double64& t_fix,
                               const double64& external_t_diff, 
                               const bool& verbose);
    ~H2ONaClThermalEquilibrator();
    
    Fluidproperties            Liquid();
    Fluidproperties            Vapor();
    Fluidproperties            Bulk();
    Fluidproperties            Salt();

    void                       ThreePhaseProperties(Fluidproperties& bulk_external,
                                                    Fluidproperties& liquid_external,
                                                    Fluidproperties& vapor_external,
                                                    Fluidproperties& salt_external);

    Fluidproperties            ReportLiquidProperties(const double64& t, 
                                                      const double64& p, 
                                                      const double64& x, 
                                                      const double64& h);
    Fluidproperties            ReportVaporProperties( const double64& t, 
                                                      const double64& p, 
                                                      const double64& x, 
                                                      const double64& h);
    Fluidproperties            ReportBulkProperties(  const double64& t, 
                                                      const double64& p, 
                                                      const double64& x, 
                                                      const double64& h);
    Fluidproperties            ReportSaltProperties(  const double64& t, 
                                                      const double64& p, 
                                                      const double64& x, 
                                                      const double64& h);

    double64                   Resid();
    int                        n_iterations();
    bool                       Equilibrated();
    bool                       Fatal();
    
  private:

    const double64&            mass_rock;
    const double64&            cp_rock;
    const double64&            rho_rock;
    const double64&            phi;
    const double64&            mass_fluid;
    const double64&            wt_current;
    const double64&            t_previous;
    const double64&            p_current;
    const double64&            H_current;
    const double64&            H_previous;
    const double64&            t_fix;
    const double64&            t_diff;
    const bool&                fixed_t;
    const bool&                verbose;
    //    const bool& lognode;
    bool                       equilibrated;
    bool                       fatal;

    const double64             convergence_criterion; // fraction of fluid enthalpy
    const double64             minimum_dt;

    double64                   t_eq;
    double64                   x_current_eq;
    double64                   h_fluid_eq;
    double64                   h_fluid_test;
    double64                   tmin;
    double64                   tmax;
    double64                   resid;
    double64                   Hmax; 
    double64                   Hmin;
    double64                   mass_rock_eq;
    double64                   cp_rock_eq;
    double64                   mass_fluid_eq;
    double64                   wt_current_eq;
    double64                   t_previous_eq;
    double64                   p_current_eq;
    double64                   H_current_eq;
    double64                   H_previous_eq;
    double64                   H_test;
    double64                   tdummy;
    double64                   hdummy;
    double64                   t_previous_ini;
    // new debug stuff
    double64                   mass_rock_ini;
    double64                   cp_rock_ini;
    double64                   rho_rock_ini;
    double64                   porosity_ini;
    double64                   mass_fluid_ini;
    double64                   wt_current_ini;
    double64                   x_current_ini;
    double64                   p_current_ini;
    double64                   H_current_ini;
    double64                   H_previous_ini;
    // end new debug stuff
				 
    int                        icrit;
    int                        icrit_max;
    //id long myid;

    //    ofstream logfile,Hfile;

    Rock                       rock;
    H2ONaClFluidProperties     fluid;
    H2OThermalEquilibrator     water_equilibrator;

    double64                   ComputeTotalEnthalpyAtTemperature( const double64& t );
    void                       ErrorCheckHmin(const int& i);
    void                       ErrorCheckHmax(const int& i);
    void                       FluidPropertiesErrorCheck(); 
    void                       PrintStatusToCerr();
    void                       ErrorConditionsToScreen();
    void                       FindInitialValues();

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
     - const double64& external_mass_rock  [kg]    : the mass of rock that shall be equilibrated with the fluid
     - const double64& external_cp_rock    [J/kg/K]: the rock's isobaric heat capacity (assumed to be independent of temperature)
     - const double64& external_rho_rock   [kg/m^3]: the rock's mass density
     - const double64& external_porosity   [ ]     : (currently not used) the rock's porosity (this is allowed to be different from the actual fluid volume); 
     - const double64& external_mass_fluid [kg]    : the mass of "fluid" - this is kg H2O + kg NaCl (solid+dissolved); including the solid is essential here as salt may be dissolved during equilibration
     - const double64& external_wt_current [ ]     : 100 * [kg NaCl (solid+dissolved)] / [kg NaCl (solid+dissolved) + kg H2O]
     - const double64& external_t_previous [C]     : temperature at thermal equilibrium from the previous time step
     - const double64& external_p_current  [bar]   : preesure ar which thermal equilibration will be performed
     - const double64& external_H_current  [J]     : total enthalpy in system (total enthalpy of rock + total enthalpy of fluid); this is the result of a transport step, including advection and thermal diffusion
     - const double64& external_H_previous [J]     : the above at  the previous time step
     - const bool&      fixed_external_t    [ ]     : a bool, mostly for use at fixed temperature boundary condition; if true, no thermal equilibration will be performed, the fluid properties and volume at p_current and t_fix (see next entry) will be calculated
     - const double64& t_fix               [C]     : the fixed temperature to be used if fixed_external_t == true
     - const double64& external_t_diff     [C]     : rock temperature after the thermal diffusion step and prior to equilibration. By default, if thermal equilibration fails, this will be assumed to be the best proxy to an equilibration result

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
