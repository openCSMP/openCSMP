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

using namespace std;


namespace csmp
{
  class H2OThermalEquilibrator
  {
  public: 
    H2OThermalEquilibrator(const double64& external_mass_rock,  // [kg]
                           const double64& external_cp_rock,    // [J/kg/K]
                           const double64& external_rho_rock,   // [kg/m^3]
                           const double64& external_porosity,   // [dimensionless]
                           const double64& external_mass_fluid, // [kg]
                           const double64& external_t_previous, // [C]
                           const double64& external_p_current,  // [Pa]
                           const double64& external_H_current,  // [J] NOT [J/kg]
                           const double64& external_H_previous, // [J]
                           const bool&      fixed_external_t,
                           const double64& t_fix,
                           const double64& external_t_diff);

    ~H2OThermalEquilibrator();
    
    Fluidproperties        Liquid();
    Fluidproperties        Vapor();
    Fluidproperties        Bulk();

    Fluidproperties        ReportLiquidProperties( const double64& t, const double64& p, const double64& h );
    Fluidproperties        ReportVaporProperties(  const double64& t, const double64& p, const double64& h );
    Fluidproperties        ReportBulkProperties(   const double64& t, const double64& p, const double64& h );

    int                    n_iterations();
    int                    EqType();
    bool                   Equilibrated();
    bool                   Fatal();
    
    void                   TemperatureDependentHeatCapacityRock( double64 cpr_min_ext, 
                                                                 double64 t_min_ext,
                                                                 double64 cpr_max_ext, 
                                                                 double64 t_max_ext );

    double64               TemperatureDependent_cpr( double64 T );
    
  private:

    const double64&        mass_rock;
    const double64&        cp_rock;
    const double64&        rho_rock;
    const double64&        phi;
    const double64&        mass_fluid;
    const double64&        t_previous;
    const double64&        p_current;
    const double64&        H_current;
    const double64&        H_previous;
    const double64&        t_fix;
    const double64&        t_diff;
    const bool&            fixed_t;
				   
    bool                   equilibrated;
    bool                   fatal;

    const double64         convergence_criterion; // fraction of fluid enthalpy

    double64               t_eq; 
    double64               h_fluid_eq;
    double64               tmin;
    double64               tmax;
    double64               resid; 
    double64               Hmax; 
    double64               Hmin;

    int                    icrit;
    int                    icrit_max;

    double64               mass_rock_eq;
    double64               cp_rock_eq;
    double64               mass_fluid_eq;
    double64               t_previous_eq;
    double64               p_current_eq;
    double64               H_current_eq;
    double64               H_previous_eq;
    double64               H_test;
    double64               tdummy;
    double64               hdummy;

    double64               cpr_min; 
    double64               cpr_max; 
    double64               t_min; 
    double64               t_max; 
    double64               cpr_t_dep;

    bool                   t_dependent_cpr;

    Rock                   rock;
    H2OFluidProperties     fluid;

    Fluidproperties        liquidprops;
    Fluidproperties        vaporprops;
    Fluidproperties        bulkprops;

    double64               ComputeTotalEnthalpyAtTemperature( const double64& t );
    void                   ErrorCheckHmin(const int& i);
    void                   ErrorCheckHmax(const int& i);

    void                   Equilibrate();
    void                   FindInitialValues();
    void                   AssignAllPropertiesViaReport();
    void                   CheckForTminTmaxOutOfRange();
    void                   FluidPropertiesErrorCheck(); 
    void                   PrintStatusToCerr();

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

     Another set of calls - ReportLiquidProperties(const double64& t, const double64& p, const double64& h), ReportVaporProperties(const double64& t, const double64& p, const double64& h), ReportBulkProperties(const double64& t, const double64& p, const double64& h) - lets you query the properties for a desired combination of temperature, pressure and specific fluid enthalpy. Use with care because getting the logic for this right in a simulation context is non-trivial.  
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
