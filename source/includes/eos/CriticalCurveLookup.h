#ifndef  CRITICALCURVELOOKUP_H
#define  CRITICALCURVELOOKUP_H

#include "CSMP_definitions.h"
#include "CriticalPointH2O.h"
#include "ErrorHandler.h"

#include <cmath>

/* Changelog
   2007 and before, Thomas Driesner: original CSP version
   May 13, 2008, Thomas Driesner: add comments, clean up interface
   May  6, 2009, Thomas Driesner: add out-of-range handling
   March 23, 2010, Thomas Driesner: removed an inconsistency in the pressure calculation ( the previous version didn't use the published formulas )
   unknown date, pre-2014: initial port to CSMP++
   February 10. 2015 : continued port to csmp++, adaptation of pressure units to Pa
*/
namespace csmp
{
  /// Creation of lookup tables temperature-pressure-compsition relations for the critical curve of the system H2O-NaCl, based on class CriticalCurve
  class CriticalCurveLookup
  {

  public:
    CriticalCurveLookup(const double& externaltemperature);
    ~CriticalCurveLookup();

    double                Temperature();
    double                Pressure();
    double                MassFractionNaCl();
    double                TfromP(const double& press); 

    // the following ones may become obsolete in the course of code clean-up    
    //    double                MoleFractionNaCl();
    double                Density();
    double                Enthalpy();
    double                HeatCapacity();
    double                Compressibility();
    double                Viscosity();
    double                ValueOf(const int& property_index);

  private:

    const double&         temperature;

    double                tcurrent;        ///< temperature [C] for internal use
    double                pcurrent;        ///< pressure [Pa] for internal use
    double                t_res; 
    double                tnorm; 
    double                pnorm;

    long                    it; 
    long                    it_min; 
    long                    it_max; 
    long                    t_dim; 
    long                    i_guess; 
    long                    i_max; 
    long                    i_min;

    std::vector<double>   storage_vector;

    CriticalPointH2O        cp_h2o;
    void                    GetTemperatureIndex(const double& t);
    void                    GetTemperatureIndexCriticalCurve(const double& t);
    ErrorHandler&           csmp_error;
  };

  /**
     @class CriticalCurveLookup CriticalCurveLookup.h "eos/h2o_nacl/CriticalCurveLookup.h"

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact 
     thomas.driesner@erdw.ethz.ch

     @section motivation Motivation
     CriticalCurveLookup generates and queries a lookup table that stores the temperature-pressure-composition [Celsius, bar, mole fraction NaCl] coordinates of the critical curve of the system H2O-NaCl, according to the paper

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.

     In addition, mass density [kg/m3], specific enthalpy [J/kg], specific isobaric heat capacity [J/kg/C], isothermal compressibility [Pa^-1], dynamic viscosity [Pa s] on the critical curve are available, and are computed by applying the formulations from the paper

     Driesner T. (2007): The system H2O-NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000oC, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4902-4919.

     at the temperature-pressure-composition coordinates of the critical curve.

     The lookup table "CriticalCurveLookupTable.bin" is written as a binary file, and stores the properties listed above in a single, 1-dimensional std::vector<double>, named storage_vector. storage_vector comprises subsequent blocks in each of which the values for the one property are stored in sequence of ascending temperature from the critical temperature of water (see file CriticalPointH2O.cpp) to 1000.0e0 degrees Celsius. Each block is of size t_dim (the number of entries, the first thing that is computed in the constructor), the sequence is defined in the file "LookupPropertyIndex.h". The temperature spacing is variable (see member GetTemperatureIndex for details), being densest near the critical temperature of water and largest at the highest temperatures.
   
     @section testing Testing
     =======
     The lookup version was tested against the original non-lookup code (May 26, 2008). Max deviations are:

     - Composition: -0.5% at 374.03 C, above 374.3 C, deviations are <0.1% and decay rapidly to negligible values.
     - Pressure: deviations are negligible over the whole temperature range (.00x % or less)
     - Density: (nota bene: this is NOT used in fluid flow simulations in the original ETH scheme as programmed by Thomas Driesner et al.) up to 10% at 374.04, rapidly decaying via 3% at 374.6 to less than 1% above 375.5. Notice that the uncertainties are always maximum between lookup grid point and are "exact" at the grid points. (Oct 2010, Thomas Driesner: as far as I can see, these deviations have been reduced after some debugging; fruther testing required)
     - Enthalpy: (again: not used in fluid flow simulations). 2% at 374.04, rapidly decaying, < 0.05% at t > 380 C
     - HeatCapacity: (again: not used in fluid flow simulations): large deviations at t < 375.5, which are solely due to the "hard-wired" value assigned at tcrit_h20 (1e10, this happens in H2OLookup). From 424 C: < 0.1%
     - Compressibility: (again: not used ...) Near tcrit_h2o suffers from the same problem as HeatCapacity, <1% is reached at about 385 C.
     - Viscosity: (again ...) 7% at 374.04, < 1% at t> 374.5
  */
}// namespace csmp
#endif
