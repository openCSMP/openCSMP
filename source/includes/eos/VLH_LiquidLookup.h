#ifndef   VLH_LIQUIDLOOKUP_H_
#define   VLH_LIQUIDLOOKUP_H_

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

#include "TriplePointNaCl.h"
#include "States.h"

namespace csmp
{

  class VLH_LiquidLookup 
  {

  public:
    VLH_LiquidLookup(const double64& externaltemperature);
    ~VLH_LiquidLookup();

    double64                ValueOf(const int& property_index);

    double64                Temperature();
    double64                Pressure();
    double64                MassFractionNaCl();
    double64                Density();
    double64                Enthalpy();
    double64                HeatCapacity();
    double64                Compressibility();
    double64                Viscosity();
    double64                Pmax();
    double64                Tmax();
    double64                TfromP(const double64& press, const double64& t_estimate); 
    double64                DPressureDT();
    double64                DEnthalpyDT();

    std::vector<double64>   properties_at_tmax;

  private:

    const double64&         temperature;

    double64                tcurrent;
    double64                pcurrent;
    double64                xcurrent;
    double64                tmax;
    double64                pmax;
    double64                pnorm;
    double64                t_res; 
    double64                tnorm; 
    double64                dp;

    long                    it;
    long                    t_dim;
    long                    i_max;
    long                    i_min;
    long                    i_guess;
    long                    it_p_max;

    States                  state;

    std::vector<double64>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point

    TriplePointNaCl         tp_nacl;

    void GetTemperatureIndex(const double64& t);

    ErrorHandler&           csmp_error;
  };

  inline double64 VLH_LiquidLookup::Tmax(){ return tmax; }
  inline double64 VLH_LiquidLookup::Pmax(){ return pmax; }

  /**
     author: Thomas Driesner, ETH Zuerich
     contact: thomas.driesner@erdw.ethz.ch
     latest modification: May 27, 2008

     VLH_LiquidLookup generates and queries a lookup table that stores the temperature-pressure-composition [Celsius, bar, mole fraction] coordinates of the halite-saturated liquid curve on the three-phase halite+liquid+vapor surface of the system H2O-NaCl, according to the paper

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.

     In addition, mass density [kg/m3], specific enthalpy [J/kg], specific isobaric heat capacity [J/kg/C], isothermal compressibility [Pa^-1], dynamic viscosity [Pa s] on the curve are available, and are computed by applying the formulations from the paper

     Driesner T. (2007): The system H2O-NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000oC, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4902-4919.

     The lookup table "VLH_LiquidLookupTable.bin" is written as a binary file, and stores the properties listed above in a single, 1-dimensional std::vector<double64>.  The vector comprises subsequent blocks in each of which the values for one property are stored in sequence of ascending temperature from 0 C to the triple point of NaCl (re-check, pleae). Each block is of size t_dim (the number of entries, the first thing that is computed in the constructor), the sequence of property blocks is defined in the file "LookupPropertyIndex.h". The temperature spacing is variable (see member GetTemperatureIndex for details), being densest near the critical temperature of water and largest at the lowest and highest temperatures.
   
     Known issues
     Pressure() still return bars, a legacy from the SoWat development.

     Testing
     The lookup version was tested against the original non-lookup code (May 26, 2007). Max deviations are:
     Composition: -0.5% at 374.03 C, above 374.3 C, deviations are <0.1% and decay rapidly to negligible values.
     Pressure: deviations are negligible over the whole temperature range (.00x % or less)
     Density: (nota bene: this is NOT used in fluid flow simulations in the original ETH scheme by Thomas Driesner) up to 10% at 374.04, rapidly decaying via 3% at 374.6 to less than 1% above 375.5. Notice that the uncertainties are always maximum between lookup grid point and are "exact" at the grid points.
     Enthalpy: (again: not used in fluid flow simulations). 2% at 374.04, rapidly decaying, < 0.05% at t > 380 C
     HeatCapacity: (again: not used in fluid flow simulations): large deviations at t < 375.5, which are solely due to the "hard-wired" value assigned at tcrit_h20 (1e10, this happens in H2OLookup). From 424 C: < 0.1%
     Compressibility: (again: not used ...) Near tcrit_h2o suffers from the same problem as HeatCapacity, <1% is reached at about 385 C.
     Viscosity: (again ...) 0.07% at 374.04, generally negligible
  */
}// namespace csmp
#endif
