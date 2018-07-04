#ifndef NACLSATURATEDVAPOR_H
#define NACLSATURATEDVAPOR_H

#include "CSMP_definitions.h"

#include "NaClSublimationCurve.h"
#include "CriticalCurve.h"
#include "CriticalPointH2O.h"
#include "HaliteLiquidus.h"

/* Changelog
   October 16, 2014, Thomas Driesner - re-ported to CSMP++:
   - removed functionality that is currently not needed (e.g., T and P derivatives)
   - changed interface to [Pa] for pressure and [mass fraction NaCl] for composition
   - pushed formatting towards CSMP++ style guide
   - added doxygen basic documentation
   - basic (manual) testing against CSP version - OK
 */

namespace csmp
{
  /// Mass fraction NaCl of halite-saturated vapor in the system H2O-NaCl as a function of temperature and pressure
  class NaClSaturatedVapor 
  {

  public:

    NaClSaturatedVapor(const double64& externaltemperature, // [C]
                       const double64& externalpressure);   // [Pa]
    ~NaClSaturatedVapor();
    
    double64 MassFractionNaCl();                            // [mass fraction NaCl]

  private:
    NaClSaturatedVapor();

    const double64& temperature_;                           ///< internal reference to temperature [C] in flow code
    const double64& pressure_;                              ///< internal reference to pressure [Pa] in flow code
    
    const double64 k0_;                                     ///< parameter k0 of Driesner & Heinrich, GCA, 2007
    const double64 k1_;                                     ///< parameter k1 of Driesner & Heinrich, GCA, 2007
    const double64 k2_;                                     ///< parameter k2 of Driesner & Heinrich, GCA, 2007
    const double64 k3_;                                     ///< parameter k3 of Driesner & Heinrich, GCA, 2007
    const double64 k4_;                                     ///< parameter k4 of Driesner & Heinrich, GCA, 2007
    const double64 k5_;                                     ///< parameter k5 of Driesner & Heinrich, GCA, 2007
    const double64 k6_;                                     ///< parameter k6 of Driesner & Heinrich, GCA, 2007
    const double64 k7_;                                     ///< parameter k7 of Driesner & Heinrich, GCA, 2007
    const double64 k8_;                                     ///< parameter k8 of Driesner & Heinrich, GCA, 2007
    const double64 k9_;                                     ///< parameter k9 of Driesner & Heinrich, GCA, 2007
    const double64 k10_;                                    ///< parameter k10 of Driesner & Heinrich, GCA, 2007
    const double64 k11_;                                    ///< parameter k11 of Driesner & Heinrich, GCA, 2007
    const double64 k12_;                                    ///< parameter k12 of Driesner & Heinrich, GCA, 2007
    const double64 k13_;                                    ///< parameter k13 of Driesner & Heinrich, GCA, 2007
    const double64 k14_;                                    ///< parameter k14 of Driesner & Heinrich, GCA, 2007
    const double64 k15_;                                    ///< parameter k15 of Driesner & Heinrich, GCA, 2007

    double64 pdummy_;                                       ///< internal pressure parameter [Pa]
    double64 tdummy_;                                       ///< internal temperature parameter [C]
    double64 tcurrent_;                                     ///< internal temperature parameter [C]
    double64 pcurrent_;                                     ///< internal pressure parameter [bar]
    double64 pnorm_;                                        ///< normalized pressure, equation 16 of Driesner & Heinrich, GCA, 2007
    double64 massfractionnacl_;                             ///< self-explanatory ;-)
    double64 xsat_;                                         ///< mole fraction NaCl on liquidus
    double64 xsatatpnacl_;                                  ///< mole fraction NaCl on liquidus at NaCl sublimation pressure
    double64 delxsat_;                                      ///< represents one of the term in equation 15 of Driesner & Heinrich, GCA, 2007
    double64 delxgsat_;                                     ///< represents one of the term in equation 15 of Driesner & Heinrich, GCA, 2007
    double64 delx_;                                         ///< variable used for computing equation 15 of Driesner & Heinrich, GCA, 2007
    double64 j0_;                                           ///< parameter j0 of Driesner & Heinrich, GCA, 2007
    double64 j1_;                                           ///< parameter j1 of Driesner & Heinrich, GCA, 2007
    double64 j2_;                                           ///< parameter j2 of Driesner & Heinrich, GCA, 2007
    double64 j3_;                                           ///< parameter j3 of Driesner & Heinrich, GCA, 2007
    double64 pnacl_;                                        ///< NaCl sublimation pressure [bar]
    double64 pcrit_;                                        ///< Critical pressure of H2O-NaCl [bar]
    
    void   CheckState();
    double64 MassFractionNaCl(const double64& t, const double64& p);
 
    NaClSublimationCurve   naclsubl;
    CriticalCurve          critcurve;
    CriticalPointH2O       cp_h2o;
    HaliteLiquidus         liquidus;
  };

  /**
     @class NaClSaturatedVapor NaClSaturatedVapor.h "eos/h2o_nacl/NaClSaturatedVapor.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "NaClSaturatedVapor" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the composition of halite-saturated vapor as a function of temperature and pressure. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "NaClSaturatedVapor" with the temperature (in C) and pressure (in Pa) as constructor arguments. "NaClSaturatedVapor" has an internal mechanism to make sure that it always uses the current value of that temperature and pressure as they exist in your code. Public member names should be self-explanatory, I hope.

     The NaCl concentration can be retrieved in units of mass fraction. If you want to convert these or need values for molar masses of NaCl and H2O, please use ConvertConcentrationUnitsNaCl.h to stay consistent.

     @code

     double64 t(somevalue), p(anothervalue); // temperature [C] and pressure [in Pa] in user's application
     ...
     NaClSaturatedVapor naclsatvap(t,p); 
     ...
     // here, t and p may be updated
     ...
     x = naclsatvap.MassFractionNaCl(); // compute and assign to x the mass fraction of NaCl in halite-saturated vapor for current t and p

     @endcode

     @section testing Testing
     Manual testing against CSP version - OK

     @attention No range check for temperature and pressure is being performed; these check are considered the responsibility of the user if using this class directly
  */

}//csmp

#endif





