// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TWOPHASELIQUID_H
#define TWOPHASELIQUID_H

#include "TriplePointNaCl.h"
#include "NaClBoilingCurve.h"
#include "ThreephaseHLV.h"
#include "CriticalCurve.h"
#include "CriticalPointH2O.h"
#include "HaliteLiquidus.h"

/* Changelog
   
   October 16, 2014, Thomas Driesner
   - re-ported to CSMP++
   - cleaned up interface, eliminating functionaility that is currently not being used anymore
   - changed units in interface to [Pa] for pressure and [mass fraction NaCl] for composition
   - added basic doxygen documentation
   - changed coding style towrads CSMP++ style guide recommendations
   - performed manual testing against CSP version - OK
 */
namespace csmp
{
  /// Composition of liquid phase on twophase vapor+liquid envelope
  class TwophaseLiquid
  {

  public:
    TwophaseLiquid(const double& externaltemperature, // [C]
                   const double& externalpressure,    // [Pa]
                   const double& ph2o);               // [Pa]
    ~TwophaseLiquid();
    
    double MassFractionNaCl();                        // [mass fraction NaCl]

  private:

    TwophaseLiquid();

    const double&        temperature_;                ///< internal reference to temperature [C] in flow code
    const double&        pressure_;                   ///< internal reference to fluid pressure [P] in flow code
    const double&        ph2o_;                       ///< internal reference to water boiling pressure (from Water object, which must have been constructed with the same temperature and fluid pressure variables as constructor arguments)

    const double         h1_;                         ///< parameter h1 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h2_;                         ///< parameter h2 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h3_;                         ///< parameter h3 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h4_;                         ///< parameter h4 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h5_;                         ///< parameter h5 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h6_;                         ///< parameter h6 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h7_;                         ///< parameter h7 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h8_;                         ///< parameter h8 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h9_;                         ///< parameter h9 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h10_;                        ///< parameter h10 of equation 11 of Driesner & Heinrich, GCA, 2007
    const double         h11_;                        ///< parameter h11 of equation 11 of Driesner & Heinrich, GCA, 2007

    double               pdummy_;                     ///< internal pressure variable [Pa] for interfacing with HaliteLiquidus
    double               tdummy_;                     ///< internal temperature variable [C] for interfacing with HaliteLiquidus
    double               tcurrent_;                   ///< internal temperature variable [C]
    double               pcurrent_;                   ///< internal pressure variable [bar]
    double               massfractionnacl_;           ///< internal variable for composition 
    double               xsat_;                       ///< internal variable for liquidus composition [mole fraction NaCl]
    double               xl_;                         ///< internal variable for liquid composition [mole fraction NaCl]
    double               pcrit_;                      ///< pressure on critical curve [bar]
    double               xcrit_;                      ///< composition on critical curve [mole fraction NaCl]
    double               pnacl_;                      ///< NaCl boiling pressure [bar]
    double               psat_;                       ///< pressure of threephase vapor+liquid+halite coexistence [bar]
    double               ph2obar_;                    ///< water boiling pressure [bar]
    double               pcritminph2o_;               ///< pcrit_-ph2o_
    double               psatsq_;                     ///< psat_*psat_
    double               ph2osq_;                     ///< ph2obar_*ph2obar_
    double               g1_;                         ///< parameter g1 of equation 11 of Driesner & Heinrich, GCA, 2007
    double               g2_;                         ///< parameter g2 of equation 11 of Driesner & Heinrich, GCA, 2007
    
    void                   CheckState();
    void                   UpdateParameters();
    double               MassFractionNaCl(const double& t, const double& p);

    TriplePointNaCl        nacl_triple;
    NaClBoilingCurve       naclboil;
    ThreephaseHLV          hlv;
    CriticalCurve          critcurve;
    CriticalPointH2O       cp_h2o;
    HaliteLiquidus         liquidus;
  };

  /**
     @class TwophaseLiquid TwophaseLiquid.h "eos/h2o_nacl/TwophaseLiquid.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "TwophaseLiquid" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the composition of liquid on the twophase liquid+vapor envelope as a function of temperature and pressure. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "TwophaseLiquid" with the temperature ("externaltemperature", in C), fluid pressure ("externalpressure", in Pa) and water boiling pressure ("ph2o", in Pa, see attention remark below) as constructor arguments.

     "TwophaseLiquid" has an internal mechanism to make sure that it always uses the current value of temperature, fluid pressure, and water boiling pressure as they exist in your code. Public member names should be self-explanatory, I hope.

     The NaCl concentration can be retrieved in units of mass fraction. If you want to convert these or need values for molar masses of NaCl and H2O, please use ConvertConcentrationUnitsNaCl.h to stay consistent.

     @attention The water boling pressure as constructor argument requires that a Water object is being instantiated before twophase liquid, with the same temperature and fluid pressure as its constructor arguments. Also, before querying TwophaseLiquid::MassFractionNaCl(), one has to makie sure that ph2o has been updated. See code example. I apologize for this odd design but for performance reasons it made sense when created originally.

     @code
     double t(somevalue);          // temperature [C] in user's application
     double p(anothervalue);       // fluid pressure [in Pa] in user's application
     double ph2o(yetanothervalue); // water boling pressure [in Pa] in user's application
     ...
     Water          water(t,p);
     TwophaseLiquid twophase_liquid(t,p,ph2o); 
     ...
     // Updating t,p,ph2o:
     t    = my_new_t_value;
     p    = my_new_p_value;
     ph2o = water.SaturationPressure();
     // And retrieving the new value of mass fraction NaCl in TwophaseLiquid:
     x = twophase_liquid.MassFractionNaCl(); // compute and assign to x the mass fraction of NaCl of TwophaseLiquid at t and p
     ...
     @endcode

     @section testing Testing
     Tested manually against CSP version OK
  */

}//csmp
#endif





