// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef BRINE_H
#define BRINE_H

#include "CSMP_definitions.h"

#include "Water.h"
#include "ThreephaseHLV.h"
#include "States.h"
#include "CriticalCurve.h"

/*   Changelog
     Nov. 27, 2012, Thomas Driesner: initial port to CSMP++
     October 2014, Thomas Driesner: 
     - redid porting to CSMP++
     - changed pressure units of interface to Pa
     - changed concentration units to mass fraction NaCl
     - changed formatting to adhere better to CSMP++ style guide
     - added doxygen documentation
     - basic, manual testing against CSP version - OK
*/

namespace csmp
{
  /// Thermodynamic properties of H2O-NaCl solutions
  class Brine
  {

  public:

    Brine(const double& externaltemperature,   // [C]
          const double& externalpressure,      // [Pa]
          const double& externalcomposition);  // [mass fraction NaCl]
    ~Brine();

    double         MolarVolume();              // [cm3 mole-1]
    double         DMolarVolumeDT();           // [cm3 mole C-1]
    double         Density();                  // [kg m-3]
    double         CompressibilityBar();       // [bar-1]
    double         Compressibility();          // [Pa-1]
    double         Expansivity();              // [C-1]
    double         IsochoreSlope();            // [Pa C-1]
    double         Enthalpy();                 // [J kg-1]
    double         HeatCapacity();             // [J kg-1 C-1]
    double         Viscosity();                // [Pa s]
 
  private:

    Brine();                                     ///< default constructor private as public makes no sense

    const  double& temperature_;               ///< const ref to temperature [C] in flow code
    const  double& pressure_;                  ///< const ref to fluid pressure [Pa] in flow code
    const  double& composition_;               ///< const ref to salinity [mass fraction NaCl] in flow code

    const  double  visc_u0_;                   ///< internal parameter of Viscosity function
    const  double  visc_u1_;                   ///< internal parameter of Viscosity function
    const  double  visc_muc800_;               ///< internal parameter of Viscosity function

    double         tcurrent_;                  ///< internal temperature variable [C]
    double         pcurrent_;                  ///< internal pressure variable [bar]
    double         xcurrent_;                  ///< internal composition variable [mole fraction NaCl]

    double         density_;                   ///< internal variable for density [kg m-3]
    double         molarvolume_;               ///< internal variable for molar volume of solution [cm3 mole-1]
    double         dmolarvolumedt_;            ///< temperature derivative of molarvolume_
    double         dmolarvolumedp_;            ///< pressure derivative of molarvolume_
    double         expansivity_;               ///< internal variable for themal expansivity [C-1]
    double         compressibility_;           ///< internal variable for isothermal compressibility [bar-1]
    double         isochore_;                  ///< internal variable for isochore slope [bar C-1]

    double         enthalpy_;                  ///< internal variable for specific enthalpy [J kg-1]
    double         heatcapacity_;              ///< internal variable for specific heat capacity [J kg-1 C-1]

    double         tdummy_;                    ///< internal temperature variable [C] for interfacing with dummywater member
    double         pdummy_;                    ///< internal pressure variable [bar]
    double         pdummyPa_;                  ///< internal pressure variable [Pa] for interfacing with dummywater member

    double         muw_;                       ///< internal variable pure water dynamic viscosity
    double         viscosity_;                 ///< internal variable dynamic viscosity of solution 
    double         pdummy_mu_;                 ///< pressure variable within viscosity computation
    double         xm_;                        ///< internal composition variable within viscosity computation

    double         nacl_a0_;                   ///< parameter n1,XNaCl=1 in eq 11 of Driesner, GCA 2007
    double         nacl_a1_;                   ///< parameter n2,XNaCl=1 in eq 11 of Driesner, GCA 2007
    double         h_nacl_a0_;                 ///< parameter q1,XNaCl=1 in eq 25 of Driesner, GCA 2007
    double         h_nacl_a1_;                 ///< parameter q2,XNaCl=1 in eq 26 of Driesner, GCA 2007
    double         dnacl_a0dp_;                ///< pressure derivative of nacl_a0_
    double         dnacl_a1dp_;                ///< pressure derivative of nacl_a1_
    double         a0_;                        ///< parameter n1 of Driesner, GCA 2007
    double         a0_1_;                      ///< parameter n11 of Driesner, GCA 2007
    double         h_a0_;                      ///< parameter q1 of Driesner, GCA 2007
    double         h_a0_1_;                    ///< parameter q11 of Driesner, GCA 2007
    double         da0dp_;                     ///< pressure derivative of a0_
    double         da0_1dp_;                   ///< pressure derivative of a0_1_
    double         a1_0_;                      ///< parameter n21 of Driesner, GCA 2007
    double         a1_1_;                      ///< parameter n22 of Driesner, GCA 2007
    double         a1_;                        ///< parameter n2 of Driesner, GCA 2007
    double         h_a1_1_;                    ///< parameter q21 of Driesner, GCA 2007
    double         h_a1_2_;                    ///< parameter q22 of Driesner, GCA 2007
    double         h_a1_;                      ///< parameter q2 of Driesner, GCA 2007
    double         da1_0dp_;                   ///< pressure derivative of a1_0_
    double         da1_1dp_;                   ///< pressure derivative of a1_1_
    double         da1dp_;                     ///< pressure derivative of a1_
    double         dev_a0_0_;                  ///< parameter n300 of Driesner, GCA 2007
    double         dev_a0_1_;                  ///< parameter n301 of Driesner, GCA 2007
    double         dev_a0_2_;                  ///< parameter n302 of Driesner, GCA 2007
    double         dev_a0_;                    ///< parameter n30 of Driesner, GCA 2007
    double         dev_a1_0_;                  ///< parameter n310 of Driesner, GCA 2007
    double         dev_a1_1_;                  ///< parameter n311 of Driesner, GCA 2007
    double         dev_a1_2_;                  ///< parameter n312 of Driesner, GCA 2007
    double         dev_a1_;                    ///< parameter n31 of Driesner, GCA 2007
    double         ddev_a0_0dp_;               ///< pressure derivative of dev_a0_0_
    double         ddev_a0_1dp_;               ///< pressure derivative of dev_a0_1_
    double         ddev_a0_2dp_;               ///< pressure derivative of dev_a0_2_
    double         ddev_a0dp_;                 ///< pressure derivative of dev_a0_
    double         ddev_a1_0dp_;               ///< pressure derivative of dev_a1_0_
    double         ddev_a1_1dp_;               ///< pressure derivative of dev_a1_1_
    double         ddev_a1_2dp_;               ///< pressure derivative of dev_a1_2_
    double         ddev_a1dp_;                 ///< pressure derivative of dev_a1_

    void             CheckStatus();
    void             UpdateParameters(  const double& t, const double& p, const double& x);
    double         MolarVolumeFromWater();
    double         DMolarVolumeDT(    const double& t_ );
    double         DMolarVolumeDP(    const double& t_ );
    double         MolarVolumeHighTX( const double& t_, const double& p_, const double& x_);
    double         MolarVolumeLowT();
    double         Enthalpy(          const double& t_, const double& p_, const double& x_);
    double         HeatCapacity(      const double& x_);
    double         Viscosity(         const double& t_, const double& x_, const double& muw_);

    Water            dummywater;
    ThreephaseHLV    hlv;
    CriticalCurve    critcurve;
    States           mystate;

    
  }; // end class declaration Brine

  /**
     @class Brine Brine.h "eos/h2o_nacl/Brine.h"

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
     Brine is legacy code from Thomas Driesner's developments for the H2O-NaCl system, but it is REQUIRED in the CSMP context. There, it is essential to build accurate lookup tables for the H2O-NaCl system, which are then actually used by CSMP. As long as we don't distribute the lookup tables per se for the various platforms, the Brine legacy will remain in CSMP.
     
     Brine computes the properties of H<sub>2</sub>O-NaCl fluids as a function of temperature-pressure-composition [Celsius, bar, mole fraction NaCl] , according to the paper
     
     Driesner T. (2007): The system H2O-NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000oC, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4902-4919.

     @section usage Usage
     Construct an instance of "Brine" with the temperature (in C), Pressure (Pa) and NaCl Concentration (mass fraction) that live in your code part as constructor arguments. "Brine" has an internal mechanism (via const&) to make sure that it always uses the current values of temperature, pressure, and composition that you have calculated externally. Public member names should be self-explanatory. Don't experiment with making private memebrs public as that will almost inevitably lead to wrong results!

     @code
     double temperature(100.0);       // temperature in your code
     double fluid_pressure(10.0e6);   // fluid pressure in your code
     double massfraction_nacl(0.1);   // mass fraction NaCl in fluid, in your code
     
     Brine    brine(temperature, fluid_pressure, massfraction_nacl); // instantiate a Brine object

     ... // simulation goes on

     // new values for temperature, fluid pressure, and composition have been computed:
     temperature      = newtemperature;
     fluid_pressure   = new_fluid_pressure;
     massfractionnacl = new_massfraction_nacl;

     // retrieve and assign new fluid properties, no need to new communicate he new temperature, fluid pressure, and composition, they are "known" automatically by Brine
     rho = brine.Density();
     mu  = brine.Viscosity();
     cp  = brine.HeatCapacity();
     
     ... // move on
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
     Performance is not optimized. Whenever one of the public members is invoked ALL properties are being computed. There is currently no intention to improve this as Brine should only be used to compute the lookup tables.
     Pressure() still returns values in bars, a legacy from the SoWat development.
     Composition() is in mole fraction NaCl, to convert to other units, the functions in "ConvertConcentrationUnitsNaCl.h" may be useful
     
     @section testing Testing
     testing was done in the period before publication in 2007, prior to writing this documentation, no details are available anymore
  */


}
#endif
