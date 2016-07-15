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

    Brine(const double64& externaltemperature,   // [C]
          const double64& externalpressure,      // [Pa]
          const double64& externalcomposition);  // [mass fraction NaCl]
    ~Brine();

    double64         MolarVolume();              // [cm3 mole-1]
    double64         DMolarVolumeDT();           // [cm3 mole C-1]
    double64         Density();                  // [kg m-3]
    double64         CompressibilityBar();       // [bar-1]
    double64         Compressibility();          // [Pa-1]
    double64         Expansivity();              // [C-1]
    double64         IsochoreSlope();            // [Pa C-1]
    double64         Enthalpy();                 // [J kg-1]
    double64         HeatCapacity();             // [J kg-1 C-1]
    double64         Viscosity();                // [Pa s]
 
  private:

    Brine();                                     ///< default constructor private as public makes no sense

    const  double64& temperature_;               ///< const ref to temperature [C] in flow code
    const  double64& pressure_;                  ///< const ref to fluid pressure [Pa] in flow code
    const  double64& composition_;               ///< const ref to salinity [mass fraction NaCl] in flow code

    const  double64  visc_u0_;                   ///< internal parameter of Viscosity function
    const  double64  visc_u1_;                   ///< internal parameter of Viscosity function
    const  double64  visc_muc800_;               ///< internal parameter of Viscosity function

    double64         tcurrent_;                  ///< internal temperature variable [C]
    double64         pcurrent_;                  ///< internal pressure variable [bar]
    double64         xcurrent_;                  ///< internal composition variable [mole fraction NaCl]

    double64         density_;                   ///< internal variable for density [kg m-3]
    double64         molarvolume_;               ///< internal variable for molar volume of solution [cm3 mole-1]
    double64         dmolarvolumedt_;            ///< temperature derivative of molarvolume_
    double64         dmolarvolumedp_;            ///< pressure derivative of molarvolume_
    double64         expansivity_;               ///< internal variable for themal expansivity [C-1]
    double64         compressibility_;           ///< internal variable for isothermal compressibility [bar-1]
    double64         isochore_;                  ///< internal variable for isochore slope [bar C-1]

    double64         enthalpy_;                  ///< internal variable for specific enthalpy [J kg-1]
    double64         heatcapacity_;              ///< internal variable for specific heat capacity [J kg-1 C-1]

    double64         tdummy_;                    ///< internal temperature variable [C] for interfacing with dummywater member
    double64         pdummy_;                    ///< internal pressure variable [bar]
    double64         pdummyPa_;                  ///< internal pressure variable [Pa] for interfacing with dummywater member

    double64         muw_;                       ///< internal variable pure water dynamic viscosity
    double64         viscosity_;                 ///< internal variable dynamic viscosity of solution 
    double64         pdummy_mu_;                 ///< pressure variable within viscosity computation
    double64         xm_;                        ///< internal composition variable within viscosity computation

    double64         nacl_a0_;                   ///< parameter n1,XNaCl=1 in eq 11 of Driesner, GCA 2007
    double64         nacl_a1_;                   ///< parameter n2,XNaCl=1 in eq 11 of Driesner, GCA 2007
    double64         h_nacl_a0_;                 ///< parameter q1,XNaCl=1 in eq 25 of Driesner, GCA 2007
    double64         h_nacl_a1_;                 ///< parameter q2,XNaCl=1 in eq 26 of Driesner, GCA 2007
    double64         dnacl_a0dp_;                ///< pressure derivative of nacl_a0_
    double64         dnacl_a1dp_;                ///< pressure derivative of nacl_a1_
    double64         a0_;                        ///< parameter n1 of Driesner, GCA 2007
    double64         a0_1_;                      ///< parameter n11 of Driesner, GCA 2007
    double64         h_a0_;                      ///< parameter q1 of Driesner, GCA 2007
    double64         h_a0_1_;                    ///< parameter q11 of Driesner, GCA 2007
    double64         da0dp_;                     ///< pressure derivative of a0_
    double64         da0_1dp_;                   ///< pressure derivative of a0_1_
    double64         a1_0_;                      ///< parameter n21 of Driesner, GCA 2007
    double64         a1_1_;                      ///< parameter n22 of Driesner, GCA 2007
    double64         a1_;                        ///< parameter n2 of Driesner, GCA 2007
    double64         h_a1_1_;                    ///< parameter q21 of Driesner, GCA 2007
    double64         h_a1_2_;                    ///< parameter q22 of Driesner, GCA 2007
    double64         h_a1_;                      ///< parameter q2 of Driesner, GCA 2007
    double64         da1_0dp_;                   ///< pressure derivative of a1_0_
    double64         da1_1dp_;                   ///< pressure derivative of a1_1_
    double64         da1dp_;                     ///< pressure derivative of a1_
    double64         dev_a0_0_;                  ///< parameter n300 of Driesner, GCA 2007
    double64         dev_a0_1_;                  ///< parameter n301 of Driesner, GCA 2007
    double64         dev_a0_2_;                  ///< parameter n302 of Driesner, GCA 2007
    double64         dev_a0_;                    ///< parameter n30 of Driesner, GCA 2007
    double64         dev_a1_0_;                  ///< parameter n310 of Driesner, GCA 2007
    double64         dev_a1_1_;                  ///< parameter n311 of Driesner, GCA 2007
    double64         dev_a1_2_;                  ///< parameter n312 of Driesner, GCA 2007
    double64         dev_a1_;                    ///< parameter n31 of Driesner, GCA 2007
    double64         ddev_a0_0dp_;               ///< pressure derivative of dev_a0_0_
    double64         ddev_a0_1dp_;               ///< pressure derivative of dev_a0_1_
    double64         ddev_a0_2dp_;               ///< pressure derivative of dev_a0_2_
    double64         ddev_a0dp_;                 ///< pressure derivative of dev_a0_
    double64         ddev_a1_0dp_;               ///< pressure derivative of dev_a1_0_
    double64         ddev_a1_1dp_;               ///< pressure derivative of dev_a1_1_
    double64         ddev_a1_2dp_;               ///< pressure derivative of dev_a1_2_
    double64         ddev_a1dp_;                 ///< pressure derivative of dev_a1_

    void             CheckStatus();
    void             UpdateParameters(  const double64& t, const double64& p, const double64& x);
    double64         MolarVolumeFromWater();
    double64         DMolarVolumeDT(    const double64& t_ );
    double64         DMolarVolumeDP(    const double64& t_ );
    double64         MolarVolumeHighTX( const double64& t_, const double64& p_, const double64& x_);
    double64         MolarVolumeLowT();
    double64         Enthalpy(          const double64& t_, const double64& p_, const double64& x_);
    double64         HeatCapacity(      const double64& x_);
    double64         Viscosity(         const double64& t_, const double64& x_, const double64& muw_);

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
     double64 temperature(100.0);       // temperature in your code
     double64 fluid_pressure(10.0e6);   // fluid pressure in your code
     double64 massfraction_nacl(0.1);   // mass fraction NaCl in fluid, in your code
     
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
