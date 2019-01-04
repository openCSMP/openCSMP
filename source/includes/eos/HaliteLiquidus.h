#ifndef HALITELIQUIDUS_H
#define HALITELIQUIDUS_H

#include "CSMP_definitions.h"

#include "NaClMeltingCurve.h"

/*   Changelog
     March 27, 2013, Thomas Driesner: ported to CSMP++, changed pressure units in interface to Pa.
     September 9, 2014, Thomas Driesner: re-visited the port, cleaned up and adopted CSMP++ style rules
     September 25, 2014, Thomas Driesner: minor revision of documentation, basic testing against CSP version - OK
*/

namespace csmp
{

  /// NaCl concentration of brine on halite saturation surface as function of temperature and pressure 
  class HaliteLiquidus
  {
    public:
      HaliteLiquidus( const double64& externaltemperature, 
                      const double64& externalpressure );
			~HaliteLiquidus();
      
      double64         MassFractionNaCl();
      double64         MoleFractionNaCl();

    private:
      HaliteLiquidus(); // disable default construction

      const double64&  temperature_;        ///< reference to temperature [C] in flow code
      const double64&  pressure_;           ///< reference to fluid pressure [Pa] in flow code

      double64         tcurrent_;           ///< temperature [C] for internal use
      double64         pcurrent_;           ///< fluid pressure [bar] for internal use
      double64         molefraction_nacl_;  ///< mole fraction NaCl, based on molar mass in ConvertConcentrationUnitsNaCl.h

      double64         tmelt_;              ///< halite melting temperature for given pressure
      double64         myt_;                ///< nomalized temperature for internal use
      double64         e0_;                 ///< parameter e0 for eq. 8 in Driesner&Heinrich (2007)
      double64         e1_;                 ///< parameter e1 for eq. 8 in Driesner&Heinrich (2007)
      double64         e2_;                 ///< parameter e2 for eq. 8 in Driesner&Heinrich (2007)
      double64         e3_;                 ///< parameter e3 for eq. 8 in Driesner&Heinrich (2007)
      double64         e4_;                 ///< parameter e4 for eq. 8 in Driesner&Heinrich (2007)
      double64         e5_;                 ///< parameter e5 for eq. 8 in Driesner&Heinrich (2007)

      NaClMeltingCurve naclmelting;

      void             CheckState();
      void             UpdateParameters();

      double64         MolefractionNaCl( const double64& myt_ ); 
  };

  /**
     @class HaliteLiquidus HaliteLiquidus.h "eos/h2o_nacl/HaliteLiquidus.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "HaliteLiquidus" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the composition of halite-saturated brine as a function of temperature and pressure. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "HaliteLiquidus" with the temperature (in C) and pressure (in Pa) as constructor arguments. "HaliteLiquidus" has an internal mechanism to make sure that it always uses the current value of that temperature and pressure as they exist in your code. Public member names should be self-explanatory, I hope.

     The NaCl concentration can be retrieved in units of mass fraction. If you want to convert these or need values for molar masses of NaCl and H2O, please use ConvertConcentrationUnitsNaCl.h to stay consistent.

     @code

     double64 t(somevalue), p(anothervalue); // temperature [C] and pressure [in Pa] in user's application
     ...
     HaliteLiquidus liquidus(t,p); 
     ...
     // here, t and p may be updated
     ...
     x = liquidus.MassFractionNaCl(); // compute and assign to x the mass fraction of NaCl in halite-saturated brine for current t and p

     @endcode

     @section testing Testing

     testing was done in the period before publication in 2007, prior to writing this documentation, no details are available anymore
  */

}//csmp
#endif













