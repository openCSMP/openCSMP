#ifndef HALITE_H
#define HALITE_H

#include "CSMP_definitions.h"

/*   Changelog
     March 27, 2013, Thomas Driesner: ported to CSMP++, changed pressure units in interface to Pa.
     September 25, 2014, Thomas Driesner: re-visited the port, cleaned up and adopted CSMP++ style rules, basic testing - OK
*/

namespace csmp{

  /// Properties of the mineral halite (crystalline NaCl) as a function of temperature and pressure
  class Halite
  {

  public:
    Halite(const double64& externaltemperature,
           const double64& externalpressure);
    ~Halite();
    
    double64         MassFractionNaCl() const;           // [XNaCl]
    double64         Density();                          // [kg/m^3]
    double64         MolarVolume();                      // [cm^3/mol]
    double64         Compressibility();                  // [Pa^-1]
    double64         Enthalpy();                         // [J/kg] referenced to zero Enthalpy of liquid water at 0 C
    double64         HeatCapacity();                     // [J/kg/K]

  private:
    Halite();

    const double64   l0_;                                ///< parameter l0 for eq. 2 in Driesner (2007)
    const double64   l1_;                                ///< parameter l1 for eq. 2 in Driesner (2007)
    const double64   l2_;                                ///< parameter l2 for eq. 2 in Driesner (2007)
    const double64   l3_;                                ///< parameter l3 for eq. 3 in Driesner (2007)
    const double64   l4_;                                ///< parameter l4 for eq. 3 in Driesner (2007)
    const double64   l5_;                                ///< parameter l5 for eq. 3 in Driesner (2007)
    const double64&  temperature_;                       ///< reference to temperature [C] in flow code
    const double64&  pressure_;                          ///< reference to fluid pressure [Pa] in flow code

    double64         tcurrent_;                          ///< temperature [C] for internal use
    double64         pcurrent_;                          ///< fluid pressure [bar] for internal use
    double64         myt_;                               ///< reduced temperature (t-ttriple) for internal use (eq. 30 in Driesner (2007))
    
    void             CheckState();
    double64         Density(const double64& t, const double64& p);
    double64         DDensityDP(const double64& t);
    double64         DDensityDT(const double64& t, const double64& p);
    double64         HeatCapacity(const double64& t, const double64& p);
    double64         ZeroBarDensity(const double64& t);
    double64         Enthalpy(const double64& t, const double64& p);
  };


  /**
     @class Halite Halite.h "eos/h2o_nacl/Halite.h"
                                                                                  
     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch

     @section motivation Motivation
     "Halite" is legacy code from Thomas Driesner's developments for the H2O-NaCl system and is mostly used to build accurate lookup tables for the H2O-NaCl system, which are then actually used by CSMP++. It can be used instead of HaliteLookup but will be slower.

     "Halite" computes various properties of the mineral halite (solid crystalline salt, NaCl) as a function of temperature and pressure [Celsius, Pa] according to the papers

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.

     Driesner T. (2007): The system H2O-NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000oC, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4902-4919.
                               
     @section usage Usage                                                                                             
     Construct an instance of "Halite" with temperature (in C) and pressure (in Pa) as they exist in the code that is supposed to use "Halite" as constructor variables. "Halite" has an internal mechanism to make sure that it always uses the current values of temperature and pressure . Public member names should be self-explanatory, I hope.

     @code                                                                                                            
     double64 t; // temperature [C] in user's application
     double64 p; // pressure [Pa] in user's application                                                             
     Halite halite(t,p);
     ...
     t = some_value;
     p = some_other_value;
     cout << halite.Density() << endl; // will return density at the new t and p conditions
     @endcode                                                                                                         
                                                                                                                      
     @section dependencies Dependencies                                                                               
     requires "ConvertConcentrationUnitsNaCl.h" to ensure consistent unit conversions   

     @section testing Testing
     testing was done in the period before publication in 2007, prior to writing this documentation, no details are available anymore
  */                                                                                                                  

}
#endif





