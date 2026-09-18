// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
    Halite(const double& externaltemperature,
           const double& externalpressure);
    ~Halite();
    
    double         MassFractionNaCl() const;           // [XNaCl]
    double         Density();                          // [kg/m^3]
    double         MolarVolume();                      // [cm^3/mol]
    double         Compressibility();                  // [Pa^-1]
    double         Enthalpy();                         // [J/kg] referenced to zero Enthalpy of liquid water at 0 C
    double         HeatCapacity();                     // [J/kg/K]

  private:
    Halite();

    const double   l0_;                                ///< parameter l0 for eq. 2 in Driesner (2007)
    const double   l1_;                                ///< parameter l1 for eq. 2 in Driesner (2007)
    const double   l2_;                                ///< parameter l2 for eq. 2 in Driesner (2007)
    const double   l3_;                                ///< parameter l3 for eq. 3 in Driesner (2007)
    const double   l4_;                                ///< parameter l4 for eq. 3 in Driesner (2007)
    const double   l5_;                                ///< parameter l5 for eq. 3 in Driesner (2007)
    const double&  temperature_;                       ///< reference to temperature [C] in flow code
    const double&  pressure_;                          ///< reference to fluid pressure [Pa] in flow code

    double         tcurrent_;                          ///< temperature [C] for internal use
    double         pcurrent_;                          ///< fluid pressure [bar] for internal use
    double         myt_;                               ///< reduced temperature (t-ttriple) for internal use (eq. 30 in Driesner (2007))
    
    void             CheckState();
    double         Density(const double& t, const double& p);
    double         DDensityDP(const double& t);
    double         DDensityDT(const double& t, const double& p);
    double         HeatCapacity(const double& t, const double& p);
    double         ZeroBarDensity(const double& t);
    double         Enthalpy(const double& t, const double& p);
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
     double t; // temperature [C] in user's application
     double p; // pressure [Pa] in user's application                                                             
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





