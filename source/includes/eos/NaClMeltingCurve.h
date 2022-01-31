#ifndef NACL_MELTING_CURVE_H
#define NACL_MELTING_CURVE_H

/*   Changelog
     September 9, 2014, Thomas Driesner: re-ported to CSMP++ and simplified class (removed unnecessary functions etc.)
     September 25, 2014, Thomas Driesner: minor revisions, tested agains CSP version - OK
*/

#include "CSMP_definitions.h"

#include "ConvertConcentrationUnitsNaCl.h"
#include "TriplePointNaCl.h"

namespace csmp
{
  /// Temperature and pressure coordinates of dry halite melting curve
  class NaClMeltingCurve
  {
    
  public:
    NaClMeltingCurve( const double& externaltemperature, 
                      const double& externalpressure );
    ~NaClMeltingCurve();
    
    double TmeltFromP();           ///< returns melting temperature [C] for given p [Pa]
    double PmeltFromT();           ///< returns melting pressure [Pa] for given T [C]
    
  private:
    NaClMeltingCurve(); // disable default construction

    const double& temperature_;   ///< reference to temperature [C] in flow code
    const double& pressure_;      ///< reference to fluid pressure [Pa] in flow code

    double        tcurrent_;      ///< temperature [C] for internal use
    double        pcurrent_;      ///< fluid pressure [bar] for internal use

    TriplePointNaCl tp_nacl;
  };


  /**
     @class NaClMeltingCurve NaClMeltingCurve.h "eos/h2o_nacl/NaClMeltingCurve.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "NaClMeltingCurve" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the temperature-pressure coordinates of dry halite melting curve. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "NaClMeltingCurve" with the temperature (in C) and pressure (in Pa) as constructor arguments. "NaClMeltingCurve" has an internal mechanism to make sure that it always uses the current value of that temperature and pressure as they exist in your code. Public member names should be self-explanatory, I hope.

     @code

     double t(somevalue), p(anothervalue); // temperature [C] and pressure [in Pa] in user's application
     ...
     NaClMeltingCurve naclmelt(t,p); 
     ...
     // here, t may be updated
     ...
     cout << "new halite melting pressure: " << naclmelt.PmeltFromT() << "Pa\n"; // compute and output melting pressure for new temperature

     @endcode

     @section issues Known issues

     The underlying equations are valid only in a restricted temperature-pressure range; range checking and error handling has not yet been implemented, see documentation notes for NaClMeltingCurve::TmeltFromP() and NaClMeltingCurve::PmeltFromT()

     @section testing Testing

     testing was done in the period before publication in 2007, prior to writing this documentation, no details are available anymore
  */

}//csmp
#endif
