// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NACLBOILINGCURVE_H
#define NACLBOILINGCURVE_H

#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"

/*   Changelog
     September 25, 2014, Thomas Driesner: re-ported to CSMP++ including changing units to Pa, reduced functionality to what is actually needed, adopted CSMP++ coding styles, tested agains CSP version - OK
*/

namespace csmp
{
  // Vapor pressure of the NaCl boiling curve
  class NaClBoilingCurve{
    
  public:

    NaClBoilingCurve( const double& externaltemperature );
    ~NaClBoilingCurve();
    
    double VaporPressure();

  private:

    NaClBoilingCurve();

    const double& temperature_;  ///< reference to temperature [C] in flow code

    double tcurrent_;            ///< temperature [C] for internal use 

    TriplePointNaCl tp_nacl;
    
  };
  /**
     @class NaClBoilingCurve NaClBoilingCurve.h "eos/h2o_nacl/NaClBoilingCurve.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "NaClBoilingCurve" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the vapor pressure of boiling NaCl liquid. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "NaClBoilingCurve" with the temperature (in C) as constructor argument. "NaClBoilingCurve" has an internal mechanism to make sure that it always uses the current value of that temperature as they exist in your code. Public member names should be self-explanatory, I hope.

     @code

     double t(somevalue); // temperature [C] in user's application
     ...
     NaClBoilingCurve naclboil(t); 
     ...
     // here, t may be updated
     ...
     cout << "new NaCl boiling pressure: " << naclboil.VaporPressure() << "Pa\n"; // compute and output vapor pressure on boiling curve [Pa] for new temperature

     @endcode

     @section issues Known issues

     The underlying equations are valid only in a restricted temperature-pressure range; range checking and error handling has not yet been implemented

     @section testing Testing

     tested against CSP version (ok), which in turn was tested in the period before publication in 2007, prior to writing this documentation, no details are available anymore
  */

}
#endif
