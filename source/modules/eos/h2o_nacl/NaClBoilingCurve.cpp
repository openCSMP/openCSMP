// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include<cmath>

#include "NaClBoilingCurve.h"
#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{
  /** custom constructor
  @attention Use only this conctructor. NaClBoilingCurve computes the vapor pressure of liquid NaCl boiling as a function of temperature. The latter is typically computed somewhere outside this class. By using const ref to that external variable as constructor argument, NaClBoilingCurve is able to decipher the current temperature value when being queried for the vapor pressure. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  NaClBoilingCurve::NaClBoilingCurve( const double& externaltemperature )
    :
    temperature_(externaltemperature)
  {
  }


  /** destructor
   */
  NaClBoilingCurve::~NaClBoilingCurve()
  {
  }


  /** vapor pressure at boiling for given temperature (eqn. ? of Driesner & Heinrich (2007))

  @attention No range check is performed, current temperature may be out of valid range if higher than ca. 1000C and less than NaCl triple point temperature (800.7 C).
  @todo Implement range checking and error handling
  */
  double NaClBoilingCurve::VaporPressure()
  {
    tcurrent_ = temperature_;
    return 1.0e5*pow( 10.0, log10(tp_nacl.Pressure()*1.0e-5) + 9418.12* (1.0/(tp_nacl.Temperature()+273.15) - 1.0/(tcurrent_+273.15)) );
  }

}
