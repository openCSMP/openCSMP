#include "NaClSublimationCurve.h"

#include<cmath>

using namespace std;

namespace csmp{

  /** custom constructor
  @attention Use only this conctructor. NaClSublimationCurve computes the vapor pressure of NaCl sublimation as a function of temperature. The latter is typically computed somewhere outside this class. By using const ref to that external variable as constructor argument, NaClSublimationCurve is able to decipher the current temperature value when being queried for the vapor pressure. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  NaClSublimationCurve::NaClSublimationCurve(const double& externaltemperature)
    :
    temperature_(externaltemperature)
  {
  }

  /** destructor */
  NaClSublimationCurve::~NaClSublimationCurve()
  {
  }

  /** vapor pressure at sublimation for given temperature (eqn. ? of Driesner & Heinrich (2007))

  @attention No range check is performed, function makes only sense below NaCl triple point temperature (800.7 C).
  @todo Implement range checking and error handling
  */

  double NaClSublimationCurve::VaporPressure()
  {
    tcurrent_ = temperature_;
    return 1.0e5*pow(10.0,log10(tp_nacl.Pressure()*1.0e-5)+11806.1*(1.0/(tp_nacl.Temperature()+273.15)-1.0/(tcurrent_+273.15)));
  }

}//csmp
