#include "NaClMeltingCurve.h"

using namespace std;

namespace csmp
{
  /** custom constructor
  
  @attention Use only this conctructor. NaClMeltingCurve computes temperature or pressure of halite melting as a function of pressure or temperature, respectively. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, NaClMeltingCurve is able to decipher their current values when being queried for the values. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  NaClMeltingCurve::NaClMeltingCurve( const double64& externaltemperature, 
                                      const double64& externalpressure )
    :
    temperature_(externaltemperature),
    pressure_(externalpressure),
    tcurrent_(0.0),
    pcurrent_(0.0)
  {
  }


  /** destructor
   */
  NaClMeltingCurve::~NaClMeltingCurve()
  {
  }


  /** melting temperature for given pressure (eqn. 1 of Driesner & Heinrich (2007))

  @attention No range check is performed, current pressure may be out of valid range if higher than 500 MPa and less than NaCl triple point pressure.
  @todo Implement range checking and error handling
  */
  double64 NaClMeltingCurve::TmeltFromP()
  { //ok
    pcurrent_ = pressure_*1.0e-5;
    return tp_nacl.Temperature() + 0.024726*pcurrent_;
  }


  /** melting pressure for given temperature (eqn. 1 of Driesner & Heinrich (2007))

  @attention No range check is performed, current temperature may be out of valid range if higher than ca. 925C and less than NaCl triple point temperature (800.7 C).
  @todo Implement range checking and error handling
  */
  double64 NaClMeltingCurve::PmeltFromT()
  { //ok
    tcurrent_ = temperature_;
    return ( tcurrent_ - tp_nacl.Temperature() ) / 0.024726e-5;
  }

}//csmp
