#include "TwophaseVapor.h"

#include<cmath>

#include "ConvertConcentrationUnitsNaCl.h"


using namespace std;

namespace csmp{
  /** custom constructor
      
@attention Use only this conctructor. TwophaseVapor computes the NaCl concentration of liquid on the twophase liquid+vapor surface as a function of temperature and pressure. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, TwophaseVapor is able to decipher their current values when being queried for the NaCl concentration. Only this constructor is allowed to ensure that functionality, default constructor has been made private. See remarks on thirs constructor variable in full class documentation and code example.
  */

  TwophaseVapor::TwophaseVapor(const double& externaltemperature, 
                               const double& externalpressure,
                               const double& ph2o)
    : 
    temperature_(externaltemperature), 
    pressure_(externalpressure),
    ph2o_(ph2o),
    tcurrent_(0.0),
    pcurrent_(0.0),
    massfractionnacl_(0.0),
    liquidus(tcurrent_,pcurrent_),
    twophase_l(tcurrent_,pcurrent_,ph2o_),
    naclsatvap(tcurrent_,pcurrent_)
  {
  }


  /** destructor
   */
  TwophaseVapor::~TwophaseVapor()
  {
  }


  /** Value of mass fraction of dissolved NaCl in vapor on twophase liquid+vapor surface
      @attention The value is based on converting the mole fraction computed from equations 12-17 in Driesner & Heinrich (2007). Conversion is based on ConvertConcentrationUnitsNaCl.h; to stay consistent, use values and functions provide therein for further computations if you need unit conversion.
   */
  double TwophaseVapor::MassFractionNaCl()                  
  { 
    CheckState(); 
    return massfractionnacl_; 
  } 



  /** Update internal variables for current temperature and pressure and trigger computation
   */
  void TwophaseVapor::CheckState(){
    tcurrent_ = temperature_;
    pcurrent_ = pressure_;
    massfractionnacl_ = XNaCl2Massfraction(MoleFractionNaCl());
  }


  /** Equation 13 of Driesner and Heinrich (2007)
   */
  double TwophaseVapor::MoleFractionNaCl(){
    double result;
    result  = Massfraction2XNaCl(naclsatvap.MassFractionNaCl());
    result *= Massfraction2XNaCl(twophase_l.MassFractionNaCl());
    result /= Massfraction2XNaCl(liquidus.MassFractionNaCl());
    return result;
  }

}//csmp
