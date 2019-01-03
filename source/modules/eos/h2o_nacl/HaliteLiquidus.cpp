#include <cmath>
#include "HaliteLiquidus.h"
#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{

  /** custom constructor
      
@attention Use only this conctructor. HaliteLiquidus computes the NaCl concentration of halite-saturated brine as a function of temperature and pressure. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, HaliteLiquidus is able to decipher their current values when being queried for the NaCl concentration. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  HaliteLiquidus::HaliteLiquidus( const double64& externaltemperature, 
                                  const double64& externalpressure )
    :	temperature_(externaltemperature),
      pressure_(externalpressure),
      tcurrent_( temperature_-1.0 ),
      pcurrent_( pressure_-1.0 ),
      molefraction_nacl_(0.0), 
      tmelt_(0.0),
      myt_(0.0),
      e0_(0.0),
      e1_(0.0),
      e2_(0.0),
      e3_(0.0),
      e4_(0.0),
      e5_(0.0),
      naclmelting( tcurrent_, pcurrent_ )
  {
  }


  /** destructor
   */
  HaliteLiquidus::~HaliteLiquidus()
  {
  }


  /** Value of mass fraction of dissolved NaCl in halite-saturated brine
      @attention The value is based on converting the mole fraction computed from equation 8 in Driesner & Heinrich (2007). Conversion is based on ConvertConcentrationUnitsNaCl.h; to stay consistent, use values and functions provide therein for further computations 
   */
  double64 HaliteLiquidus::MassFractionNaCl()
  { 
    CheckState(); 
    return XNaCl2Massfraction( molefraction_nacl_ );
  }

  double64 HaliteLiquidus::MoleFractionNaCl()
  {
      CheckState();
      return molefraction_nacl_;
  }

  /** Update internal variables for current temperature and pressure
   */
  void HaliteLiquidus::CheckState()
  {
    tcurrent_               = temperature_;
    pcurrent_               = pressure_;
    
    UpdateParameters();
    
    molefraction_nacl_      = MolefractionNaCl(myt_);
    return;
  }


  /** Equation 8 of Driesner and Heinrich (2007)
   */
  double64 HaliteLiquidus::MolefractionNaCl( const double64& myt_ )
  {
    return e0_+myt_ * ( e1_+myt_ * ( e2_+myt_ * ( e3_+myt_ * ( e4_+e5_*myt_ ) ) ) );
  }


  /** Update parameters for MoleFractionNaCl( const double64& myt_ )
   */
  void HaliteLiquidus::UpdateParameters()
  {
    tmelt_     =  naclmelting.TmeltFromP();
    myt_       =  tcurrent_/tmelt_;

    pcurrent_ *= 1.0e-5; // following equation is for bar ...

    e0_ =  0.0989944  + pcurrent_ * (  3.30796e-06 - 4.71759e-10*pcurrent_ );
    e1_ =  0.00947257 + pcurrent_ * ( -8.66460e-06 + 1.69417e-09*pcurrent_ );
    e2_ =  0.610863   + pcurrent_ * ( -1.51716e-05 + 1.19290e-08*pcurrent_ );
    e3_ = -1.64994    + pcurrent_ * (  0.000203441 - 6.46015e-08*pcurrent_ );
    e4_ =  3.36474    + pcurrent_ * ( -0.000154023 + 8.17048e-08*pcurrent_ );
    e5_ =  1.0-e0_-e1_-e2_-e3_-e4_;
  }

}
