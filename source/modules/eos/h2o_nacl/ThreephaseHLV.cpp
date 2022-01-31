#include "ThreephaseHLV.h"

using namespace std;

namespace csmp
{
  ThreephaseHLV::ThreephaseHLV(const double& externaltemperature)
    :	temperature_(externaltemperature),
      f0 (      4.64e-3),
      f1 (      5.0e-07),
      f2 (     16.9078),
      f3 (   -269.148),
      f4 (   7632.04),
      f5 ( -49563.6),
      f6 ( 233119.0),
      f7 (-513556.0),
      f8 ( 549708.0),
      f9 (-284628.0),
      f10(  5.0e-4-(f0+f1+f2+f3+f4+f5+f6+f7+f8+f9)),
      tcurrent_(0.0),
      myt_(0.0),
      pressure_(0.0),
      dpressuredt_(0.0)
  {
  }


  ThreephaseHLV::~ThreephaseHLV()
  {
  }


  double ThreephaseHLV::Temperature() const { return temperature_; }

  double ThreephaseHLV::Pressure()          { CheckStatus(); return pressure_*1.0e5; }

  double ThreephaseHLV::DPressureDT()       { CheckStatus(); return dpressuredt_*1.0e5; }


  void ThreephaseHLV::CheckStatus()
  {
    tcurrent_           = temperature_;
    myt_                = tcurrent_/nacl_triple.Temperature();
    pressure_           = Pressure(myt_);
    dpressuredt_        = DPressureDT( myt_, tcurrent_ );

    return;
  }

  double ThreephaseHLV::Pressure(const double& myt_)
  {
    // Calculates the Pressure on the Halite-Lquid-Vapor coexistence surface
    // Current version is in bars as a function of temperature t [Celsius]
    return f0+myt_*(f1+myt_*(f2+myt_*(f3+myt_*(f4+myt_*(f5+myt_*(f6+myt_*(f7+myt_*(f8+myt_*(f9+myt_*f10)))))))));
  }


  double ThreephaseHLV::DPressureDT(const double& myt_, const double& t_)
  {
    return myt_*(f1+myt_*(2.0*f2+myt_*(3.0*f3+myt_*(4.0*f4+myt_*(5.0*f5
                                                                 +myt_*(6.0*f6+myt_*(7.0*f7+myt_*(8.0*f8+myt_*(9.0*f9+myt_*10.0*f10)))))))))/t_;
  }

  double ThreephaseHLV::Tmax()
  {
    // this was found by iteration and is exact within 1.e-13 K
    return 5.946324429298147e+02;
  }

  double ThreephaseHLV::Pmax()
  {
    // this was found by iteration and is exact within DPressureDT(Tmax)*1.e-15 K,
    return 3.90147444337959825588768580928445e+07;
  }

}//csmp
