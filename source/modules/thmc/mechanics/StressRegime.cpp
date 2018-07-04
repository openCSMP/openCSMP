#include "StressRegime.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp{

StressRegime::StressRegime( double64 Sv, double64 SH, double64 Sh, double64 trend )
  : Sv_( Sv ),
    SH_( SH ),
    Sh_( Sh ),
    trend_( trend ),
    PI_( 3.14159265358979324F ){

  CheckTrend();
  EstablishPrincipalStressUnitVectors();

}


void StressRegime::EstablishPrincipalStressUnitVectors()
{
  //establishing stress vectors for given regime
  double64 alpha;
  //vertical
  Svv_(0) = 0.;
  Svv_(1) = 0.;
  Svv_(2) = 1.;
  if( trend_ <= 90. ){
    alpha = 90. - trend_;
    //maximum horizontal
    SHv_(0) = cos( alpha*PI_ / 180 );
    SHv_(1) = -sin( alpha*PI_ / 180 );
    SHv_(2) = 0.;
  } else if( trend_ > 90. ) {
    alpha = trend_ - 90.;
    //maximum horizontal
    SHv_(0) = cos( alpha*PI_ / 180 );
    SHv_(1) = sin( alpha*PI_ / 180 );
    SHv_(2) = 0.;
  }

  //minimum horizontal
  EstablishMinimumHorizontalStressVector();

}

void StressRegime::EstablishMinimumHorizontalStressVector()
{
  double64 alpha = trend_ + 90;
  if( alpha > 180. )
    alpha -= 180.;
  if( alpha >= 90. ){
    alpha -= 90.;
    Shv_(0) = cos( alpha*PI_ / 180 );
    Shv_(1) = sin( alpha*PI_ / 180 );
    Shv_(2) = 0.;
  } else {
    alpha = 90. - alpha;
    Shv_(0) = cos( alpha*PI_ / 180 );
    Shv_(1) = -sin( alpha*PI_ / 180 );
    Shv_(2) = 0.;
  }

}

double64 StressRegime::VerticalStressMagnitude() const{

  return Sv_;

}

double64 StressRegime::MaximumHorizontalStressMagnitude() const{

  return SH_;

}

double64 StressRegime::MinimumHorizontalStressMagnitude() const{

  return Sh_;

}

double64 StressRegime::MaximumHorizontalStressTrend() const{

  return trend_;

}

VectorVariable<3U>  StressRegime::VerticalStressUnitVector() const{

  return Svv_;

}

VectorVariable<3U>  StressRegime::MinimumHorizontalStressUnitVector() const{

  return Shv_;

}

VectorVariable<3U>  StressRegime::MaximumHorizontalStressUnitVector() const{

  return SHv_;

}


void StressRegime::CheckTrend()
{
  if( trend_ > 180. && trend_ < 270. )
    trend_ -= 90.;
  if( trend_ > 270. && trend_ < 360. )
    trend_ -= 180.;
  if( trend_ > 360. || trend_ < 0. )
    throw csmp::Exception( FATAL_ERROR, "NormalFromStressRegime::CheckTrend",
                                 "Range of stress trend exceeded (0-180)",
                                 "specify valid trend value(Maximum Horizontal Stress");
}


} //csmp
