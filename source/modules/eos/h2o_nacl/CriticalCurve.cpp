// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <cmath>

#include "CriticalCurve.h"
#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{

  /** custom constructor
      @attention Use only this conctructor. CriticalCurve computes results a function of temperature, which istypically computed somewhere outside this class. By using const ref to that external temperature as constructor argument, CriticalCurve is able to decipher the current value and report correct results when being queried. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  CriticalCurve::CriticalCurve(const double& externaltemperature)
    : temperature_(externaltemperature), 
      c0_(cp_h2o.Pressure()*1.0e-5),
      c1_(-2.36),
      c2_(0.128534), 
      c3_(-0.023707), 
      c4_(0.00320089),
      c5_(-0.000138917), 
      c6_(1.02789e-07),
      c7_(-4.8376e-11), 
      c8_(2.36), 
      c9_(-0.0131417), 
      c10_(0.00298491), 
      c11_(-0.000130114),
      c12_(0.0), /* don't use, will be calculated */
      c13_(0.0), /* don't use, will be calculated */
      c14_(-0.000488336),
      c1A_(1.0), 
      c2A_(1.5), 
      c3A_(2.0), 
      c4A_(2.5), 
      c5A_(3.0), 
      c6A_(4.0), 
      c7A_(5.0), 
      c8A_(1.0),
      c9A_(2.0), 
      c10A_(2.5), 
      c11A_(3.0),
      d1_(8.0e-5), 
      d2_(1.0e-5), 
      d3_(-1.37125e-7), 
      d4_(9.46822e-10), 
      d5_(-3.50549e-12), 
      d6_(6.57369e-15),
      d7_(-4.89423e-18), 
      d8_(0.0777761), 
      d9_(2.7042e-4), 
      d10_(-4.244821e-7), 
      d11_(2.580872e-10),
      tcurrent_(0.0),
      molefraction_(0.0),
      myt_(0.0), 
      t_(0.0), 
      delta_(0.0),
      dpdt_(0.0)
  {
  }

  /** destructor
   */
  CriticalCurve::~CriticalCurve()
  {
  }

  
  /** Pressure [Pa] on critical curve for given temperature
   */
  double CriticalCurve::Pressure()
  { 
    CheckState(); 
    return Pressure(tcurrent_)*1.0e5; 
  }


  /** Mole fraction NaCl on critical curve for given temperature
   */
  double CriticalCurve::MoleFractionNaCl()
  { 
    CheckState(); 
    return Molefraction(tcurrent_);   
  }


  /** Mass fraction NaCl on critical curve for given temperature
   */
  double CriticalCurve::MassFractionNaCl()
  { 
    return XNaCl2Massfraction(MoleFractionNaCl());
  }


  /** Update internal variables for current temperature
   */
  void CriticalCurve::CheckState()
  {
    tcurrent_       = temperature_;
    return;
  }


  /** Equations 5 of Driesner and Heinrich (2007)
   */
  double CriticalCurve::Pressure( const double& t_ )
  {
    myt_ = t_-cp_h2o.Temperature();
    if(t_ >= cp_h2o.Temperature() && t_ <= 500.0e0)
      {
        return c0_ + c8_*myt_ + c9_*myt_*myt_ + c10_*pow( myt_,c10A_ ) + c11_*pow( myt_,c11A_ );
      }
    else if(t_ > 500.0)
      {
        myt_          = 500.0-cp_h2o.Temperature();
        c12_          = c0_ + c8_*myt_ + c9_*myt_*myt_ +c10_*pow( myt_,c10A_ ) + c11_*pow( myt_,c11A_ );
        c13_          = c8_ + c9_*2.0*myt_ + c10_*c10A_*pow( myt_,c10A_-1.0 ) + c11_*c11A_*pow( myt_,c11A_-1.0 );
        myt_          = t_-500.0;
        return c12_ + c13_*myt_ + c14_*myt_*myt_;
      }
    else
      {
        myt_ = cp_h2o.Temperature()-t_;
        return 
          c0_
          + c1_*myt_
          + c2_*pow( myt_,c2A_ )
          + c3_*pow( myt_,c3A_ )
          + c4_*pow( myt_,c4A_ )
          + c5_*pow( myt_,c5A_ )
          + c6_*pow( myt_,c6A_ )
          + c7_*pow( myt_,c7A_ );
      }    
  }


  /** Temperature derivative of equations 5 of Driesner and Heinrich (2007)
   */
  double CriticalCurve::DPressureDT(const double& t_)
  {
    myt_ = t_-cp_h2o.Temperature();

    if(t_ >= cp_h2o.Temperature() && t_ <= 500.0e0){ 
      return c8_ + 2.0*c9_*myt_ + c10_*c10A_*pow( myt_,(c10A_-1.0) ) + c11_*c11A_*pow( myt_,(c11A_-1.0) );
    }
    else if(t_ > 500.0)
      {
        myt_ = 500.0-cp_h2o.Temperature();
        c13_ = c8_ + c9_*2.0*myt_ + c10_*c10A_*pow( myt_,(c10A_-1.0) ) + c11_*c11A_*pow( myt_,(c11A_-1.0) );
        myt_ = t_-500.0;
        return c13_ + 2.0*c14_*myt_;
      }
    else
      {
        myt_ = cp_h2o.Temperature()-t_; // -> dmyt_/dt = -1
        return 
          -c1_
          -c2_*c2A_*pow( myt_,c2A_-1.0 )
          -c3_*c3A_*pow( myt_,c3A_-1.0 )
          -c4_*c4A_*pow( myt_,c4A_-1.0 )
          -c5_*c5A_*pow( myt_,c5A_-1.0 )
          -c6_*c6A_*pow( myt_,c6A_-1.0 )
          -c7_*c7A_*pow( myt_,c7A_-1.0 );
      }    
  }


  /** Equations 7 of Driesner and Heinrich (2007)
   */
  double CriticalCurve::Molefraction(const double& t_)
  {
    if(t_ <= cp_h2o.Temperature()) return 0.0e0;
    else
      {
        myt_ = t_-cp_h2o.Temperature();
        if( myt_ > 0.0e0 && t_ < 600.0e0){
          return(myt_* (d1_+ myt_*(d2_+ myt_*(d3_+ myt_*(d4_+ myt_*(d5_+ myt_*(d6_+myt_*d7_) ) ) ) ) ) );
        }
        if( t_ >= 600.0e0)
          {
            myt_ = t_-600.0e0;
            return( d8_+ myt_*(d9_+myt_*(d10_+myt_*d11_)) );
          }
      }
    return 0.0;
  }


  /** Temperature derivative of equations 7 of Driesner and Heinrich (2007)
   */
  double CriticalCurve::DMolefractionDT( const double& t_ )
  {
    if( t_ <= cp_h2o.Temperature() ) return 0.0;
    else
      {
        myt_ = t_-cp_h2o.Temperature();
        if( myt_ > 0.0e0 && t_ < 600.0e0 )
          {
            return d1_+ myt_*(2.0*d2_+ myt_*(3.0*d3_+ myt_*(4.0*d4_+ myt_*(5.0*d5_+ myt_*(6.0*d6_+myt_*d7_*7.0)))));
          }
        if( t_ >= 600.0e0)
          {
            myt_ = t_-600.0e0;
            return d9_+ myt_*(2.0*d10_+myt_*d11_*3.0);
          }
      }
    return 0.0;
  }



  /** Find temperature [C] on critical curve for given pressure [Pa]
      @attention NO RANGE CHECK is being performed
   */
  double CriticalCurve::TfromP( const double& press_ )
  {
    // Critical temperature for a given pressure
    // using Newton iteration
    // WARNING: this  code is not fail-save
    double mypress_ = press_ * 1.0e-5;
    int i_ = 0;

    delta_ = 1.0;
    
    if( mypress_ <= 100.)  t_ = 200.;
    if( mypress_ <= 1200.) t_ = 499.;
    else t_ = 950.;
    while( fabs(delta_) > 1.e-10 )
      {
        i_++;
        delta_ = Pressure(t_)-mypress_;
        if(fabs(delta_) <= 1.0e-10) break;
        if(i_ > 50) break;
        dpdt_  = DPressureDT(t_);
        t_    -= delta_/dpdt_;
      }
    return t_;    
  }


  /** Find mole fractions NaCl on critical curve for given pressure [Pa]
      @attention NO RANGE CHECK is being performed
   */
  double CriticalCurve::TfromMolefractionNaCl(const double& molefraction_)
  {
    // Critical temperature for a given pressure
    // using Newton iteration
    // WARNING: this  code is not fail-save

    int i_ = 0;

    delta_ = 0.1;
    
    if(      molefraction_ <=0.    ) t_ =200.;
    else if( molefraction_ <= 0.05 ) t_ = 450.;
    else if( molefraction_ <= 0.1  ) t_ = 600.; 
    else                             t_ = 950.;

    while( fabs(delta_) > 1.e-10 )
      {
        i_++;
        delta_ = Molefraction(t_)-molefraction_;
        if(fabs(delta_) <= 1.0e-10) break;
        if(i_ > 50) break;
        dpdt_  = DMolefractionDT(t_);
        t_    -= delta_/dpdt_;
      }
    return t_;
  }


}
