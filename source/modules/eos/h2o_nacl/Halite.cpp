// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include<cmath>
#include<iostream>

#include "Halite.h"
#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{
 
  /** custom constructor
      
      @attention Use only this conctructor. Halite computes halite properties as a function of temperature and pressure. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, HaliteLiquidus is able to decipher their current values when being queried for halite properties. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  Halite::Halite(const double& externaltemperature, 
                 const double& externalpressure)
    : l0_(2.1704e3),
      l1_(-2.4599e-1),
      l2_(-9.5797e-05),
      l3_(5.727e-3),
      l4_(2.715e-3),
      l5_(733.4e0),
      temperature_(externaltemperature),
      pressure_(externalpressure),
      tcurrent_(temperature_ - 1.0e0),
      pcurrent_(pressure_-1.0e0),
      myt_(0.0)
  {
  }


  /** destructor
   */
  Halite::~Halite()
  {
  }

  /** Value of mass fraction of NaCl in halite (trivial as that is always 1.0)
   */
  double Halite::MassFractionNaCl() const { return 1.0e0; }


  /** Density [kg/m3] of halite as a function of temperature and pressure
   */
  double Halite::Density()        { CheckState(); return Density(tcurrent_,pcurrent_); } 


  /** Molar volume [cm3/mol] of halite as a function of temperature and pressure; used by other classes of the H2O-NaCl eos module
   */
  double Halite::MolarVolume()    { CheckState(); return mnacl/Density(tcurrent_,pcurrent_)*1.0e3; }


  /** Isothermal compressibility [Pa-1] of halite as a function of temperature and pressure
   */
  double Halite::Compressibility(){ CheckState(); return 1.0e-5/Density(tcurrent_,pcurrent_)*DDensityDP(tcurrent_); } // to be tested: bar vs. Pa 


  /** Specific enthalpy [J/kg] of halite as a function of temperature and pressure
      @attention Enthalpy NEVER has an absolute value but is computed relative to a reference value. For the conventions adopted here, please check Driesner (2007). 
   */
  double Halite::Enthalpy()       { CheckState(); return Enthalpy(tcurrent_,pcurrent_); } 


  /** Specific heat capacity [J/kg/C] of halite as a function of temperature and pressure
   */
  double Halite::HeatCapacity()   { CheckState(); return HeatCapacity(tcurrent_,pcurrent_); }


  /** Update internal temperature and pressure values
   */
  void Halite::CheckState()
  {
    tcurrent_ = temperature_;
    pcurrent_ = pressure_*1.0e-5; // convert to bar as needed by the respective formulas
    return;
  }


  /** Equation 1 of Driesner (2007)
   */
  double Halite::Density( const double& t_, const double& p_ )
  {
    
    return( ZeroBarDensity(t_) + DDensityDP(t_)*p_ );
  }


  /** Equation 2 of Driesner (2007)
   */
  double Halite::ZeroBarDensity(const double& t_)
  { 
    return( l0_ + t_*(l1_ + l2_*t_) ); 
  }


  /** Equation 3 of Driesner (2007)
   */
  double Halite::DDensityDP(const double& t_)
  {
    return l3_+l4_*exp(t_/l5_);
  }


  /** Equation 30 of Driesner (2007)
   */
  double Halite::HeatCapacity(const double& t_, const double& p_)
  {
    myt_ = t_-800.7;
    return 1148.81+2.0*0.275774*myt_+3.0*8.8103e-05*myt_*myt_
      +( -0.0017099 - 2.0*1.91367e-06*t_ - 3.0*2.88485e-09*t_*t_ )*p_
      +( 5.29063e-08 - 2.0*4.81542e-11*t_ + 3.0*2.16915e-13*t_*t_ )*p_*p_;
  } 


  /** Combined Equations 29 and 30 of Driesner (2007) and integrated to get specific enthalpy.
      @attention The integration constant has been optimized manually to get best possible consistency between volumetric and enthalpic data. Minor inconsistencies may still arise due to accumulation of numerical roundoff errors. These should not affect a normal simulation.
   */
  double Halite::Enthalpy(const double& t_, const double& p_)
  {
    myt_ = t_-800.7;
    // although this formula looks odd, please keep it as is;
    return 
      -7.96280801296234130859375e-08
      +6.8762601870112121105+
      562881.+
      -7.13875011336131137796e+00+
      126106.0+582410.0-481803.0
      +1148.81*myt_
      +0.275774*myt_*myt_
      +8.8103e-05*myt_*myt_*myt_
      +(44.6652-0.0017099*t_-1.91367e-06*t_*t_-2.88485e-09*t_*t_*t_)*p_
      +(-7.41999e-05+5.29063e-08*t_-4.81542e-11*t_*t_+2.16915e-13*t_*t_*t_)*p_*p_;
  }

}
