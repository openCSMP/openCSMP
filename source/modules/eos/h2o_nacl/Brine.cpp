#include "Brine.h"

#include <cmath>

#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{

  /** custom constructor
      
@attention Use only this conctructor. Brine computes thermodynamic properties of H2O-NaCl solutions as a function of temperature, pressure, and composition. These governing parameters are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, Brine is able to decipher their current values when being queried for fluid properties. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  Brine::Brine(const double64&    externaltemperature,
               const double64&    externalpressure,
               const double64&    externalcomposition) 
    : 
    temperature_(externaltemperature),
    pressure_(externalpressure),
    composition_(externalcomposition),
    visc_u0_(  4.71624e-3 ),
    visc_u1_( -4.02030e-6 ),
    visc_muc800_( visc_u0_ + visc_u1_*800.0 ),
    tcurrent_( temperature_+1.0 ),
    pcurrent_( pressure_+1.0 ),
    xcurrent_( composition_-0.1 ),
    density_(0.0),
    molarvolume_(0.0),
    dmolarvolumedt_(0.0),
    dmolarvolumedp_(0.0),
    expansivity_(0.0),
    compressibility_(0.0),
    isochore_(0.0),
    enthalpy_(0.0),
    heatcapacity_(0.0),
    tdummy_(0.0),
    pdummy_(0.0),
    pdummyPa_(0.0),
    muw_(0.0),
    viscosity_(0.0),
    pdummy_mu_(0.0),
    xm_(0.0),
    nacl_a0_(0.0),
    nacl_a1_(0.0),
    h_nacl_a0_(0.0),
    h_nacl_a1_(0.0),
    dnacl_a0dp_(0.0),
    dnacl_a1dp_(0.0),
    a0_(0.0),
    a0_1_(0.0),
    h_a0_(0.0),
    h_a0_1_(0.0),
    da0dp_(0.0),
    da0_1dp_(0.0),
    a1_0_(0.0),
    a1_1_(0.0),
    a1_(0.0),
    h_a1_1_(0.0),
    h_a1_2_(0.0),
    h_a1_(0.0),
    da1_0dp_(0.0),
    da1_1dp_(0.0),
    da1dp_(0.0),
    dev_a0_0_(0.0),
    dev_a0_1_(0.0),
    dev_a0_2_(0.0),
    dev_a0_(0.0),
    dev_a1_0_(0.0),
    dev_a1_1_(0.0),
    dev_a1_2_(0.0),
    dev_a1_(0.0),
    ddev_a0_0dp_(0.0),
    ddev_a0_1dp_(0.0),
    ddev_a0_2dp_(0.0),
    ddev_a0dp_(0.0),
    ddev_a1_0dp_(0.0),
    ddev_a1_1dp_(0.0),
    ddev_a1_2dp_(0.0),
    ddev_a1dp_(0.0),
    dummywater(tdummy_,pdummyPa_),
    hlv(tcurrent_),
    critcurve(tcurrent_)
  {
  }

  /** destructor
   */
  Brine::~Brine()
  {
  }


  /** Molar volume of H2O-NaCl solution in [cm3 mol-1] */
  double64 Brine::MolarVolume()
  { 
    CheckStatus(); 
    return  molarvolume_;
  }  

  /** Temperature derivative of molar volume of H2O-NaCl solution in [cm3 mol-1 C-1] */
  double64 Brine::DMolarVolumeDT()
  { 
    CheckStatus(); 
    return  dmolarvolumedt_; 
  }    

  /** Isothermal compressibility of H2O-NaCl solution in [bar-1] */
  double64 Brine::CompressibilityBar()
  { 
    CheckStatus(); 
    return -dmolarvolumedp_/molarvolume_; 
  }

  /** Isothermal compressibility of H2O-NaCl solution in [Pa-1] */
  double64 Brine::Compressibility()
  { 
    CheckStatus(); 
    return -dmolarvolumedp_/molarvolume_*1.0e-5; 
  }

  /** Thermal expansivity of H2O-NaCl solution in [C-1] */
  double64 Brine::Expansivity()
  { 
    CheckStatus(); 
    return  dmolarvolumedt_/molarvolume_; 
  }

  /** Density of H2O-NaCl solution in [kg m-3] */
  double64 Brine::Density()
  { 
    CheckStatus(); 
    return  MolarVolumeAndX2Density(molarvolume_,xcurrent_); 
  }

  /** Isochore slope (dP/dT)_V of H2O-NaCl solution in [Pa C-1] */
  double64 Brine::IsochoreSlope()
  { 
    CheckStatus(); 
    return -dmolarvolumedt_/(dmolarvolumedp_*1.0e-5);
  }

  /** Specific enthalpy of H2O-NaCl solution in [J kg-1] */
  double64 Brine::Enthalpy()
  { 
    CheckStatus(); 
    return  enthalpy_;    
  }

  /** Specific heat capacity of H2O-NaCl solution in [J kg-1 C-1] */
  double64 Brine::HeatCapacity()
  { 
    CheckStatus(); 
    return  heatcapacity_;    
  }    

  /** Dynamic viscosity of H2O-NaCl solution in [Pa s] */
  double64 Brine::Viscosity()
  { 
    CheckStatus(); 
    return  viscosity_;
  } 



  /** Update internal pressure, temperature, composition variables and decide which formulations to use for computations of properties */
  void Brine::CheckStatus()
  {
    tcurrent_        = temperature_;
    pcurrent_        = pressure_*1.0e-5;
    xcurrent_        = Weight2XNaCl( 100.0*composition_ );
    tdummy_          = tcurrent_;
    pdummy_          = pcurrent_;
    pdummyPa_        = pcurrent_*1.0e5;


    // ***********************************
    // sequence matters, do not change !!!
    // ***********************************
    UpdateParameters( tcurrent_, pcurrent_, xcurrent_ );
	
    if(tcurrent_ > hlv.Tmax() && pcurrent_ < hlv.Pmax()*1.0e-5 && xcurrent_ > 0.1)
      {
        molarvolume_     =   MolarVolumeHighTX( tcurrent_, pcurrent_, xcurrent_ );
        dmolarvolumedt_  = ( MolarVolumeHighTX( tcurrent_+0.01, pcurrent_, xcurrent_ )   - molarvolume_ )/0.01;
        dmolarvolumedp_  = ( MolarVolumeHighTX( tcurrent_, pcurrent_+0.01, xcurrent_ )   - molarvolume_ )/0.01;
      }
    else if(pcurrent_ < 16.0 && tdummy_ >= dummywater.SaturationTemperature()-1.0e-2 && xcurrent_ > 1.0e-4)
      {
        molarvolume_     =   MolarVolumeLowT();
        
        tdummy_         += 0.01;
        dmolarvolumedt_  = ( MolarVolumeLowT() - molarvolume_ ) / 0.01;
        tdummy_         -= 0.01;
        
        pdummy_         += 0.01; pdummyPa_ = pdummy_*1.0e5; pcurrent_ += 0.01;
        dmolarvolumedp_  = ( MolarVolumeLowT() - molarvolume_ ) / 0.01;
        pdummy_         -= 0.01; pdummyPa_ = pdummy_*1.0e5; pcurrent_ -= 0.01;
      }
    else
      {
        molarvolume_     = MolarVolumeFromWater();
        dmolarvolumedt_  = DMolarVolumeDT( tcurrent_ );
        dmolarvolumedp_  = DMolarVolumeDP( tcurrent_ );
      }
    // *** crude workaround for bad compressibilities at low t ***
    if( pcurrent_ < 200.0 && tcurrent_ < 300.0 )
      {
        pdummy_          = 200.0;
        UpdateParameters( tcurrent_, pdummy_, xcurrent_ );
        dmolarvolumedp_  = DMolarVolumeDP( tcurrent_ );
        pdummy_          = pcurrent_;
        UpdateParameters( tcurrent_, pcurrent_, xcurrent_ );
      }
      
    enthalpy_     = Enthalpy( tcurrent_, pcurrent_, xcurrent_ );
    heatcapacity_ = HeatCapacity( xcurrent_ );
    // new
    viscosity_    = Viscosity( tcurrent_, xcurrent_, muw_ );
    // ***
    return;
  }

  /** Update equation parameters

  @attention Notice that many expressions in the following memeber functions have been copy/pasted "as is" from Maple I have never checked them for possible numerical pitfalls such as adding/subtracting small and large numbers etc. or for any performance issues
  */
  void Brine::UpdateParameters(const double64& t_, const double64& p_, const double64& x_)
  {
    // linear parameters
    nacl_a0_    =  330.47 + 0.942876*sqrt(p_) + 0.0817193*p_ - 2.47556e-08*p_*p_ + 3.45052e-10*p_*p_*p_;
    nacl_a1_    = -0.0370751 + 0.00237723*sqrt(p_) + 5.42049e-05*p_ + 5.84709e-09*p_*p_ - 5.99373e-13*p_*p_*p_;
    dnacl_a0dp_ =  0.5*0.942876/sqrt(p_) + 0.0817193 - 2.47556e-08*2.0e0*p_ + 3.45052e-10*3.0e0*p_*p_;
    dnacl_a1dp_ =  0.5*0.00237723/sqrt(p_) + 5.42049e-05 + 5.84709e-09*2.0e0*p_ - 5.99373e-13*3.0e0*p_*p_;
    
    a0_1_       = -54.2958 - 45.7623*exp( -0.000944785*p_ );
    da0_1dp_    = (-45.7623) * (-0.000944785) * exp( -0.000944785*p_ );
    
    a0_         = nacl_a0_ + a0_1_*(1.0-x_) + (-nacl_a0_-a0_1_)*(1.0-x_)*(1.0-x_);
    // n1 = (n1,x_nacl=1) + n11*(1-x_) + (-n1,x_nacl=1 - n11) * (1-x_)^2
    //    cout << "n1,x_nacl=1 = " << nacl_a0 << endl;
    //    cout << "n11        = " << a0_1 << endl; 
    da0dp_      = dnacl_a0dp_ + da0_1dp_*(1.0e0-x_) + (-dnacl_a0dp_-da0_1dp_)*(1.0e0-x_)*(1.0e0-x_);
    
    a1_0_       = -2.6142 - 0.000239092*p_;
    da1_0dp_    = -0.000239092;
    
    a1_1_       = 0.0356828 + 4.37235e-06*p_ + 2.0566e-09*p_*p_; 
    da1_1dp_    = 4.37235e-06 + 2.0e0*2.0566e-09*p_;
    
    a1_         = -a1_0_*sqrt(a1_1_) + 1.0 + a1_0_*sqrt(x_+a1_1_) + x_*a1_0_*sqrt(a1_1_) - x_ - x_*a1_0_*sqrt(1.0+a1_1_) + x_*nacl_a1_;
    da1dp_      = -da1_0dp_*sqrt(a1_1_)
      - a1_0_ / sqrt(a1_1_) * da1_1dp_/2.0
      + da1_0dp_ * sqrt( x_+a1_1_ )
      + a1_0_ / sqrt(x_+a1_1_) * da1_1dp_/2.0
      + x_ * da1_0dp_ * sqrt(a1_1_)
      + x_ * a1_0_ / sqrt(a1_1_) * da1_1dp_/2.0
      - x_ * da1_0dp_ * sqrt(1.0+a1_1_)
      - x_ * a1_0_ / sqrt(1.0+a1_1_) * da1_1dp_ / 2.0
      + x_ * dnacl_a1dp_;
    
    // deviation function stuff
    dev_a0_0_    =  7.60664e6/(p_+472.051)/(p_+472.051);
    ddev_a0_0dp_ = -2.0*7.60664e6/pow(p_+472.051,3);
    dev_a0_1_    = -50.0-86.1446*exp(-0.000621128*p_);
    ddev_a0_1dp_ = (-86.1446)*(-0.000621128)*exp(-0.000621128*p_);
    dev_a0_2_    = 294.318*exp(-0.00566735*p_);
    ddev_a0_2dp_ = 294.318*(-0.00566735)*exp(-0.00566735*p_);
    dev_a0_      = dev_a0_0_*exp(dev_a0_1_*x_)-dev_a0_0_+dev_a0_2_*x_;
    ddev_a0dp_   = ddev_a0_0dp_* exp(dev_a0_1_*x_)
      +dev_a0_0_*ddev_a0_1dp_*x_*exp(dev_a0_1_*x_)
      -ddev_a0_0dp_
      +ddev_a0_2dp_*x_;
    
    dev_a1_0_    = -0.0732761*exp(-0.0023772*p_)-5.2948e-5*p_;
    ddev_a1_0dp_ = (-0.0732761)*(-0.0023772)*exp(-0.0023772*p_)-5.2948e-5;
    dev_a1_1_    = -47.2747+24.3653*exp(-0.00125533*p_);
    ddev_a1_1dp_ = 24.3653*(-0.00125533)*exp(-0.00125533*p_);
    dev_a1_2_    = -0.278529-0.00081381*p_;
    ddev_a1_2dp_ = -0.00081381;
    dev_a1_      = dev_a1_0_*exp(dev_a1_1_*x_)+dev_a1_2_*x_;
    ddev_a1dp_   = ddev_a1_0dp_*exp(dev_a1_1_*x_)
      +dev_a1_0_*ddev_a1_1dp_*exp(dev_a1_1_*x_)
      +ddev_a1_2dp_*x_;
    
    pdummy_   = p_;
    pdummyPa_ = pdummy_*1.0e5;
    tdummy_   = a0_+a1_*t_+dev_a0_*exp(dev_a1_*t_);
    //	cout << "n1 = " << a0 << ", n2 = " << a1 << endl;    
    
    // new
    if(   tcurrent_ < 373.976e0 && pcurrent_ <= dummywater.SaturationPressure()*1.01e-5 // 1.01 is quickfix for lookup-table problems
          && ( mystate == L || mystate == F ) )
      {
        muw_       = dummywater.LiquidViscosityFromT(tcurrent_);
      }
    else if( tcurrent_ >= 373.976e0 && pcurrent_ <= critcurve.Pressure() 
             && ( mystate == L || mystate == F ) )
      {
        pdummy_mu_ = critcurve.Pressure();
        muw_       = dummywater.ViscosityFromTandP(tcurrent_,pdummy_mu_);
      }
    else
      {
        muw_       = dummywater.ViscosityFromTandP(tcurrent_,pcurrent_);
      }

    return;
  }


  /** Molar volume from equation 7 of Driesner, GCA, 2007 */
  double64 Brine::MolarVolumeFromWater()
  { 
    return dummywater.MolarVolume();
  }


  /** Molar volume from equation 18 of Driesner, GCA, 2007 */
  double64 Brine::MolarVolumeHighTX(const double64& t_, const double64& p_, const double64& x_)
  {
    double64 pmax_    = hlv.Pmax()*1.0e-5;
    UpdateParameters( t_, pmax_ , x_ );
    double64 vpmax_   = MolarVolumeFromWater(); //check
    double64 dvdpmax_ = DMolarVolumeDP( t_ ); //check
    UpdateParameters( t_, 1000.0e0, x_ );
    double64 v1000_   = MolarVolumeFromWater(); //check
    double64 a2_      = 1000.0;//pow(10.0e0,(794.405+x*(-3909.85+x*(6262.01-3116.54*x)))+(-0.554022+x*(3.13968+x*(-5.58252+2.96633*x)))*t);
    double64 dummy_ = 
      ((- log(pmax_+a2_)*dvdpmax_*pmax_ 
        - log(pmax_+a2_)*dvdpmax_*a2_ 
        + vpmax_
        + log(1000.0+a2_)*dvdpmax_*pmax_
        + log(1000.0+a2_)*dvdpmax_*a2_-v1000_ 
        )
       /(- log(pmax_+a2_)*pmax_
         - log(pmax_+a2_)*a2_
         + pmax_ 
         + log(1000.0+a2_)*pmax_
         + log(1000.0+a2_)*a2_
         - 1000.0)
       - dvdpmax_
       )
      * (pmax_+a2_)*log(pmax_+a2_)
      - (- log(pmax_+a2_)*dvdpmax_*pmax_
         - log(pmax_+a2_)*dvdpmax_*a2_
         + vpmax_
         + log(1000.0+a2_)*dvdpmax_*pmax_
         + log(1000.0+a2_)*dvdpmax_*a2_
         - v1000_
         )
      / (- log(pmax_+a2_)*pmax_
         - log(pmax_+a2_)*a2_
         + pmax_
         + log(1000.0+a2_)*pmax_
         + log(1000.0+a2_)*a2_
         - 1000.0
         )
      * pmax_;
    double64 v_ = 
      dummy_ 
      + vpmax_
      - ((- log(pmax_+a2_)*dvdpmax_*pmax_
          - log(pmax_+a2_)*dvdpmax_*a2_
          + vpmax_
          + log(1000.0+a2_)*dvdpmax_*pmax_
          + log(1000.0+a2_)*dvdpmax_*a2_
          - v1000_
          )
         /
         ( - log(pmax_+a2_)*pmax_
           - log(pmax_+a2_)*a2_
           + pmax_
           + log(1000.0+a2_)*pmax_
           + log(1000.0+a2_)*a2_
           - 1000.0
           )
         - dvdpmax_
         )
      * (pmax_+a2_)*log(p_+a2_)
      + (- log(pmax_+a2_)*dvdpmax_*pmax_
         - log(pmax_+a2_)*dvdpmax_*a2_
         + vpmax_
         + log(1000.0+a2_)*dvdpmax_*pmax_
         + log(1000.0+a2_)*dvdpmax_*a2_
         - v1000_
         )
      /(- log(pmax_+a2_)*pmax_
        - log(pmax_+a2_)*a2_
        + pmax_
        + log(1000.0+a2_)*pmax_
        + log(1000.0+a2_)*a2_
        - 1000.0
        )
      *p_;
    return v_;
  }
  

  /** Molar volume from equation 17 of Driesner, GCA, 2007 */
  double64 Brine::MolarVolumeLowT()
  {
    double64 tsat_    = dummywater.SaturationTemperature();
    double64 v0_      = dummywater.LiquidSaturationMolarVolumeForP();
    double64 rho0_    = dummywater.LiquidSaturationDensityForP();
    //        in the next line, the expression in brackets is an empirical fit 
    //        to PROST data for drhol/dt at water Saturation curve
    double64 dvdtsat_ = -1000*mh2o*( 151.142-151.848* pow( pcurrent_+2.13768e-07,0.000792594 )
                                     -0.00992859*pcurrent_)/rho0_/rho0_;
    double64 logp_    = log10(pcurrent_);
    double64 lowT_a2_ = 2.0125e-07+3.29977e-09*exp(-4.31279*logp_)-1.17748e-07*logp_+7.58009e-08*logp_*logp_;
    return  2.0*lowT_a2_*tsat_*tsat_*tsat_
      -tsat_*dvdtsat_+v0_
      -3.0*tdummy_*lowT_a2_*tsat_*tsat_
      +tdummy_*dvdtsat_+lowT_a2_*tdummy_*tdummy_*tdummy_;
  }


  /** Temperature derivative of eq. 7 of Driesner, GCA, 2007 */
  double64 Brine::DMolarVolumeDT(const double64& t_)
  { 
    return dummywater.Dvdt()*(a1_+dev_a0_*dev_a1_*exp(dev_a1_*t_));
  }
  

  /** Temperature derivative of eq. 7 of Driesner, GCA, 2007 */
  double64 Brine::DMolarVolumeDP(const double64& t_)
  {
    return dummywater.Dvdpbar()
      + dummywater.Dvdt()*(da0dp_+da1dp_*t_+ddev_a0dp_*exp(dev_a1_*t_)
                           +dev_a0_*ddev_a1dp_*t_*exp(dev_a1_*t_));
  }
  

  /** Specific enthalpy, equation 21 of Driesner, GCA, 2007 */
  double64 Brine::Enthalpy(const double64& t_, const double64& p_, const double64& x_)
  { 
    h_nacl_a0_ = 47.9048 - 0.00936994*p_ + 6.51059e-6*p_*p_;
    h_a0_1_    = -32.1724 + 0.0621255*p_;
    h_a0_      = h_nacl_a0_ 
      + h_a0_1_*(1.0e0-x_) 
      + (-h_nacl_a0_-h_a0_1_)*(1.0e0-x_)*(1.0e0-x_);
    h_nacl_a1_ =  0.241022 + 3.45087e-5*p_ -4.28356e-9*p_*p_;
    h_a1_1_    = -1.69513 - 0.000452781*p_ -6.04279e-08*p_*p_;
    h_a1_2_    =  0.0612567 + 1.88082e-05*p_;
    h_a1_      = -h_a1_1_*sqrt(h_a1_2_) 
      + 1.0e0 
      + h_a1_1_*sqrt(x_+h_a1_2_)
      + x_*h_a1_1_*sqrt(h_a1_2_)
      - x_
      - x_*h_a1_1_*sqrt(1.0e0+h_a1_2_)
      + x_*h_nacl_a1_;
    tdummy_    = h_a0_ + h_a1_*t_;

    return dummywater.Enthalpy();
  }


  /** Specific isobaric heat capacity, temperature derivative of Driesner, GCA, 2007 */
  double64 Brine::HeatCapacity(const double64& x_){ 
    //TD: adjusted to correct version!!!
    return dummywater.HeatCapacity()*h_a1_/(1.0+x_); //CAUTION!!!!!, didn't I remove the /(1.0+x)?
  }


  /** Dynamic viscosity, using Palliser & McKibbin, 1998 */
  double64 Brine::Viscosity(const double64& t_, const double64& x_, const double64& muw_){
    // pure water, use shortcut
    if(x_ == 0.0e0){
      return muw_;
    }
    // with salt, use Palliser&McKibbin's linear interpolation
    xm_ = XNaCl2Weight(x_)/100.0e0;
    if(t_ < 800.0){
      return
        (  muw_ * (1.0 + 3.0*xm_) * pow(((800.0-t_)/800.0),9.0)
           + pow((t_/800.0),9.0)  * (muw_ * (1.0-xm_) + visc_muc800_*xm_))
        / (pow(((800.0-t_)/800.0),9.0)+pow((t_/800.0),9.0));
    }
    else{
      return muw_*(1.0-xm_)+(visc_u0_+visc_u1_*t_)*xm_;
    }
  }
    

}//csmp
