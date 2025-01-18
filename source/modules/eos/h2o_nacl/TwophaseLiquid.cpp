#include "TwophaseLiquid.h"

#include<cmath>

#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp{

  /** custom constructor
      
@attention Use only this conctructor. TwophaseLiquid computes the NaCl concentration of liquid on the twophase liquid+vapor surface as a function of temperature and pressure. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, TwophaseLiquid is able to decipher their current values when being queried for the NaCl concentration. Only this constructor is allowed to ensure that functionality, default constructor has been made private. See reamrks on thirs constructor variable in full class documentation and code example.
  */
  TwophaseLiquid::TwophaseLiquid(const double& externaltemperature, 
                                 const double& externalpressure,
                                 const double& ph2o)
    : 
    temperature_(externaltemperature), 
    pressure_(externalpressure), 
    ph2o_(ph2o),
    h1_(0.00168486),
    h2_(0.000219379),
    h3_(438.58),
    h4_(18.4508),
    h5_(-5.6765e-10),
    h6_(6.73704e-06),
    h7_(1.44951e-07),
    h8_(384.904),
    h9_(7.07477),
    h10_(6.06896e-05),
    h11_(0.00762859),
    pdummy_(0.0), 
    tdummy_(0.0), 
    tcurrent_(0.0), 
    pcurrent_(0.0), 
    massfractionnacl_(0.0), 
    xsat_(0.0), 
    xl_(0.0), 
    pcrit_(0.0), 
    xcrit_(0.0), 
    pnacl_(0.0), 
    psat_(0.0),
    pcritminph2o_(0.0),
    psatsq_(0.0), 
    ph2osq_(0.0),
    g1_(0.0), 
    g2_(0.0),
    naclboil(tcurrent_),
    hlv(tcurrent_), 
    critcurve(tcurrent_), 
    liquidus(tdummy_,pdummy_)
  {
  }


  /** destructor
   */
  TwophaseLiquid::~TwophaseLiquid()
  {
  }


  /** Value of mass fraction of dissolved NaCl in liquid on twophase liquid+vapor surface
      @attention The value is based on converting the mole fraction computed from equation 11 in Driesner & Heinrich (2007). Conversion is based on ConvertConcentrationUnitsNaCl.h; to stay consistent, use values and functions provide therein for further computations if you need unit conversion.
   */
  double TwophaseLiquid::MassFractionNaCl()
  { 
    CheckState(); 
    return massfractionnacl_; 
  } 


  /** Update internal variables for current temperature and pressure and trigger computation
   */
  void TwophaseLiquid::CheckState()
  {
    tcurrent_ = temperature_;
    pcurrent_ = pressure_*1.0e-5;
    tdummy_   = tcurrent_;
    pdummy_   = pcurrent_;
 
    UpdateParameters();

    massfractionnacl_ = MassFractionNaCl(tcurrent_,pcurrent_);
    return;
  }


  /** Update parameters for equation 11 of Driesner and Heinrich (2007)
   */
   void TwophaseLiquid::UpdateParameters()
  {
    xcrit_   = Massfraction2XNaCl(critcurve.MassFractionNaCl());
    pcrit_   = critcurve.Pressure()*1.0e-5;
    ph2obar_ = ph2o_*1.0e-5;

    g1_     = h2_+(h1_-h2_) / (1.0+exp((tcurrent_-h3_)/h4_))+h5_*tcurrent_*tcurrent_;
    g2_     = h7_+(h6_-h7_) / (1.0+exp((tcurrent_-h8_)/h9_))+h10_*exp(-h11_*tcurrent_);

    if(tcurrent_ <= nacl_triple.Temperature())
      {
        psat_    = hlv.Pressure()*1.0e-5;
        psatsq_  = psat_*psat_;
        pdummy_  = psat_*1.0e5;
        xsat_    = Massfraction2XNaCl(liquidus.MassFractionNaCl());
        if(tcurrent_ < cp_h2o.Temperature())
          {
            ph2osq_       = ph2obar_*ph2obar_; 
            pcritminph2o_ = pcrit_-ph2obar_;
            if(pcritminph2o_ < 0.0e0) pcritminph2o_ = 0.0e0; // may be relevant VERY close to Tcrit !!!
          }
      }
    else
      {
        pnacl_ = naclboil.VaporPressure()*1.0e-5;
        xsat_  = 1.0;
      }
  }


  /** Equation 11 of Driesner and Heinrich (2007), modified to return mass fraction rather than mole fraction NaCl
   */
  double TwophaseLiquid::MassFractionNaCl(const double& t_, const double& p_)
  {
    if(t_ < cp_h2o.Temperature())
      {
        xl_ = -( sqrt(pcrit_-p_)*(g2_*(psatsq_-ph2osq_+2.0*pcrit_*(ph2obar_-psat_))+g1_*(ph2obar_-psat_)-xsat_)
                 +sqrt(pcritminph2o_)*(g1_*(psat_-p_)+g2_*(p_*p_-psatsq_+2.0*pcrit_*(psat_-p_))+xsat_)
                 +sqrt(pcrit_-psat_)*(g2_*(ph2osq_-p_*p_+2.0*pcrit_*(p_-ph2obar_))+g1_*(p_-ph2obar_)))
          /(sqrt(pcrit_-psat_)-sqrt(pcritminph2o_));
      }
    else if(t_ <= nacl_triple.Temperature() )
      {
        xl_ = xcrit_-(xcrit_+g1_*(pcrit_-psat_)+g2_*(pcrit_*pcrit_-2.0*pcrit_*psat_+psatsq_)-xsat_)/sqrt(pcrit_-psat_)*sqrt(pcrit_-p_)+g1_*(pcrit_-p_)+g2_*(pcrit_-p_)*(pcrit_-p_);
      }
    else
      {
        xl_ = xcrit_-(xcrit_+g1_*(pcrit_-pnacl_)+g2_*(pcrit_*pcrit_-2.0*pcrit_*pnacl_+pnacl_*pnacl_)-1.0)/sqrt(pcrit_-pnacl_)*sqrt(pcrit_-p_)+g1_*(pcrit_-p_)+g2_*(pcrit_-p_)*(pcrit_-p_);
      }
    return XNaCl2Massfraction(xl_) ;
  }


}//csmp
