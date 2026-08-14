#include "NaClSaturatedVapor.h"

#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{

  /** custom constructor
      
@attention Use only this conctructor. NaClSaturatedVapor computes the NaCl concentration of halite-saturated brine as a function of temperature and pressure. These are typically computed somewhere outside this class. By using const refs to those external variables as constructor argument, NaClSaturatedVapor is able to decipher their current values when being queried for the NaCl concentration. Only this constructor is allowed to ensure that functionality, default constructor has been made private.
  */
  NaClSaturatedVapor::NaClSaturatedVapor(const double& externaltemperature, 
                                         const double& externalpressure)
    : 
    temperature_(externaltemperature), 
    pressure_(externalpressure), 
    k0_(-0.235694),
    k1_(-0.188838),
    k2_(0.004),
    k3_(0.0552466),
    k4_(0.66918),
    k5_(3.96848e2),
    k6_(4.50e1),
    k7_(-3.2719e-07),
    k8_(1.41699e2),
    k9_(-0.292631),
    k10_(-0.00139991),
    k11_(1.95965e-06),
    k12_(-7.3653e-10),
    k13_(0.904411),
    k14_(0.000769766),
    k15_(-1.18658e-06),
    pdummy_(0.0), 
    tdummy_(0.0), 
    tcurrent_(0.0), 
    pcurrent_(0.0), 
    pnorm_(0.0), 
    massfractionnacl_(0.0), 
    xsat_(0.0), 
    xsatatpnacl_(0.0), 
    delxsat_(0.0), 
    delxgsat_(0.0), 
    delx_(0.0),
    j0_(0.0), 
    j1_(0.0), 
    j2_(0.0), 
    j3_(0.0), 
    pnacl_(0.0), 
    pcrit_(0.0),
    naclsubl(tcurrent_), 
    critcurve(tcurrent_), 
    liquidus(tdummy_,pdummy_)
  {
  }


  /** destructor
   */
  NaClSaturatedVapor::~NaClSaturatedVapor()
  {
  }


  /** Value of mass fraction of dissolved NaCl in halite-saturated brine
      @attention The value is based on converting the mole fraction computed from equation 9 in Driesner & Heinrich (2007). Conversion is based on ConvertConcentrationUnitsNaCl.h; to stay consistent, use values and functions provide therein if you need other units 
   */
  double NaClSaturatedVapor::MassFractionNaCl()
  {
    CheckState();  
    return massfractionnacl_; 
  } 


  /** Update internal variables and trigger calculation of mass fraction
   */
  void NaClSaturatedVapor::CheckState()
  {
    tcurrent_    = temperature_;
    pcurrent_    = pressure_*1.0e-5;

    tdummy_      = tcurrent_;
    pdummy_      = pcurrent_*1.0e5;

    pcrit_       = critcurve.Pressure()*1.0e-5;
    pnacl_       = naclsubl.VaporPressure()*1.0e-5;
    pnorm_       = (pcurrent_-pnacl_)/(pcrit_-pnacl_);
	
    xsat_        = Weight2XNaCl(liquidus.MassFractionNaCl()*100.);

    j0_          = k0_+k1_*exp(-k2_*tcurrent_);
    j1_          = k4_+(k3_-k4_)/(1.0+exp((tcurrent_-k5_)/k6_))+k7_*pow(tcurrent_+k8_,2);
    j2_          = k9_+tcurrent_*(k10_+tcurrent_*(k11_+k12_*tcurrent_));
    j3_          = k13_+tcurrent_*(k14_+k15_*tcurrent_);

    massfractionnacl_ = MassFractionNaCl(tcurrent_,pcurrent_);

    return;
  }


   /** Equations 9 and 12-17 of Driesner and Heinrich (2007)
   */
 double NaClSaturatedVapor::MassFractionNaCl(const double& t_, const double& p_)
  {
    pdummy_         = pnacl_*1.0e5;
    xsatatpnacl_    = Weight2XNaCl(liquidus.MassFractionNaCl()*100.);
    delxgsat_       = log10(xsatatpnacl_);
    delxsat_        = log10(pnacl_/pcrit_);
    delx_           = pnorm_;
    delx_           = 1.0
      +j0_*pow(1.0-delx_,j1_)
      +j2_*(1.0-delx_)
      +j3_*pow(1.0-delx_,2)
      +(-1.0-j0_-j2_-j3_)*pow(1.0-delx_,3);
    delx_          *= (delxsat_-delxgsat_);
    delx_          += delxgsat_;
    delx_          -= log10(pnacl_/p_);
    delx_           = pow(10.0,delx_);
    pdummy_         = p_;
    delx_           = xsat_/delx_;
    return  XNaCl2Massfraction(delx_);
  }

}//csmp    
