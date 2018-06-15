#include <cmath>

#include "BrooksCoreySaturationFunctionswithHysteresis.h"
#include "FlowFunctions.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

// 1) destructor of Brooks Corey Pressure feilds
template<size_t dim, template<size_t> class USER> BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::BrooksCoreySaturationFunctionswithHysteresis()
{
}

// 3) a  constructor of Brooks Corey Pressure feilds
template<size_t dim, template<size_t> class USER>
BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::BrooksCoreySaturationFunctionswithHysteresis(double64 aw, double64 ao,  double64 cw, double64 co) :
aw_(aw),ao_(ao),cw_(cw),co_(co)
{
}

// 4) Brooks & Corey capillary pressure formula
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::Pc ( TARGET_PLACEMENT& p ) {
   
    auto sH2O = p.Interpolate(User()->key_sH2O);
    auto sCO2 = p.Interpolate(User()->key_sCO2);
    auto srH2O = p.Interpolate(User()->key_srH2O);
    auto srCO2 = p.Interpolate(User()->key_srCO2);
    
    return cw_*pow((1.0 - srH2O)/(sH2O - srH2O),aw_) +
           co_*pow((1.0 - srCO2)/(sCO2 - srCO2),ao_);
}

//
// Eq. 3 from Skjaeveland et al. 2000, we calculate the oil entry pressure "co".
// double64 BrooksCoreyCapillaryPressure::EstimateOilEntryPressure(double64 S) {
//    return -cw_*pow((1.-S-Sor_)/(1.-Sor_),ao_)*pow((1.-Swr_)/(S-Swr_),aw_);
//}
//

// 5) The first derivative of Brooks & Corey capillary pressure formula
//  by assuming the dSCO2 = -dSH2O
//
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::dpcds( TARGET_PLACEMENT& p) {
   
   auto sH2O = p.Interpolate(User()->key_sH2O);
   auto srH2O = p.Interpolate(User()->key_srH2O);
   auto sCO2 = p.Interpolate(User()->key_sCO2);
   auto srCO2 = p.Interpolate(User()->key_srCO2);
    
   return -aw_*cw_*pow((1.0 - srH2O) / (sH2O - srH2O), aw_) / (sH2O - srH2O) +
           ao_*co_*pow((1.0 - srCO2) / (sCO2 - srCO2), ao_) / (sCO2 - srCO2);
}

// 7) Oil residual saturation has been estimated from Land's formula based on the initial oil saturations.

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::OilResidualSaturation( TARGET_PLACEMENT& p) {
  return 1./(C_land_+1./p.Interpolate(User()->key_sCO2)) ;
}

//
// Water residual saturation has been estimated from the intersection of imbibition curve with saturation axis.
// We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
// the water residual saturation and oil residual saturation.

//double64 BrooksCoreyCapillaryPressure::WaterResidualSaturation( double64 S, double64 Sor) {
//  return 1.-(1.-S)/(1.-(pow((-cw_/co_),1./aw_)*pow(((1.-Sor)/(1.-S-Sor)),(ao_/aw_)))) ;
//}
//


// 8) For calculating the relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
// The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
// the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
//

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::krw( TARGET_PLACEMENT& p ) {
  
  auto sH2O = p.Interpolate(User()->key_sH2O);
  auto srH2O = p.Interpolate(User()->key_srH2O);
  auto sCO2 = p.Interpolate(User()->key_sCO2);
  auto srCO2 = p.Interpolate(User()->key_srCO2);

  double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
  double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000

  double64 krwww = pow(Snw, 3.+2.*aw_) ; // refer to Eq. 10 from Skjaeveland et al. 2000
  double64 krwow = (1. - pow(Sno, 2.*ao_+1.)) * square(1. - Sno);

  return (cw_*krwww - co_*krwow) / (cw_ - co_) ;
}

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::krn( TARGET_PLACEMENT& p) {
  
  auto sH2O = p.Interpolate(User()->key_sH2O);
  auto srH2O = p.Interpolate(User()->key_srH2O);
  auto sCO2 = p.Interpolate(User()->key_sCO2);
  auto srCO2 = p.Interpolate(User()->key_srCO2);

  double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
  double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000

  double64 kroww = (1. - pow(Snw,2.*aw_+1.)) * square(1. - Snw);
  double64 kroow = pow(Sno, 3.+2.*ao_) ; // refer to Eq. 12 from Skjaeveland et al. 2000

  return (cw_*kroww - co_*kroow) / (cw_ - co_) ;
}

// 9) For calculating the 1st derivative of relative permeability

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::dkrwds( TARGET_PLACEMENT& p ) {
    
    auto sH2O  = p.Interpolate(User()->key_sH2O);
    auto srH2O = p.Interpolate(User()->key_srH2O);
    auto sCO2  = p.Interpolate(User()->key_sCO2);
    auto srCO2 = p.Interpolate(User()->key_srCO2);

    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000

    double64 dSo = -1.; // if change of water saturation will change co2 saturations.

    double64 dSnw = 1./(1.-srH2O-srCO2) ;

    double64 dSno = dSo * dSnw; // dSo/(1.-p.Interpolate(User()->key_srH2O)-p.Interpolate(User()->key_srCO2)) ;

    double64 dkrwww = (3.+2.*aw_)*pow(Snw,(2.+2.*aw_))*dSnw ;
    double64 dkrwow = (-(2*ao_+1)*pow(Sno,(2.*ao_))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*ao_+1.)))*(1.-Sno))*dSno ;

    return (cw_*dkrwww-co_*dkrwow)/(cw_-co_) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::dkrnds( TARGET_PLACEMENT& p) {

    auto sH2O  = p.Interpolate(User()->key_sH2O);
    auto srH2O = p.Interpolate(User()->key_srH2O);
    auto sCO2  = p.Interpolate(User()->key_sCO2);
    auto srCO2 = p.Interpolate(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-srH2O-srCO2) ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Interpolate(User()->key_srH2O)-p.Interpolate(User()->key_srCO2)) ;
    
    double64 dkroww = (-(2*aw_+1)*pow(Snw,(2.*aw_))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*aw_+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*ao_)*pow(Sno,(2.+2.*ao_))*dSno  ;

    return (cw_*dkroww-co_*dkroow)/(cw_-co_) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000

}
  
  
  // 13) the water fractional flow is calculated from the water relative permeability and oil relative permeability.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::WaterFractionalFlow( TARGET_PLACEMENT& p, double64 S) {
    
    p.Interpolate(User()->key_sH2O) = S ;  // I am not sure this is the correct way doing in C++...
    
    auto mu_w  = p.Interpolate(User()->key_muH2O);
    auto mu_o  = p.Interpolate(User()->key_muCO2);

    return (krw(p)/mu_w)/(krw(p)/mu_w+krn(p)/mu_o) ;
    
  }
  
  // 14) the 1s derivative of water fractional flow is calculated from the water relative permeability and oil relative permeability.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::FirstDerivativeWaterFractionalFlow( TARGET_PLACEMENT& p, double64 S){
    
    p.Interpolate(User()->key_sH2O) = S ;  // I am not sure this is the correct way doing in C++...
    auto mu_w  = p.Interpolate(User()->key_muH2O);
    auto mu_o  = p.Interpolate(User()->key_muCO2);
    
    auto Numerator_fw = krw(p)/mu_w ;
    auto Denumerator_fw = (krw(p)/mu_w+krn(p)/mu_o) ;
    
    auto dNumerator_fw = dkrwds(p)/mu_w ;
    auto dDenumerator_fw = dkrwds(p)/mu_w+dkrnds(p)/mu_o ;
    
    
    return ((dNumerator_fw*Denumerator_fw-dDenumerator_fw*Numerator_fw)/(Denumerator_fw*Denumerator_fw));
    
  }
  
  // 15) the oil fractional flow is calculated from the water relative permeability and oil relative permeability.
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::OilFractionalFlow( TARGET_PLACEMENT& p, double64 S){
    
    p.Interpolate(User()->key_sH2O) = S ;  // I am not sure this is the correct way doing in C++...
    auto mu_w  = p.Interpolate(User()->key_muH2O);
    auto mu_o  = p.Interpolate(User()->key_muCO2);
    
    return (krn(p)/mu_o)/(krw(p)/mu_w+krn(p)/mu_o) ;
    
  }
  
  // 15) The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book. This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::InflectionPointSaturation(TARGET_PLACEMENT& p){
    
    auto S = 1-p.Interpolate(User()->key_srCO2) ;
    
    double64 DS ;
    double64 F1,Fold ;
    
    DS = 0.001;
    Fold = -10000;
    
    auto Maxiter(4) ;
    auto it(1) ;
    
    while (it < Maxiter){
      
      F1 = FirstDerivativeWaterFractionalFlow (p, S);
      
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = FirstDerivativeWaterFractionalFlow (p, S);
      }
      
      S = S + 2*DS ;
      DS = DS/10. ;
      it++ ;
    }
    
    S = S-10*DS ;
    
    return S ;
  }
  
  // 15) The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::TangentPointSaturation(TARGET_PLACEMENT& p){
    
    auto Si = InflectionPointSaturation();
    auto Sf = 1-p.Interpolate(User()->key_srCO2) ;
    
    return TheSecantMethod(p,Si,Sf) ;
    
  }
  
  // 16) The Shock front wave calculated after estimation of tangent Saturation point.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::ShockFrontVelocity(TARGET_PLACEMENT& p){
    
    auto S = TangentPointSaturation(p) ;
    
    return FirstDerivativeWaterFractionalFlow(S) ;
    
  }
  
  // 16) The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::BuckleyLeverettFunction( TARGET_PLACEMENT& p, double64 S) {
    
    auto srH2O = p.Interpolate(User()->key_srH2O);

    return FirstDerivativeWaterFractionalFlow(p, S)-(WaterFractionalFlow(p, S)-WaterFractionalFlow(p, srH2O))/(S-srH2O) ;
  }
  
  
  // 17) Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionswithHysteresis<dim,USER>::TheSecantMethod( TARGET_PLACEMENT& p, double64 S1, double64 S2 ){
    
    double it(1) ;
    
    double NewPoint ;
    
    double F1, F2 ;
    
    F1 = 10000 ;
    NewPoint = 0 ;
    
    while (abs(F1)>1e-10) {
      
      F1 = BuckleyLeverettFunction (p, S1);
      F2 = BuckleyLeverettFunction (p, S2);
      
      NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
      
      S2 = S1 ;
      S1 = NewPoint  ;
      
      it++ ;
      
    }
    
    return NewPoint ;
    
  }


}
