#include <cmath>

#include "BrooksCoreySaturationFunctionsWithHysteresis.h"
#include "FlowFunctionsDraft2.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

// 1) destructor of Brooks Corey Pressure feilds
template<size_t dim, template<size_t> class USER> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::BrooksCoreySaturationFunctionsWithHysteresis()
{
}

// 3) a  constructor of Brooks Corey Pressure feilds
template<size_t dim, template<size_t> class USER>
BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::BrooksCoreySaturationFunctionsWithHysteresis(double64 aw, double64 ao,  double64 cw, double64 co) :
aw_(aw),ao_(ao),cw_(cw),co_(co)
{
}
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& p )
  {
    
    this->Swmin = p.Obtain(User()->key_srH2O) ;
    this->Swmax = 1.0-p.Obtain(User()->key_srCO2);
    
    this->aw_ = p.Obtain(User()->key_aw) ;  // Mahyar: please change key_bcp to key_aw
    this->ao_ = p.Obtain(User()->key_ao) ;  // Mahyar: please change key_bcp to key_ao
    
    this->cw_ = p.Obtain(User()->key_cw) ;   // Mahyar: please change key_bcp to key_cw
    this->co_ = p.Obtain(User()->key_co) ;   // Mahyar: please change key_bcp to key_co
    
  }
  
  //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteElementPlacement<1U,ELEMENT>& ) ;
  //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteVolumePlacement<1U,NODE>& ) ;

  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& p, double64 aw, double64 ao,  double64 cw, double64 co)
  {
    
    this->Swmin = p.Obtain(User()->key_srH2O) ;
    this->Swmax = 1.0-p.Obtain(User()->key_srCO2);
    
    this->aw_ = aw ;
    this->ao_ = ao ;
    
    this->cw_ = cw ;
    this->co_ = co ;
    
  }
    
  //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteElementPlacement<1U,ELEMENT>&, double64, double64,  double64, double64 ) ;
  //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteVolumePlacement<1U,NODE>&, double64, double64,  double64, double64 ) ;
    
    
    
    template<size_t dim, template<size_t> class USER>
    template<class TARGET_PLACEMENT>
    void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& p, double64 a[2][2], double64 c[2][2], process ProcessPath)
    {
        
        this->Swmin = p.Obtain(User()->key_srH2O) ;
        this->Swmax = 1.0-p.Obtain(User()->key_srCO2);
        
        this->aw_ = a[H2O][ProcessPath] ;
        this->ao_ = a[CO2][ProcessPath] ;
        
        this->cw_ = c[H2O][ProcessPath] ;
        this->co_ = c[CO2][ProcessPath] ;
        
    }
    
    //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteElementPlacement<1U,ELEMENT>&, double64 a[2][2] , double64 c[2][2], process ProcessPath) ;
    //template void FlowFunctionsDraft2<1U>::InitialiseBrooksCoreyParameters( FiniteVolumePlacement<1U,NODE>&, double64 a[2][2], double64 c[2][2], process ProcessPath) ;
    
    
    template<size_t dim, template<size_t> class USER>
    template<class TARGET_PLACEMENT>
    double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::EffectiveSaturation( TARGET_PLACEMENT& p ) const
    {
        // seff = (sw - swr) / (1 - swr - snr)
        double64 seff =  (p.Obtain(User()->key_sH2O) - p.Obtain(User()->key_srH2O)) /
        (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
        
        return std::min( std::max( seff, 0. ), 1. );
    }
    
    //template double64 FlowFunctionsDraft2<1U>::EffectiveSaturation( FiniteElementPlacement<1U,ELEMENT>& ) const ;
    //template double64 FlowFunctionsDraft2<1U>::EffectiveSaturation( FiniteVolumePlacement<1U,NODE>& ) const ;


// 4) Brooks & Corey capillary pressure formula
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::pc( TARGET_PLACEMENT& p ) const
  {
   
    auto sH2O = p.Obtain(User()->key_sH2O);
    auto sCO2 = p.Obtain(User()->key_sCO2);
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    return cw_*pow((1.0 - srH2O)/(sH2O - srH2O),aw_) +
           co_*pow((1.0 - srCO2)/(sCO2 - srCO2),ao_);
}

template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::pc( FiniteVolumePlacement<1U,NODE>& ) const ;
template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::pc( FiniteElementPlacement<1U,ELEMENT>& ) const ;
  
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
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dpcds( TARGET_PLACEMENT& p) const
  {
   
   auto sH2O = p.Obtain(User()->key_sH2O);
   auto srH2O = p.Obtain(User()->key_srH2O);
   auto sCO2 = p.Obtain(User()->key_sCO2);
   auto srCO2 = p.Obtain(User()->key_srCO2);
    
   return -aw_*cw_*pow((1.0 - srH2O) / (sH2O - srH2O), aw_) / (sH2O - srH2O) +
           ao_*co_*pow((1.0 - srCO2) / (sCO2 - srCO2), ao_) / (sCO2 - srCO2);
}
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dpcds( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dpcds( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dpcds( FiniteElementPlacement<3U,ELEMENT>& ) const;


// 7) Oil residual saturation has been estimated from Land's formula based on the initial oil saturations.

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::OilResidualSaturation( TARGET_PLACEMENT& p) const
  {
  return 1./(C_land_+1./p.Obtain(User()->key_sCO2)) ;
}
    
template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::OilResidualSaturation( FiniteVolumePlacement<1U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::OilResidualSaturation( FiniteElementPlacement<1U,ELEMENT>& ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::OilResidualSaturation( FiniteElementPlacement<1U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::OilResidualSaturation( FiniteElementPlacement<2U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::OilResidualSaturation( FiniteElementPlacement<3U,NODE>& ) const;
  */
    
    // 24/July/2018
    // Mahyar: This is residual water saturation for the transit process from drainage to imbibition
    //
// Water residual saturation has been estimated from the intersection of imbibition curve with drianage curve
// We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
// the water residual saturation.

    template<size_t dim, template<size_t> class USER>
    template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterResidualSaturation( TARGET_PLACEMENT& p, double64 a[2][2], double64 c[2][2]) const
    {
        
        auto sH2O = p.Obtain(User()->key_sH2O);
        auto srH2O = p.Obtain(User()->key_srH2O);
        auto sCO2 = p.Obtain(User()->key_sCO2);
        auto srCO2 = p.Obtain(User()->key_srCO2);
        
        double64 awi_ = a[H2O][Imbibition] ;
        double64 aoi_ = a[CO2][Imbibition] ;
        
        double64 awd_ = a[H2O][Drainage] ;
        double64 aod_ = a[CO2][Drainage] ;
        
        double64 cwi_ = c[H2O][Imbibition] ;
        double64 coi_ = c[CO2][Imbibition] ;
        
        double64 cwd_ = c[H2O][Drainage] ;
        double64 cod_ = c[CO2][Drainage] ;
        
        double64 RHS = cwd_*pow((1.0 - srH2O)/(sH2O - srH2O),awd_) + cod_*pow((1.0 - srCO2)/(sCO2 - srCO2),aod_)- coi_*pow((1.0 - srCO2)/(sCO2 - srCO2),aoi_) ;
        
        RHS = cwi_/RHS ;
        RHS = pow(RHS,1./awi_) ;
        
        return  1.-(1.-sH2O)/(1.-RHS) ;
    }
    
    template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::WaterResidualSaturation( FiniteVolumePlacement<1U,NODE>& , double64 a[2][2] , double64 c[2][2]) const;
    template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::WaterResidualSaturation( FiniteElementPlacement<1U,ELEMENT>&, double64 a[2][2] , double64 c[2][2] ) const;

    
    // 24/July/2018
    // Mahyar: Here , I calculate the residual water saturation  and residual oil saturation for the transit process from  imbibition to drainage
    //
    // These  residual saturation has been estimated from the intersection of imbibition curve with drianage curve
    // We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
    // the water residual saturation.
    
    template<size_t dim, template<size_t> class USER>
    template<class TARGET_PLACEMENT> void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterAndOilResidualSaturationImbibitionToDrainage( TARGET_PLACEMENT& p, double64 Swr_, double64 Sor_)
    {
        
        auto sH2O = p.Obtain(User()->key_sH2O);
        auto srH2O = p.Obtain(User()->key_srH2O);
        auto sCO2 = p.Obtain(User()->key_sCO2);
        auto srCO2 = p.Obtain(User()->key_srCO2);
        
        double64 awi_, awd_, aoi_, aod_;
        double64 cwi_, cwd_, coi_, cod_ ;
        
        Swr_ = srH2O ;    // This means the algorithm use the current residual water saturation for the next steps
        Sor_ = srCO2 ;    // This means the algorithm use the current residual oil saturation for the next steps
        
        double64 Swt_,  Sot_ ;  // target saturation of oil and water
        
        double64 F1, F2;
        double64 dF1dSwr, dF1dSor;
        double64 dF2dSwr, dF2dSor ;
        
        double64 deltaJ(1e+3) ;
        
        
        while (abs(deltaJ)>1e-10) {   // Newton Raphson for finding residual water and oil trapped in the model.
            
            F1 = cwd_*pow((1.0 - Swr_)/(sH2O - Swr_),awd_) + cod_*pow((1.0 - Sor_)/(sCO2 - Sor_),aod_) -
                 cwi_*pow((1.0 - srH2O)/(sH2O - srH2O),awi_) - coi_*pow((1.0 - srCO2)/(sCO2 - srCO2),aoi_) ;  // the cost function whihc has two root (Sw_)
        
            F2 = cwd_*pow((1.0 - Swr_)/(Swt_ - Swr_),awd_) + cod_*pow((1.0 - Sor_)/(Sot_ - Sor_),aod_) -
                 cwi_*pow((1.0 - srH2O)/(Swt_ - srH2O),awi_) - coi_*pow((1.0 - srCO2)/(Sot_ - srCO2),aoi_) ;  // the cost function whihc has two root (Sw_)
        
            dF1dSwr =  awd_*cwd_*pow((1.-Swr_)/(sH2O-Swr_),awd_)/((sH2O-Swr_)*(1.-Swr_))*(1.-sH2O) ;   // derivative cost function to Swr
            dF1dSor =  aod_*cod_*pow((1.-Sor_)/(sCO2-Sor_),aod_)/((sCO2-Sor_)*(1.-Sor_))*(1.-sCO2) ;
        
            dF2dSwr =  awd_*cwd_*pow((1.-Swr_)/(Swt_-Swr_),awd_)/((Swt_-Swr_)*(1.-Swr_))*(1.-Swt_) ;
            dF2dSor =  aod_*cod_*pow((1.-Sor_)/(Sot_-Sor_),aod_)/((Sot_-Sor_)*(1.-Sor_))*(1.-Sot_) ;
        
            deltaJ  = dF1dSwr*dF2dSor - dF1dSor*dF2dSwr ;
            
            Swr_ -= (dF2dSor*Swr_ - dF2dSwr*Sor_)/deltaJ ;
            Sor_ -= (-dF1dSor*Swr_ + dF1dSwr*Sor_)/deltaJ ;
         
         
         }
    }
    
    
    
    // 24/July/2018
    // Mahyar: Here , I calculate the residual water saturation  and residual oil saturation for the transit process from drainage to imbibitions
    //
    // These  residual saturation has been estimated from the intersection of imbibition curve with drianage curve
    // We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
    // the water residual saturation.
    
    template<size_t dim, template<size_t> class USER>
    template<class TARGET_PLACEMENT> void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterAndOilResidualSaturationDrainageToImbibition( TARGET_PLACEMENT& p, double64 Swr_, double64 Sor_)
    {
        
        auto sH2O = p.Obtain(User()->key_sH2O);
        auto srH2O = p.Obtain(User()->key_srH2O);
        auto sCO2 = p.Obtain(User()->key_sCO2);
        auto srCO2 = p.Obtain(User()->key_srCO2);
        
        double64 awi_, awd_, aoi_, aod_;
        double64 cwi_, cwd_, coi_, cod_ ;
        
        Swr_ = srH2O ;    // This means the algorithm use the current residual water saturation for the next steps
        Sor_ = srCO2 ;    // This means the algorithm use the current residual oil saturation for the next steps
        
        double64 Swt_,  Sot_ ;  // target saturation of oil and water
        
        double64 F1, F2;
        double64 dF1dSwr, dF1dSor;
        double64 dF2dSwr, dF2dSor ;
        
        double64 deltaJ(1e+3) ;
        
        
        while (abs(deltaJ)>1e-10) {   // Newton Raphson for finding residual water and oil trapped in the model.
            
            F1 = cwi_*pow((1.0 - Swr_)/(sH2O - Swr_),awi_) + coi_*pow((1.0 - Sor_)/(sCO2 - Sor_),aoi_) -
                 cwd_*pow((1.0 - srH2O)/(sH2O - srH2O),awd_) - cod_*pow((1.0 - srCO2)/(sCO2 - srCO2),aod_) ;  // the cost function whihc has two root (Sw_)
            
            F2 = cwi_*pow((1.0 - Swr_)/(Swt_ - Swr_),awi_) + coi_*pow((1.0 - Sor_)/(Sot_ - Sor_),aoi_) -
                 cwd_*pow((1.0 - srH2O)/(Swt_ - srH2O),awd_) - cod_*pow((1.0 - srCO2)/(Sot_ - srCO2),aod_) ;  // the cost function whihc has two root (Sw_)
            
            dF1dSwr =  awi_*cwi_*pow((1.-Swr_)/(sH2O-Swr_),awi_)/((sH2O-Swr_)*(1.-Swr_))*(1.-sH2O) ;   // derivative cost function to Swr
            dF1dSor =  aoi_*coi_*pow((1.-Sor_)/(sCO2-Sor_),aoi_)/((sCO2-Sor_)*(1.-Sor_))*(1.-sCO2) ;
            
            dF2dSwr =  awi_*cwi_*pow((1.-Swr_)/(Swt_-Swr_),awi_)/((Swt_-Swr_)*(1.-Swr_))*(1.-Swt_) ;
            dF2dSor =  aoi_*coi_*pow((1.-Sor_)/(Sot_-Sor_),aoi_)/((Sot_-Sor_)*(1.-Sor_))*(1.-Sot_) ;
            
            deltaJ  = dF1dSwr*dF2dSor - dF1dSor*dF2dSwr ;
            
            Swr_ -= (dF2dSor*Swr_ - dF2dSwr*Sor_)/deltaJ ;
            Sor_ -= (-dF1dSor*Swr_ + dF1dSwr*Sor_)/deltaJ ;
            
            
        }
    }
    

// 8) For calculating the relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
// The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
// the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
//

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw( TARGET_PLACEMENT& p ) const
  {
  
  auto sH2O = p.Obtain(User()->key_sH2O);
  auto srH2O = p.Obtain(User()->key_srH2O);
  auto sCO2 = p.Obtain(User()->key_sCO2);
  auto srCO2 = p.Obtain(User()->key_srCO2);

  double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
  double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
  Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
  Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    

  double64 krwww = pow(Snw, 3.+2.*aw_) ; // refer to Eq. 10 from Skjaeveland et al. 2000
  double64 krwow = (1. - pow(Sno, 2.*ao_+1.)) * square(1. - Sno);

  return (cw_*krwww - co_*krwow) / (cw_ - co_) ;
}
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<3U,ELEMENT>& ) const;


  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<1U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<2U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<3U,NODE>& ) const;
  */
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw( TARGET_PLACEMENT& p , double64 S) const
  {
    
    auto sH2O = S ;
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2 = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 krwww = pow(Snw, 3.+2.*aw_) ; // refer to Eq. 10 from Skjaeveland et al. 2000
    double64 krwow = (1. - pow(Sno, 2.*ao_+1.)) * square(1. - Sno);
    
    return (cw_*krwww - co_*krwow) / (cw_ - co_) ;
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<1U,NODE>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<2U,NODE>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krw( FiniteElementPlacement<3U,NODE>&, double64 ) const;
  */
  

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn( TARGET_PLACEMENT& p) const
  {
  
  auto sH2O = p.Obtain(User()->key_sH2O);
  auto srH2O = p.Obtain(User()->key_srH2O);
  auto sCO2 = p.Obtain(User()->key_sCO2);
  auto srCO2 = p.Obtain(User()->key_srCO2);

  double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
  double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
  Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
  Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
  double64 kroww = (1. - pow(Snw,2.*aw_+1.)) * square(1. - Snw);
  double64 kroow = pow(Sno, 3.+2.*ao_) ; // refer to Eq. 12 from Skjaeveland et al. 2000

  return (cw_*kroww - co_*kroow) / (cw_ - co_) ;
}
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<3U,ELEMENT>& ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<1U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<2U,NODE>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<3U,NODE>& ) const;
  */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn( TARGET_PLACEMENT& p, double64 S) const
  {
    
    auto sH2O = S ;
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2 = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 kroww = (1. - pow(Snw,2.*aw_+1.)) * square(1. - Snw);
    double64 kroow = pow(Sno, 3.+2.*ao_) ; // refer to Eq. 12 from Skjaeveland et al. 2000
    
    return (cw_*kroww - co_*kroow) / (cw_ - co_) ;
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

  
  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<1U,NODE>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<2U,NODE>&, double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::krn( FiniteElementPlacement<3U,NODE>&, double64 ) const;
  */

// 9) For calculating the 1st derivative of relative permeability

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds( TARGET_PLACEMENT& p ) const
  {
    
    auto sH2O  = p.Obtain(User()->key_sH2O);
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2  = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);

    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1

    double64 dSo = -1.; // if change of water saturation will change co2 saturations.

    double64 dSnw = 1./(1.-srH2O-srCO2) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;

    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;

    double64 dkrwww = (3.+2.*aw_)*pow(Snw,(2.+2.*aw_))*dSnw ;
    double64 dkrwow = (-(2*ao_+1)*pow(Sno,(2.*ao_))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*ao_+1.)))*(1.-Sno))*dSno ;

    return (cw_*dkrwww-co_*dkrwow)/(cw_-co_) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<3U,ELEMENT>& ) const;

  
  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<3U,ELEMENT>& ) const;
  */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds( TARGET_PLACEMENT& p , double64 S) const
  {
    
    auto sH2O  = S ;
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2  = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-srH2O-srCO2) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkrwww = (3.+2.*aw_)*pow(Snw,(2.+2.*aw_))*dSnw ;
    double64 dkrwow = (-(2*ao_+1)*pow(Sno,(2.*ao_))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*ao_+1.)))*(1.-Sno))*dSno ;
    
    return (cw_*dkrwww-co_*dkrwow)/(cw_-co_) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrwds( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  */

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds( TARGET_PLACEMENT& p) const
  {

    auto sH2O  = p.Obtain(User()->key_sH2O);
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2  = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-srH2O-srCO2) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;

    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkroww = (-(2*aw_+1)*pow(Snw,(2.*aw_))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*aw_+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*ao_)*pow(Sno,(2.+2.*ao_))*dSno  ;

    return (cw_*dkroww-co_*dkroow)/(cw_-co_) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000

}
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<3U,ELEMENT>& ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<3U,ELEMENT>& ) const;
  */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds( TARGET_PLACEMENT& p, double64 S ) const
  {
    
    auto sH2O  = S ;
    auto srH2O = p.Obtain(User()->key_srH2O);
    auto sCO2  = p.Obtain(User()->key_sCO2);
    auto srCO2 = p.Obtain(User()->key_srCO2);
    
    double64 Snw = (sH2O - srH2O) / (1. - srH2O - srCO2);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - srCO2) / (1. - srH2O - srCO2);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-srH2O-srCO2) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkroww = (-(2*aw_+1)*pow(Snw,(2.*aw_))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*aw_+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*ao_)*pow(Sno,(2.+2.*ao_))*dSno  ;
    
    return (cw_*dkroww-co_*dkroow)/(cw_-co_) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::dkrnds( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  */
  /*
  // 13) the water fractional flow is calculated from the water relative permeability and oil relative permeability.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterFractionalFlow( TARGET_PLACEMENT& p, double64 S) const
  {
    
    auto mu_w  = p.Obtain(User()->key_muH2O);
    auto mu_o  = p.Obtain(User()->key_muCO2);

    return (krw(p, S)/mu_w)/(krw(p, S)/mu_w+krn(p, S)/mu_o) ;
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::WaterFractionalFlow( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::WaterFractionalFlow( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::WaterFractionalFlow( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  
  
  // 14) the 1s derivative of water fractional flow is calculated from the water relative permeability and oil relative permeability.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::FirstDerivativeWaterFractionalFlow( TARGET_PLACEMENT& p, double64 S) const
  {
    
    auto mu_w  = p.Obtain(User()->key_muH2O);
    auto mu_o  = p.Obtain(User()->key_muCO2);
    
    auto Numerator_fw = krw(p, S)/mu_w ;
    auto Denumerator_fw = (krw(p, S)/mu_w+krn(p, S)/mu_o) ;
    
    auto dNumerator_fw = dkrwds(p, S)/mu_w ;
    auto dDenumerator_fw = dkrwds(p, S)/mu_w+dkrnds(p, S)/mu_o ;
    
    
    return ((dNumerator_fw*Denumerator_fw-dDenumerator_fw*Numerator_fw)/(Denumerator_fw*Denumerator_fw));
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::FirstDerivativeWaterFractionalFlow( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::FirstDerivativeWaterFractionalFlow( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::FirstDerivativeWaterFractionalFlow( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  
  
  // 15) the oil fractional flow is calculated from the water relative permeability and oil relative permeability.
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::OilFractionalFlow( TARGET_PLACEMENT& p, double64 S) const
  {
    
    auto mu_w  = p.Obtain(User()->key_muH2O);
    auto mu_o  = p.Obtain(User()->key_muCO2);
    
    return (krn(p, S)/mu_o)/(krw(p, S)/mu_w+krn(p, S)/mu_o) ;
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::OilFractionalFlow( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::OilFractionalFlow( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::OilFractionalFlow( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  
  */
  // 15) The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book. This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InflectionPointSaturation(TARGET_PLACEMENT& p) const
  {
    
    auto S = 1.-p.Obtain(User()->key_srCO2) ;
    
    double64 DS ;
    double64 F1,Fold ;
    
    DS = 0.001;
    Fold = -10000;
    
    auto Maxiter(4) ;
    auto it(1) ;
    
    while (it < Maxiter){
      
      F1 = User()->dfds(p, 0U, S);
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = User()->dfds(p, 0U, S);
        
        if((S<Swmin)||(S>Swmax)) return S=-1;
        
      }
      
      S = S + 2*DS ;
      DS = DS/10. ;
      it++ ;
      
      if((S<Swmin)||(S>Swmax)) return S=-1;

    }
    
    S = S-10*DS ;
    
    if((S<Swmin)||(S>Swmax)) return S=-1;

    
    return S ;
  }
  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::InflectionPointSaturation( FiniteElementPlacement<1U,NODE>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::InflectionPointSaturation( FiniteElementPlacement<2U,NODE>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::InflectionPointSaturation( FiniteElementPlacement<3U,NODE>&  ) const;
  */
  
  // 15) The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::TangentPointSaturation(TARGET_PLACEMENT& p) const
  {
    
    auto Si = this->InflectionPointSaturation(p);
    auto Sf = 1-p.Obtain(User()->key_srCO2) ;
    
    auto St=TheSecantMethod(p,Si,Sf) ;
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::TangentPointSaturation( FiniteElementPlacement<1U,ELEMENT>&  ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::TangentPointSaturation( FiniteElementPlacement<1U,NODE>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::TangentPointSaturation( FiniteElementPlacement<2U,NODE>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::TangentPointSaturation( FiniteElementPlacement<3U,NODE>&  ) const;
  */
  
  // 16) The Shock front wave calculated after estimation of tangent Saturation point.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::ShockFrontVelocity(TARGET_PLACEMENT& p) const
  {
    
    auto S = this->TangentPointSaturation(p) ;
    
    S = std::min( std::max( S, 0. ), 1. );
    
   // return FirstDerivativeWaterFractionalFlow(p, S) ;
    return User()->dfds(p, 0U, S);
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<1U,ELEMENT>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<2U,ELEMENT>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<3U,ELEMENT>&  ) const;

  
  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<1U,ELEMENT>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<2U,ELEMENT>&  ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::ShockFrontVelocity( FiniteElementPlacement<3U,ELEMENT>&  ) const;
  */
  
  
  // 16) The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::BuckleyLeverettFunction( TARGET_PLACEMENT& p, double64 S) const
  {
    

    auto srH2O = p.Obtain(User()->key_srH2O);

   // return FirstDerivativeWaterFractionalFlow(p, S)-(WaterFractionalFlow(p, S)-WaterFractionalFlow(p, srH2O))/(S-srH2O) ;
    return User()->dfds(p, 0U, S)-(User()->f(p, 0U, S)-User()->f(p, 0U, srH2O))/(S-srH2O);
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::BuckleyLeverettFunction( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::BuckleyLeverettFunction( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::BuckleyLeverettFunction( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::BuckleyLeverettFunction( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  */
  
  
  // 17) Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::TheSecantMethod( TARGET_PLACEMENT& p, double64 S1, double64 S2 ) const
  {
    
    double64 it(1.0) ;
    
    double64 NewPoint ;
    
    double64 F1, F2 ;
    
    F1 = 10000 ;
    NewPoint = 0 ;
    
    while (abs(F1)>1e-10) {
      
      F1 = this->BuckleyLeverettFunction(p, S1);
      F2 = this->BuckleyLeverettFunction(p, S2);
      
      NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
      
      S2 = S1 ;
      S1 = NewPoint  ;
      
     // if((S1<Swmin)||(S1>Swmax)) return NewPoint=-1;

      
      it++ ;
      
    }
    return  NewPoint;
    
  }
  
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::TheSecantMethod( FiniteElementPlacement<1U,ELEMENT>& , double64 , double64 ) const;

  /*
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>::TheSecantMethod( FiniteElementPlacement<1U,NODE>& , double64 , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>::TheSecantMethod( FiniteElementPlacement<2U,NODE>& , double64 , double64 ) const;
  template double64 BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>::TheSecantMethod( FiniteElementPlacement<3U,NODE>& , double64 , double64 ) const;
  */
  template class BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsDraft2>;
  template class BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsDraft2>;
  template class BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsDraft2>;


}
