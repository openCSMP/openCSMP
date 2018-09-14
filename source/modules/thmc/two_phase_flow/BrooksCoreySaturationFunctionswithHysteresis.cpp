#include "CSMP_mathUtilities.h"
#include "BrooksCoreySaturationFunctionsWithHysteresis.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
/**
   The default constructor of Brooks Corey Saturation Functions With Hysteresis class
*/
template<size_t dim, template<size_t> class USER> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::BrooksCoreySaturationFunctionsWithHysteresis()
  : ac_params_(8)
  {
  }
  
 
 

// INLINE METHODS

/// get seff at the element barycentre
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::EffectiveSaturation( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
  
    double64 seff =  (sH2O - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }





/// get seff from the supplied saturation value
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::EffectiveSaturation_at( const Element<dim>* const e, double64 sw ) const
  {
    double64 seff =  (sw - e->Read(User()->key_srH2O)) /
    (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }





/**
   
  Oil residual saturation has been estimated from Land's formula based on the initial oil saturations.
  
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::OilResidualSaturation( const Element<dim>* const e) const
  {
    const double64 sCO2 = e->PropertyValueAtBaryCenter( User()->key_sCO2 );

    return 1. / ( C_land_ + 1. / sCO2 );
}





/**
 
    Check whether the process is imbibition or drainage. check the change of saturation before and after at barycenter.
 
*/
template<size_t dim, template<size_t> class USER>
TWO_PHASE_FLOW_PROCESS BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::FlowProcess( const Element<dim>* const e )
  {
    const double64 S_old = e->PropertyValueAtBaryCenter( User()->key_sCO2 );
    const double64 S_new = e->PropertyValueAtBaryCenter( User()->key_sCO2_1 );
    
    if (S_old > S_new) return IMBIBITION;
    
    if (S_old < S_new) return DRAINAGE;
   
   return DRAINAGE;
    
} // end FlowProcess


  
  
  
  
  
  
/**
 The main Initialisation of functions...Read the parameters from the input variables file...
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( const Element<dim>* const e )
  {
    e->Read(User()->key_kri_param, ac_params_);
    
    awd_ = ac_params_[AWD];
    aod_ = ac_params_[AOD];
    
    cwd_ = ac_params_[CWD];
    cod_ = ac_params_[COD];
    
    awi_ = ac_params_[AWI];
    aoi_ = ac_params_[AOI];
    
    cwi_ = ac_params_[CWI];
    coi_ = ac_params_[COI];
    
    TWO_PHASE_FLOW_PROCESS ProcessPath = this->FlowProcess(e);
    
    switch (ProcessPath) {
      case DRAINAGE:
        aw_ =  ac_params_[AWD];
        ao_ =  ac_params_[AOD];
        cw_ =  ac_params_[CWD];
        co_ =  ac_params_[COD];
        break;
        
      case IMBIBITION:
        aw_ =  ac_params_[AWI];
        ao_ =  ac_params_[AOI];
        cw_ =  ac_params_[CWI];
        co_ =  ac_params_[COI];
        break;
    }
    
  }
 
  
 
  
  
  
  
  
/**
 
 Second option for Initialisation  of parameters...Read the parameters on the fly and check that Process is Drainage or Imbibitions
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( const Element<dim>* const e,
                                                                                              std::array<std::array<double64, 2>,2> a_,
                                                                                              std::array<std::array<double64, 2>,2> c_)
  {
    TWO_PHASE_FLOW_PROCESS ProcessPath = this->FlowProcess(e);  // Here , the process is checked ... either is Imbibitions or Drainage
    
    aw_ = a_[H2O][ProcessPath] ;
    ao_ = a_[CO2][ProcessPath] ;
    cw_ = c_[H2O][ProcessPath] ;
    co_ = c_[CO2][ProcessPath] ;
    
    awd_ = a_[H2O][DRAINAGE] ;
    aod_ = a_[CO2][DRAINAGE] ;
    
    cwd_ = c_[H2O][DRAINAGE] ;
    cod_ = c_[CO2][DRAINAGE] ;
    
    awi_ = a_[H2O][IMBIBITION] ;
    aoi_ = a_[CO2][IMBIBITION] ;
    
    cwi_ = c_[H2O][IMBIBITION] ;
    coi_ = c_[CO2][IMBIBITION] ;
    
  }


  

  
  
  
 
/**
   
   Third option for Initialisation  of parameters...Read the parameters on the fly and force the Process to be Drainage or Imbibitions
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::InitialiseBrooksCoreyParameters( const Element<dim>* const e,
                    std::array<std::array<double64, 2>,2> a_, std::array<std::array<double64, 2>,2> c_, TWO_PHASE_FLOW_PROCESS ProcessPath)
  {
    aw_ = a_[H2O][ProcessPath] ;
    ao_ = a_[CO2][ProcessPath] ;
    cw_ = c_[H2O][ProcessPath] ;
    co_ = c_[CO2][ProcessPath] ;
    
    awd_ = a_[H2O][DRAINAGE] ;
    aod_ = a_[CO2][DRAINAGE] ;
    
    cwd_ = c_[H2O][DRAINAGE] ;
    cod_ = c_[CO2][DRAINAGE] ;
    
    awi_ = a_[H2O][IMBIBITION] ;
    aoi_ = a_[CO2][IMBIBITION] ;
    
    cwi_ = c_[H2O][IMBIBITION] ;
    coi_ = c_[CO2][IMBIBITION] ;
   }
 

  
  
  
  

  
  

  

/**
   
  Brooks & Corey capillary pressure formula, here on spot the capilary pressure is bound to main drainage and imbibition curves.
  If it is out of bound will be replaced with the Drainage (Or Imbibition) capillary pressure.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::pc( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(User()->key_psrCO2)  ;
    double64 PseudoSwr_ = e->Read(User()->key_psrH2O)  ;  // To have primary Drianage and Imibition parameters
    
    double64 pc(0);
    
    if ((sH2O > PseudoSwr_)&&(sCO2 > PseudoSor_))
      pc = cw_*pow((1. - PseudoSwr_)/(sH2O - PseudoSwr_),aw_) + co_*pow((1.0 - PseudoSor_)/(sCO2 - PseudoSor_),ao_ );
    
    if (pc > MaxCapillaryPressure)  return MaxCapillaryPressure ; // Check if it is not larger than the Maximum Capillary pressure
    if (pc < -MaxCapillaryPressure) return -MaxCapillaryPressure ;// Check if it is not larger than the Minimum Capillary pressure

    
    std::pair<double64,double64> Pc_limits = this->CapillaryPressureLimits(e); // Check if it is in the range of Main Drainage and Imbibition limits
    
    if (pc <= Pc_limits.second) return Pc_limits.second;
    if (pc >= Pc_limits.first)  return Pc_limits.first;
    
    return pc;
}
  
  
  
  
  
/**
   
  The first derivative of Brooks & Corey capillary pressure formula by assuming the variation of CO2 saturation is opposite of H2O saturations
  Here we also check the capillary function bound, and then if it ois out of bound the derivative.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dpcds( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(User()->key_psrCO2)  ;
    double64 PseudoSwr_ = e->Read(User()->key_psrH2O)  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_ );
  
    double64 Dpc(0);
    
    // first calculate the upper limits and lower limits of the capilary curve.
    if (sH2O < PseudoSwr_ )     Dpc = -MaxCapillaryPressureDerivative;
    else if (sCO2 < PseudoSor_) Dpc = -MaxCapillaryPressureDerivative ;
    else
      Dpc = -aw_*cw_*pow((1.0 - PseudoSwr_) / (sH2O - PseudoSwr_), aw_) / (sH2O - PseudoSwr_) + ao_*co_*pow((1.0 - PseudoSor_) / (sCO2 - PseudoSor_), ao_) / (sCO2 - PseudoSor_);
    
    return Dpc ;
}
  
 
  
  
  
  
  
  
  
  
/**
   
   This is residual water saturation for the transit process from drainage to imbibition Water residual saturation has been estimated
   from the intersection of imbibition curve with drainage curve We assume the rock properties are constant, i.e. aw, ao, co, and cw
   are constant. We need to estimate the residual water saturation.
  
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterResidualSaturation( const Element<dim>* const e, double64 Sro) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(User()->key_psrCO2)  ;
    double64 PseudoSwr_ = e->Read(User()->key_psrH2O)  ;  // To have primary Drianage and Imibition parameters
  
    double64 srH2O = PseudoSwr_ ;
    
    
    double64 RHS = cwd_*pow((1.0 - PseudoSwr_)/(sH2O - PseudoSwr_),awd_) + cod_*pow((1.0 - PseudoSor_)/(sCO2 - PseudoSor_),aod_)- coi_*pow((1.0 - Sro)/(sCO2 - Sro),aoi_) ;
    
    RHS = cwi_/RHS ;
    RHS = pow(RHS,1./awi_) ;
    
    if ( (sH2O > PseudoSwr_) && (sH2O < 1.-PseudoSor_)) srH2O = 1.-(1.-sH2O)/(1.-RHS) ;
    
    return  srH2O ;
  }
  
 
  
  

  
  
  
  
/**
  Here , I calculate the residual water saturation and residual oil saturation for the transit process from  imbibition to drainage.
 
  These  residual saturation has been estimated from the intersection of imbibition curve with drianage curve
  We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
  the water residual saturation.
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterAndOilResidualSaturationImbibitionToDrainage( const Element<dim>* const e, double64& Swr_, double64& Sor_) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
  
    double64 Sw1(e->Read(User()->key_SwDrToImb)) ;
    double64 Sw2(e->Read(User()->key_SwImbToDr)) ;
    
    double64 So2(1.-Sw2);
    double64 So1(1.-Sw1);
    
    double64 srH2O(e->Read(User()->key_srH2O)) ;    // This means the algorithm use the current residual water saturation for the next steps
    double64 srCO2(e->Read(User()->key_srCO2)) ;    // This means the algorithm use the current residual oil saturation for the next steps
    
    
    
    Swr_ = e->Read(User()->key_srH2O) ;    // This means the algorithm use the current residual water saturation for the next steps
    Sor_ = e->Read(User()->key_srCO2) ;    // This means the algorithm use the current residual oil saturation for the next steps
    
    
    double64 RHS = coi_*pow((1.0 - srCO2)/(So2 - srCO2),aoi_) + cwi_*pow((1.0 - srH2O)/(Sw2 - srH2O),awi_)- cwd_*pow((1.0 - srH2O)/(Sw2 - srH2O),awd_) ;
    
    RHS = cod_/RHS ;
    RHS = pow(RHS,1./aod_) ;
  
    if ((sH2O > srH2O) && (sH2O < 1.-srCO2))  Sor_ = 1.-(1.-So2)/(1.-RHS) ;   // Mahyar: this way with assumeing the highest saturation of water is more effected by the assumptatic behaviour of oil part of curve.
    
    
    Swr_ = srH2O ;    // This means the algorithm use the current residual water saturation for the next steps
    
    RHS = coi_*pow((1.0 - srCO2)/(So1 - srCO2),aoi_) + cwi_*pow((1.0 - srH2O)/(Sw1 - srH2O),awi_)- cod_*pow((1.0 - srCO2)/(So1 - srCO2),aod_) ;
    
    RHS = cwd_/RHS ;
    RHS = pow(RHS,1./awd_) ;
    
    if ((sH2O > srH2O) && (sH2O < 1.-srCO2))  Swr_ = 1.-(1.-Sw1)/(1.-RHS) ;  // Mahyar: this way with assuming the lowest saturation of water is more effected by the assumptotic behaviour of water part of curve.
}

  
  
  
  
  
  
  
/**
  
  For calculating the water relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
  The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
  the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
 
*/
  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    // Mahyar: This function will check if the capillary pressure stay in the bounds, and if not replace the correct process
    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_ );

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 krwww = pow(Snw, 3.+2.*aw_) ; // refer to Eq. 10 from Skjaeveland et al. 2000
    double64 krwow = (1. - pow(Sno, 2.* ao_+1.)) * square(1. - Sno);
    
    return (cw_*krwww - co_*krwow) / (cw_ - co_) ;  // Mahyar: If sH2O< srH2O or sH2O > 1-srCO2 .. The krw will be between 0 or 1
}
  
  
  
  
  
  
  
  
  
  
  
/**
   
   Second form of the water relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw_at( const Element<dim>* const e, double64 S) const
  {
    double64 sH2O(S) ;
    double64 sCO2(1.-S);
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    // Mahyar: This function will check if the capillary pressure stay in the bounds, and if not replace the correct process
    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_);

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 krwww = pow(Snw, 3.+2.*aw_) ; // refer to Eq. 10 from Skjaeveland et al. 2000
    double64 krwow = (1. - pow(Sno, 2.*ao_+1.)) * square(1. - Sno);
    
    return (cw_*krwww - co_*krwow) / (cw_ - co_) ;  // Mahyar: If sH2O< srH2O or sH2O > 1-srCO2 .. The krw will be between 0 or 1
}



  
  
  
  
  
  
  
/**
   
   For calculating the oil relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
   The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
   the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    
    // Mahyar: This function will check if the capillary pressure stay in the bounds, and if not replace the correct process
    CheckPcLimitsAndResetResiduals( e, PseudoSwr_, PseudoSor_ );

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 kroww = (1. - pow(Snw,2.*aw_+1.)) * square(1. - Snw);
    double64 kroow = pow(Sno, 3.+2.*ao_) ; // refer to Eq. 12 from Skjaeveland et al. 2000
    
    return (cw_*kroww - co_*kroow) / (cw_ - co_) ; // Mahyar: If sCO2< srCO2 or sCO2 > 1-srH2O .. The krn will be between 0 or 1
}





  
/**
   
   Second form of the oil relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn_at( const Element<dim>* const e, double64 S) const
  {
    double64 sH2O = S ;
    double64 sCO2 = 1.-S ;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_);

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 kroww = (1. - pow(Snw,2.*aw_+1.)) * square(1. - Snw);
    double64 kroow = pow(Sno, 3.+2.*ao_) ; // refer to Eq. 12 from Skjaeveland et al. 2000
    
    return (cw_*kroww - co_*kroow) / (cw_ - co_) ;
}




  
  
  
  
  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_);

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkrwww = (3.+2.*aw_)*pow(Snw,(2.+2.*aw_))*dSnw ;
    double64 dkrwow = (-(2*ao_+1)*pow(Sno,(2.*ao_))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*ao_+1.)))*(1.-Sno))*dSno ;
    
    return (cw_*dkrwww-co_*dkrwow)/(cw_-co_) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

  
  
  
  
  
  
  
  
/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_at( const Element<dim>* const e, double64 S) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals( e, PseudoSwr_, PseudoSor_ );

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkrwww = (3.+2.*aw_)*pow(Snw,(2.+2.*aw_))*dSnw ;
    double64 dkrwow = (-(2*ao_+1)*pow(Sno,(2.*ao_))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*ao_+1.)))*(1.-Sno))*dSno ;
    
    return (cw_*dkrwww-co_*dkrwow)/(cw_-co_) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

  
  
  
  
  
  
  
  
  
  
  
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals( e, PseudoSwr_, PseudoSor_ );

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkroww = (-(2*aw_+1)*pow(Snw,(2.*aw_))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*aw_+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*ao_)*pow(Sno,(2.+2.*ao_))*dSno  ;
    
    return (cw_*dkroww-co_*dkroow)/(cw_-co_) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000
    
}
  
  
  
  
  
  
  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_at( const Element<dim>* const e, double64 S ) const
  {
    const double64 sH2O = S;
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(User()->key_psrCO2))  ;
    double64 PseudoSwr_(e->Read(User()->key_psrH2O))  ;  // To have primary Drianage and Imibition parameters
    
    CheckPcLimitsAndResetResiduals( e, PseudoSwr_, PseudoSor_ );

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw; // dSo/(1.-p.Obtain(User()->key_srH2O)-p.Obtain(User()->key_srCO2)) ;
    
    double64 dkroww = (-(2*aw_+1)*pow(Snw,(2.*aw_))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*aw_+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*ao_)*pow(Sno,(2.+2.*ao_))*dSno  ;
    
    return (cw_*dkroww-co_*dkroow)/(cw_-co_) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000
    
}

  
  
  
  
  
  
  
  
  
  

/**
 
 Check The Limits of Capillary presure and reset the residuals
 
 TODO: a function that checks should not change anything
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::CheckPcLimitsAndResetResiduals( const Element<dim>* const e, double64& newSrH2O , double64& newSrCO2 ) const
  {
    const double64 srH2O(e->Read(User()->key_srH2O));  // To have primary Drianage and Imibition parameters
    const double64 srCO2(e->Read(User()->key_srCO2));
    
    std::pair<double64,double64> Pc_limits = CapillaryPressureLimits(e);
    
    double64 local_pc = pc(e);
    
    //check if hit prmary drainage
    if (local_pc >= Pc_limits.second ) {
      
      newSrCO2 = srCO2 ;
      newSrH2O = srH2O ;
      
      aw_ = awd_ ;
      cw_ = cwd_ ;
      ao_ = aod_ ;
      co_ = cod_ ;
    }
    
    //check if hit prmary Imbibition
    if (local_pc <= Pc_limits.first ) {
      
      newSrCO2 = srCO2 ;
      newSrH2O = srH2O ;
      
      aw_ = awi_ ;
      cw_ = cwi_ ;
      ao_ = aoi_ ;
      co_ = coi_ ;
    }
}
  
  
  
  
  
  
  
  
  
  
  
/**
   
   calculate the upper and lowwer limits of Capillary presure
   
*/
template<size_t dim, template<size_t> class USER>
  std::pair<double64,double64> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::CapillaryPressureLimits( const Element<dim>* const e ) const
  {
    std::pair<double64,double64> pc_PrimaryImbibitionDrainage;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 srCO2(e->Read(User()->key_srCO2))  ;
    double64 srH2O(e->Read(User()->key_srH2O))  ;  // To have primary Drianage and Imibition parameters
    
    // first calculate the upper limits and lower limits of the capilary curve.
    if (sH2O < srH2O){
      pc_PrimaryImbibitionDrainage.first  = MaxCapillaryPressure ;
      pc_PrimaryImbibitionDrainage.second = MaxCapillaryPressure ;
    }else if (sCO2 < srCO2) {
      pc_PrimaryImbibitionDrainage.first  = -MaxCapillaryPressure ;
      pc_PrimaryImbibitionDrainage.second = -MaxCapillaryPressure ;
    }else{
      
      pc_PrimaryImbibitionDrainage.first  = cwd_*pow((1.0 - srH2O)/(sH2O - srH2O),awd_) + cod_*pow((1.0 - srCO2)/(sCO2 - srCO2),aod_);
      pc_PrimaryImbibitionDrainage.second = cwi_*pow((1.0 - srH2O)/(sH2O - srH2O),awi_) + coi_*pow((1.0 - srCO2)/(sCO2 - srCO2),aoi_);
      
    }
    
    return pc_PrimaryImbibitionDrainage;
  }

  
  
  
  
  
  
  
  
  
  // Numerical derivative added
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }



 
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  

template class BrooksCoreySaturationFunctionsWithHysteresis<1U,CO2H2O_FunctionsModule1>;
template class BrooksCoreySaturationFunctionsWithHysteresis<2U,CO2H2O_FunctionsModule1>;
template class BrooksCoreySaturationFunctionsWithHysteresis<3U,CO2H2O_FunctionsModule1>;

} // csmp
