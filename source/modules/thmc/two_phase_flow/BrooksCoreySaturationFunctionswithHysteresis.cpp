#include "CSMP_mathUtilities.h"
#include "BrooksCoreySaturationFunctionswithHysteresis.h"
#include "FlowFunctionsModule.h"
#include "Element.h"
#include "Model.h"


using namespace std;

namespace csmp {
  
template<size_t dim, template<size_t> class USER>
BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::BrooksCoreySaturationFunctionsWithHysteresis( PropertyDatabase<dim>& pref )

  : key_SwImbToDr_(pref.StorageKey("previous imbibition endpoint")),
    key_SwDrToImb_(pref.StorageKey("previous drainage endpoint")),
    key_prsH2O_(pref.StorageKey("pseudo residual saturation aqueous phase")),
    key_prsCO2_(pref.StorageKey("pseudo residual saturation carbonic phase"))
 {
   /* in case these will have to be created dynamically
   
     model.CreateProperty( "previous imbibition endpoint", "none", SCALAR, ELEMENT);
     model.CreateProperty( "previous drainage endpoint", "none", SCALAR, ELEMENT);
     model.CreateProperty( "pseudo residual water saturation", "none", SCALAR, ELEMENT);
     model.CreateProperty( "pseudo residual CO2 saturation", "none", SCALAR, ELEMENT);
   */
 }
  
  
  
  
  
  
  
/// get seff at the element barycentre
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::EffectiveSaturation( Element<dim>* const e ) const
  {

    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;
    
    double64 seff =  (sH2O - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }





  
  
  
  
  
  
  
  
/// get seff from the supplied saturation value
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::EffectiveSaturation_at( Element<dim>* const e, double64 sw ) const
  {
    double64 seff =  (sw - e->Read(User()->key_srH2O)) /
    (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    
    return std::min( std::max( seff, 0. ), 1. );
  }





  
  
  
  
  
/**
   
  Oil residual saturation has been estimated from Land's formula based on the initial oil saturations.
  
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::OilResidualSaturation( Element<dim>* const e) const
  {

    const double64 sCO2 = e->PropertyValueAtBaryCenter( User()->key_sCO2 );
    assert(sCO2 >= 0 and sCO2 <=1) ;

    return 1. / ( C_land_ + 1. / sCO2 );
}




  
  
  
  
  

/**
 
    Check whether the process is imbibition or drainage. check the change of saturation before and after at barycenter.
 
*/
template<size_t dim, template<size_t> class USER>
typename BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::TWO_PHASE_FLOW_PROCESS
BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::FlowProcess( const Element<dim>* const e ) const
  {
    const double64 S_old = e->PropertyValueAtBaryCenter( User()->key_sCO2_0 );
    assert(S_old >= 0 and S_old <=1) ;

    const double64 S_new = e->PropertyValueAtBaryCenter( User()->key_sCO2 );
    assert(S_new >= 0 and S_new <=1) ;

    
    if (S_old > S_new) return IMBIBITION;
    if (S_old < S_new) return DRAINAGE;
   
   return DRAINAGE;
    
} // end FlowProcess


  
  
  
  
  
  
  
  /**
   
   This function set the parameters for the Imbibitions of Drainage curves
   
   */
  template<size_t dim, template<size_t> class USER>
  void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::SetBrooksCoreyCurvesParameters( const Element<dim>* const e, std::array<double64, 2>& a,std::array<double64, 2>& c ) const {
  
  
  e->Read(User()->key_kri_param,ac_params_) ;
  
  TWO_PHASE_FLOW_PROCESS ProcessPath = FlowProcess(e);
  
  switch (ProcessPath) {
    case DRAINAGE:
      a[H2O] =  ac_params_[AWD];
      a[CO2] =  ac_params_[AOD];
      c[H2O] =  ac_params_[CWD];
      c[CO2] =  ac_params_[COD];
      break;
      
    case IMBIBITION:
      a[H2O] =  ac_params_[AWI];
      a[CO2] =  ac_params_[AOI];
      c[H2O] =  ac_params_[CWI];
      c[CO2] =  ac_params_[COI];
      break;

  }
  
  }
  
  
  
  
 
  
  
  
  
  
  
/**
   
   This function update the pseudo residual and end point saturations
   
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::UpdatePseudoResidualAndEndpointSaturations( Element<dim>* e ) const
{ 
      const double64 dsCO2_old = e->Read(User()->key_dsCO2);
  
      const double64 S_old = e->PropertyValueAtBaryCenter( User()->key_sCO2_0 );
      assert(S_old >= 0 and S_old <=1) ;
  
      const double64 S_new = e->PropertyValueAtBaryCenter( User()->key_sCO2 );
      assert(S_new >= 0 and S_new <=1) ;

      const double64 dsCO2_new = S_new - S_old;
      e->Store( User()->key_dsCO2, makeScalar( e->Status( User()->key_dsCO2), dsCO2_new ));
    
      const double64 oldSro =  e->Read(User()->key_srCO2) ;
      const double64 oldSrw =  e->Read(User()->key_srH2O) ;
  
      if (isnan(e->Read(key_prsCO2_)))    e->Store(key_prsCO2_, makeScalar( e->Status(key_prsCO2_), oldSro ));
      if (isnan(e->Read(key_prsH2O_)))    e->Store(key_prsH2O_, makeScalar( e->Status(key_prsH2O_), oldSrw ));
  
      if (isnan(e->Read(key_SwDrToImb_))) e->Store(key_SwDrToImb_, makeScalar( ANY , 0 ));
      if (isnan(e->Read(key_SwImbToDr_))) e->Store(key_SwImbToDr_, makeScalar( ANY , 1 ));
  
      if (dsCO2_old > 0. && dsCO2_new < 0.)   {
        e->Store(key_SwDrToImb_, makeScalar( e->Status(key_SwDrToImb_), e->PropertyValueAtBaryCenter( User()->key_sH2O ) ));
      }
      else if ( dsCO2_old < 0. && dsCO2_new > 0. ) {
        e->Store(key_SwImbToDr_, makeScalar( e->Status(key_SwImbToDr_), e->PropertyValueAtBaryCenter( User()->key_sH2O ) ));
      }
  
      //do nothing if saturation has not changed
      if (dsCO2_new == 0.) return;
      
      //update parameter values if required
      double64 Sw     = e->PropertyValueAtBaryCenter( User()->key_sH2O );      //water saturation
      assert(Sw >= 0     and Sw <=1) ;

      double64 Sw_min = e->Read( key_SwDrToImb_ ); //previous drainage endpoint
      assert(Sw_min >= 0 and Sw_min <=1) ;

      double64 Sw_max = e->Read( key_SwImbToDr_ ); //previous imbibition endpoint
      assert(Sw_max >= 0 and Sw_max <=1) ;

      double64 newSro(0); //new pseudo residual saturation carbonic phase
      double64 newSrw(0); //new pseudo residual saturation aqueous phase      
      double64 tol_(0.000001) ;
  
      bool flag_(false) ;      
      
      if (dsCO2_old > 0. && dsCO2_new < 0.) {//turning from DRAINAGE to IMBIBITION 
          if(Sw_min==0.)
            flag_ = (Sw < tol_*0.1);
          else          
            flag_ = (abs(Sw-Sw_min)/Sw_min < tol_);
          if ( flag_ ) {
            //compute and store new pseudo residual saturations
            newSro = OilResidualSaturation(e); 
            newSrw = WaterResidualSaturation(e, newSro) ;
            
            assert(newSro >= 0 and newSro <=1) ;
            assert(newSrw >= 0 and newSrw <=1) ;
            
            e->Store(key_prsCO2_, makeScalar( e->Status(key_prsCO2_), newSro ));
            e->Store(key_prsH2O_, makeScalar( e->Status(key_prsH2O_), newSrw ));

          }  
            
      } else if ( dsCO2_old < 0. && dsCO2_new > 0. ) { //turning from IMBIBITION to DRAINAGE
          flag_ = (abs(Sw-Sw_max)/Sw_max < tol_) and (Sw_max > Sw_min) ;
          if ( flag_ ) {
            //compute and store new pseudo residual saturations
            
            newSro = e->Read(key_prsCO2_) ;
            newSrw = e->Read(key_prsH2O_) ;
            
            WaterAndOilResidualSaturationImbibitionToDrainage(e, newSrw, newSro);
            
            assert(newSro >= 0 and newSro <=1) ;
            assert(newSrw >= 0 and newSrw <=1) ;
            
            e->Store( key_prsCO2_, makeScalar( e->Status( key_prsCO2_), newSro ));
            e->Store( key_prsH2O_, makeScalar( e->Status( key_prsH2O_), newSrw ));
            
          }  
      }       
}  

  
  

  
  
  
  
  
  
  
  
  
  


/**
   
  Brooks & Corey capillary pressure formula, here on spot the capilary pressure is bound to main drainage and imbibition curves.
  If it is out of bound will be replaced with the Drainage (Or Imbibition) capillary pressure.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::pc( Element<dim>* const e ) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    
    // update the pesudo residuals if it is necessory
    this->UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(key_prsCO2_)  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_ = e->Read(key_prsH2O_)  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    
    double64 pc(0);
    
    if ((sH2O > PseudoSwr_)&&(sCO2 > PseudoSor_))
      {
      pc = c[H2O]*pow((1. - PseudoSwr_)/(sH2O - PseudoSwr_),a[H2O]) + c[CO2]*pow((1.0 - PseudoSor_)/(sCO2 - PseudoSor_),a[CO2] );
      }
    else if (sH2O <= PseudoSwr_) pc =  MaxCapillaryPressure ;
    else if (sCO2 <= PseudoSor_) pc = -MaxCapillaryPressure ;
    
    std::pair<double64,double64> Pc_limits = this->CapillaryPressureLimits(e); // Check if it is in the range of Main Drainage and Imbibition limits
    
    if (pc <= Pc_limits.second) return Pc_limits.second;
    if (pc >= Pc_limits.first)  return Pc_limits.first;
    
    if (pc >  MaxCapillaryPressure) return  MaxCapillaryPressure ; // Check if it is not larger than the Maximum Capillary pressure
    if (pc < -MaxCapillaryPressure) return -MaxCapillaryPressure ; // Check if it is not larger than the Minimum Capillary pressure
    
    return pc;
}
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // Capilary pressure at the specific saturations
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::pc_at( Element<dim>* const e, double64 s ) const
  {
    
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    
    // update the pesudo residuals if it is necessory
    this->UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = s;
    assert(sH2O >= 0 and sH2O <=1) ;
    
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(key_prsCO2_)  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;
    
    double64 PseudoSwr_ = e->Read(key_prsH2O_)  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;
    
    
    double64 pc(0);
    
    if ((sH2O > PseudoSwr_)&&(sCO2 > PseudoSor_))
    {
      pc = c[H2O]*pow((1. - PseudoSwr_)/(sH2O - PseudoSwr_),a[H2O]) + c[CO2]*pow((1.0 - PseudoSor_)/(sCO2 - PseudoSor_),a[CO2] );
    }
    else if (sH2O <= PseudoSwr_) pc = MaxCapillaryPressure ;
    else if (sCO2 <= PseudoSor_) pc = -MaxCapillaryPressure ;
    
    std::pair<double64,double64> Pc_limits = this->CapillaryPressureLimits(e); // Check if it is in the range of Main Drainage and Imbibition limits
    
    if (pc <= Pc_limits.second) return Pc_limits.second;
    if (pc >= Pc_limits.first)  return Pc_limits.first;
    
    if (pc > MaxCapillaryPressure)  return MaxCapillaryPressure ; // Check if it is not larger than the Maximum Capillary pressure
    if (pc < -MaxCapillaryPressure) return -MaxCapillaryPressure ;// Check if it is not larger than the Minimum Capillary pressure
    
    return pc;
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
/**
   
  The first derivative of Brooks & Corey capillary pressure formula by assuming the variation of CO2 saturation is opposite of H2O saturations
  Here we also check the capillary function bound, and then if it ois out of bound the derivative.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dpcds( Element<dim>* const e ) const
  {
  
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    const double64 sCO2 = 1. - sH2O;
    

    double64 PseudoSor_ = e->Read(key_prsCO2_)  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_ = e->Read(key_prsH2O_)  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;
  
    double64 Dpc(0);
    
    // first calculate the upper limits and lower limits of the capilary curve.

    if (sH2O < PseudoSwr_ )     Dpc = -MaxCapillaryPressureDerivative;
    else if (sCO2 < PseudoSor_) Dpc = -MaxCapillaryPressureDerivative ;
    else
      Dpc = -a[H2O]*c[H2O]*pow((1.0 - PseudoSwr_) / (sH2O - PseudoSwr_), a[H2O]) / (sH2O - PseudoSwr_) + a[CO2]*c[CO2]*pow((1.0 - PseudoSor_) / (sCO2 - PseudoSor_), a[CO2]) / (sCO2 - PseudoSor_);
    
    return Dpc ;
}
  
 

  
  
  
  
  
  
  
  
// The first derivative of Capilary pressure at the specific saturations
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dpcds_at( Element<dim>* const e, double64 sH2O ) const
  {
    assert( sH2O >= 0. );
    assert( sH2O <= 1. );
  
    const double64 sCO2 = 1. - sH2O;
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    double64 PseudoSor_ = e->Read(key_prsCO2_)  ;
    double64 PseudoSwr_ = e->Read(key_prsH2O_)  ;  // To have primary Drianage and Imibition parameters
    
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;
    
    double64 Dpc(0);
    
    // first calculate the upper limits and lower limits of the capilary curve.
    if (sH2O < PseudoSwr_ )     Dpc = -MaxCapillaryPressureDerivative;
    else if (sCO2 < PseudoSor_) Dpc = -MaxCapillaryPressureDerivative ;
    else
      Dpc = -a[H2O]*c[H2O]*pow((1.0 - PseudoSwr_) / (sH2O - PseudoSwr_), a[H2O]) / (sH2O - PseudoSwr_) + a[CO2]*c[CO2]*pow((1.0 - PseudoSor_) / (sCO2 - PseudoSor_), a[CO2]) / (sCO2 - PseudoSor_);
    
    return Dpc ;
}    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
   
   This is residual water saturation for the transit process from drainage to imbibition Water residual saturation has been estimated
   from the intersection of imbibition curve with drainage curve We assume the rock properties are constant, i.e. aw, ao, co, and cw
   are constant. We need to estimate the residual water saturation.
  
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterResidualSaturation( Element<dim>* const e, double64 Sro) const
  {
    
    const double64 awd_ =  ac_params_[AWD];
    const double64 aod_ =  ac_params_[AOD];
    const double64 cwd_ =  ac_params_[CWD];
    const double64 cod_ =  ac_params_[COD];
    
    const double64 awi_ =  ac_params_[AWI];
    const double64 aoi_ =  ac_params_[AOI];
    const double64 cwi_ =  ac_params_[CWI];
    const double64 coi_ =  ac_params_[COI];
   
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_ = e->Read(key_prsCO2_)  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_ = e->Read(key_prsH2O_)  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

  
    double64 srH2O = PseudoSwr_ ;
    
    
    double64 RHS = cwd_*pow((1.0 - PseudoSwr_)/(sH2O - PseudoSwr_),awd_) + cod_*pow((1.0 - PseudoSor_)/(sCO2 - PseudoSor_),aod_)- coi_*pow((1.0 - Sro)/(sCO2 - Sro),aoi_) ;
    
    RHS = cwi_/RHS ;
    RHS = pow(RHS,1./awi_) ;
    
    if ( (sH2O > PseudoSwr_) && (sH2O < 1.-PseudoSor_)) srH2O = 1.-(1.-sH2O)/(1.-RHS) ;
    
    return  srH2O ;
  }
  
 
  
  

  
  
  
  
  
  
  
  
  
  
/**
  Here , I calculate the residual water saturation and residual oil saturation for the transit process from  imbibition to drainage.
 
  These  residual saturation has been estimated from the intersection of imbibition curve with drianage curve.
  We assume the rock properties are constant, i.e. aw, ao, co, and cw are constant. Just we need to estimate
  the water and oil residual saturations.
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::WaterAndOilResidualSaturationImbibitionToDrainage( Element<dim>* const e, double64& Swr_, double64& Sor_) const
  {
    const size_t    Nr(3);   // number of iteration of solving Non-linear system of equations to get the psedo-resduals.
    const double64 awd_ =  ac_params_[AWD];
    const double64 aod_ =  ac_params_[AOD];
    const double64 cwd_ =  ac_params_[CWD];
    const double64 cod_ =  ac_params_[COD];
    
    const double64 awi_ =  ac_params_[AWI];
    const double64 aoi_ =  ac_params_[AOI];
    const double64 cwi_ =  ac_params_[CWI];
    const double64 coi_ =  ac_params_[COI];
    
    const double64 iaod_(1.0/aod_) ;
    const double64 iawd_(1.0/awd_) ;
    
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    double64 Sw1(e->Read(key_SwDrToImb_)) ;
    assert(Sw1 >= 0 and Sw1 <=1) ;

    double64 Sw2(e->Read(key_SwImbToDr_)) ;
    assert(Sw2 >= 0 and Sw2 <=1) ;

    double64 So2(1.-Sw2);
    double64 So1(1.-Sw1);
    
    double64 psrH2O(e->Read(key_prsH2O_)) ;    // This means the algorithm use the current residual water saturation for the next steps
    assert(psrH2O >= 0 and psrH2O <=1) ;

    double64 psrCO2(e->Read(key_prsCO2_)) ;    // The same for residual oil saturation for the next steps
    assert(psrCO2 >= 0 and psrCO2 <=1) ;
    
    double64 PcS2 = coi_*pow((1.0 - psrCO2)/(So2 - psrCO2),aoi_) + cwi_*pow((1.0 - psrH2O)/(Sw2 - psrH2O),awi_) ;
    double64 PcS1 = coi_*pow((1.0 - psrCO2)/(So1 - psrCO2),aoi_) + cwi_*pow((1.0 - psrH2O)/(Sw1 - psrH2O),awi_) ;
    
    Swr_= psrH2O ;
    
    for (size_t iter=0; iter<Nr ; iter++) { // do estimate for pesdu residual saturation
      
      // step one estimates the oil residual saturation from upper turning point
      Sor_ = 1.0 - (1.0-So2)/(1.0 - pow(cod_/(PcS2 - cwd_*pow((1.0 - Swr_)/(Sw2 - Swr_),awd_)) ,iaod_)) ;
      
      // step two estimates the water residual saturation from lowwer turning point
      Swr_ = 1.0 - (1.0-Sw1)/(1.0 - pow(cwd_/(PcS1 - cod_*pow((1.0 - Sor_)/(So1 - Sor_),aod_)) ,iawd_)) ;
      
    } ;
    
    if (isnan(Swr_)) Swr_ = psrH2O ;  // in case the iterative solver give us nan, we use the default value which is previous pesdu saturations
    if (isnan(Sor_)) Sor_ = psrCO2 ;  // in case the iterative solver give us nan, we use the default value which is previous pesdu saturations
    
    //if (Swr_ >= 0 and Swr_ <= 1) Swr_ = psrH2O ;
    //if (Sor_ >= 0 and Sor_ <= 1) Sor_ = psrCO2 ;

}

  
  
  
  
  
  
  
  
  
  
  
  
/**
  
  For calculating the water relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
  The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
  the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
 
*/
  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw( Element<dim>* const e ) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    const double64 sCO2 = 1. - sH2O;

    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    
    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 krwww = pow(Snw, 3.+2.*a[H2O]) ; // refer to Eq. 10 from Skjaeveland et al. 2000
    double64 krwow = (1. - pow(Sno, 2.* a[CO2]+1.)) * square(1. - Sno);
    
    double64 krw = (c[H2O]*krwww - c[CO2]*krwow) / (c[H2O] - c[CO2]) ;
    
    std::pair<double64, double64> limits  = this->krwLimits_at(e, sH2O) ;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    if ((PseudoSor_ != Sor_) || (PseudoSwr_ != Swr_))
    {
      krw = krw *(limits.second- limits.first) + limits.first ;
    }
    
    return krw ;  // Mahyar: If sH2O< srH2O or sH2O > 1-srCO2 .. The krw will be between 0 or 1
}
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
   
   Second form of the water relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krw_at( Element<dim>* const e, double64 S) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    double64 sH2O(S) ;
    double64 sCO2(1.-S);
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 krwww = pow(Snw, 3.+2.*a[H2O]) ; // refer to Eq. 10 from Skjaeveland et al. 2000
    double64 krwow = (1. - pow(Sno, 2.*a[CO2]+1.)) * square(1. - Sno);
    
    double64 krw = (c[H2O]*krwww - c[CO2]*krwow) / (c[H2O] - c[CO2]) ;
    
    std::pair<double64, double64> limits  = this->krwLimits_at(e, sH2O) ;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    if ((PseudoSor_ != Sor_) || (PseudoSwr_ != Swr_))
    {
      krw = krw *(limits.second- limits.first) + limits.first ;
    }
    
    return krw ;
    
  }



  
  
  
  

  
  
  
  
  
/**
   
   For calculating the oil relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
   The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
   the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn( Element<dim>* const e ) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;

    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;
    
    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 kroww = (1. - pow(Snw,2.*a[H2O]+1.)) * square(1. - Snw);
    double64 kroow = pow(Sno, 3.+2.*a[CO2]) ; // refer to Eq. 12 from Skjaeveland et al. 2000
    
    double64 krn = (c[H2O]*kroww - c[CO2]*kroow) / (c[H2O] - c[CO2]) ;
    
    std::pair<double64, double64> limits  = this->krnLimits_at(e, sH2O) ;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    if ((PseudoSor_ != Sor_) || (PseudoSwr_ != Swr_))
    {
      krn = krn *(limits.second- limits.first) + limits.first ;
    }
    
    return krn ;
    
    }




  
  
  
  
  

  
/**
   
   Second form of the oil relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krn_at( Element<dim>* const e, double64 S) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    double64 sH2O = S ;
    double64 sCO2 = 1.-S ;
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    CheckPcLimitsAndResetResiduals(e, PseudoSwr_, PseudoSor_, a, c);

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 kroww = (1. - pow(Snw,2.*a[H2O]+1.)) * square(1. - Snw);
    double64 kroow = pow(Sno, 3.+2.*a[CO2]) ; // refer to Eq. 12 from Skjaeveland et al. 2000
    
    double64 krn = (c[H2O]*kroww - c[CO2]*kroow) / (c[H2O] - c[CO2]) ;
    
    std::pair<double64, double64> limits  = this->krnLimits_at(e, sH2O) ;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    if ((PseudoSor_ != Sor_) || (PseudoSwr_ != Swr_))
    {
      krn = krn *(limits.second- limits.first) + limits.first ;
    }
    
    return krn ;
  }




  
  
  
  
  
  
  
  
  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds( Element<dim>* const e ) const
  {
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000

    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw;
    
    double64 dkrwww = (3.+2.*a[H2O])*pow(Snw,(2.+2.*a[H2O]))*dSnw ;
    double64 dkrwow = (-(2*a[CO2]+1)*pow(Sno,(2.*a[CO2]))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*a[CO2]+1.)))*(1.-Sno))*dSno ;
    
    return (c[H2O]*dkrwww-c[CO2]*dkrwow)/(c[H2O]-c[CO2]) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

  
  
  
  
  
  
  
  
  
  
  
  
/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_at( Element<dim>* const e, double64 S) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O (S);
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(key_prsCO2_));
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw;
    
    double64 dkrwww = (3.+2.*a[H2O])*pow(Snw,(2.+2.*a[H2O]))*dSnw ;
    double64 dkrwow = (-(2*a[CO2]+1)*pow(Sno,(2.*a[CO2]))*square(1.-Sno) - 2.*(1.-pow(Sno,(2.*a[CO2]+1.)))*(1.-Sno))*dSno ;
    
    return (c[H2O]*dkrwww-c[CO2]*dkrwow)/(c[H2O]-c[CO2]) ; // 1st derivative of  Eq. 14a from Skjaeveland et al. 2000
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds( Element<dim>* const e ) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert(sH2O >= 0 and sH2O <=1) ;

    
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw;
    
    double64 dkroww = (-(2*a[H2O]+1)*pow(Snw,(2.*a[H2O]))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*a[H2O]+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*a[CO2])*pow(Sno,(2.+2.*a[CO2]))*dSno  ;
    
    return (c[H2O]*dkroww-c[CO2]*dkroow)/(c[H2O]-c[CO2]) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000
    
}
  
  
  

  
  
  
  
  
  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_at( Element<dim>* const e, double64 S ) const
  {
    
    std::array<double64, 2> a ;
    std::array<double64, 2> c ;
    
    SetBrooksCoreyCurvesParameters(e , a, c) ;
    
    // update the pesudo residuals if it is necessory
    UpdatePseudoResidualAndEndpointSaturations(e) ;
    
    const double64 sH2O = S;
    const double64 sCO2 = 1. - sH2O;
    
    double64 PseudoSor_(e->Read(key_prsCO2_))  ;
    assert(PseudoSor_ >= 0 and PseudoSor_ <=1) ;

    double64 PseudoSwr_(e->Read(key_prsH2O_))  ;  // To have primary Drianage and Imibition parameters
    assert(PseudoSwr_ >= 0 and PseudoSwr_ <=1) ;

    double64 Snw = (sH2O - PseudoSwr_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - PseudoSor_) / (1. - PseudoSwr_ - PseudoSor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    
    double64 dSo = -1.; // if change of water saturation will change co2 saturations.
    
    double64 dSnw = 1./(1.-PseudoSwr_-PseudoSor_) ;
    
    if (Snw*(1.-Snw)*Sno*(1.-Sno)==0.) dSnw = 0.0 ;
    
    double64 dSno = dSo * dSnw;
    
    double64 dkroww = (-(2*a[H2O]+1)*pow(Snw,(2.*a[H2O]))*square(1.-Snw) - 2.*(1.-pow(Snw,(2.*a[H2O]+1.)))*(1.-Snw))*dSnw ;
    double64 dkroow = (3.+2.*a[CO2])*pow(Sno,(2.+2.*a[CO2]))*dSno  ;
    
    return (c[H2O]*dkroww-c[CO2]*dkroow)/(c[H2O]-c[CO2]) ; // 1st derivative of  Eq. 14b from Skjaeveland et al. 2000
    
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
 
 Check The Limits of Capillary presure and reset the residuals
 
*/
template<size_t dim, template<size_t> class USER>
void BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::CheckPcLimitsAndResetResiduals( Element<dim>* const e , double64& newSrH2O, double64& newSrCO2, std::array<double64, 2>& a,std::array<double64, 2>& c) const {
  
    const double64 srH2O(e->Read(User()->key_srH2O));  // To have primary Drianage and Imibition parameters
    assert(srH2O >= 0 and srH2O <=1) ;

    const double64 srCO2(e->Read(User()->key_srCO2));
    assert(srCO2 >= 0 and srCO2 <=1) ;

    std::pair<double64,double64> Pc_limits = CapillaryPressureLimits(e);
    
    double64 local_pc = pc(e);
    
    //check if hit primary drainage
    if (local_pc >= Pc_limits.second ) {
      
      newSrCO2 = srCO2 ;
      newSrH2O = srH2O ;
      
      a[H2O] =  ac_params_[AWD];
      a[CO2] =  ac_params_[AOD];
      c[H2O] =  ac_params_[CWD];
      c[CO2] =  ac_params_[COD];
    }
    
    //check if hit primary Imbibition
    if (local_pc <= Pc_limits.first ) {
      
      newSrCO2 = srCO2 ;
      newSrH2O = srH2O ;
      
      a[H2O] =  ac_params_[AWI];
      a[CO2] =  ac_params_[AOI];
      c[H2O] =  ac_params_[CWI];
      c[CO2] =  ac_params_[COI];
    }
}
  

  
  
  
  
  
  
  
/**
   
   calculate the upper and lowwer limits of Capillary presure
   
*/
template<size_t dim, template<size_t> class USER>
  std::pair<double64,double64> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::CapillaryPressureLimits( Element<dim>* const e ) const
  {
    std::pair<double64,double64> pc_PrimaryImbibitionDrainage;
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double64 sCO2 = 1. - sH2O;
    
    double64 srCO2(e->Read(User()->key_srCO2))  ;
    double64 srH2O(e->Read(User()->key_srH2O))  ;  // To have primary Drianage and Imibition parameters
    
    const double64 awd_ =  ac_params_[AWD];
    const double64 aod_ =  ac_params_[AOD];
    const double64 cwd_ =  ac_params_[CWD];
    const double64 cod_ =  ac_params_[COD];
    
    const double64 awi_ =  ac_params_[AWI];
    const double64 aoi_ =  ac_params_[AOI];
    const double64 cwi_ =  ac_params_[CWI];
    const double64 coi_ =  ac_params_[COI];
    
    
    // first calculate the upper limits and lower limits of the capilary curve.
    
    pc_PrimaryImbibitionDrainage.first  = cwd_*pow((1.0 - srH2O)/(sH2O - srH2O),awd_) + cod_*pow((1.0 - srCO2)/(sCO2 - srCO2),aod_);
    pc_PrimaryImbibitionDrainage.second = cwi_*pow((1.0 - srH2O)/(sH2O - srH2O),awi_) + coi_*pow((1.0 - srCO2)/(sCO2 - srCO2),aoi_);
    
    if (pc_PrimaryImbibitionDrainage.first  > MaxCapillaryPressure)   pc_PrimaryImbibitionDrainage.first  = MaxCapillaryPressure ;
    if (pc_PrimaryImbibitionDrainage.second > MaxCapillaryPressure)   pc_PrimaryImbibitionDrainage.second = MaxCapillaryPressure ;
    
    if (pc_PrimaryImbibitionDrainage.first  < -MaxCapillaryPressure)  pc_PrimaryImbibitionDrainage.first  = -MaxCapillaryPressure ;
    if (pc_PrimaryImbibitionDrainage.second < -MaxCapillaryPressure)  pc_PrimaryImbibitionDrainage.second = -MaxCapillaryPressure ;
    
    if (isinf(pc_PrimaryImbibitionDrainage.first))  pc_PrimaryImbibitionDrainage.first  = MaxCapillaryPressure ;
    if (isinf(pc_PrimaryImbibitionDrainage.second)) pc_PrimaryImbibitionDrainage.second = MaxCapillaryPressure ;
    
    return pc_PrimaryImbibitionDrainage;
  }

  
  
  
  
  
  
  
  /**
   
   calculate the upper and lowwer limits of relative permeability of water
   
   */
  template<size_t dim, template<size_t> class USER>
  std::pair<double64,double64> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krwLimits_at( Element<dim>* const e, double64 const S ) const
  {
    std::pair<double64,double64> krw_PrimaryImbibitionDrainage;
    const double64 sH2O = S ;
    const double64 sCO2 = 1. - sH2O;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    const double64 awd_ =  ac_params_[AWD];
    const double64 aod_ =  ac_params_[AOD];
    const double64 cwd_ =  ac_params_[CWD];
    const double64 cod_ =  ac_params_[COD];
    
    const double64 awi_ =  ac_params_[AWI];
    const double64 aoi_ =  ac_params_[AOI];
    const double64 cwi_ =  ac_params_[CWI];
    const double64 coi_ =  ac_params_[COI];
    
    
    // first calculate the upper limits and lower limits of the capilary curve.
    
  
    double64 Snw = (sH2O - Swr_) / (1. - Swr_ - Sor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - Sor_) / (1. - Swr_ - Sor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    krw_PrimaryImbibitionDrainage.first  = (cwd_*pow(Snw, 3.+2.*awd_)  - cod_*(1. - pow(Sno, 2.* aod_+1.)) * square(1. - Sno)) / (cwd_ - cod_) ;
    krw_PrimaryImbibitionDrainage.second = (cwi_*pow(Snw, 3.+2.*awi_)  - coi_*(1. - pow(Sno, 2.* aoi_+1.)) * square(1. - Sno)) / (cwi_ - coi_) ;
    
    return krw_PrimaryImbibitionDrainage;
  }

  
  
  
  
  
  
  
  
  /**
   
   calculate the upper and lowwer limits of relative permeability of non-wet phase
   
   */
  template<size_t dim, template<size_t> class USER>
  std::pair<double64,double64> BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::krnLimits_at( Element<dim>* const e, double64 const S ) const
  {
    std::pair<double64,double64> krn_PrimaryImbibitionDrainage;
    const double64 sH2O = S ;
    const double64 sCO2 = 1. - sH2O;
    
    double64 Sor_(e->Read(User()->key_srCO2))  ;
    double64 Swr_(e->Read(User()->key_srH2O))  ;
    
    const double64 awd_ =  ac_params_[AWD];
    const double64 aod_ =  ac_params_[AOD];
    const double64 cwd_ =  ac_params_[CWD];
    const double64 cod_ =  ac_params_[COD];
    
    const double64 awi_ =  ac_params_[AWI];
    const double64 aoi_ =  ac_params_[AOI];
    const double64 cwi_ =  ac_params_[CWI];
    const double64 coi_ =  ac_params_[COI];
    
    
    // first calculate the upper limits and lower limits of the capilary curve.
    
    
    double64 Snw = (sH2O - Swr_) / (1. - Swr_ - Sor_);   // refer to Eq. 11 from Skjaeveland et al. 2000
    double64 Sno = (sCO2 - Sor_) / (1. - Swr_ - Sor_);   // refer to Eq. 13 from Skjaeveland et al. 2000
    
    Snw = std::min( std::max( Snw, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    Sno = std::min( std::max( Sno, 0. ), 1. ); // be sure the normalize saturation varies between 0 and 1
    
    krn_PrimaryImbibitionDrainage.first   = (cwi_*(1. - pow(Snw,2.*awi_+1.)) * square(1. - Snw)  - coi_*pow(Sno, 3.+2.*aoi_)) / (cwi_ - coi_) ;
    krn_PrimaryImbibitionDrainage.second  = (cwd_*(1. - pow(Snw,2.*awd_+1.)) * square(1. - Snw)  - cod_*pow(Sno, 3.+2.*aod_)) / (cwd_ - cod_) ;

    return krn_PrimaryImbibitionDrainage;
  }

  
  
  
  
  
  
  
  
  /**
   
  Numerical derivative of relative water permeability
  
  */
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_Numerical( Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  
  
  
  
  
  /**
   
  Numerical derivative of relative water permeability at specific saturation
   
  */
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrwds_at_Numerical( Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }



 
  
  
  
  
  
  
  
  
  /**
   
   Numerical derivative of relative co2 permeability
   
   */
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_Numerical( Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  
  
  
  
  
  /**
   
   Numerical derivative of relative co2 permeability at specific saturation
   
   */
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctionsWithHysteresis<dim,USER>::dkrnds_at_Numerical( Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  

template class BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsModule2>;
template class BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsModule2>;
template class BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsModule2>;

template class BrooksCoreySaturationFunctionsWithHysteresis<1U,FlowFunctionsModule5>;
template class BrooksCoreySaturationFunctionsWithHysteresis<2U,FlowFunctionsModule5>;
template class BrooksCoreySaturationFunctionsWithHysteresis<3U,FlowFunctionsModule5>;

} // csmp
