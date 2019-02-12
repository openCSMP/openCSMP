#include "FlowFunctionsModule.h"
#include "H2O_CO2_NaCl_FlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "Element.h"

/*
    SKM testing log - full three-phase case:  brine - CO2 - halite
*/

using namespace std;

namespace csmp {

// helper functions

/**
    Reports the existence of aqueous, carbonic and halite phases on all nodes of the element.
    Returns true if both fluid phases co-exist and halite is absent.
    Else, returns false.
*/
template<size_t dim, template<size_t> class USER>
bool H2O_CO2_NaCl_FlowFunctions<dim,USER>::ContinuousPhases( Element<dim>* const e, bool& aqueous, bool& carbonic, bool& halite ) const
 {
    assert( e != nullptr );
    const size_t nodes(e->Nodes());
    aqueous = carbonic = halite = true;

    for ( size_t i=0; i<nodes; ++i ) {
         double64 s = e->N(i)->Read( User()->key_sH2O );
         if ( isnan(s) && s <= numeric_limits<double64>::epsilon()*2 ) aqueous = false;
         s = e->N(i)->Read( User()->key_sCO2 );
         if ( isnan(s) && s <= numeric_limits<double64>::epsilon()*2 ) carbonic = false;
         s = e->N(i)->Read( User()->key_NaCl );
         if ( isnan(s) && s <= numeric_limits<double64>::epsilon()*2 ) halite = false;
      }
   
    if ( aqueous == true && carbonic == true && halite == false ) return true;
   
    return false;
 }
 
 
 
 

template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::InterpolateAqueousPhaseSaturation( const Element<dim>* const e ) const
 {
    assert( e != nullptr );
    e->N_AtBaryCenter( e->FE()->NRST );
    const size_t nodes(e->Nodes());
    double sw(0.);

    for ( size_t i=0; i<nodes; ++i ) {
         double64 s = e->N(i)->Read( User()->key_sH2O );
         sw += e->FE()->NRST[i] * s;
      }
   
    return sw;
 }

 
 
 
 
/**
    Interpolates saturations of all mobile phases to barycenter.
 
    @return returns sum of sw + snw, which should be 1 - halite_saturation
*/
template<size_t dim, template<size_t> class USER>
double64  H2O_CO2_NaCl_FlowFunctions<dim,USER>::InterpolateSaturations( Element<dim>* const e, double64& sw, double64& sn ) const
 {
    assert( e != nullptr );
    e->N_AtBaryCenter( e->FE()->NRST );
    const size_t nodes(e->Nodes());
   
    sw = sn = 0.;

    for ( size_t i=0; i<nodes; ++i ) {
         double64 s = e->N(i)->Read( User()->key_sH2O );
         sw += e->FE()->NRST[i] * s;
         s = e->N(i)->Read( User()->key_sCO2 );
         sn += e->FE()->NRST[i] * s;
      }

    return sw + sn;
 }



 

/**
    Interpolates the properties of the aqueous phase, if there is any, else returns false.

    @return true if method was able to interpolate a density and false if it was not.
 
    @todo TODO: check whether this could be improved by reading the element phase state.
 
    @note compared to NEITHER as a state where no phase is continous across the nodes,
    looking at the saturations offers more information ?

*/
template<size_t dim, template<size_t> class USER>
bool H2O_CO2_NaCl_FlowFunctions<dim,USER>::InterpolateAqueousPhase( Element<dim>* const e, double64& sw, double64& rhow, double64& muw ) const
 {
    assert( e != nullptr );
    e->N_AtBaryCenter( e->FE()->NRST );
    const size_t nodes(e->Nodes());
    double64     ipol_sum(0.);
   
    sw = rhow = muw = 0.;

    for ( size_t i=0; i<nodes; ++i ) {
         double64 s = e->N(i)->Read( User()->key_sH2O );
         if ( s > 0. ) {
             sw   += e->FE()->NRST[i] * s;
             rhow += e->FE()->NRST[i] * e->N(i)->Read( User()->key_rhoH2O );
             muw  += e->FE()->NRST[i] * e->N(i)->Read( User()->key_muH2O );
             ipol_sum += e->FE()->NRST[i];
          }
      }
    // correcting for the changed weighting factor when nodes have no density or viscosity values
    if ( ipol_sum <= numeric_limits<double64>::epsilon() ) return false;
    rhow *= 1. / ipol_sum;
    muw  *= 1. / ipol_sum;
   
    return true;
   
 } // end InterpolateAqueousPhase





/**
 
     @return true if method was able to interpolate a density and false if it was not.

*/
template<size_t dim, template<size_t> class USER>
bool H2O_CO2_NaCl_FlowFunctions<dim,USER>::InterpolateCarbonicPhase( Element<dim>* const e,
                                                                     double64& snw, double64& rhon, double64& mun ) const
 {
    assert( e != nullptr );
    e->N_AtBaryCenter( e->FE()->NRST );
    const size_t nodes(e->Nodes());
    double64     ipol_sum(0.);
   
    snw = rhon = mun = 0.;

    for ( size_t i=0; i<nodes; ++i ) {
        // carbonic phase
        double64 sn = e->N(i)->Read( User()->key_sCO2 );
        if ( sn > 0. ) {
             snw  += e->FE()->NRST[i] * sn;
             rhon += e->FE()->NRST[i] * e->N(i)->Read( User()->key_rhoCO2 );
             mun  += e->FE()->NRST[i] * e->N(i)->Read( User()->key_muCO2);
             ipol_sum += e->FE()->NRST[i];
          }
       }
    // correcting for the changed weighting factor when nodes have no density or viscosity values
    if ( ipol_sum <= numeric_limits<double64>::epsilon() ) return false;
    rhon *= 1. / ipol_sum;
    mun  *= 1. / ipol_sum;
   
    return true;

 } // end InterpolateCarbonicPhase






template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::InterpolateSystem( Element<dim>* const e,
                                                                  double64& sw, double64& rhow, double64& muw,
                                                                  double64& snw, double64& rhon, double64& mun ) const
 {
    assert( e != nullptr );
    e->N_AtBaryCenter( e->FE()->NRST );
    const size_t nodes(e->Nodes());
    double64 halite_saturation(0.);
    double64 ipol_sum1(0.), ipol_sum2(0.);

    sw = rhow = muw = snw = rhon = mun = 0.;

    for ( size_t i=0; i<nodes; ++i ) {
         double64 s = e->N(i)->Read( User()->key_sH2O );
         if ( s > 0. ) {
              sw   += e->FE()->NRST[i] * s;
              rhow += e->FE()->NRST[i] * e->N(i)->Read( User()->key_rhoH2O );
              muw  += e->FE()->NRST[i] * e->N(i)->Read( User()->key_muH2O );
              ipol_sum1 += e->FE()->NRST[i];
           }
         s = e->N(i)->Read( User()->key_sCO2 );
         if ( s > 0. ) {
              snw  += e->FE()->NRST[i] * s;
              rhon += e->FE()->NRST[i] * e->N(i)->Read( User()->key_rhoCO2 );
              mun  += e->FE()->NRST[i] * e->N(i)->Read( User()->key_muCO2);
              ipol_sum2 += e->FE()->NRST[i];
           }
         halite_saturation += e->FE()->NRST[i] * e->N(i)->Read( User()->key_NaCl );
      }
    // correcting for the changed weighting factor when nodes have no density or viscosity values
    if ( ipol_sum1 > 0. ) {
         rhow *= 1. / ipol_sum1;
         muw  *= 1. / ipol_sum1;
      }
    if ( ipol_sum2 > 0. ) {
         rhon *= 1. / ipol_sum2;
         mun  *= 1. / ipol_sum2;
      }

    assert( fabs(sw + snw + halite_saturation) <= 1 + numeric_limits<double64>::epsilon()*2. );
    assert( fabs(sw + snw + halite_saturation) >= 1 - numeric_limits<double64>::epsilon()*2. );
   
    return halite_saturation;
    
} // InterpolateSystem
  









// FLOW FUNCTIONS
  
/**
    Saturation is always expected to have a value between 0..1.
*/
// SKM OK
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Sw( Element<dim>* const e ) const
  {
     assert( e != nullptr );
     return e->PropertyValueAtBaryCenter( User()->key_sH2O );
  }
  
  


// SKM OK
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Snw( Element<dim>* const e ) const
  {
     assert( e != nullptr );
     return e->PropertyValueAtBaryCenter( User()->key_sCO2 );
  }

  
  
  
/**
      Mobility of phase i, lambda_i = rho_i * kri(sw) / mu_i.
 
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Mobility( Element<dim>* const e, size_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    // aqueous phase
    const bool with_aqueous_phase_properties = InterpolateAqueousPhase( e, sw_, rhow_, muw_ );

    if ( phase == 0U ) {
         if ( !with_aqueous_phase_properties || sw_ <= e->Read( User()->key_srH2O ) ) return 0.;
         assert( User()->krw_at(e,sw_) >= 0. );
         assert( User()->krw_at(e,sw_) <= 1. );
         return rhow_ * (User()->krw_at(e,sw_) / muw_);
      }
   
    // carbonic phase
    if ( sw_ >= 1. - numeric_limits<double64>::epsilon() ) return 0.;
    const bool with_carbonic_phase_properties = InterpolateCarbonicPhase( e, sn_, rhon_, mun_ );
    if ( !with_carbonic_phase_properties || sn_ <= e->Read( User()->key_srCO2 ) ) return 0.;
    assert( User()->krn_at(e,sw_) >= 0. );
    assert( User()->krn_at(e,sw_) <= 1. );
    return rhon_ * (User()->krn_at(e,sw_) / mun_);
 }


  
   
/**
    Mobility of phase i, lambda_i = rho_i * kri(sw) / mu_i.
 
    Using prescribed sw value, instead of value intepolated to element barycentre.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Mobility_at( Element<dim>* const e, size_t phase, double64 sw ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    const double64 srw  = e->Read( User()->key_srH2O );
    const double64 srn  = e->Read( User()->key_srCO2 );

    if ( phase == 0U ) {
         if ( sw <= srw  || !InterpolateAqueousPhase( e, sw_, rhow_, muw_ ) ) return 0.;
         // note: the input sw must be used here
         assert( User()->krw_at(e,sw) >= 0. );
         assert( User()->krw_at(e,sw) <= 1. );
         return rhow_ * (User()->krw_at(e,sw) / muw_);
      }
   
    // carbonic phase
    if ( sw >= 1. - srn || !InterpolateCarbonicPhase( e, sn_, rhon_, mun_ ) ) return 0.;
    assert( User()->krn_at(e,sw) >= 0. );
    assert( User()->krn_at(e,sw) <= 1. );
    return rhon_ * (User()->krn_at(e,sw) / mun_);
 }
  



/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityDerivative( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U ) {
          if ( !InterpolateAqueousPhase( e, sw_, rhow_, muw_ ) || sw_ <= e->Read( User()->key_srH2O ) ) return 0.;
          return rhow_ * (User()->dkrwds_at(e,sw_) / muw_);
       }
    // carbonic phase
    if ( sw_ >= 1. - e->Read( User()->key_srCO2 ) || !InterpolateCarbonicPhase( e, sn_, rhon_, mun_ ) ) return 0.;
    assert( User()->krn_at(e,sw_) >= 0. );
    assert( User()->krn_at(e,sw_) <= 1. );

    return rhon_ * (User()->dkrnds_at(e,sw_) / mun_);
 }
  



/**
    Mobility saturation derivative for aqueous (i=0) or carbonic (i=1) phase.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityDerivative_at( Element<dim>* const e, size_t phase, double64 sw ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U ) {
           if ( sw <= e->Read( User()->key_srH2O ) || !InterpolateAqueousPhase( e, sw_, rhow_, muw_ ) ) return 0.;
          return rhow_ * (User()->dkrwds_at(e,sw) / muw_);
       }
    // carbonic phase
    if ( sw >= 1. - e->Read( User()->key_srCO2 ) || !InterpolateCarbonicPhase( e, sn_, rhon_, mun_ ) ) return 0.;
    assert( User()->krn_at(e,sw) >= 0. );
    assert( User()->krn_at(e,sw) <= 1. );

    return rhon_ * (User()->dkrnds_at(e,sw) / mun_);
  }
  

  
  
  
/**
    Sum of mass mobilities (not multiplied with permeability).
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TotalMobility(Element<dim>* const e ) const
  {
    assert( e != nullptr );
    
    const bool with_aq_phase_props = InterpolateAqueousPhase( e, sw_, rhow_, muw_ );
    assert( sw_ >= 0. && sw_ <= 1. );
   
    // if there is a mobile aqueous phase and the carbonic phase is immobile
    assert( !isnan(User()->krw_at(e,sw_)) );
    double64 mob_t = (with_aq_phase_props==true) ? rhow_ * User()->krw_at(e,sw_) / muw_ : 0.;
    const double64 srn = e->Read( User()->key_srCO2 );
    if ( sw_ >= (1. - srn) ) return mob_t;

    // if there is a mobile carbonic phase as well
    if ( !InterpolateCarbonicPhase( e, sn_, rhon_, mun_ ) ) return mob_t;
    assert( !isnan(User()->krn_at(e,sw_)) );
    if ( sn_ > srn ) mob_t += rhon_ * User()->krn_at(e,sw_) / mun_;
    
    return mob_t;
 }




/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TotalMobility_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    assert( sw >= 0. && sw <= 1. );

    // if there is a mobile aqueous phase and the carbonic phase is immobile
    assert( !isnan(User()->krw_at(e,sw)) );
    double64 mob_t = (InterpolateAqueousPhase( e, sw_, rhow_, muw_ )) ? rhow_ * User()->krw_at(e,sw) / muw_ : 0.;
    const double64 srn = e->Read( User()->key_srCO2 );
    if ( sw >= (1. - srn) ) return mob_t;

    // if there is a mobile carbonic phase as well
    InterpolateCarbonicPhase( e, sn_, rhon_, mun_ );
    assert( !isnan(User()->krn_at(e,sw)) );
    if ( sn_ > srn ) mob_t += rhon_ * User()->krn_at(e,sw) / mun_;
     
    return mob_t;
 }


  
  



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProduct( Element<dim>* const e ) const
 {
    assert( e != nullptr );
    const double64 lw_rhow = Mobility(e,0U);
    const double64 ln_rhon = Mobility(e,1U);
    return  (lw_rhow * ln_rhon) / (lw_rhow + ln_rhon);
 }




/**
   Mobility product evaluated at the users supplied water saturation.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProduct_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    const double64 lw_rhow = Mobility_at(e,0U,sw);
    const double64 ln_rhon = Mobility_at(e,1U,sw);
    return  (lw_rhow * ln_rhon) / (lw_rhow + ln_rhon);
 }

  
  
  


/** 
    Saturation derivative of mobility product.
 
    rhow*rhon*((diff(lw(sw), sw))*ln(sw)^2*rhon*sw+lw(sw)^2*rhow)/((rhow*lw(sw)+rhon*ln(sw))^2*sw);
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative( Element<dim>* const e, bool evaluate_numerically ) const
 {
    assert( e != nullptr );
    // product is zero at endmember saturations
    if ( User()->EffectiveSaturation(e) <= 0. || User()->EffectiveSaturation(e) >= 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(e);
    // mobilities already contain densities
    const double64 lw  = Mobility( e, 0U );
    const double64 ln  = Mobility( e, 1U );
    // NB: if the mobilities are greated than 0, the phase properties will be available
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = (lw > 0.) ? rhow_ * User()->dkrwds(e) / muw_ : 0.;
    const double64 dlnds = (ln > 0.) ? rhon_ * User()->dkrnds(e) / mun_ : 0.;
 
    // TODO: check whether this is mathematically correct
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative", "method not properly implemented yet.");
    return ( dlwds*ln2 + dlnds*lw2 ) / lt2;
  }
  

  
  
  


  
  /**
      Saturation derivative of mobility product.
   */
  template<size_t dim, template<size_t> class USER>
  double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative_at( Element<dim>* const e, double64 sw ) const
  {
    assert( e != nullptr );
    const double64 srw  = e->Read( User()->key_srH2O );
    const double64 srn  = e->Read( User()->key_srCO2 );
    // mobility product is zero at endmember saturations
    if ( sw < srw ) return 0.;
    if ( sw > 1. - srn ) return 0.;
    
    const double64 lw  = Mobility_at(e,0U,sw);
    const double64 ln  = Mobility_at(e,1U,sw);
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    // TODO: check whether this is mathematically correct
    const double64 dlwds = (lw > 0.) ? rhow_ * User()->dkrwds_at(e, sw) / muw_ : 0.;
    const double64 dlnds = (ln > 0.) ? rhon_ * User()->dkrnds_at(e, sw) / mun_ : 0.;
    
    return ( dlwds * ln2 + dlnds * lw2 ) / lt2;
  }
  

  
  
  



/**
    Computes the fractional mass flow of the wetting (phase=0) or non-wetting
    (phase=1) phases using the relative k's. and viscosities. Note that
    Initialize() must be called first.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::f( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( Mobility( e, phase ) > 0. );
    assert( TotalMobility(e) > 0. );
    
    return Mobility( e, phase ) / TotalMobility(e);
  }
  

  
  
  
  

/**

Computes the fractional flow of the wetting (phase=0) and non-wetting
(phase=1) phases using the prescribed water saturation value, in stead of the intepolated value.

*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::f_at( Element<dim>* const e, size_t phase, double64 sw ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( sw >= 0. && sw <= 1. );

    const double64 srw = e->Read( User()->key_srH2O );
    const double64 srn = e->Read( User()->key_srCO2 );

    if ( phase == 0 ) {
         if ( sw <= srw ) return 0.;
         else if ( sw  >= 1. - srn ) return 1.;
      }
    else if ( phase == 1 ) {
         if ( sw <= srw ) return 1.;
         else if ( sw  >= 1. - srn ) return 0.;
      }

    assert( Mobility_at( e, phase, sw ) > 0. );
    assert( TotalMobility_at(e,sw) > 0. );

    return Mobility_at( e, phase, sw) / TotalMobility_at(e, sw);
 }
  

  
  
  
  

/**
    Derivative of fractional flow function with regard to the water saturation.
    Positivie if phase=aqueous, negative otherwise.
 
    @todo TODO: test whether this produces plausible results
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    InterpolateSaturations( e, sw_, sn_ );
    const double64 srw = e->Read( User()->key_srH2O );
    const double64 srn = e->Read( User()->key_srCO2 );
    // the fractional flow derivative is zero beyond the end-point saturations
    if ( sw_ <= srw ) return 0.;
    if ( sn_ <= srn ) return 0.;
 
    // if all phase are present
    double64 halite_sat = InterpolateSystem( e, sw_, rhow_, muw_, sn_, rhon_, mun_ );
    assert( fabs(sw_ + sn_ + halite_sat - 1.) <numeric_limits<double64>::epsilon() );
   
    const double64 lw  = Mobility(e,0U);
    const double64 ln  = Mobility(e,1U);
    double64       lt  = lw + ln;
    const double64 lt2 = lt * lt;
    
    const double64 dlwds = (lw > 0.) ? rhow_ * (User()->dkrwds_at(e,sw_) / muw_) : 0.;
    const double64 dlnds = (ln > 0.) ? rhon_ * (User()->dkrnds_at(e,sw_) / mun_) : 0.;
    
    return (dlwds * ln - dlnds * lw) / lt2;
 }
  

  
  
  
  

/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
// TODO: refactor because it falls down in single phase case where the physical properties of the phase are not known
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    assert( sw >= 0. and sw <= 1. );

    // the fractional flow derivative is zero beyond the end-point saturations
    const double64 srw = e->Read( User()->key_srH2O );
    const double64 srn = e->Read( User()->key_srCO2 );
    if      ( sw <= srw )       return 0.;
    else if ( sw  >= 1. - srn ) return 0.;
   
    // if all phase are present
    InterpolateSystem( e, sw_, rhow_, muw_, sn_, rhon_, mun_ );
   
    const double64 lw  = Mobility_at(e,0U,sw);
    const double64 ln  = Mobility_at(e,1U,sw);
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
   
    const double64 dlwds = (sw  > 0.) ? rhow_ * (User()->dkrwds_at(e,sw) / muw_) : 0.;
    const double64 dlnds = (sn_ > 0.) ? rhon_ * (User()->dkrnds_at(e,sw) / mun_) : 0.;
   
    return (dlwds * ln - dlnds * lw) / lt2;
 }
  

  
  
  
  
  

/**
     fractional flow derivative for wetting phase = 0.
     
     If not overloaded, this returns the derivative of the fractional flow
     function at the current saturation of the wetting phase (see Helmig, 1997,
     p. 108, eqn. 3.74, term 2 (first part).
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::AdvectionMultiplier( Element<dim>* const e ) const
{
   assert( e != nullptr );
   return dfds( e, 0U );
}


  
  
  
  
/**
     Mass-based formulation for pressure equation:  k * (kri/mui * rho_i^2 + kri/muj * rhoj^2) g grad_Y
 
     @note Since the gravity term is a vector variable in the case of lower dimensional elements
     it is written back to the model.
 
     @param dip_vc the resulting gravity term is returned into the supplied vector.
 */
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::GravityTerm( Element<dim>* const e,
                                                        VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
    // TODO: use the specific acceleration of gravity that is stored on the actual model.
    e->Read( User()->key_dip, dip_vc );
    assert( !isnan(e->Read(User()->key_kV)) );

    // aqueous phase
    double64 gravity_term = (InterpolateAqueousPhase( e, sw_, rhow_, muw_ )) ? rhow_ * rhow_ * (User()->krw_at(e,sw_) / muw_) : 0.;
    // if there is only aqueous phase
    if ( sw_ >= 1. - numeric_limits<double64>::epsilon() ) {
         dip_vc *= e->Read(User()->key_kV) * User()->acceleration_of_gravity_ * gravity_term;
         return;
      }
   
    // carbonic phase
    InterpolateCarbonicPhase( e, sn_, rhon_, mun_ );
    if ( sn_ > 0. ) gravity_term += rhon_ * rhon_ * (User()->krn_at(e,sw_) / mun_);

    dip_vc *= e->Read(User()->key_kV) * User()->acceleration_of_gravity_ * gravity_term;
}




/**
     Mass-based formulation for transport equation:  k * (kri/mui * rho_i^2) g grad_Y
 
     @note Since the gravity term is a vector variable in the case of lower dimensional elements
     it is written back to the model.
 
     @param dip_vc the resulting gravity term is returned into the supplied vector.
 */
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::GravityMultiplier_phase( Element<dim>* const e,
                                                                    VectorVariable<dim>& dip_vc, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    e->Read( User()->key_dip, dip_vc );
   
    // aqueous phase
    if ( phase == 0U ) {
        if ( InterpolateAqueousPhase( e, sw_, rhow_, muw_ ) )
          dip_vc *= e->Read(User()->key_kV) * (User()->krw_at(e,sw_) / muw_) * rhow_ * rhow_ * User()->acceleration_of_gravity_;
        else
          dip_vc *= 0.;
        return;
     }
   
    // carbonic phase
    sw_ = InterpolateAqueousPhaseSaturation( e );
    if ( InterpolateCarbonicPhase( e, sn_, rhon_, mun_ ) )
      dip_vc *= e->Read(User()->key_kV) * (User()->krn_at(e,sw_) / mun_) * rhon_ * rhon_ * User()->acceleration_of_gravity_;
    else dip_vc *= 0.;
}  
    

  
  
  
/**
     Gravity term multiplied with lambda_overbar for the mass-based transport scheme.
 */
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::GravityMultiplier_G( Element<dim>* const e,
                                                                VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    GravityTerm(e,dip_vc);
    dip_vc *= MobilityProduct(e);
 }


  
  

  
  
/**
 Computes multiplier for advection gravity coefficient. The divergence
 lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
 p. 108, eqn. 3.74, term 2 (second part).
 */
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::GravityMultiplier_dGds( Element<dim>* const e,
                                                                   VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    GravityTerm(e,dip_vc);
    dip_vc *= MobilityProductDerivative(e);
 }
  

  
  

  
  
/**
     Returns the diffusivity of CO2 in the porous medium saturated with water,
     using the diffusivity calculated by the PhaseStateFinder as an input.
     The calculation of the diffusivity of water in CO2 (phase == 0) is not implemented yet.
 
     The scaling of the diffusivity is done with the porosity, not the permeability.
 
     @note For phase=1 this coefficient is defined only for aq, aq_carb and aq_carb_salt phase states,
     else zero is returned.
 
     @discussion Diffusion caclulations make sense only for a higher order accurate transport scheme
     because the first-order scheme already is more diffusive than is physically realistic.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
   
    if ( phase == 0U )
      throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier",
                            "diffusion coefficient of water in CO2-saturated porous medium not implemented yet");
    assert( phase == 1U );
    SYSTEM_STATE state = static_cast<SYSTEM_STATE>(e->Read( User()->key_ePHS ));
    if (  state == SYSTEM_STATE::aq || state == SYSTEM_STATE::aq_carb || state == SYSTEM_STATE::aq_carb_salt )
      {
         // TODO: we need a meaningful expression here for p,T dependent diffusion in a mass-based framework
         return e->Read( User()->key_phi ) * e->Read( User()->key_diff ); //  rho, muw etc.
      }
   
    return 0.;
 }
  

  
  
  

  
  
/**
    When capillary spreading is simulated using the FE method, it can be modelled as a nonlinear diffusion process.
    Nonlinear because the driving gradient varies with saturation.
 
    qm_capillary = rhon * k * lambda_overbar * dpc/dsw * grad_sw
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        diffusion multiplier

    See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
    permeability in direction of flow  x  lambda_overbar  x pc-gradient.
 
    @note lambda overbar contains the fluid densities and viscosities.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( Element<dim>* const e ) const
{
   assert( e != nullptr );
  
   // dealing with the endpoints where lambda_overbar is zero
   const double64 lambda_overbar = MobilityProduct(e);
   if ( lambda_overbar <= numeric_limits<double64>::epsilon() ) return lambda_overbar;
   assert( User()->key_k.type == SCALAR );
   const double64 k     = e->Read(User()->key_k);
   assert( !isnan(k) );
   const double64 dpcds = User()->dpcds(e);
   assert( !isnan(dpcds) );
   return k * lambda_overbar * dpcds;
}
 
 
 
 
 
// TODO: create versions for tensor permeability
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier_Phase( Element<dim>* const e, size_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( User()->key_k.type == SCALAR );
    
    return  e->Read(User()->key_k) * ( (phase==0U) ? Mobility( e, 0U ) : Mobility( e, 1U ) ) * User()->dpcds(e);
  }
    



  // =============================================================================================================
  
  // NODE-BASED COMPUTATIONS using element parameters, but saturations and fluid properties from the current node

  // =============================================================================================================
  
    /// density * lambda = kri(sw)/mu_i  of the phase i: 0 for water, 1 for the non-wetting phase
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Mobility( Element<dim>* const eptr, size_t n, size_t phase ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    assert( phase == 0U || phase == 1U );
   
    // aqueous case
    if ( phase == 0U )
      return eptr->N(n)->Read(User()->key_rhoH2O) * (User()->krw_at(eptr,eptr->N(n)->Read(User()->key_sH2O)) / eptr->N(n)->Read(User()->key_muH2O));
   
    // carbonic phase
    return eptr->N(n)->Read(User()->key_rhoCO2) * (User()->krn_at(eptr,eptr->N(n)->Read(User()->key_sH2O)) / eptr->N(n)->Read(User()->key_muCO2));
   
 } // end Mobility




    /// d lambda_i / dsw
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityDerivative( Element<dim>* const eptr, size_t n, size_t phase ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    assert( phase == 0U || phase == 1U );

    // aqueous phase
    if ( phase == 0U ) {
          const double64 sw(eptr->N(n)->Read(User()->key_sH2O)),
                         rhow(eptr->N(n)->Read(User()->key_rhoH2O)),
                         muw(eptr->N(n)->Read(User()->key_muH2O));
      
          if ( sw <= eptr->Read( User()->key_srH2O ) ) return 0.;
          return rhow * (User()->dkrwds_at(eptr,sw) / muw);
       }
   
    // carbonic phase
    const double64 sw(eptr->N(n)->Read(User()->key_sH2O)),
                   snw(eptr->N(n)->Read(User()->key_sCO2)),
                   rhon(eptr->N(n)->Read(User()->key_rhoCO2)), mun(eptr->N(n)->Read(User()->key_muCO2));

    if ( sw >= 1. - eptr->Read( User()->key_srCO2 ) ) return 0.;
    assert( User()->krn_at(eptr,sw) >= 0. );
    assert( User()->krn_at(eptr,sw) <= 1. );

    return rhon * (User()->dkrnds_at(eptr,sw) / mun);

 } // end MobilityDerivative
 
 
 
 
    /// lambda_t: sum of phase-mobility * density products
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TotalMobility( Element<dim>* const eptr, size_t n ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );

    // if there is a mobile aqueous phase and the carbonic phase is immobile
    const double64 sw(eptr->N(n)->Read(User()->key_sH2O));
    double64 mob_t = ( sw > 0. ) ? eptr->N(n)->Read(User()->key_rhoH2O) * User()->krw_at(eptr,sw) / eptr->N(n)->Read(User()->key_muH2O) : 0.;
    const double64 srn = eptr->Read( User()->key_srCO2 );
    if ( sw >= (1. - srn) ) return mob_t;

    // if there is a mobile carbonic phase as well
    const double64 snw(eptr->N(n)->Read(User()->key_sCO2));
    if ( snw > 0. ) mob_t += eptr->N(n)->Read(User()->key_rhoCO2) * User()->krn_at(eptr,sw) / eptr->N(n)->Read(User()->key_muCO2);
   
    return mob_t;

 } // end TotalMobility
 
 
  
    /// lambda overbar: mobility product l_overbar = (li * rhow * lj * rhonw) / (li*rhow + lj*rhonw),  also known as G
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProduct( Element<dim>* const eptr, size_t node ) const
 {
    assert( eptr != nullptr );
    assert( node < eptr->Nodes() );
   
    const double64 lambda_w = Mobility( eptr, node, 0U );
    const double64 lambda_n = Mobility( eptr, node, 1U );

    return (lambda_w * lambda_n) / (lambda_w + lambda_n);
   
 } // end MobilityProduct
 
 

    /// d lambda overbar / dsw also known as dGds
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative( Element<dim>* const e, size_t n, bool evaluate_numerically ) const
 {
    assert( e != nullptr );
    assert( n < e->Nodes() );

    // product is zero at endmember saturations
    const double64 sw(e->N(n)->Read(User()->key_sH2O));
    if ( User()->EffectiveSaturation_at(e,sw) <= 0. || User()->EffectiveSaturation_at(e,sw) >= 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically )
      throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative", "numerical version not implemented yet");

    const double64 lw   = Mobility( e, n, 0U );
    const double64 ln   = Mobility( e, n, 1U );
    const double64 rhow = (sw > 0.) ? e->N(n)->Read(User()->key_rhoH2O) : 0.;
    const double64 srn  = e->Read( User()->key_srCO2 );
    const double64 rhon = (srn > 0.) ? e->N(n)->Read(User()->key_rhoCO2) : 0.;
    const double64 lt   = lw * rhow + ln * rhon;
    const double64 lt2  = lt * lt;
    const double64 ln2  = ln * ln;
    const double64 lw2  = lw * lw;
   
    const double64 dlwds = User()->dkrwds(e) / e->N(n)->Read(User()->key_muH2O);
    const double64 dlnds = User()->dkrnds(e) / e->N(n)->Read(User()->key_muCO2);
 
    return ( rhow * dlwds*ln2 + rhon * dlnds*lw2 ) / lt2;

 } // end MobilityProductDerivative
 
 
 

     /// fractional mass flow; 0=water, 1=non-wetting phase
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::f( Element<dim>* const eptr, size_t node, size_t phase ) const
 {
    assert( eptr != nullptr );
    assert( node < eptr->Nodes() );
    assert( phase == 0U || phase == 1U );
   
    const double64 lw = Mobility( eptr, node, 0U );
    const double64 ln = Mobility( eptr, node, 1U );

    if ( phase == 0U ) return lw / (lw + ln);
      
    return ln / (lw + ln);

 } // end f
 
 
 

/**
    Derivative of fractional flow function w.r.t water saturation; note that fw+fn=1, dfw_dsw = -dfn_dsn.
 
    @attention In the absence of gravitational forces dfds is positive for the derivative w.r.t. the fractional flow of the aqueous phase
    and negative with regard to the fractional flow of the carbonic phase.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds( Element<dim>* const e, size_t n, size_t phase ) const
 {
    assert( e != nullptr );
    assert( n < e->Nodes() );
    assert( phase == 0U || phase == 1U );

    // the fractional flow derivative is zero beyond the end-point saturations
    const double64 sw  = e->N(n)->Read(User()->key_sH2O);
    const double64 sn  = e->N(n)->Read(User()->key_sCO2);
    const double64 srw = e->Read( User()->key_srH2O );
    const double64 srn = e->Read( User()->key_srCO2 );
   
    // if there is no change of the fractional flow with saturation
    if ( sw <= srw ) return 0.;
    if ( sn <= srn ) return 0.;
 
    // if both phases are mobile
    const double64 rhow = (sw > 0.) ? e->N(n)->Read(User()->key_rhoH2O) : 0.;
    const double64 rhon = (sn > 0.) ? e->N(n)->Read(User()->key_rhoCO2) : 0.;
    const double64 muw = e->N(n)->Read(User()->key_muH2O);
    const double64 mun = e->N(n)->Read(User()->key_muCO2);
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    const double64 lw  = rhow * (User()->krw_at(e,sw) / muw);
    const double64 ln  = rhon * (User()->krn_at(e,sw) / mun);
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
   
    const double64 dlwds = rhow * (User()->dkrwds_at(e,sw) / muw);
    const double64 dlnds = rhon * (User()->dkrnds_at(e,sw) / mun);
   
    const double dfds_value = (dlwds * ln - dlnds * lw) / lt2;
    return (phase == 0U) ? dfds_value : -dfds_value;

 } // end dfds
 
 
 

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::AdvectionMultiplier( Element<dim>* const eptr, size_t n ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::AdvectionMultiplier", "numerical version not implemented yet");

 } // end AdvectionMultiplier
 
 
 

    /// k * kri(sw)/mi * rho_i^2 projected onto the dip vector of the current element; writes result to dip vector
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::GravityTerm(Element<dim>* const eptr, size_t n, VectorVariable<dim>& dip_vc ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
    // TODO: use the specific acceleration of gravity that is stored on the actual model.
    eptr->Read( User()->key_dip, dip_vc );
    assert( !isnan(eptr->Read(User()->key_kV)) );

    // aqueous phase
    const double64 sw(eptr->N(n)->Read(User()->key_sH2O)),
                   rhow(eptr->N(n)->Read(User()->key_rhoH2O)),
                   muw(eptr->N(n)->Read(User()->key_muH2O));
    double64 gravity_term = (sw > 0.) ? rhow * rhow * (User()->krw_at(eptr,sw) / muw) : 0.;
     // if there is only aqueous phase
    if ( sw >= 1. - numeric_limits<double64>::epsilon()*2. ) {
         dip_vc *= eptr->Read(User()->key_kV) * User()->acceleration_of_gravity_ * gravity_term;
         return;
      }

    // carbonic phase
    const double64 snw(eptr->N(n)->Read(User()->key_sCO2)),
                   rhon(eptr->N(n)->Read(User()->key_rhoCO2)),
                   mun(eptr->N(n)->Read(User()->key_muCO2));
    if ( snw > 0. ) gravity_term += rhon * rhon * (User()->krn_at(eptr,sw) / mun);

    dip_vc *= eptr->Read(User()->key_kV) * User()->acceleration_of_gravity_ * gravity_term;

 } // end GravityTerm
 
 
 

    /// average mass diffusion coefficient for CO2 in the aqueous phase
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier( Element<dim>* const eptr, size_t n, size_t phase ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    assert( phase == 0U || phase == 1U );
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier", "numerical version not implemented yet");

 } // end DiffusionMultiplier
 
 
 
  
    /// mass diffusion coefficient for the non-linear diffusion of saturation due to the saturation dependent dpc/ds
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( Element<dim>* const eptr, size_t n ) const
 {
    assert( eptr != nullptr );
    assert( n < eptr->Nodes() );
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier", "numerical version not implemented yet");

 } // end CapillaryDiffusionMultiplier

  
  
  // ===============================================================================================
  
  // Numerical derivatives (not cross-checked yet)
  
  // ===============================================================================================
  
  template<size_t dim, template<size_t> class USER>
  double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds_Numerical( Element<dim>* const e, size_t phase, double64 h ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );

     const double64 dSedSw( 1.0/ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
     const double64 seff(User()->EffectiveSaturation(e));
     
     if ( seff < 0.+h )
       return ( f_at( e, phase, seff + h ) - f_at( e, phase, seff ) ) / h * dSedSw;
    
     if ( seff > 1.-h )
       return ( f_at( e, phase, seff ) - f_at( e, phase, seff - h ) ) / h * dSedSw;
     
     return ( f_at( e, phase, seff+h ) - f_at( e, phase, seff-h ) )/ (2.0*h)* dSedSw;
     

  }
  

  
  
  
  
  
  

  /**
   Derivative of fractional flow function at specific sw, for a particular phase (0 = wetting, 1 = non-wetting)
   
   @todo check whether code for end-member cases has to be reinstated.
   */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds_at_Numerical( Element<dim>* const e, double64 sw, double64 h) const
  {
    assert( e != nullptr );
    double64 Denumerator = TotalMobility_at(e,sw);
    double64 dDenumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U) + User()->dkrnds_at_Numerical(e,sw,h)/User()->Viscosity(e, 1U);
    
    double64 Numerator = Mobility_at(e,0U,sw);
    double64 dNumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U);
    
    return (dNumerator*Denumerator-dDenumerator*Numerator)/(Denumerator*Denumerator);
  }


  
  
  

template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dGds_Numerical( Element<dim>* const e, double64 h ) const
  {
    assert( e != nullptr );
    const double64 lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;

    const double64 dlwds = User()->dkrwds_Numerical(e, h ) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds_Numerical(e, h ) / User()->Viscosity( e, 1U );

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;

  }


  



/// derivative of wetting phase mobility
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dlwds_Numerical( Element<dim>* const e, double64 h ) const
 {
    assert( e != nullptr );
    return User()->dkrwds_Numerical( e, h ) / User()->Viscosity( e, 0U );
 }


  
  
  
  


/// derivative of non-wetting phase mobility
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dlnds_Numerical( Element<dim>* const e, double64 h ) const
 {
    assert( e != nullptr );
    return User()->dkrnds_Numerical( e, h ) / User()->Viscosity( e, 1U );

 }


  
  
  



/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::InflectionPointSaturation( Element<dim>* const e ) const
  {
    assert( e != nullptr );
    double64 S = 1.-e->Read(User()->key_srCO2) ;
    
    double64 Swmin = e->Read(User()->key_srH2O) ;
    double64 Swmax = 1.0-e->Read(User()->key_srCO2);
    
    double64 DS(0.001);
    double64 Fold(-10000.);
    
    int Maxiter(4) ;
    int it(1) ;
    
    while (it < Maxiter){
      
      double64 F1 = dfds_at(e, S);
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = dfds_at(e, S);
        
        if((S<Swmin)||(S>Swmax)) return S=0;
        
      }
      
      S = S + 2*DS ;
      DS = DS/10. ;
      it++ ;
      
      if((S<Swmin)||(S>Swmax)) return S=0;
      
    }
    
    S = S-10*DS ;
    
    if((S<Swmin)||(S>Swmax)) return S=0;
    
    
    return S ;
}
  

  
  
  
  
  
  
  
  
/**
 
  The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.

 */

template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TangentPointSaturation( Element<dim>* const e ) const
  {
    assert( e != nullptr );
    const double64 Si = InflectionPointSaturation(e);
    const double64 Sf = 1-e->Read(User()->key_srCO2) ;
    double64 St = FindRootSecantMethod(e,Si,Sf) ;
    
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
}


  
  
  
  
  
  
/**

 The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TangentOfFractionalFlowFunction( Element<dim>* const e, double64 S ) const
  {
     assert( e != nullptr );
     const double64 srH2O =e->Read(User()->key_srH2O);
     return dfds_at(e, S) - (f_at(e, 0U, S) - f_at(e, 0U, srH2O)) / (S - srH2O);
  }
  

  
  
  
  

  
  
/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::FindRootSecantMethod( Element<dim>* const e, double64 S1, double64 S2 ) const
  {
    assert( e != nullptr );
    double64 F1 = 10000.;

    while ( abs(F1)>1e-10 ) {
      
        F1 = TangentOfFractionalFlowFunction(e, S1);
        double64 F2 = TangentOfFractionalFlowFunction(e, S2);
        
        double64 NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
        
        S2 = S1;
        S1 = NewPoint;
      }
    
    return  S1;
}


  
  
  
  
  
  
/**
 @attension: Generally "MaxFractionalqFlowDerivative" is not speed. It is after Buckley-Leverett problem
 */

template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MaxFractionalFlowDerivative( Element<dim>* const e ) const
{
    assert( e != nullptr );
  double64 S = InflectionPointSaturation(e);   // This is correct maximum fractional flow derivative of water phase.
  
  return dfds_at(e, S );
}
  

  
  
  


  
/// shock speed base on Buckley Leverett theory
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::ShockSpeed( Element<dim>* const e ) const
{
    assert( e != nullptr );
  return ShockFrontVelocity(e) ;
}
  

  
  

  
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::ShockHeight( Element<dim>* const e ) const
{
    assert( e != nullptr );
  //return TangentPointSaturation(e);
  
    double64  dfds_max(0.), dfds_s_max(0.), s_shock(1.), dfds; 
    double64 swr = e->Read(User()->key_srH2O);
    double64 snr = e->Read(User()->key_srCO2);    
       
    // loop over the saturation interval finding the maximum value of the fractional flow derivative
    // note the bounds! - only within these dfds is actually defined
    for ( double64 sw=swr; sw<=(1.-snr); sw+=0.005 ) {
        dfds = dfds_at(e,sw);
        // max fractional flow derivative
        dfds_max = std::max( dfds_max, dfds );
        // dfds at shock front and shock height
        double64 dfds_s = dfds * sw;
        if ( dfds_s > dfds_s_max ) {
            dfds_s_max = dfds_s;
            s_shock = sw;
        }
    }
    return s_shock;    
}
  

  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
void H2O_CO2_NaCl_FlowFunctions<dim,USER>::ShockSpeedAndHeight( Element<dim>* const e, double64& speed, double64& height) const
{
    assert( e != nullptr );
  height = ShockHeight(e) ;
  speed  = ShockSpeed(e)  ;
  
}


  
  
  


/**

  The Shock front wave calculated after estimation of tangent Saturation point.
  
*/

template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::ShockFrontVelocity( Element<dim>* const e ) const
  {
     assert( e != nullptr );
    double64 S = TangentPointSaturation(e);
    S = std::min( std::max( S, 0. ), 1. );
    
    return dfds_at(e, S);
}
  

  
  
  
template class H2O_CO2_NaCl_FlowFunctions<1U,FlowFunctionsModule4>;
template class H2O_CO2_NaCl_FlowFunctions<2U,FlowFunctionsModule4>;
template class H2O_CO2_NaCl_FlowFunctions<3U,FlowFunctionsModule4>;
  

template class H2O_CO2_NaCl_FlowFunctions<1U,FlowFunctionsModule5>;
template class H2O_CO2_NaCl_FlowFunctions<2U,FlowFunctionsModule5>;
template class H2O_CO2_NaCl_FlowFunctions<3U,FlowFunctionsModule5>;


template class H2O_CO2_NaCl_FlowFunctions<1U,FlowFunctionsModule6>;
template class H2O_CO2_NaCl_FlowFunctions<2U,FlowFunctionsModule6>;
template class H2O_CO2_NaCl_FlowFunctions<3U,FlowFunctionsModule6>;


} // end namespace csmp





