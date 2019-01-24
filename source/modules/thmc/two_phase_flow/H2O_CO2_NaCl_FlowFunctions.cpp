#include "FlowFunctionsModule.h"
#include "H2O_CO2_NaCl_FlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
/**
    Saturation is always expected to have a value between 0..1.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Sw( Element<dim>* const e ) const
  {
     assert( e != nullptr );
     return e->PropertyValueAtBaryCenter( User()->key_sH2O );
  }
  

  
  
/**
      Mobility of phase i, lambda_i = rho_i * kri(sw) / mu_i.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::Mobility( Element<dim>* const e, size_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( User()->Viscosity( e, phase ) > 0. );
    assert( User()->Viscosity( e, phase ) < 1. );
 
    if ( phase == 0U ) {
         assert( User()->krw(e) >= 0. );
         assert( User()->krw(e) <= 1. );
         return User()->Density( e, 0U ) * (User()->krw(e) / User()->Viscosity( e, 0U ));
      }
   
    assert( User()->krn(e) >= 0. );
    assert( User()->krn(e) <= 1. );
    return User()->Density( e, 1U ) * (User()->krn(e) / User()->Viscosity( e, 1U ));
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
    assert( User()->Viscosity( e, phase ) > 0. );
    assert( User()->Viscosity( e, phase ) < 1. );

    if ( phase == 0U ) {
         assert( User()->krw_at(e,sw) >= 0. );
         assert( User()->krw_at(e,sw) <= 1. );
         return User()->Density( e, 0U ) * (User()->krw_at(e,sw) / User()->Viscosity( e, 0U ));
      }
   
    assert( User()->krn_at(e,sw) >= 0. );
    assert( User()->krn_at(e,sw) <= 1. );
    return User()->Density( e, 1U ) * (User()->krn_at(e,sw) / User()->Viscosity( e, 1U ));
 }
  



/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityDerivative( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->Density( e, 0U ) * (User()->dkrwds(e) / User()->Viscosity( e, 0U ));
    
    return User()->Density( e, 1U ) * (User()->dkrnds(e) / User()->Viscosity( e, 1U ));
  }
  



  /**
   Mobility saturation derivative for phase i.
   */
  template<size_t dim, template<size_t> class USER>
  double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityDerivative_at( Element<dim>* const e, size_t phase, double64 sw ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->Density( e, 0U ) * (User()->dkrwds_at(e,sw) / User()->Viscosity( e, 0U ));
    
    return User()->Density( e, 1U ) * (User()->dkrnds_at(e,sw) / User()->Viscosity( e, 1U ));
  }
  

  
  
  
  /**
      Sum of mass mobilities (not multiplied with permeability).
   */
  template<size_t dim, template<size_t> class USER>
  double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TotalMobility(Element<dim>* const e ) const
  {
    assert( e != nullptr );
    return User()->Density( e, 1U ) * (User()->krn(e) / User()->Viscosity( e, 1U )) +
           User()->Density( e, 0U ) * (User()->krw(e) / User()->Viscosity( e, 0U ));
 }




/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::TotalMobility_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    return User()->Density( e, 1U ) * (User()->krn_at(e,sw) / User()->Viscosity( e, 1U )) +
           User()->Density( e, 0U ) * (User()->krw_at(e,sw) / User()->Viscosity( e, 0U ));
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
    const double64 lw = Mobility(e,0U);
    const double64 ln = Mobility(e,1U);
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    return  (lw * rhow * ln * rhon) / (lw * rhow + ln * rhon);
 }




/**
   Mobility product evaluated at the users supplied water saturation.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProduct_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    const double64 lw = Mobility_at(e,0U,sw);
    const double64 ln = Mobility_at(e,1U,sw);
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    return  (lw * rhow * ln * rhon) / (lw * rhow + ln * rhon);
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

    const double64 lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    const double64 lt  = lw * rhow + ln * rhon;
    const double64 lt2 = lt * lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = User()->dkrwds(e) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds(e) / User()->Viscosity( e, 1U );
 
    // TODO: check whether this is mathematically correct
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative", "method not properly implemented yet.");
    return ( rhow * dlwds*ln2 + rhon * dlnds*lw2 ) / lt2;
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
    
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    const double64 lw  = User()->krw_at(e, sw) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn_at(e, sw) / User()->Viscosity( e, 1U );
    const double64 lt  = lw * rhow + ln * rhon;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    // TODO: check whether this is mathematically correct
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::MobilityProductDerivative_at", "method not properly implemented yet.");
    const double64 dlwds = User()->dkrwds_at(e, sw) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds_at(e, sw) / User()->Viscosity( e, 1U );
    
    return ( rhow*dlwds*ln2 + rhon*dlnds*lw2 )/lt2;
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
    assert( Mobility_at( e, phase, sw ) > 0. );
    assert( TotalMobility_at(e,sw) > 0. );

    return Mobility_at( e, phase, sw) / TotalMobility_at(e, sw);
 }
  

  
  
  
  

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo TODO: test whether this produces plausible results
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );

    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    const double64 muw = User()->Viscosity( e, 0U );
    const double64 mun = User()->Viscosity( e, 1U );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    const double64 lw  = rhow * (User()->krw(e) / muw);
    const double64 ln  = rhon * (User()->krn(e) / mun);
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
    
    const double64 dlwds = rhow * (User()->dkrwds(e) / muw);
    const double64 dlnds = rhon * (User()->dkrnds(e) / mun);
    
    return (dlwds * ln - dlnds * lw) / lt2;
 }
  

  
  
  
  

/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::dfds_at( Element<dim>* const e, double64 sw ) const
 {
    assert( e != nullptr );
    assert( sw >= 0. and sw <= 1. );
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    const double64 muw = User()->Viscosity( e, 0U );
    const double64 mun = User()->Viscosity( e, 1U );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    const double64 lw  = rhow * (User()->krw_at(e, sw) / muw);
    const double64 ln  = rhon * (User()->krn_at(e, sw) / mun);
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
   
    const double64 dlwds = User()->dkrwds_at(e, sw) / muw;
    const double64 dlnds = User()->dkrnds_at(e, sw) / mun;
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
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
    const double64 rhow = User()->Density(e,0U);
    const double64 rhon = User()->Density(e,1U);
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    const double64 muw = User()->Viscosity( e, 0U );
    const double64 mun = User()->Viscosity( e, 1U );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    const double64 termw = rhow * rhow * (User()->krw(e) / muw);
    const double64 termn = rhon * rhon * (User()->krn(e) / mun);

    // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
    // TODO: use the specific acceleration of gravity that is stored on the actual model.
    e->Read( User()->key_dip, dip_vc );
    assert( !isnan(e->Read(User()->key_kV)) );
    dip_vc *= e->Read(User()->key_kV) * -ACC_GRAVITY * (termw + termn);
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
    
    if(phase == 0U) {
        const double64 rhow = User()->Density(e,0U);
        assert( !isnan(rhow) );
        const double64 muw = User()->Viscosity( e, 0U );
        assert( !isnan(muw) );
        e->Read( User()->key_dip, dip_vc );
        dip_vc *= e->Read(User()->key_kV) * (User()->krw(e) / muw) * rhow * rhow * -ACC_GRAVITY;
    } else if (phase == 1U){
        const double64 rhon = User()->Density(e,1U);
        assert( !isnan(rhon) );
        const double64 mun = User()->Viscosity( e, 1U );
        assert( !isnan(mun) );
        e->Read( User()->key_dip, dip_vc );
        dip_vc *= e->Read(User()->key_kV) * (User()->krn(e) / mun) * rhon * rhon * -ACC_GRAVITY;    
    }  
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
     Returns the self diffusion coefficient (m2/s) for the phase of interest.
 
     Rethink: diffusion makes sense only for a higher order scheme,
     but in multicompenent transport the diffusivities of the components in the phase are needed,
     which is a vector of quantities.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 1U or phase == 0U );
 
    throw csmp::Exception( ERROR, "H2O_CO2_NaCl_FlowFunctions<dim,USER>::DiffusionMultiplier", "this method must still be implemented");
  
    // TODO: Make sure that this uses the permeability in the direction of the facet normal
    return e->Read(User()->key_k) / ( (phase==1U) ? User()->Viscosity( e, 0U ) : User()->Viscosity( e, 1U ) );
 }
  

  
  
  

  
  
/**
    When capillary spreading is simulated using the FE method, it can be modelled as a nonlinear diffusion process.
    Nonlinear because the driving gradient varies with saturation.
 
          diffusion multiplier
    Dpc = k * lambda_overbar * dpc/dsw * grad_sw
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
    See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
    permeability in direction of flow  x  lambda_overbar  x pc-gradient.
 */
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( Element<dim>* const e ) const
{
   assert( e != nullptr );
   assert( User()->key_k.type == SCALAR );
   const double64 k     = e->Read(User()->key_k);
   const double64 dpcds = User()->dpcds(e);
   assert( !isnan(k) );
   assert( !isnan(dpcds) );
   return k * MobilityProduct(e) * dpcds;
}


//by supplying mobility product
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( Element<dim>* const e, double64 lambda_overbar ) const
{
   assert( e != nullptr );
   assert( User()->key_k.type == SCALAR );
   const double64 k     = e->Read(User()->key_k);
   const double64 dpcds = User()->dpcds(e);
   assert( !isnan(k) );
   assert( !isnan(dpcds) );
   return k * lambda_overbar * dpcds;
}
 
 
  
template<size_t dim, template<size_t> class USER>
double64 H2O_CO2_NaCl_FlowFunctions<dim,USER>::CapillaryDiffusionMultiplier_Phase( Element<dim>* const e, size_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return  e->Read(User()->key_k) * ( (phase==0U) ? Mobility( e, 0U ) : Mobility( e, 1U ) )* User()->dpcds(e);
  }
    

  
  
  
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





