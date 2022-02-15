#include "FlowFunctionsModule.h"
#include "TwoPhaseFlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::Sw( Element<dim>* const e ) const
  {
     assert( e != nullptr );
     return e->PropertyValueAtBaryCenter( User()->key_sH2O );
  }
  

  
  
/**
      Mobility of phase i, kri(sw) / mu_i.
 */
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::Mobility( Element<dim>* const e, uint32_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( User()->Viscosity( e, phase ) > 0. );
    assert( User()->Viscosity( e, phase ) < 1. );
 
    if ( phase == 0U ) {
         assert( User()->krw(e) >= 0. );
         assert( User()->krw(e) <= 1. );
         return User()->krw(e) / User()->Viscosity( e, 0U );
      }
   
    assert( User()->krn(e) >= 0. );
    assert( User()->krn(e) <= 1. );
    return User()->krn(e) / User()->Viscosity( e, 1U );
 }


  
   
/**
    Mobility of phase i, kri(sw) / mu_i.
    Using prescribed sw value, instead of intepolated value.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::Mobility_at( Element<dim>* const e, uint32_t phase, double sw ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( User()->Viscosity( e, phase ) > 0. );
    assert( User()->Viscosity( e, phase ) < 1. );

    if ( phase == 0U ) {
         assert( User()->krw_at(e,sw) >= 0. );
         assert( User()->krw_at(e,sw) <= 1. );
         return User()->krw_at(e,sw) / User()->Viscosity( e, 0U );
      }
   
    assert( User()->krn_at(e,sw) >= 0. );
    assert( User()->krn_at(e,sw) <= 1. );
    return User()->krn_at(e,sw) / User()->Viscosity( e, 1U );
 }
  



/**
    Mobility saturation derivative for phase i.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative( Element<dim>* const e, uint32_t phase, bool evaluate_numerically ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds(e) / User()->Viscosity( e, 0U );
    
    return User()->dkrnds(e) / User()->Viscosity( e, 1U );
  }
  



  /**
   Mobility saturation derivative for phase i.
   */
  template<uint32_t dim, template<uint32_t> class USER>
  double TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative_at( Element<dim>* const e, uint32_t phase, double sw ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds_at(e,sw) / User()->Viscosity( e, 0U );
    
    return User()->dkrnds_at(e,sw) / User()->Viscosity( e, 1U );
  }
  

  
  
  
  /**
      Sum of mobilities (not multiplied with permeability).
   */
  template<uint32_t dim, template<uint32_t> class USER>
  double TwoPhaseFlowFunctions<dim,USER>::TotalMobility(Element<dim>* const e ) const
  {
    assert( e != nullptr );
    return User()->krn(e) / User()->Viscosity( e, 1U )
         + User()->krw(e) / User()->Viscosity( e, 0U );
 }




/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::TotalMobility_at( Element<dim>* const e, double sw ) const
 {
    assert( e != nullptr );
    return User()->krn_at(e,sw) / User()->Viscosity( e, 1U )
         + User()->krw_at(e,sw) / User()->Viscosity( e, 0U );
 }


  
  



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::MobilityProduct( Element<dim>* const e ) const
 {
    assert( e != nullptr );
    return Mobility(e,0U) * Mobility(e,1U) / TotalMobility(e);
 }


  
  
  
  


/** 
    Saturation derivative of mobility product.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative( Element<dim>* const e, bool evaluate_numerically ) const
 {
    assert( e != nullptr );
    // product is zero at endmember saturations
    if ( User()->EffectiveSaturation(e) <= 0. || User()->EffectiveSaturation(e) >= 1. )
      return static_cast<double>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(e);

    const double lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double lt  = lw + ln;
    const double lt2 = lt*lt;
    const double ln2 = ln*ln;
    const double lw2 = lw*lw;
    
    const double dlwds = User()->dkrwds(e) / User()->Viscosity( e, 0U );
    const double dlnds = User()->dkrnds(e) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  

  
  
  


  
  /**
   Saturation derivative of mobility product.
   */
  template<uint32_t dim, template<uint32_t> class USER>
  double TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative_at( Element<dim>* const e , double sw) const
  {
    assert( e != nullptr );
    // product is zero at endmember saturations
    //    if ( User()->EffectiveSaturation(e) <= 0. || User()->EffectiveSaturation(e) >= 1. )
    //    return static_cast<double>(0.);
    
    const double lw  = User()->krw_at(e, sw) / User()->Viscosity( e, 0U );
    const double ln  = User()->krn_at(e, sw) / User()->Viscosity( e, 1U );
    const double lt  = lw + ln;
    const double lt2 = lt*lt;
    const double ln2 = ln*ln;
    const double lw2 = lw*lw;
    
    const double dlwds = User()->dkrwds_at(e, sw) / User()->Viscosity( e, 0U );
    const double dlnds = User()->dkrnds_at(e, sw) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  

  
  
  
  
  




/**
 
Computes the fractional flow of the wetting (ehase=1) and non-wetting
(ehase=2) phases using the relative k's. and viscosities. Note that
Initialize() must be called first.  
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::f( Element<dim>* const e, uint32_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    assert( Mobility( e, phase ) > 0. );
    assert( TotalMobility(e) > 0. );
    
    return Mobility( e, phase ) / TotalMobility(e);
  }
  

  
  
  
  
  
  

/**

Computes the fractional flow of the wetting (ehase=1) and non-wetting
(ehase=2) phases using prescribed sw value, instead of intepolated value.

*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::f_at( Element<dim>* const e, uint32_t phase, double sw ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    const double srw = e->Read( User()->key_srH2O );
    const double srn = e->Read( User()->key_srCO2 );

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
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dfds( Element<dim>* const e, uint32_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );

    const double lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double lt  = lw + ln;
    const double lt2 = lt * lt;
    
    const double dlwds = User()->dkrwds(e) / User()->Viscosity( e, 0U );
    const double dlnds = User()->dkrnds(e) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
 }
  

  
  
  
  

/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dfds_at( Element<dim>* const e, double sw ) const
 {
   assert( e != nullptr );
   const double lw  = User()->krw_at(e, sw) / User()->Viscosity( e, 0U );
   const double ln  = User()->krn_at(e, sw) / User()->Viscosity( e, 1U );
   const double lt  = lw + ln;
   const double lt2 = lt * lt;
   
   const double dlwds = User()->dkrwds_at(e, sw) / User()->Viscosity( e, 0U );
   const double dlnds = User()->dkrnds_at(e, sw) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
  }
  

  
  
  
  
  

/**
     fractional flow derivative for wetting phase = 0.
     
     If not overloaded, this returns the derivative of the fractional flow
     function at the current saturation of the wetting phase (see Helmig, 1997,
     p. 108, eqn. 3.74, term 2 (first part).
 */
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::AdvectionMultiplier( Element<dim>* const e ) const
{
   assert( e != nullptr );
   return dfds( e, 0U );
}


  
  
  
  
/**
     k * delta_rho * g
 */
template<uint32_t dim, template<uint32_t> class USER>
void TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier( Element<dim>* const e,
                                                        VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    const double rhow = User()->Density(e,0U);
    const double rhon = User()->Density(e,1U);
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    
    // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
    const double delta_rho = rhow- rhon;
    
    // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
    // TODO: use the specific acceleration of gravity that is stored on the actual model.
    e->Read( User()->key_dip, dip_vc );
    const size_t v( (dim==1u) ? 0u : 1u );
    if(isnan(dip_vc(v))) { //dip vector has not been initialised
        if(dim==1u) {dip_vc(0u) = -1.;}
        else if(dim==2u) {dip_vc(0u) = 0.; dip_vc(1u) = -1.;}
        else {dip_vc(0u) = 0.; dip_vc(1u) = -1.; dip_vc(2u) = 0.;}
    }
    const double kV = (User()->key_k.type==SCALAR) ? e->Read(User()->key_k) : e->Read(User()->key_kV);
    assert( !isnan(kV) );
    dip_vc *= kV * -ACC_GRAVITY * delta_rho;  
}   
 

  

  
  
  
/**
     See Sebastian Geiger's thesis (2004), closed form, i.e.
 
     G = lambda_overbar * k * delta_rho * g
 */
template<uint32_t dim, template<uint32_t> class USER>
void TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_G( Element<dim>* const e,
                                                           VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    GravityMultiplier(e,dip_vc);
    dip_vc *= MobilityProduct(e);
 } 


  
  

  
  
/**
 Computes multiplier for advection gravity coefficient. The divergence
 lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
 p. 108, eqn. 3.74, term 2 (second part).
 */
template<uint32_t dim, template<uint32_t> class USER>
void TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_dGds( Element<dim>* const e,
                                                                   VectorVariable<dim>& dip_vc ) const
 {
    assert( e != nullptr );
    GravityMultiplier(e,dip_vc);
    dip_vc *= MobilityProductDerivative(e);
 }   
  
  

  
  
/**
 Returns the diffusion coefficient for the phase of interest. If not
 overloaeded, the hydraulic conductivity is returned.
 */
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::DiffusionMultiplier( Element<dim>* const e, size_t phase ) const
 {
    assert( e != nullptr );
    assert( phase == 1U or phase == 0U );
  
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return e->Read(User()->key_k) / ( (phase==1U) ? User()->Viscosity( e, 0U ) : User()->Viscosity( e, 1U ) );
 }
  

  
  
  

  
  
/**
 See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
 permeability in direction of flow  x  lambda_overbar  x pc-gradient.
 */
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( Element<dim>* const e ) const
{
    assert( e != nullptr );
  // TODO: Make sure that this is the permeability in the direction of the facet normal
  return e->Read(User()->key_k) * MobilityProduct(e) * User()->dpcds(e);
}
  
  
  
  
  
  
  
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier_Phase( Element<dim>* const e, size_t phase ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return  e->Read(User()->key_k) * ( (phase==0U) ? Mobility( e, 0U ) : Mobility( e, 1U ) )* User()->dpcds(e);
  }
  
 
  
  

  
  
  
  // ===============================================================================================
  
  // Numerical derivatives
  
  // ===============================================================================================
  
  template<uint32_t dim, template<uint32_t> class USER>
  double TwoPhaseFlowFunctions<dim,USER>::dfds_Numerical( Element<dim>* const e, size_t phase, double h ) const
  {
    assert( e != nullptr );
    assert( phase == 0U or phase == 1U );

     const double dSedSw( 1.0/ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
     const double seff(User()->EffectiveSaturation(e));
     
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
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dfds_at_Numerical( Element<dim>* const e, double sw, double h) const
  {
    assert( e != nullptr );
    double Denumerator = TotalMobility_at(e,sw);
    double dDenumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U) + User()->dkrnds_at_Numerical(e,sw,h)/User()->Viscosity(e, 1U);
    
    double Numerator = Mobility_at(e,0U,sw);
    double dNumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U);
    
    return (dNumerator*Denumerator-dDenumerator*Numerator)/(Denumerator*Denumerator);
  }


  
  
  

template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dGds_Numerical( Element<dim>* const e, double h ) const
  {
    assert( e != nullptr );
    const double lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double lt  = lw + ln;
    const double lt2 = lt*lt;
    const double ln2 = ln*ln;
    const double lw2 = lw*lw;

    const double dlwds = User()->dkrwds_Numerical(e, h ) / User()->Viscosity( e, 0U );
    const double dlnds = User()->dkrnds_Numerical(e, h ) / User()->Viscosity( e, 1U );

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;

  }


  



/// derivative of wetting phase mobility
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dlwds_Numerical( Element<dim>* const e, double h ) const
 {
    assert( e != nullptr );
    return User()->dkrwds_Numerical( e, h ) / User()->Viscosity( e, 0U );
 }


  
  
  
  


/// derivative of non-wetting phase mobility
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::dlnds_Numerical( Element<dim>* const e, double h ) const
 {
    assert( e != nullptr );
    return User()->dkrnds_Numerical( e, h ) / User()->Viscosity( e, 1U );

 }


  
  
  



/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::InflectionPointSaturation( Element<dim>* const e ) const
  {
    assert( e != nullptr );
    double S = 1.-e->Read(User()->key_srCO2) ;
    
    double Swmin = e->Read(User()->key_srH2O) ;
    double Swmax = 1.0-e->Read(User()->key_srCO2);
    
    double DS(0.001);
    double Fold(-10000.);
    
    int Maxiter(4) ;
    int it(1) ;
    
    while (it < Maxiter){
      
      double F1 = dfds_at(e, S);
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

template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::TangentPointSaturation( Element<dim>* const e ) const
  {
    assert( e != nullptr );
    const double Si = InflectionPointSaturation(e);
    const double Sf = 1-e->Read(User()->key_srCO2) ;
    double St = FindRootSecantMethod(e,Si,Sf) ;
    
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
}


  
  
  
  
  
  
/**

 The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::TangentOfFractionalFlowFunction( Element<dim>* const e, double S ) const
  {
     assert( e != nullptr );
     const double srH2O =e->Read(User()->key_srH2O);
     return dfds_at(e, S) - (f_at(e, 0U, S) - f_at(e, 0U, srH2O)) / (S - srH2O);
  }
  

  
  
  
  

  
  
/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::FindRootSecantMethod( Element<dim>* const e, double S1, double S2 ) const
  {
    assert( e != nullptr );
    double F1 = 10000.;

    while ( abs(F1)>1e-10 ) {
      
        F1 = TangentOfFractionalFlowFunction(e, S1);
        double F2 = TangentOfFractionalFlowFunction(e, S2);
        
        double NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
        
        S2 = S1;
        S1 = NewPoint;
      }
    
    return  S1;
}


  
  
  
  
  
  
/**
 @attension: Generally "MaxFractionalqFlowDerivative" is not speed. It is after Buckley-Leverett problem
 */

template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::MaxFractionalFlowDerivative( Element<dim>* const e ) const
{
    assert( e != nullptr );
  double S = InflectionPointSaturation(e);   // This is correct maximum fractional flow derivative of water phase.
  
  return dfds_at(e, S );
}
  

  
  
  


  
/// shock speed base on Buckley Leverett theory
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::ShockSpeed( Element<dim>* const e ) const
{
    assert( e != nullptr );
  return ShockFrontVelocity(e) ;
}
  

  
  

  
template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::ShockHeight( Element<dim>* const e ) const
{
    assert( e != nullptr );
  //return TangentPointSaturation(e);
  
    double  dfds_max(0.), dfds_s_max(0.), s_shock(1.), dfds; 
    double swr = e->Read(User()->key_srH2O);
    double snr = e->Read(User()->key_srCO2);    
       
    // loop over the saturation interval finding the maximum value of the fractional flow derivative
    // note the bounds! - only within these dfds is actually defined
    for ( double sw=swr; sw<=(1.-snr); sw+=0.005 ) {
        dfds = dfds_at(e,sw);
        // max fractional flow derivative
        dfds_max = std::max( dfds_max, dfds );
        // dfds at shock front and shock height
        double dfds_s = dfds * sw;
        if ( dfds_s > dfds_s_max ) {
            dfds_s_max = dfds_s;
            s_shock = sw;
        }
    }
    return s_shock;    
}
  

  
  
  
  
  
  
template<uint32_t dim, template<uint32_t> class USER>
void TwoPhaseFlowFunctions<dim,USER>::ShockSpeedAndHeight( Element<dim>* const e, double& speed, double& height) const
{
    assert( e != nullptr );
  height = ShockHeight(e) ;
  speed  = ShockSpeed(e)  ;
  
}


  
  
  


/**

  The Shock front wave calculated after estimation of tangent Saturation point.
  
*/

template<uint32_t dim, template<uint32_t> class USER>
double TwoPhaseFlowFunctions<dim,USER>::ShockFrontVelocity( Element<dim>* const e ) const
  {
     assert( e != nullptr );
    double S = TangentPointSaturation(e);
    S = std::min( std::max( S, 0. ), 1. );
    
    return dfds_at(e, S);
}
  

  
  
template class TwoPhaseFlowFunctions<1U,FlowFunctionsModule1>;
template class TwoPhaseFlowFunctions<2U,FlowFunctionsModule1>;
template class TwoPhaseFlowFunctions<3U,FlowFunctionsModule1>;
  

template class TwoPhaseFlowFunctions<1U,FlowFunctionsModule2>;
template class TwoPhaseFlowFunctions<2U,FlowFunctionsModule2>;
template class TwoPhaseFlowFunctions<3U,FlowFunctionsModule2>;


template class TwoPhaseFlowFunctions<1U,FlowFunctionsModule3>;
template class TwoPhaseFlowFunctions<2U,FlowFunctionsModule3>;
template class TwoPhaseFlowFunctions<3U,FlowFunctionsModule3>;


} // end namespace csmp





