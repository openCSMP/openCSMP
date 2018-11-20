#include "CO2H2O_FunctionsModule1.h"
#include "TwoPhaseFlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::Sw( const Element<dim>* const e ) const
  {
      return e->PropertyValueAtBaryCenter( User()->key_sH2O );
  }
  

  
  
/**
      Mobility of phase i, kri(sw) / mu_i.
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::Mobility( const Element<dim>* const e, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    
    // salinity=0
    if ( phase == 0U ) return User()->krw(e) / User()->Viscosity( e, 0U );
   
    return User()->krn(e) / User()->Viscosity( e, 1U );
 }


  
   
/**
    Mobility of phase i, kri(sw) / mu_i.
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::Mobility_at( const Element<dim>* const e, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      
      return User()->krw_at(e,sw) / User()->Viscosity( e, 0U );
   
    return User()->krn_at(e,sw) / User()->Viscosity( e, 1U );
 }
  



/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative( const Element<dim>* const e, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds(e) / User()->Viscosity( e, 0U );
    
    return User()->dkrnds(e) / User()->Viscosity( e, 1U );
  }
  



  /**
   Mobility saturation derivative for phase i.
   */
  template<size_t dim, template<size_t> class USER>
  double64 TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative_at( const Element<dim>* const e, size_t phase, double64 sw ) const
  {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds_at(e,sw) / User()->Viscosity( e, 0U );
    
    return User()->dkrnds_at(e,sw) / User()->Viscosity( e, 1U );
  }
  

  
  
  
  /**
      Sum of mobilities (not multiplied with permeability).
   */
  template<size_t dim, template<size_t> class USER>
  double64 TwoPhaseFlowFunctions<dim,USER>::TotalMobility(const Element<dim>* const e ) const
  {
    return User()->krn(e) / User()->Viscosity( e, 1U )
         + User()->krw(e) / User()->Viscosity( e, 0U );
 }




/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::TotalMobility_at( const Element<dim>* const e, double64 sw ) const
 {
    return User()->krn_at(e,sw) / User()->Viscosity( e, 1U )
         + User()->krw_at(e,sw) / User()->Viscosity( e, 0U );
 }


  
  



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProduct( const Element<dim>* const e ) const
 {
    //assert( key_k.type == SCALAR );
    return Mobility(e,0U) * Mobility(e,1U) / TotalMobility(e);
 }


  
  
  
  


/** 
    Saturation derivative of mobility product.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative( const Element<dim>* const e, bool evaluate_numerically ) const
 {
    // product is zero at endmember saturations
    if ( User()->EffectiveSaturation(e) <= 0. || User()->EffectiveSaturation(e) >= 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(e);

    const double64 lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = User()->dkrwds(e) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds(e) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  

  
  
  


  
  /**
   Saturation derivative of mobility product.
   */
  template<size_t dim, template<size_t> class USER>
  double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative_at( const Element<dim>* const e , double64 sw) const
  {
    // product is zero at endmember saturations
    //    if ( User()->EffectiveSaturation(e) <= 0. || User()->EffectiveSaturation(e) >= 1. )
    //    return static_cast<double64>(0.);
    
    const double64 lw  = User()->krw_at(e, sw) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn_at(e, sw) / User()->Viscosity( e, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = User()->dkrwds_at(e, sw) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds_at(e, sw) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  

  
  
  
  
  




/**
 
Computes the fractional flow of the wetting (ehase=1) and non-wetting
(ehase=2) phases using the relative k's. and viscosities. Note that
Initialize() must be called first.  
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::f( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return Mobility( e, phase ) / TotalMobility(e);
  }
  

  
  
  
  
  
  

/**

Computes the fractional flow of the wetting (ehase=1) and non-wetting
(ehase=2) phases using prescribed sw value, instead of intepolated value.

*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::f_at( const Element<dim>* const e, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return Mobility_at( e, phase, sw) / TotalMobility_at(e, sw);
 }
  

  
  
  
  
  


  

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::dfds( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );

    const double64 lw  = User()->krw(e) / User()->Viscosity( e, 0U );
    const double64 ln  = User()->krn(e) / User()->Viscosity( e, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
    
    const double64 dlwds = User()->dkrwds(e) / User()->Viscosity( e, 0U );
    const double64 dlnds = User()->dkrnds(e) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
 }
  

  
  
  
  



/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::dfds_at( const Element<dim>* const e, double64 sw ) const
 {

   const double64 lw  = User()->krw_at(e, sw) / User()->Viscosity( e, 0U );
   const double64 ln  = User()->krn_at(e, sw) / User()->Viscosity( e, 1U );
   const double64 lt  = lw + ln;
   const double64 lt2 = lt * lt;
   
   const double64 dlwds = User()->dkrwds_at(e, sw) / User()->Viscosity( e, 0U );
   const double64 dlnds = User()->dkrnds_at(e, sw) / User()->Viscosity( e, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
  }
  

  
  
  
  
  
  


/**
     fractional flow derivative for wetting phase = 0.
     
     If not overloaded, this returns the derivative of the fractional flow
     function at the current saturation of the wetting phase (see Helmig, 1997,
     p. 108, eqn. 3.74, term 2 (first part).
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::AdvectionMultiplier( const Element<dim>* const e ) const
{
  return dfds( e, 0U );
}


  
  
  

  
  
/**
     k * delta_rho * g
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityTerm( const Element<dim>* const e) const
{
  // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
  const double64 delta_rho = User()->Density( e, 0U ) - User()->Density( e, 1U );
  
  // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
  // TODO: use the specific acceleration of gravity that is stored on the actual model.
   const double64 k_g_drho = e->Read(User()->key_kV) * -ACC_GRAVITY * delta_rho;
  
  // else compute result using G saturation derivative
  return k_g_drho;
}


  
  
  
  
  


  
  
  
/**
     See Sebastian Geiger's thesis (2004), closed form, i.e.
 
     G = lambda_overbar * k * delta_rho * g
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_G( const Element<dim>* const e ) const
{
  return GravityTerm(e) * MobilityProduct(e);
}


  
  
  
  
  
  
  
  
  
/**
 Computes multiplier for advection gravity coefficient. The divergence
 lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
 p. 108, eqn. 3.74, term 2 (second part).
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_dGds( const Element<dim>* const e ) const
{
  return GravityTerm(e) * MobilityProductDerivative(e);
}
  

  
  
  
  



  
  
  
  
/**
 Returns the diffusion coefficient for the phase of interest. If not
 overloaeded, the hydraulic conductivity is returned.
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::DiffusionMultiplier( const Element<dim>* const e, size_t phase ) const
{
  assert( phase == 1U or phase == 0U );
  
  // TODO: Make sure that this is the permeability in the direction of the facet normal
  return e->Read(User()->key_k) / ( (phase==1U) ? User()->Viscosity( e, 0U ) : User()->Viscosity( e, 1U ) );
}
  

  
  
  
  
  
  
  
  
  
  
  
  
  
/**
 See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
 permeability in direction of flow  x  lambda_overbar  x pc-gradient.
 */
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( const Element<dim>* const e ) const
{
  // TODO: Make sure that this is the permeability in the direction of the facet normal
  return e->Read(User()->key_k) * MobilityProduct(e) * User()->dpcds(e);
}
  
  
  
  
  
  
  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier_Phase( const Element<dim>* const e, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return  e->Read(User()->key_k) * ( (phase==0U) ? Mobility( e, 0U ) : Mobility( e, 1U ) )* User()->dpcds(e);
  }
  
 
  
  
  
  
  
  
  
  
  
  
  // ===============================================================================================
  
  // Numerical derivatives
  
  // ===============================================================================================
  
  template<size_t dim, template<size_t> class USER>
  double64 TwoPhaseFlowFunctions<dim,USER>::dfds_Numerical( const Element<dim>* const e, size_t phase, double64 h ) const
  {
    assert( phase == 0U or phase == 1U );

     // first version: direct differentiation
    
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
double64 TwoPhaseFlowFunctions<dim,USER>::dfds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h) const
  {
    double64 Denumerator = TotalMobility_at(e,sw);
    double64 dDenumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U) + User()->dkrnds_at_Numerical(e,sw,h)/User()->Viscosity(e, 1U);
    
    double64 Numerator = Mobility_at(e,0U,sw);
    double64 dNumerator = User()->dkrwds_at_Numerical(e,sw,h)/User()->Viscosity(e, 0U);
    
    return (dNumerator*Denumerator-dDenumerator*Numerator)/(Denumerator*Denumerator);
  }


  
  
  
  
  






template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::dGds_Numerical( const Element<dim>* const e, double64 h ) const
  {
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
double64 TwoPhaseFlowFunctions<dim,USER>::dlwds_Numerical( const Element<dim>* const e, double64 h ) const
 {
    //const double64 seff(User()->EffectiveSaturation(e));
    return User()->dkrwds_Numerical( e, h ) / User()->Viscosity( e, 0U );
 }


  
  
  
  
  
  
  
  


/// derivative of non-wetting phase mobility
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::dlnds_Numerical( const Element<dim>* const e, double64 h ) const
 {
    //const double64 seff(User()->EffectiveSaturation(e));
    return User()->dkrnds_Numerical( e, h ) / User()->Viscosity( e, 1U );

 }


  
  
  
  
  
  
  
  



/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::InflectionPointSaturation( const Element<dim>* const e ) const
  {
    
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
double64 TwoPhaseFlowFunctions<dim,USER>::TangentPointSaturation( const Element<dim>* const e ) const
  {
    
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
double64 TwoPhaseFlowFunctions<dim,USER>::TangentOfFractionalFlowFunction( const Element<dim>* const e, double64 S ) const
  {
     const double64 srH2O =e->Read(User()->key_srH2O);
     return dfds_at(e, S) - (f_at(e, 0U, S) - f_at(e, 0U, srH2O)) / (S - srH2O);
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::FindRootSecantMethod( const Element<dim>* const e, double64 S1, double64 S2 ) const
  {
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
double64 TwoPhaseFlowFunctions<dim,USER>::MaxFractionalFlowDerivative( const Element<dim>* const e ) const
{
  double64 S = InflectionPointSaturation(e);   // This is correct maximum fractional flow derivative of water phase.
  
  return dfds_at(e, S );
}
  

  
  
  
  
  
  


  
/// shock speed base on Buckley Leverett theory
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockSpeed( const Element<dim>* const e ) const
{
  return ShockFrontVelocity(e) ;
}
  

  
  
  
  
  
  


  
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockHeight( const Element<dim>* const e ) const
{
  //return TangentPointSaturation(e);
  
    double64  dfds_max(0.), dfds_s_max(0.), s_shock, dfds; 
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
void TwoPhaseFlowFunctions<dim,USER>::ShockSpeedAndHeight( const Element<dim>* const e, double64& speed, double64& height) const
{
  height = ShockHeight(e) ;
  speed  = ShockSpeed(e)  ;
  
}


  
  
  
  
  



/**

  The Shock front wave calculated after estimation of tangent Saturation point.
  
*/

template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockFrontVelocity( const Element<dim>* const e ) const
  {
    
    double64 S = TangentPointSaturation(e);
    S = std::min( std::max( S, 0. ), 1. );
    
    return dfds_at(e, S);
}
  

  
  
  
template class TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule0>;
template class TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule0>;
template class TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule0>;
  

template class TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule1>;
template class TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule1>;
template class TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule1>;


template class TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>;


} // end namespace csmp





