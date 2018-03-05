#include "FlowFunctions.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

/**
    Mobility of phase i, kri(sw) / mu_i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::Mobility( TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      return this->krw(p) / this->Viscosity( p, 0U );
   
    return this->krn(p) / this->Viscosity( p, 1U );
 }

template double64 FlowFunctions<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctions<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::Mobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::Mobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;



/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::MobilityDerivative( TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );
   
    if ( evaluate_numerically )
      return (phase == 0U) ? dlwds_Numerical(p) : dlnds_Numerical(p);

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( ERROR, "FlowFunctions<dim>::MobilityDerivative:", "analytic version not implemented yet." );

    // TODO: deal with the non-numerical case
    return numeric_limits<double64>::quiet_NaN();
 }
 
template double64 FlowFunctions<3U>::MobilityDerivative( FiniteElementPlacement<3U,NODE>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::MobilityDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::MobilityDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, bool ) const;




/**
    Sum of mobilities (not multiplied with permeability).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::TotalMobility( TARGET_PLACEMENT& p ) const
 {
    return this->krn(p) / this->Viscosity( p, 1U )
         + this->krw(p) / this->Viscosity( p, 0U );
 }

template double64 FlowFunctions<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::TotalMobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::TotalMobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::MobilityProduct( TARGET_PLACEMENT& p ) const
 {
    assert( key_k.type == SCALAR );
    return Mobility(p,0U) * Mobility(p,1U) / TotalMobility(p);
 }

template double64 FlowFunctions<3U>::MobilityProduct( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::MobilityProduct( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::MobilityProduct( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::MobilityProduct( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



/** 
    Saturation derivative of mobility product.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::MobilityProductDerivative( TARGET_PLACEMENT& p, bool evaluate_numerically ) const
 {
    // product is zero at endmember saturations
    if ( this->EffectiveSaturation(p) <= 0. || this->EffectiveSaturation(p) >= 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(p);

    const double64 lw  = this->krw(p) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;

    const double64 dlwds = this->dkrwds(p) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds(p) / this->Viscosity( p, 1U );

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
 }

template double64 FlowFunctions<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT>&, bool ) const;
template double64 FlowFunctions<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;



/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::f( TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
   
    return Mobility( p, phase ) / TotalMobility(p);
 }

template double64 FlowFunctions<3U>::f( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctions<3U>::f( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::f( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::f( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

 

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::dfds( TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );

    // end-member derivatives do not require relperms
//    if ( this->EffectiveSaturation(sw) <= 0. || this->EffectiveSaturation(sw) >= 1. )
//      return static_cast<double64>(0.);
   
    if ( evaluate_numerically ) return dfds_Numerical(p,phase,true);

    const double64 lw  = this->krw(p) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;

    const double64 dlwds = this->dkrwds(p) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds(p) / this->Viscosity( p, 1U );

    return ( dlwds*ln - dlnds*lw ) / lt2;
 }

template double64 FlowFunctions<3U>::dfds( FiniteElementPlacement<3U,ELEMENT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::dfds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::dfds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 FlowFunctions<3U>::dfds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, bool ) const;




template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::MaxFractionalFlowDerivative( TARGET_PLACEMENT& p ) const
 {
    double64 speed, height;
    ShockSpeedHeight( p, speed, height );
    return speed;
 }

template double64 FlowFunctions<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



/// linearized fractional flow derivative
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::ShockSpeed( TARGET_PLACEMENT& p ) const
 {
    double64 speed, height;
    ShockSpeedHeight( p, speed, height );
    return speed;
 }

template double64 FlowFunctions<3U>::ShockSpeed( FiniteElementPlacement<3U,ELEMENT>& ) const;



template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::ShockHeight( TARGET_PLACEMENT& p ) const
 {
    double64 speed, height;
    ShockSpeedHeight( p, speed, height );
    return height;
 }

template double64 FlowFunctions<3U>::ShockHeight( FiniteElementPlacement<3U,ELEMENT>& ) const;



template<size_t dim>
template<class TARGET_PLACEMENT>
void FlowFunctions<dim>::ShockSpeedHeight( TARGET_PLACEMENT& p,
                                           double64& speed, double64& height ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( ERROR, "FlowFunctions<dim>::ShockSpeedHeight:", "analytic version not implemented yet." );

  double64 se(0.), dfds_s;
  speed = 0.;

  // TODO: improve this funky implementation
  while ( speed < (dfds_s=se) ) // must be non-wetting phase saturation
    {
       speed  = dfds_s;
       height = dfds_s * se;
       se += 0.01;
    }
}

template void FlowFunctions<3U>::ShockSpeedHeight( FiniteElementPlacement<3U,ELEMENT>&, double64&, double64& ) const;



/**
    fractional flow derivative for wetting phase = 0.
 
    If not overloaded, this returns the derivative of the fractional flow
    function at the current saturation of the wetting phase (see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (first part).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::AdvectionMultiplier( TARGET_PLACEMENT& p, bool evaluate_numerically ) const
 {
    if ( evaluate_numerically ) return dfds_Numerical( p, 0U );
    return dfds(p,0U,evaluate_numerically);
 }

template double64 FlowFunctions<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,ELEMENT>&, bool ) const;
template double64 FlowFunctions<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;



/** 
    k * delta_rho * g
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::GravityTerm( TARGET_PLACEMENT& p ) const
{
   // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
   const double64 delta_rho = this->Density( p, 0U ) - this->Density( p, 1U );
  
   // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
   // TODO: use the specific acceleration of gravity that is stored on the actual model.
   const double64 k_g_drho = p.Interpolate(key_kV) * -ACC_GRAVITY * delta_rho;

   // else compute result using G saturation derivative
   return k_g_drho;
}

template double64 FlowFunctions<3U>::GravityTerm( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::GravityTerm( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::GravityTerm( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::GravityTerm( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



/**
    See Sebastian Geiger's thesis (2004), closed form, i.e.
 
    G = lambda_overbar * k * delta_rho * g
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::GravityMultiplier_G( TARGET_PLACEMENT& p ) const
{
   return GravityTerm(p) * MobilityProduct(p);
}

template double64 FlowFunctions<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


 
/**
    Computes multiplier for advection gravity coefficient. The divergence
    lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (second part).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::GravityMultiplier_dGds( TARGET_PLACEMENT& p, bool evaluate_numerically ) const
{
   return GravityTerm(p) * MobilityProductDerivative(p,evaluate_numerically);
}

template double64 FlowFunctions<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,ELEMENT>&, bool ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 FlowFunctions<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;



/**
    Returns the diffusion coefficient for the phase of interest. If not
    overloaeded, the hydraulic conductivity is returned.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::DiffusionMultiplier( TARGET_PLACEMENT& p, size_t phase ) const
{
     assert( phase == 1U or phase == 2U );
  
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Interpolate(key_kfn) / ( (phase==1u) ? this->Viscosity( p, 0U ) : this->Viscosity( p, 1U ) );
} 

template double64 FlowFunctions<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctions<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctions<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;



/**
    See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
    permeability in direction of flow  x  lambda_overbar  x pc-gradient.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::CapillaryDiffusionMultiplier( TARGET_PLACEMENT& p ) const
{
   // TODO: Make sure that this is the permeability in the direction of the facet normal
   return p.Interpolate(this->key_kfn) * MobilityProduct(p) * this->dpcds(p);
} 

template double64 FlowFunctions<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctions<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctions<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;





// linearized diffusion multiplier for large-timestep calculations
/*
template<size_t dim>
double64 FlowFunctions<dim>::DiffusionCharacteristic( size_t ) const
{
cout <<"\nFlowFunctions<"<<  dim <<">::DiffusionCharacteristic (base class): ";
cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
return std::numeric_limits<double64>::quiet_NaN();
}
*/



     



// ===============================================================================================

// Numerical derivatives

// ===============================================================================================

template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::dfds_Numerical( TARGET_PLACEMENT& p, size_t phase, double64 h ) const
{
  assert( phase == 0U or phase == 1U );
  /*
  // first version: direct differentiation
  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( fw_at( seff_ + h ) - fw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( fw_at( seff_ ) - fw_at( seff_ - h ) ) / h * dSedSw;

  return ( fw_at( seff_+h ) - fw_at( seff_-h ) )/ (2.0*h)* dSedSw;

  */

  const double64 seff(this->EffectiveSaturation(p));
  // second version: mixed analytical and numerical differentiation
  const double64 lw  = this->krw(seff) / this->Viscosity( p, 0U );
  const double64 ln  = this->krn(seff) / this->Viscosity( p, 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;

  const double64 dlwds = this->dkrwds_Numerical( seff,h ) / this->Viscosity( p, 0U );
  const double64 dlnds = this->dkrnds_Numerical( seff,h ) / this->Viscosity( p, 1U );

  return ( dlwds * ln - dlnds * lw ) / lt2;
  //*/

}

template double64 FlowFunctions<3U>::dfds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctions<3U>::dfds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctions<3U>::dfds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctions<3U>::dfds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;




template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::dGds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
{
  /*
  // first version: direct differentiation
  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );
  h = 0.0000001;
  if( seff_ < 0.+h )
      return ( G_at( seff_ + h ) - G_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( G_at( seff_ ) - G_at( seff_ - h ) ) / h * dSedSw;


  return ( G_at( seff_+h ) - G_at( seff_-h ) )/ (2.0*h) * dSedSw;
  */

  // second version: mixed analytical and numerical differentiation
  const double64 seff(this->EffectiveSaturation(p));
  const double64 lw  = this->krw(seff) / this->Viscosity( p, 0U );
  const double64 ln  = this->krn(seff) / this->Viscosity( p, 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;
  const double64 ln2 = ln*ln;
  const double64 lw2 = lw*lw;

  const double64 dlwds = this->dkrwds_Numerical( h ) / this->Viscosity( p, 0U );
  const double64 dlnds = this->dkrnds_Numerical( h ) / this->Viscosity( p, 1U );

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}

template double64 FlowFunctions<3U>::dGds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 FlowFunctions<3U>::dGds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dGds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dGds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;




/// derivative of wetting phase mobility
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::dlwds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
 {
    const double64 seff(this->EffectiveSaturation(p));
    return this->dkrwds_Numerical( seff, h ) / this->Viscosity( p, 0U );
 }

template double64 FlowFunctions<3U>::dlwds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlwds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlwds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlwds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



/// derivative of non-wetting phase mobility
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctions<dim>::dlnds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
 {
    const double64 seff(this->EffectiveSaturation(p));
    return this->dkrnds_Numerical( seff, h ) / this->Viscosity( p, 1U );

 }

template double64 FlowFunctions<3U>::dlnds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlnds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlnds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 FlowFunctions<3U>::dlnds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;




// ====================================================================================

//    NON-MEMBER FUNCTIONS

// ====================================================================================


/// Interpolations
double64 spline_value( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (1. - t)*y1 + t*y2 + t*(1.-t)*( a*(1.-t) + b*t);

}


double64 spline_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (y2-y1)/( x2-x1 ) + (1.-2.*t)*( a*(1.-t)+b*t)/(x2-x1) + t*(1.-t)*(b-a)/(x2-x1);

}


double64 spline_second_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return 2.*( b-2.*a +(a-b)*3.*t)/(x2-x1)/(x2-x1);

}


template class FlowFunctions<1U>;
template class FlowFunctions<2U>;
template class FlowFunctions<3U>;


} // end namespace csmp





