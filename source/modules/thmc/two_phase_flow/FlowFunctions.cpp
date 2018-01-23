#include "FlowFunctions.h"

using namespace std;

namespace csmp {

template<size_t dim>
double64 FlowFunctions<dim>::Sw() const
 {
    return this->sw; // TODO: use Andrew's scheme to get current saturation
 }



template<size_t dim>
double64 FlowFunctions<dim>::Mobility( double64 sw, size_t phase ) const
 {
    assert( phase == 1U or phase == 2U );
    // salinity=0
    if ( phase == 1U ) return this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
    return this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1 );
 }



/**
    Mobility saturation derivative
*/
template<size_t dim>
double64 FlowFunctions<dim>::MobilityDerivative( double64 sw, size_t phase, bool evaluate_numerically ) const
 {
    if ( evaluate_numerically ) return (phase == 0U ) ? dlwds_Numerical(sw) : dlnds_Numerical(sw);
    // TODO: deal with the non-numerical case
    return numeric_limits<double64>::quiet_NaN();
 }
 
 

template<size_t dim>
double64 FlowFunctions<dim>::TotalMobility( double64 sw ) const
 {
    return this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U )
         + this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
 }



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim>
double64 FlowFunctions<dim>::MobilityProduct( double64 sw ) const
 {
    assert( key_k.type == SCALAR );
    return Mobility(sw,0U) * Mobility(sw,1U) / TotalMobility(sw);
 }



/** 
    saturation derivative of mobility product
*/
template<size_t dim>
double64 FlowFunctions<dim>::MobilityProductDerivative( double64 sw, bool evaluate_numerically ) const
 {
    // product is zero at endmember saturations
    if ( this->EffectiveSaturation(sw) < 0. || this->EffectiveSaturation(sw) > 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(sw);

    const double64 lw  = this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
    const double64 ln  = this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;

    const double64 dlwds = this->dkrwds(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
    const double64 dlnds = this->dkrnds(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
 }



/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim>
double64 FlowFunctions<dim>::f( double64 sw, size_t phase ) const
 {
    assert( phase == 1U or phase == 2U );
    return Mobility( sw, phase ) / TotalMobility(sw);
 }


 

/// derivative of fractional flow function (advection multipliers)
template<size_t dim>
double64 FlowFunctions<dim>::dfds( double64 sw, size_t phase, bool evaluate_numerically ) const
 {
    // end-member derivatives do not require relperms
    if ( this->EffectiveSaturation(sw) < 0. || this->EffectiveSaturation(sw) > 1. )
      return static_cast<double64>(0.);
   
    if ( evaluate_numerically ) return dfds_Numerical(sw,phase);

    const double64 lw  = this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
    const double64 ln  = this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;

    const double64 dlwds = this->dkrwds(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
    const double64 dlnds = this->dkrnds(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );

    return ( dlwds*ln - dlnds*lw )/lt2;

 }
 


template<size_t dim>
double64 FlowFunctions<dim>::MaxFractionalFlowDerivative() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }



/**
    fractional flow derivative for water
 
    If not overloaded, this returns the derivative of the fractional flow
    function at the current saturation of the wetting phase (see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (first part).
*/
template<size_t dim>
double64 FlowFunctions<dim>::AdvectionMultiplier( double64 sw, bool evaluate_numerically ) const
 {
    if ( evaluate_numerically ) return dfds_Numerical( sw, 0U );
    return dfds(sw,0U);
 } 



/// linearized fractional flow derivative
template<size_t dim>
double64 FlowFunctions<dim>::ShockSpeed() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }



template<size_t dim>
double64 FlowFunctions<dim>::ShockHeight() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return height;
 }



template<size_t dim>
void FlowFunctions<dim>::ShockSpeedHeight( double64& speed, double64& height )const
{
  double64 se(0.), dfds_s;
  speed = 0.;

  // TODO: improve this funky implementation
  while ( speed < (dfds_s=dfds(se,1U)) ) // must be non-wetting phase saturation
    {
       speed  = dfds_s;
       height = dfds_s * se;
       se += 0.01;
    }
}



/** 
    k * delta_rho * g
*/
template<size_t dim>
double64 FlowFunctions<dim>::GravityTerm( double64 sw ) const
{
   // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
   const double64 delta_rho = this->Density( this->Pressure(), this->Temperature() ) -
                              this->Density( this->Pressure(), this->Temperature(), 1U );
  
   const double64 k_g_drho = this->k * -this->acc_grav * delta_rho;

   // else compute result using G saturation derivative
   return k_g_drho;
}



/** See Sebastian Geiger's thesis (2004), closed form.
*/
template<size_t dim>
double64 FlowFunctions<dim>::GravityMultiplier_G( double64 sw ) const
{
   return GravityTerm(sw) * MobilityProduct(sw);
}


 
/**
    Computes multiplier for advection gravity coefficient. The divergence
    lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (second part).
*/
template<size_t dim>
double64 FlowFunctions<dim>::GravityMultiplier_dGds( double64 sw, bool evaluate_numerically ) const
{
   return GravityTerm(sw) * MobilityProductDerivative(sw,evaluate_numerically);
}



/**
    Returns the diffusion coefficient for the phase of interest. If not
    overloaeded, the hydraulic conductivity is returned.
*/
template<size_t dim>
double64 FlowFunctions<dim>::DiffusionMultiplier( double64 sw, size_t phase ) const
{
     assert( phase == 1U or phase == 2U );
     return this->k / ( (phase==1u) ?
                          this->Viscosity( this->Pressure(), this->Temperature() ) :
                          this->Viscosity( this->Pressure(), this->Temperature() ), 1U );
} 
 


/**
    See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
    permeability in direction of flow  x  lambda_overbar  x pc-gradient.
*/
template<size_t dim>
double64 FlowFunctions<dim>::CapillaryDiffusionMultiplier( double64 sw ) const
{
   return this->k * MobilityProduct(sw) * this->dpcds(sw);
} 






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
double64 FlowFunctions<dim>::dfds_Numerical( double64 sw, size_t phase, double64 h ) const
{
  /*
  // first version: direct differentiation
  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( fw_at( seff_ + h ) - fw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( fw_at( seff_ ) - fw_at( seff_ - h ) ) / h * dSedSw;

  return ( fw_at( seff_+h ) - fw_at( seff_-h ) )/ (2.0*h)* dSedSw;

  */

  // second version: mixed analytical and numerical differentiation
  const double64 lw  = this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
  const double64 ln  = this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;

  const double64 dlwds = this->dkrwds_Numerical( sw,h ) / this->Viscosity( this->Pressure(), this->Temperature() );
  const double64 dlnds = this->dkrnds_Numerical( sw,h ) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );

  return ( dlwds * ln - dlnds * lw ) / lt2;
  //*/

}

template<size_t dim>
double64 FlowFunctions<dim>::dGds_Numerical( double64 sw, double64 h ) const
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
  const double64 lw  = this->krw(sw) / this->Viscosity( this->Pressure(), this->Temperature() );
  const double64 ln  = this->krn(sw) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;
  const double64 ln2 = ln*ln;
  const double64 lw2 = lw*lw;


  const double64 dlwds = this->dkrwds_Numerical( h ) / this->Viscosity( this->Pressure(), this->Temperature() );
  const double64 dlnds = this->dkrnds_Numerical( h ) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}





/// derivative of wetting phase mobility
template<size_t dim>
double64 FlowFunctions<dim>::dlwds_Numerical( double64 sw, double64 h ) const
 {
    return this->dkrwds_Numerical( sw, h ) / this->Viscosity( this->Pressure(), this->Temperature() );

 }



/// derivative of non-wetting phase mobility
template<size_t dim>
double64 FlowFunctions<dim>::dlnds_Numerical( double64 sw, double64 h) const
 {
    return this->dkrnds_Numerical( sw, h ) / this->Viscosity( this->Pressure(), this->Temperature(), 1U );

 }





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





