#include "FlowFunctionsBC_Hysteretic.h"
#include "FlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {
  
  
  template<size_t dim>
  FlowFunctionsBC_Hysteretic<dim>::FlowFunctionsBC_Hysteretic( const PropertyDatabase<dim>& db)
  : variables::VariableSet_CO2GeoSequestration(db)
  {
  }
  
  /**
   Mobility of phase i, kri(sw) / mu_i.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::Mobility( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U ) return this->krw(p) / this->Viscosity( p, 0U );
   
    return this->krn(p) / this->Viscosity( p, 1U );
 }


  
 
template double64 FlowFunctionsBC_Hysteretic<1U>::Mobility(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::Mobility(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::Mobility(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
  
  
  
   
/**
    Mobility of phase i, kri(sw) / mu_i.
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::Mobility_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      
      return this->krw_at(p,sw)/this->Viscosity( p, 0U );
   
    return this->krn_at(p,sw)/this->Viscosity( p, 1U );
 }
  
  template double64 FlowFunctionsBC_Hysteretic<1U>::Mobility_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsBC_Hysteretic<2U>::Mobility_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsBC_Hysteretic<3U>::Mobility_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  


/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::MobilityDerivative( const  TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return this->dkrwds(p) / this->Viscosity( p, 0U );
    
    return this->dkrnds(p) / this->Viscosity( p, 1U );
  }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::MobilityDerivative(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, bool ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::MobilityDerivative(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t , bool) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::MobilityDerivative(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t , bool) const;


  /**
   Mobility saturation derivative for phase i.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::MobilityDerivative_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
  {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return this->dkrwds_at(p,sw) / this->Viscosity( p, 0U );
    
    return this->dkrnds_at(p,sw) / this->Viscosity( p, 1U );
  }
  
   template double64 FlowFunctionsBC_Hysteretic<1U>::MobilityDerivative_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
   template double64 FlowFunctionsBC_Hysteretic<2U>::MobilityDerivative_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
   template double64 FlowFunctionsBC_Hysteretic<3U>::MobilityDerivative_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  
  
  
  
  
  /**
   Sum of mobilities (not multiplied with permeability).
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::TotalMobility( const  TARGET_PLACEMENT& p ) const
  {
    return this->krn(p) / this->Viscosity( p, 1U )
         + this->krw(p) / this->Viscosity( p, 0U );
 }



template double64 FlowFunctionsBC_Hysteretic<1U>::TotalMobility(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::TotalMobility(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::TotalMobility(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;




/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::TotalMobility_at( const  TARGET_PLACEMENT& p, double64 sw ) const
 {
    return this->krn_at(p,sw) / this->Viscosity( p, 1U )
         + this->krw_at(p,sw) / this->Viscosity( p, 0U );
 }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::TotalMobility_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::TotalMobility_at(  const  FiniteElementPlacement<2U,ELEMENT>&,  double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::TotalMobility_at(  const  FiniteElementPlacement<3U,ELEMENT>&,  double64 ) const;



/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::MobilityProduct( const  TARGET_PLACEMENT& p ) const
 {
    //assert( key_k.type == SCALAR );
    return Mobility(p,0U) * Mobility(p,1U) / TotalMobility(p);
 }


template double64 FlowFunctionsBC_Hysteretic<1U>::MobilityProduct(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::MobilityProduct(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::MobilityProduct(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;



/** 
    Saturation derivative of mobility product.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::MobilityProductDerivative( const  TARGET_PLACEMENT& p, bool evaluate_numerically ) const
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
  
  


template double64 FlowFunctionsBC_Hysteretic<1U>::MobilityProductDerivative(  const  FiniteElementPlacement<1U,ELEMENT>&, bool ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::MobilityProductDerivative(  const  FiniteElementPlacement<2U,ELEMENT>&, bool ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::MobilityProductDerivative(  const  FiniteElementPlacement<3U,ELEMENT>&, bool ) const;



  
  /**
   Saturation derivative of mobility product.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::MobilityProductDerivative_at( const  TARGET_PLACEMENT& p , double64 sw) const
  {
    // product is zero at endmember saturations
    //    if ( this->EffectiveSaturation(p) <= 0. || this->EffectiveSaturation(p) >= 1. )
    //    return static_cast<double64>(0.);
    
    const double64 lw  = this->krw_at(p, sw) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn_at(p, sw) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = this->dkrwds_at(p, sw) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds_at(p, sw) / this->Viscosity( p, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  
  
  template double64 FlowFunctionsBC_Hysteretic<1U>::MobilityProductDerivative_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 FlowFunctionsBC_Hysteretic<2U>::MobilityProductDerivative_at(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 FlowFunctionsBC_Hysteretic<3U>::MobilityProductDerivative_at(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
  
  
  /*
   Permeability helper
   */
  namespace {
    
    template<class VARIABLE_TYPE> double64 scalarPermeability( const VARIABLE_TYPE& );
    
    double64 scalarPermeability( const ScalarVariable& var )
    {
      return var();
    }
    
    template<size_t dim>
    double64 scalarPermeability( const TensorVariable<dim>& var )
    {
      return var.Trace() / static_cast<double64>(dim);
    }
    
  }
  
  
  /**
   permeability
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::Permeability( const  TARGET_PLACEMENT& p ) const
  {
    assert( this->key_k.type == TENSOR or this->key_k.type == SCALAR );
    
    typename VariableTypeTraits<dim, decltype(this->key_k)::VariableType>::VariableType K;
    
    p.Obtain(this->key_k, K);
    double64 k = scalarPermeability(K);
    return k;
}

template double64 FlowFunctionsBC_Hysteretic<1U>::Permeability(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::Permeability(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::Permeability(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;




/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::f( const  TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return Mobility( p, phase ) / TotalMobility(p);
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::f(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::f(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::f(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;


  
  /**

Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using prescribed sw value, instead of intepolated value.

*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::f_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return this->Mobility_at( p, phase, sw) / this->TotalMobility_at(p, sw);
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::f_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::f_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::f_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;

  

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::dfds( const  TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
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

  
template double64 FlowFunctionsBC_Hysteretic<1U>::dfds(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t,  bool  ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dfds(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t,  bool  ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dfds(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t,  bool  ) const;


/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::dfds_at( const  TARGET_PLACEMENT& p, double64 sw ) const
 {

   const double64 lw  = this->krw_at(p, sw) / this->Viscosity( p, 0U );
   const double64 ln  = this->krn_at(p, sw) / this->Viscosity( p, 1U );
   const double64 lt  = lw + ln;
   const double64 lt2 = lt * lt;
   
   const double64 dlwds = this->dkrwds_at(p, sw) / this->Viscosity( p, 0U );
   const double64 dlnds = this->dkrnds_at(p, sw) / this->Viscosity( p, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
  }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::dfds_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;

template double64 FlowFunctionsBC_Hysteretic<2U>::dfds_at(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;

template double64 FlowFunctionsBC_Hysteretic<3U>::dfds_at(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;


  
  
  
  /*
   @attension: Generally "MaxFractionalFlowDerivative" is not speed. It is after Buckley-Leverett problem
   */
  
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::MaxFractionalFlowDerivative( const  TARGET_PLACEMENT& p ) const
  {
    double64 S = this->InflectionPointSaturation(p);   // This is correct maximum fractional flow derivative of water phase.
    
    return dfds_at(p, S );
  }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;



  
  // shock speed base on Buckley Leverett theory
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::ShockSpeed( const  TARGET_PLACEMENT& p ) const
  {
    return this->ShockFrontVelocity(p) ;
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::ShockSpeed(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::ShockSpeed(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::ShockSpeed(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;




  
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::ShockHeight( const  TARGET_PLACEMENT& p ) const
  {
    return this->TangentPointSaturation(p);
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::ShockHeight(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::ShockHeight(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::ShockHeight(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
  
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  void FlowFunctionsBC_Hysteretic<dim>::ShockSpeedHeight( const TARGET_PLACEMENT& p, double64& speed, double64& height) const
  {
    height = ShockHeight(p) ;
    speed  = ShockSpeed(p)  ;
    
  }
  
  template void FlowFunctionsBC_Hysteretic<1U>::ShockSpeedHeight(  const  FiniteElementPlacement<1U,ELEMENT>& , double64& , double64& ) const;
  template void FlowFunctionsBC_Hysteretic<2U>::ShockSpeedHeight(  const  FiniteElementPlacement<2U,ELEMENT>& , double64& , double64& ) const;
  template void FlowFunctionsBC_Hysteretic<3U>::ShockSpeedHeight(  const  FiniteElementPlacement<3U,ELEMENT>& , double64& , double64& ) const;
  
  
  /**
   fractional flow derivative for wetting phase = 0.
   
   If not overloaded, this returns the derivative of the fractional flow
   function at the current saturation of the wetting phase (see Helmig, 1997,
   p. 108, eqn. 3.74, term 2 (first part).
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::AdvectionMultiplier( const  TARGET_PLACEMENT& p, bool evaluate_numerically ) const
  {
    return dfds(p,0U, evaluate_numerically );
  }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::AdvectionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>&, bool ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::AdvectionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>&, bool ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::AdvectionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>&, bool ) const;



  
  
  /**
   k * delta_rho * g
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::GravityTerm( const  TARGET_PLACEMENT& p ) const
  {
    // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
    const double64 delta_rho = this->Density( p, 0U ) - this->Density( p, 1U );
    
    // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
    // TODO: use the specific acceleration of gravity that is stored on the actual model.
    //const double64 k_g_drho = p.Interpolate(key_kV) * -ACC_GRAVITY * delta_rho;
    const double64 k_g_drho = p.Obtain(key_kV) * -ACC_GRAVITY * delta_rho;
    
    // else compute result using G saturation derivative
    return k_g_drho;
  }
  

  
  template double64 FlowFunctionsBC_Hysteretic<1U>::GravityTerm(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 FlowFunctionsBC_Hysteretic<2U>::GravityTerm(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 FlowFunctionsBC_Hysteretic<3U>::GravityTerm(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
  /**
   See Sebastian Geiger's thesis (2004), closed form, i.e.
   
   G = lambda_overbar * k * delta_rho * g
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::GravityMultiplier_G( const  TARGET_PLACEMENT& p ) const
  {
    return GravityTerm(p) * MobilityProduct(p);
  }
  


template double64 FlowFunctionsBC_Hysteretic<1U>::GravityMultiplier_G(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::GravityMultiplier_G(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::GravityMultiplier_G(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;



  
  
  /**
   Computes multiplier for advection gravity coefficient. The divergence
   lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
   p. 108, eqn. 3.74, term 2 (second part).
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::GravityMultiplier_dGds( const  TARGET_PLACEMENT& p , bool evaluate_numerically ) const
  {
    return GravityTerm(p) * MobilityProductDerivative(p);
  }
  
  
template double64 FlowFunctionsBC_Hysteretic<1U>::GravityMultiplier_dGds(  const  FiniteElementPlacement<1U,ELEMENT>& , bool ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::GravityMultiplier_dGds(  const  FiniteElementPlacement<2U,ELEMENT>& , bool ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::GravityMultiplier_dGds(  const  FiniteElementPlacement<3U,ELEMENT>& , bool ) const;



  
  
  
  
  /**
   Returns the diffusion coefficient for the phase of interest. If not
   overloaeded, the hydraulic conductivity is returned.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::DiffusionMultiplier( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 1U or phase == 0U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(key_kfn) / ( (phase==1U) ? this->Viscosity( p, 0U ) : this->Viscosity( p, 1U ) );
  }
  

template double64 FlowFunctionsBC_Hysteretic<1U>::DiffusionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::DiffusionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::DiffusionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;


//template double64 FlowFunctionsBC_Hysteretic<1U>::DiffusionMultiplier( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
//template double64 FlowFunctionsBC_Hysteretic<2U>::DiffusionMultiplier( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
//template double64 FlowFunctionsBC_Hysteretic<3U>::DiffusionMultiplier( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;



  
  
  
  
  
  
  
  
  /**
   See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
   permeability in direction of flow  x  lambda_overbar  x pc-gradient.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::CapillaryDiffusionMultiplier( const  TARGET_PLACEMENT& p ) const
  {
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(this->key_kfn) * MobilityProduct(p) * this->dpcds(p);
  }
  

//template double64 FlowFunctionsBC_Hysteretic<1U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
//template double64 FlowFunctionsBC_Hysteretic<2U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
//template double64 FlowFunctionsBC_Hysteretic<3U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsBC_Hysteretic<1U>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;



  
  
  
  
  
  
template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::CapillaryDiffusionMultiplier_Phase( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(this->key_kfn) * ( (phase==0U) ? Mobility( p, 1U ) : Mobility( p, 0U ) )* this->dpcds(p);
  }
  
 
  template double64 FlowFunctionsBC_Hysteretic<1U>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
  template double64 FlowFunctionsBC_Hysteretic<2U>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
  template double64 FlowFunctionsBC_Hysteretic<3U>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;

  // linearized diffusion multiplier for large-timestep calculations
  /*
   template<size_t dim>
   double64 FlowFunctionsBC_Hysteretic<dim>::DiffusionCharacteristic( size_t ) const
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
  double64 FlowFunctionsBC_Hysteretic<dim>::dfds_Numerical( const  TARGET_PLACEMENT& p, size_t phase, double64 h ) const
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
    
    //const double64 seff(this->EffectiveSaturation(p));
    // second version: mixed analytical and numerical differentiation
    const double64 lw  = this->krw(p) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    
    const double64 dlwds = this->dkrwds_Numerical( p,h ) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds_Numerical( p,h ) / this->Viscosity( p, 1U );
    
    return ( dlwds * ln - dlnds * lw ) / lt2;
    //*/
    
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::dfds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dfds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dfds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;



  /**
   Derivative of fractional flow function at specific sw, for a particular phase (0 = wetting, 1 = non-wetting)
   
   @todo check whether code for end-member cases has to be reinstated.
   */
template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::dfds_at_Numerical( const  TARGET_PLACEMENT& p, double64 sw, double64 h) const
  {
    
    double64 Numerator, Denumerator, dNumerator, dDenumerator;
    
    Denumerator = TotalMobility_at(p,sw);
    dDenumerator = this->dkrwds_at_Numerical(p,sw,h)/this->Viscosity(p, 0U) + this->dkrnds_at_Numerical(p,sw,h)/this->Viscosity(p, 1U);
    
    Numerator = Mobility_at(p,0U,sw);
    dNumerator = this->dkrwds_at_Numerical(p,sw,h)/this->Viscosity(p, 0U);
    return (dNumerator*Denumerator-dDenumerator*Numerator)/(Denumerator*Denumerator);
  }

template double64 FlowFunctionsBC_Hysteretic<1U>::dfds_at_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dfds_at_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dfds_at_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64, double64 ) const;


template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsBC_Hysteretic<dim>::dGds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
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
    //const double64 seff(this->EffectiveSaturation(p));
  
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
  //const double64 seff(this->EffectiveSaturation(p));
  const double64 lw  = this->krw(p) / this->Viscosity( p, 0U );
  const double64 ln  = this->krn(p) / this->Viscosity( p, 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;
  const double64 ln2 = ln*ln;
  const double64 lw2 = lw*lw;

  const double64 dlwds = this->dkrwds_Numerical(p, h ) / this->Viscosity( p, 0U );
  const double64 dlnds = this->dkrnds_Numerical(p, h ) / this->Viscosity( p, 1U );

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}

template double64 FlowFunctionsBC_Hysteretic<1U>::dGds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dGds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dGds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;




/// derivative of wetting phase mobility
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::dlwds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
 {
    //const double64 seff(this->EffectiveSaturation(p));
    return this->dkrwds_Numerical( p, h ) / this->Viscosity( p, 0U );
 }

template double64 FlowFunctionsBC_Hysteretic<1U>::dlwds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dlwds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dlwds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;



/// derivative of non-wetting phase mobility
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::dlnds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
 {
    //const double64 seff(this->EffectiveSaturation(p));
    return this->dkrnds_Numerical( p, h ) / this->Viscosity( p, 1U );

 }

template double64 FlowFunctionsBC_Hysteretic<1U>::dlnds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<2U>::dlnds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 FlowFunctionsBC_Hysteretic<3U>::dlnds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;




/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::InflectionPointSaturation( const TARGET_PLACEMENT& p ) const
  {
    
    double64 S = 1.-p.Obtain(key_srCO2) ;
    
    double64 Swmin = p.Obtain(key_srH2O) ;
    double64 Swmax = 1.0-p.Obtain(key_srCO2);
    
    double64 DS(0.001);
    double64 Fold(-10000.);
    
    int Maxiter(4) ;
    int it(1) ;
    
    while (it < Maxiter){
      
      double64 F1 = dfds_at(p, S);
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = dfds_at(p, S);
        
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

  
template double64 FlowFunctionsBC_Hysteretic<1U>::InflectionPointSaturation( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<2U>::InflectionPointSaturation( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<3U>::InflectionPointSaturation( const FiniteElementPlacement<3U,ELEMENT>& ) const ;


  
  
  
  
  
  
  
  
  
  
/**
 
  The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.

 */

template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::TangentPointSaturation( const TARGET_PLACEMENT& p ) const
  {
    
    double64 Si = this->InflectionPointSaturation(p);
    double64 Sf = 1-p.Obtain(key_srCO2) ;
    
    double64 St=FindRootSecantMethod(p,Si,Sf) ;
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
}
template double64 FlowFunctionsBC_Hysteretic<1U>::TangentPointSaturation( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<2U>::TangentPointSaturation( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<3U>::TangentPointSaturation( const FiniteElementPlacement<3U,ELEMENT>& ) const ;


  
  
  
  
  
  
  
/**

  The Shock front wave calculated after estimation of tangent Saturation point.
  
*/

template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::ShockFrontVelocity( const TARGET_PLACEMENT& p ) const
  {
    
    double64 S = this->TangentPointSaturation(p) ;
    S = std::min( std::max( S, 0. ), 1. );
    
    return dfds_at(p, S);
}
  
template double64 FlowFunctionsBC_Hysteretic<1U>::ShockFrontVelocity( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<2U>::ShockFrontVelocity( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 FlowFunctionsBC_Hysteretic<3U>::ShockFrontVelocity( const FiniteElementPlacement<3U,ELEMENT>& ) const ;


  
  
  
  
  
  
  
  
  
  
/**

 The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
  
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::BuckleyLeverettFunction( const TARGET_PLACEMENT& p, double64 S ) const
  {
     const double64 srH2O = p.Obtain(this->key_srH2O);
     return dfds_at(p, S) - (f_at(p, 0U, S) - f_at(p, 0U, srH2O)) / (S - srH2O);
  }
  
template double64 FlowFunctionsBC_Hysteretic<1U>::BuckleyLeverettFunction( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const ;
template double64 FlowFunctionsBC_Hysteretic<2U>::BuckleyLeverettFunction( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const ;
template double64 FlowFunctionsBC_Hysteretic<3U>::BuckleyLeverettFunction( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const ;


  
  
  
  
  
  
  
  
  
  
/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsBC_Hysteretic<dim>::FindRootSecantMethod( const TARGET_PLACEMENT& p, double64 S1, double64 S2 ) const
  {
    double64 F1 = 10000.;

    while ( abs(F1)>1e-10 ) {
      
        F1 = BuckleyLeverettFunction(p, S1);
        double64 F2 = BuckleyLeverettFunction(p, S2);
        
        double64 NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
        
        S2 = S1;
        S1 = NewPoint;
      }
    
    return  S1;
}

template double64 FlowFunctionsBC_Hysteretic<1U>::FindRootSecantMethod( const FiniteElementPlacement<1U,ELEMENT>&, double64 , double64 )const ;
template double64 FlowFunctionsBC_Hysteretic<2U>::FindRootSecantMethod( const FiniteElementPlacement<2U,ELEMENT>&, double64 , double64 ) const ;
template double64 FlowFunctionsBC_Hysteretic<3U>::FindRootSecantMethod( const FiniteElementPlacement<3U,ELEMENT>&, double64 , double64 )const ;

  
  
  
  

  




template class FlowFunctionsBC_Hysteretic<1U>;
template class FlowFunctionsBC_Hysteretic<2U>;
template class FlowFunctionsBC_Hysteretic<3U>;


} // end namespace csmp





