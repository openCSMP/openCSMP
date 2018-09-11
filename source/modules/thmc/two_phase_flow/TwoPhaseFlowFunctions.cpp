#include "FlowFunctions1.h"
#include "TwoPhaseFlowFunctions.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {
  
/**
      Mobility of phase i, kri(sw) / mu_i.
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::Mobility( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U ) return User()->krw(p) / User()->Viscosity( p, 0U );
   
    return User()->krn(p) / User()->Viscosity( p, 1U );
 }


  
 
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<1U,NODE>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<2U,NODE>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteVolumePlacement<3U,NODE>&, size_t ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;  
  
   
/**
    Mobility of phase i, kri(sw) / mu_i.
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::Mobility_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      
      return User()->krw_at(p,sw) / User()->Viscosity( p, 0U );
   
    return User()->krn_at(p,sw) / User()->Viscosity( p, 1U );
 }
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<1U,NODE>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<2U,NODE>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::Mobility_at(  const  FiniteVolumePlacement<3U,NODE>&, size_t, double64 ) const;

/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative( const  TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds(p) / User()->Viscosity( p, 0U );
    
    return User()->dkrnds(p) / User()->Viscosity( p, 1U );
  }
  

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t , bool) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t , bool) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityDerivative(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, bool ) const;


  /**
   Mobility saturation derivative for phase i.
   */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::MobilityDerivative_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
  {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return User()->dkrwds_at(p,sw) / User()->Viscosity( p, 0U );
    
    return User()->dkrnds_at(p,sw) / User()->Viscosity( p, 1U );
  }
  
   template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityDerivative_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
   template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityDerivative_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
   template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityDerivative_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  
  
  
  
  
  /**
      Sum of mobilities (not multiplied with permeability).
   */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::TotalMobility( const  TARGET_PLACEMENT& p ) const
  {
    return User()->krn(p) / User()->Viscosity( p, 1U )
         + User()->krw(p) / User()->Viscosity( p, 0U );
 }



template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::TotalMobility_at( const  TARGET_PLACEMENT& p, double64 sw ) const
 {
    return User()->krn_at(p,sw) / User()->Viscosity( p, 1U )
         + User()->krw_at(p,sw) / User()->Viscosity( p, 0U );
 }
  

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility_at(  const  FiniteElementPlacement<2U,ELEMENT>&,  double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility_at(  const  FiniteElementPlacement<3U,ELEMENT>&,  double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<1U,NODE>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<2U,NODE>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::TotalMobility_at(  const  FiniteVolumePlacement<3U,NODE>&, double64 ) const;


/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProduct( const  TARGET_PLACEMENT& p ) const
 {
    //assert( key_k.type == SCALAR );
    return Mobility(p,0U) * Mobility(p,1U) / TotalMobility(p);
 }


template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityProduct(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityProduct(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProduct(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProduct(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


/** 
    Saturation derivative of mobility product.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative( const  TARGET_PLACEMENT& p, bool evaluate_numerically ) const
 {
    // product is zero at endmember saturations
    if ( User()->EffectiveSaturation(p) <= 0. || User()->EffectiveSaturation(p) >= 1. )
      return static_cast<double64>(0.);

    if ( evaluate_numerically ) return dGds_Numerical(p);

    const double64 lw  = User()->krw(p) / User()->Viscosity( p, 0U );
    const double64 ln  = User()->krn(p) / User()->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = User()->dkrwds(p) / User()->Viscosity( p, 0U );
    const double64 dlnds = User()->dkrnds(p) / User()->Viscosity( p, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  
  


template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<1U,ELEMENT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<2U,ELEMENT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<3U,ELEMENT>&, bool ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProductDerivative(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;

  
  /**
   Saturation derivative of mobility product.
   */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::MobilityProductDerivative_at( const  TARGET_PLACEMENT& p , double64 sw) const
  {
    // product is zero at endmember saturations
    //    if ( User()->EffectiveSaturation(p) <= 0. || User()->EffectiveSaturation(p) >= 1. )
    //    return static_cast<double64>(0.);
    
    const double64 lw  = User()->krw_at(p, sw) / User()->Viscosity( p, 0U );
    const double64 ln  = User()->krn_at(p, sw) / User()->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = User()->dkrwds_at(p, sw) / User()->Viscosity( p, 0U );
    const double64 dlnds = User()->dkrnds_at(p, sw) / User()->Viscosity( p, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  
  
  template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MobilityProductDerivative_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MobilityProductDerivative_at(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MobilityProductDerivative_at(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
  
  




/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::f( const  TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return Mobility( p, phase ) / TotalMobility(p);
  }
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::f(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::f(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::f(  const  FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::f(  const  FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f(  const  FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;


  
  /**

Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using prescribed sw value, instead of intepolated value.

*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::f_at( const  TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    
    return Mobility_at( p, phase, sw) / TotalMobility_at(p, sw);
  }
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::f_at(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::f_at(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f_at(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::f_at(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::f_at(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::f_at(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const; 
  

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::dfds( const  TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
 {
    assert( phase == 0U or phase == 1U );

    if ( evaluate_numerically ) return dfds_Numerical(p,phase,true);

    const double64 lw  = User()->krw(p) / User()->Viscosity( p, 0U );
    const double64 ln  = User()->krn(p) / User()->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
    
    const double64 dlwds = User()->dkrwds(p) / User()->Viscosity( p, 0U );
    const double64 dlnds = User()->dkrnds(p) / User()->Viscosity( p, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
 }

  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t,  bool  ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t,  bool  ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t,  bool  ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, bool ) const;



/**
    Derivative of fractional flow function at specific sw, for wetting phase 0
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::dfds_at( const  TARGET_PLACEMENT& p, double64 sw ) const
 {

   const double64 lw  = User()->krw_at(p, sw) / User()->Viscosity( p, 0U );
   const double64 ln  = User()->krn_at(p, sw) / User()->Viscosity( p, 1U );
   const double64 lt  = lw + ln;
   const double64 lt2 = lt * lt;
   
   const double64 dlwds = User()->dkrwds_at(p, sw) / User()->Viscosity( p, 0U );
   const double64 dlnds = User()->dkrnds_at(p, sw) / User()->Viscosity( p, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
  }
  

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds_at(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds_at(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_at(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds_at(  const  FiniteVolumePlacement<1U,NODE>&,  double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds_at(  const  FiniteVolumePlacement<2U,NODE>&,  double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_at(  const  FiniteVolumePlacement<3U,NODE>&,  double64 ) const;
  

/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::InflectionPointSaturation( const TARGET_PLACEMENT& p ) const
  {
    double64 S = 1.-p.Obtain(User()->key_srCO2) ;
    
    double64 Swmin = p.Obtain(User()->key_srH2O) ;
    double64 Swmax = 1.0-p.Obtain(User()->key_srCO2);
    
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

  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::InflectionPointSaturation( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::InflectionPointSaturation( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::InflectionPointSaturation( const FiniteElementPlacement<3U,ELEMENT>& ) const ;

  
  
  /*
   @attension: Generally "MaxFractionalFlowDerivative" is not speed. It is after Buckley-Leverett problem
   */
  
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::MaxFractionalFlowDerivative( const TARGET_PLACEMENT& p ) const
{
  double64 S = InflectionPointSaturation(p);   // This is correct maximum fractional flow derivative of water phase.
  
  return dfds_at(p, S );
}
  

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::MaxFractionalFlowDerivative(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;

  
// shock speed base on Buckley Leverett theory
/*
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockSpeed( const  TARGET_PLACEMENT& p ) const
{
  return ShockFrontVelocity(p) ;
}
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::ShockSpeed(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::ShockSpeed(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::ShockSpeed(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;




  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::ShockHeight( const  TARGET_PLACEMENT& p ) const
  {
    return TangentPointSaturation(p);
  }
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::ShockHeight(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::ShockHeight(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::ShockHeight(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  void TwoPhaseFlowFunctions<dim,USER>::ShockSpeedHeight( const TARGET_PLACEMENT& p, double64& speed, double64& height) const
  {
    height = ShockHeight(p) ;
    speed  = ShockSpeed(p)  ;
    
  }
  
  template void TwoPhaseFlowFunctions<1U,FlowFunctions2>::ShockSpeedHeight(  const  FiniteElementPlacement<1U,ELEMENT>& , double64& , double64& ) const;
  template void TwoPhaseFlowFunctions<2U,FlowFunctions2>::ShockSpeedHeight(  const  FiniteElementPlacement<2U,ELEMENT>& , double64& , double64& ) const;
  template void TwoPhaseFlowFunctions<3U,FlowFunctions2>::ShockSpeedHeight(  const  FiniteElementPlacement<3U,ELEMENT>& , double64& , double64& ) const;
*/
  
  
  
/**
 fractional flow derivative for wetting phase = 0.
 
 If not overloaded, this returns the derivative of the fractional flow
 function at the current saturation of the wetting phase (see Helmig, 1997,
 p. 108, eqn. 3.74, term 2 (first part).
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::AdvectionMultiplier( const  TARGET_PLACEMENT& p, bool evaluate_numerically ) const
{
  return dfds(p,0U, evaluate_numerically );
}

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>&, bool ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::AdvectionMultiplier(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;

  
  
/**
 k * delta_rho * g
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityTerm( const TARGET_PLACEMENT& p ) const
{
  // note that the projected gravity acts opposite the y-axis, term rhow - rhoo
  const double64 delta_rho = User()->Density( p, 0U ) - User()->Density( p, 1U );
  
  // here the vertical permeability (key_kV) must be used since this is the direction in which gravity acts
  // TODO: use the specific acceleration of gravity that is stored on the actual model.
  //const double64 k_g_drho = p.Interpolate(key_kV) * -ACC_GRAVITY * delta_rho;
  const double64 k_g_drho = p.Obtain(User()->key_kV) * -ACC_GRAVITY * delta_rho;
  
  // else compute result using G saturation derivative
  return k_g_drho;
}

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<1U,NODE>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<2U,NODE>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<3U,NODE>& ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityTerm(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityTerm(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityTerm(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

  
  
/**
 See Sebastian Geiger's thesis (2004), closed form, i.e.
 
 G = lambda_overbar * k * delta_rho * g
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_G( const TARGET_PLACEMENT& p ) const
  {
    return GravityTerm(p) * MobilityProduct(p);
  }

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_G(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;

  
  
/**
 Computes multiplier for advection gravity coefficient. The divergence
 lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
 p. 108, eqn. 3.74, term 2 (second part).
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::GravityMultiplier_dGds( const  TARGET_PLACEMENT& p , bool evaluate_numerically ) const
  {
    return GravityTerm(p) * MobilityProductDerivative(p);
  }
  
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<1U,ELEMENT>& , bool ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<2U,ELEMENT>& , bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<3U,ELEMENT>& , bool ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::GravityMultiplier_dGds(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, bool ) const;

  
  
  
  
/**
 Returns the diffusion coefficient for the phase of interest. If not
 overloaeded, the hydraulic conductivity is returned.
 */
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::DiffusionMultiplier( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 1U or phase == 0U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(User()->key_k) / ( (phase==0U) ? User()->Viscosity( p, 0U ) : User()->Viscosity( p, 1U ) );
  }
  

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::DiffusionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::DiffusionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::DiffusionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::DiffusionMultiplier(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
//template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::DiffusionMultiplier( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
//template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::DiffusionMultiplier( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
//template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::DiffusionMultiplier( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;



  
  
  
  
  
  
  
  
  /**
   See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
   permeability in direction of flow  x  lambda_overbar  x pc-gradient.
   */
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier( const  TARGET_PLACEMENT& p ) const
  {
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(User()->key_k) * MobilityProduct(p) * User()->dpcds(p);
  }
  

//template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
//template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
//template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;

  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::CapillaryDiffusionMultiplier_Phase( const  TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(User()->key_k) * ( (phase==0U) ? Mobility( p, 0U ) : Mobility( p, 1U ) )* User()->dpcds(p);
  }
  
 
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ,size_t ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::CapillaryDiffusionMultiplier_Phase(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;  

  // linearized diffusion multiplier for large-timestep calculations
  /*
   template<size_t dim, template<size_t> class USER>
   double64 TwoPhaseFlowFunctions<dim,USER>::DiffusionCharacteristic( size_t ) const
   {
   cout <<"\nFlowFunctions<"<<  dim <<">::DiffusionCharacteristic (base class): ";
   cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
   return std::numeric_limits<double64>::quiet_NaN();
   }
   */
  
  
  // same as above
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockSaturation( const TARGET_PLACEMENT& p, size_t phase, bool evaluate_numerically ) const
{
   /*
   cerr <<"\nTwoPhaseFlowFunctions<dim,USER>::ShockSaturation: not implemented yet.\n";
 // if ( phase == 0 ) return TangentPointSaturation(p);
 // return 1. - TangentPointSaturation(p);
   return -1.;
   */
    assert( phase == 0U or phase == 1U );
    
    double64  dfds_max(0.), dfds_s_max(0.), s_shock, speed, height; 
    
    double64 swr = p.Obtain(User()->key_srH2O);
    double64 snr = p.Obtain(User()->key_srCO2);    
    
    if(phase==0U){    
        // loop over the saturation interval finding the maximum value of the fractional flow derivative
        // note the bounds! - only within these dfds is actually defined
        for ( double64 sw=swr; sw<=(1.-snr); sw+=0.005 ) {
            double64 dfds;
            if (evaluate_numerically) dfds = dfds_at_Numerical(p,sw);
            else dfds = dfds_at(p,sw);
            // max fractional flow derivative
            dfds_max = std::max( dfds_max, dfds );
            // dfds at shock front and shock height
            double64 dfds_s = dfds * sw;
            if ( dfds_s > dfds_s_max ) {
                dfds_s_max = dfds_s;
                speed = dfds;
                height = dfds_s;
                s_shock = sw;
            }
        }
        return s_shock;
    }
        
    for ( double64 sn=snr; sn<=(1.-swr); sn+=0.005 ) {
         double64 dfds;
         if (evaluate_numerically) dfds = dfds_at_Numerical(p,1.-sn); //dfndsn = dfwdsw
         else dfds = dfds_at(p,1.-sn);             
         // max fractional flow derivative
         dfds_max = std::max( dfds_max, dfds );
         // dfds at shock front and shock height
         double64 dfds_s = dfds * sn;
         if ( dfds_s > dfds_s_max ) {
              dfds_s_max = dfds_s;
              speed = dfds;
              height = dfds_s;
              s_shock = sn;
         }
    }
    return s_shock;   
}
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::ShockSaturation(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::ShockSaturation(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, bool ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::ShockSaturation(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, bool ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::ShockSaturation(  const  FiniteVolumePlacement<1U,NODE>&, size_t, bool) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::ShockSaturation(  const  FiniteVolumePlacement<2U,NODE>&, size_t, bool) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::ShockSaturation(  const  FiniteVolumePlacement<3U,NODE>&, size_t, bool) const;  
  
  
  
  // ===============================================================================================
  
  // Numerical derivatives
  
  // ===============================================================================================
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::dfds_Numerical( const  TARGET_PLACEMENT& p, size_t phase, double64 h ) const
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
    
    //const double64 seff(User()->EffectiveSaturation(p));
    // second version: mixed analytical and numerical differentiation
    const double64 lw  = User()->krw(p) / User()->Viscosity( p, 0U );
    const double64 ln  = User()->krn(p) / User()->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    
    const double64 dlwds = User()->dkrwds_Numerical( p,h ) / User()->Viscosity( p, 0U );
    const double64 dlnds = User()->dkrnds_Numerical( p,h ) / User()->Viscosity( p, 1U );
    
    return ( dlwds * ln - dlnds * lw ) / lt2;
    //*/
    
  }
  
template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_Numerical(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;



  /**
   Derivative of fractional flow function at specific sw, for a particular phase (0 = wetting, 1 = non-wetting)
   
   @todo check whether code for end-member cases has to be reinstated.
   */
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::dfds_at_Numerical( const  TARGET_PLACEMENT& p, double64 sw, double64 h) const
  {
    double64 Denumerator = TotalMobility_at(p,sw);
    double64 dDenumerator = User()->dkrwds_at_Numerical(p,sw,h)/User()->Viscosity(p, 0U) + User()->dkrnds_at_Numerical(p,sw,h)/User()->Viscosity(p, 1U);
    
    double64 Numerator = Mobility_at(p,0U,sw);
    double64 dNumerator = User()->dkrwds_at_Numerical(p,sw,h)/User()->Viscosity(p, 0U);
    
    return (dNumerator*Denumerator-dDenumerator*Numerator)/(Denumerator*Denumerator);
  }

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64, double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteVolumePlacement<1U,NODE>&, double64, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteVolumePlacement<2U,NODE>&, double64, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dfds_at_Numerical(  const  FiniteVolumePlacement<3U,NODE>&, double64, double64 ) const;



template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 TwoPhaseFlowFunctions<dim,USER>::dGds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
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
    //const double64 seff(User()->EffectiveSaturation(p));
  
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
  //const double64 seff(User()->EffectiveSaturation(p));
  const double64 lw  = User()->krw(p) / User()->Viscosity( p, 0U );
  const double64 ln  = User()->krn(p) / User()->Viscosity( p, 1U );
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;
  const double64 ln2 = ln*ln;
  const double64 lw2 = lw*lw;

  const double64 dlwds = User()->dkrwds_Numerical(p, h ) / User()->Viscosity( p, 0U );
  const double64 dlnds = User()->dkrnds_Numerical(p, h ) / User()->Viscosity( p, 1U );

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dGds_Numerical(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



/// derivative of wetting phase mobility
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::dlwds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
 {
    //const double64 seff(User()->EffectiveSaturation(p));
    return User()->dkrwds_Numerical( p, h ) / User()->Viscosity( p, 0U );
 }

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlwds_Numerical(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



/// derivative of non-wetting phase mobility
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::dlnds_Numerical( const  TARGET_PLACEMENT& p, double64 h ) const
 {
    //const double64 seff(User()->EffectiveSaturation(p));
    return User()->dkrnds_Numerical( p, h ) / User()->Viscosity( p, 1U );

 }

template double64 TwoPhaseFlowFunctions<1U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,FlowFunctions2>::dlnds_Numerical(  const  FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;  
  
  
  
  

  

template class TwoPhaseFlowFunctions<1U,FlowFunctions2>;
template class TwoPhaseFlowFunctions<2U,FlowFunctions2>;
template class TwoPhaseFlowFunctions<3U,FlowFunctions2>;


} // end namespace csmp





