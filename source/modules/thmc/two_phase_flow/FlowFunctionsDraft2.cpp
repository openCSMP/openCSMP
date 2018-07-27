#include "FlowFunctionsDraft2.h"
#include "Fluid.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {

  
template<size_t dim>
FlowFunctionsDraft2<dim>::FlowFunctionsDraft2( const PropertyDatabase<dim>& db)
    : variables::Variables_TwoPhaseFlow(db) 
{
}

/**
    Mobility of phase i, kri(sw) / mu_i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::Mobility( TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      return this->krw(p) / this->Viscosity( p, 0U );
   
    return this->krn(p) / this->Viscosity( p, 1U );
 }
  
template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;

/*
template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
  
template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
  
template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,NODE>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,NODE>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,NODE>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
*/
/**
    Mobility of phase i, kri(sw) / mu_i.
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::Mobility( TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
    // salinity=0
    if ( phase == 0U )
      return this->krw(p,sw) / this->Viscosity( p, 0U );
   
    return this->krn(p,sw) / this->Viscosity( p, 1U );
 }
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;

/*
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteVolumePlacement<1U,NODE>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteVolumePlacement<2U,NODE>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteVolumePlacement<3U,NODE>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::Mobility( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::Mobility( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::Mobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  */
/**
  Mobility saturation derivative for phase i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::MobilityDerivative( TARGET_PLACEMENT& p, size_t phase ) const
  {
    assert( phase == 0U or phase == 1U );
    
    if ( phase == 0U )
      return this->dkrwds(p) / this->Viscosity( p, 0U );
    
    return this->dkrnds(p) / this->Viscosity( p, 1U );
  }
  
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;

  /*
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
  
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
  
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
  */

/**
    Mobility saturation derivative for phase i.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::MobilityDerivative( TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
   
   if ( phase == 0U )
     return this->dkrwds(p,sw) / this->Viscosity( p, 0U );
   
     return this->dkrnds(p,sw) / this->Viscosity( p, 1U );
 }
  
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;

 /*
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64) const;
template double64 FlowFunctionsDraft2<1U>::MobilityDerivative( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64) const;

template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64) const;
template double64 FlowFunctionsDraft2<2U>::MobilityDerivative( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64) const;
  
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64) const;
template double64 FlowFunctionsDraft2<3U>::MobilityDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64) const;
*/



/**
    Sum of mobilities (not multiplied with permeability).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::TotalMobility( TARGET_PLACEMENT& p ) const
 {
    return this->krn(p) / this->Viscosity( p, 1U )
         + this->krw(p) / this->Viscosity( p, 0U );
 }


template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,ELEMENT>& ) const;
  
/*
template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;

  
template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
*/


/**
    Sum of mobilities (not multiplied with permeability).
    Using prescribed sw value, instead of intepolated value.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::TotalMobility( TARGET_PLACEMENT& p, double64 sw ) const
 {
    return this->krn(p,sw) / this->Viscosity( p, 1U )
         + this->krw(p,sw) / this->Viscosity( p, 0U );
 }
 /*
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& , double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,ELEMENT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,ELEMENT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT>& , double64 ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& , double64 ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& , double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::TotalMobility( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::TotalMobility( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::TotalMobility( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& , double64 ) const;
*/

/**
    G - parameter known as mobility product, lambda overbar.
 
    Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn
    and de Neef (1998). Note that Initialize() must be called first.
 */
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::MobilityProduct( TARGET_PLACEMENT& p ) const
 {
    //assert( key_k.type == SCALAR );
    return Mobility(p,0U) * Mobility(p,1U) / TotalMobility(p);
 }
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteElementPlacement<1U,ELEMENT>& ) const;

/*
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProduct( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProduct( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProduct( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
*/


/** 
    Saturation derivative of mobility product.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::MobilityProductDerivative( TARGET_PLACEMENT& p ) const
 {
    // product is zero at endmember saturations
    //if ( this->EffectiveSaturation(p) <= 0. || this->EffectiveSaturation(p) >= 1. )
     // return static_cast<double64>(0.);

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
  
template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT>& ) const;

/*
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,ELEMENT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT>& ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
*/
  /**
   Saturation derivative of mobility product.
   */
  template<size_t dim>
  template<class TARGET_PLACEMENT>
  double64 FlowFunctionsDraft2<dim>::MobilityProductDerivative( TARGET_PLACEMENT& p , double64 sw) const
  {
    // product is zero at endmember saturations
//    if ( this->EffectiveSaturation(p) <= 0. || this->EffectiveSaturation(p) >= 1. )
  //    return static_cast<double64>(0.);
    
    const double64 lw  = this->krw(p, sw) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p, sw) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;
    
    const double64 dlwds = this->dkrwds(p, sw) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds(p, sw) / this->Viscosity( p, 1U );
    
    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
  }
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;

  /*
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
  
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::MobilityProductDerivative( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::MobilityProductDerivative( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::MobilityProductDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
*/
/*
   Permeability helper
*/
namespace {

template<class VARIABLE_TYPE>
double64
scalarPermeability(VARIABLE_TYPE& var);

inline double64
scalarPermeability(ScalarVariable& var)
{
    return var();
}

template<size_t dim>
inline double64
scalarPermeability(TensorVariable<dim>& var)
{
    return var.Trace() / (double64)dim;
}

}


/**
    permeability
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::Permeability( TARGET_PLACEMENT& p ) const
{
    assert( this->key_k.type == TENSOR or this->key_k.type == SCALAR );

    typename VariableTypeTraits<dim, decltype(this->key_k)::VariableType>::VariableType K;

    p.Obtain(this->key_k, K);
    double64 k = scalarPermeability(K);
    return k;
}

  /*
template double64 FlowFunctionsDraft2<1U>::Permeability( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::Permeability( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::Permeability( FiniteElementPlacement<3U,ELEMENT>& ) const;
*/


/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::f( TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
   
    return Mobility( p, phase ) / TotalMobility(p);
 }
  
template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;


  /*
template double64 FlowFunctionsDraft2<1U>::f( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::f( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::f( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
*/

/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using prescribed sw value, instead of intepolated value.

*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::f( TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
 {
    assert( phase == 0U or phase == 1U );
   
    return Mobility( p, phase, sw) / TotalMobility(p, sw);
 }
  
  template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;


  /*
  template double64 FlowFunctionsDraft2<1U>::f( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64  ) const;
  template double64 FlowFunctionsDraft2<2U>::f( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::f( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<1U>::f( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<2U>::f( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
  template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
  template double64 FlowFunctionsDraft2<3U>::f( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
*/

/**
    Derivative of fractional flow function (advection multipliers).
 
    @todo check whether code for end-member cases has to be reinstated.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::dfds( TARGET_PLACEMENT& p, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );

    const double64 lw  = this->krw(p) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;

    const double64 dlwds = this->dkrwds(p) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds(p) / this->Viscosity( p, 1U );

    return ( dlwds*ln - dlnds*lw ) / lt2;
 }
  
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;

/*
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t ) const;
 
 */

template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::dfds( TARGET_PLACEMENT& p, size_t phase, double64 sw ) const
{
    assert( phase == 0U or phase == 1U );
    
    const double64 lw  = this->krw(p, sw) / this->Viscosity( p, 0U );
    const double64 ln  = this->krn(p, sw) / this->Viscosity( p, 1U );
    const double64 lt  = lw + ln;
    const double64 lt2 = lt * lt;
    
    const double64 dlwds = this->dkrwds(p, sw) / this->Viscosity( p, 0U );
    const double64 dlnds = this->dkrnds(p, sw) / this->Viscosity( p, 1U );
    
    return ( dlwds*ln - dlnds*lw ) / lt2;
}
  
  
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;

  
  /*
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<1U>::dfds( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<2U>::dfds( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t, double64 ) const;
template double64 FlowFunctionsDraft2<3U>::dfds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, size_t, double64 ) const;
  */

/*
 @attension: Generally "MaxFractionalFlowDerivative" is not speed. It is after Buckley-Leverett problem
*/
  
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::MaxFractionalFlowDerivative( TARGET_PLACEMENT& p ) const
 {
   double64 S ;
   S = this->InflectionPointSaturation(p) ;   // This is correct maximum fractional flow derivative of water phase.
   
   
   return dfds(p, 0U, S);
 }
  
template double64 FlowFunctionsDraft2<1U>::MaxFractionalFlowDerivative( FiniteElementPlacement<1U,ELEMENT>& ) const;
  
/*
template double64 FlowFunctionsDraft2<1U>::MaxFractionalFlowDerivative( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::MaxFractionalFlowDerivative( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::MaxFractionalFlowDerivative( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::MaxFractionalFlowDerivative( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;


template double64 FlowFunctionsDraft2<2U>::MaxFractionalFlowDerivative( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::MaxFractionalFlowDerivative( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::MaxFractionalFlowDerivative( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::MaxFractionalFlowDerivative( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  

template double64 FlowFunctionsDraft2<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::MaxFractionalFlowDerivative( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  */

// shock speed base on Buckley Leverett theory
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::ShockSpeed( TARGET_PLACEMENT& p ) const
 {
  return this->ShockFrontVelocity(p) ;
 }

template double64 FlowFunctionsDraft2<1U>::ShockSpeed( FiniteElementPlacement<1U,ELEMENT>& ) const;



template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::ShockHeight( TARGET_PLACEMENT& p ) const
 {
    return this->TangentPointSaturation(p);
 }

template double64 FlowFunctionsDraft2<1U>::ShockHeight( FiniteElementPlacement<1U,ELEMENT>& ) const;


/**
    fractional flow derivative for wetting phase = 0.
 
    If not overloaded, this returns the derivative of the fractional flow
    function at the current saturation of the wetting phase (see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (first part).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::AdvectionMultiplier( TARGET_PLACEMENT& p ) const
 {
    return dfds(p,0U);
 }
  
/*
template double64 FlowFunctionsDraft2<1U>::AdvectionMultiplier( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::AdvectionMultiplier( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::AdvectionMultiplier( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::AdvectionMultiplier( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<2U>::AdvectionMultiplier( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::AdvectionMultiplier( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::AdvectionMultiplier( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::AdvectionMultiplier( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  
template double64 FlowFunctionsDraft2<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::AdvectionMultiplier( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  */
  


/** 
    k * delta_rho * g
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::GravityTerm( TARGET_PLACEMENT& p ) const
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

  /*
template double64 FlowFunctionsDraft2<1U>::GravityTerm( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityTerm( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityTerm( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<1U>::GravityTerm( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityTerm( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityTerm( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityTerm( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<2U>::GravityTerm( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityTerm( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityTerm( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityTerm( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<3U>::GravityTerm( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityTerm( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityTerm( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityTerm( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;

*/

/**
    See Sebastian Geiger's thesis (2004), closed form, i.e.
 
    G = lambda_overbar * k * delta_rho * g
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::GravityMultiplier_G( TARGET_PLACEMENT& p ) const
{
   return GravityTerm(p) * MobilityProduct(p);
}

  /*
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_G( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_G( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_G( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_G( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_G( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_G( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_G( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_G( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_G( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_G( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_G( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_G( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  */
  
 
/**
    Computes multiplier for advection gravity coefficient. The divergence
    lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997,
    p. 108, eqn. 3.74, term 2 (second part).
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::GravityMultiplier_dGds( TARGET_PLACEMENT& p ) const
{
   return GravityTerm(p) * MobilityProductDerivative(p);
}

  /*
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_dGds( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_dGds( FiniteElementPlacement<1U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_dGds( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<1U>::GravityMultiplier_dGds( FiniteElementPlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_dGds( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_dGds( FiniteElementPlacement<2U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_dGds( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::GravityMultiplier_dGds( FiniteElementPlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
  
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::GravityMultiplier_dGds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  
*/

/**
    Returns the diffusion coefficient for the phase of interest. If not
    overloaeded, the hydraulic conductivity is returned.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::DiffusionMultiplier( TARGET_PLACEMENT& p, size_t phase ) const
{
     assert( phase == 1U or phase == 2U );
  
    // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(key_kfn) / ( (phase==1U) ? this->Viscosity( p, 0U ) : this->Viscosity( p, 1U ) );
} 

  /*
template double64 FlowFunctionsDraft2<1U>::DiffusionMultiplier( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::DiffusionMultiplier( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::DiffusionMultiplier( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::DiffusionMultiplier( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::DiffusionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
*/


/**
    See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
    permeability in direction of flow  x  lambda_overbar  x pc-gradient.
*/
template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::CapillaryDiffusionMultiplier( TARGET_PLACEMENT& p ) const
{
   // TODO: Make sure that this is the permeability in the direction of the facet normal
   return p.Obtain(this->key_kfn) * MobilityProduct(p) * this->dpcds(p);
} 

  /*
template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const;

template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
*/

template<size_t dim>
template<class TARGET_PLACEMENT>
double64 FlowFunctionsDraft2<dim>::CapillaryDiffusionMultiplier_Phase( TARGET_PLACEMENT& p, size_t phase ) const
{
    assert( phase == 0U or phase == 1U );
    
   // TODO: Make sure that this is the permeability in the direction of the facet normal
    return p.Obtain(this->key_kfn) * ( (phase==0U) ? Mobility( p, 1U ) : Mobility( p, 0U ) ) * this->dpcds(p);
} 

  /*
template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier_Phase( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier_Phase( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ,size_t ) const;
template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier_Phase( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<1U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<1U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, size_t ) const;

template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<2U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<2U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, size_t ) const;
  
template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<3U,ELEMENT>&, size_t ) const;
template double64 FlowFunctionsDraft2<3U>::CapillaryDiffusionMultiplier_Phase( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, size_t ) const;
  */
  


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

template class FlowFunctionsDraft2<1U>;
template class FlowFunctionsDraft2<2U>;
template class FlowFunctionsDraft2<3U>;


} // end namespace csmp





