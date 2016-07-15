#ifndef CSMP_TWO_PHASE_MODEL_H
#define CSMP_TWO_PHASE_MODEL_H

#include "Node.h"
#include "Element.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "FiniteVolumeStencil.h"
#include "TensorVariable.h"
#include "TensorVariable1.h"
#include "TensorVariable2.h"

namespace csmp {

/** base class for 2-phase flow models

base class for 2-phase flow models (Linear, Brooks Corey, Van Genuchten, Richards)
now you can inherit BrooksCoreyWetting, BrooksCoreyNonWetting... 

@todo (3) To improve calculation speed, change relperm models from dynamic to static polymorphism (C)
 
 */
template<size_t dim>
class TwoPhaseModel {
  public:
    /// for fixed values of fluid viscosities and densities
    TwoPhaseModel( const PropertyDatabase<dim>& database,
                   double64 viscosity_nw, double64 viscosity_w,
                   double64 density_nw, double64 density_w,
                   const char* kkk,   // permeability
                   const char* sat,   // saturation wetting phase
                   const char* snr,   // residual saturation non-wetting phase
                   const char* swr,   // residual saturation wetting phase
                   const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    /// for interpolated values
    TwoPhaseModel( const PropertyDatabase<dim>& database,
                   const char* kkk,   // permeability
                   const char* mun,
                   const char* muw,
                   const char* rhn,
                   const char* rhw,
                   const char* sat,   // saturation wetting phase
                   const char* snr,   // residual saturation non-wetting phase
                   const char* swr,   // residual saturation wetting phase
                   const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false
                                     
    virtual ~TwoPhaseModel();
   
    /// reads element properties snr, swr, k, mtrl_param, fluid_param and interpolates to target placement
    virtual void Initialize( const Element<dim>& e );
    
    /// interpolates node properties to element barycenter
    virtual void InitializeForBaryCenter( const Element<dim>& e );

    virtual void InitializeForNode( const Element<dim>& e,
                                    size_t node );
                                    
    /// finite element integration points
    virtual void InitializeForIntegrationPoint( size_t ip, 
                                                const Element<dim>& e );
    /// finite volume facets
    virtual void InitializeForFacetIntegrationPoint( size_t facet, size_t ip,
                                                     const Element<dim>& e );
    /// finite volume sectors
    virtual void InitializeForSectorIntegrationPoint( size_t sector, size_t ip,
                                                      const Element<dim>& e );
     
    virtual void InterpolateNodeProperties( const Element<dim>& e );

    /// to set properties on the fly for compressible flow modelling
    void SaturationWettingPhase( double64 s_wetting );
    void ViscosityNonWettingPhase( double64 visc );
    void ViscosityWettingPhase( double64 visc );
    void DensityNonWettingPhase( double64 dens );
    void DensityWettingPhase( double64 dens );

    /// Return of keys for calculation of two phase mobilities in transport
    csmp::Index WettingPhaseSaturationKey() const;
    csmp::Index NonWettingPhaseSaturationKey() const;
    csmp::Index TotalMobilityKey() const;

    /// accessors
    double64   Permeability()  const;
    void       Permeability( double64 permeability );
    TensorVariable<dim>   TensorPermeability()  const;
    void       TensorPermeability( TensorVariable<dim> permeability );
    double64   Swr() const;
    double64   Snr() const;
    void       Swr( double64 swr );
    void       Snr( double64 snr );    

    /// fluid properties
    double64   ViscosityNonWettingPhase() const;
    double64   ViscosityWettingPhase() const;
    double64   DensityNonWettingPhase() const;
    double64   DensityWettingPhase() const;
    double64   Saturation( size_t phase ) const;
    double64   MobilityPhase( size_t phase ) const; // kr_i/mu_i
    double64   ViscosityRatio() const;

    /// sum of all phase mobilities * k
    double64 TotalMobility() const;
    
    /// sum of all phase mobilities
    double64 TotalMobilityMultiplier() const;

    // multiplier for diffusion coefficient in the case of non-linear diffusion
    // uses viscosity(phase) 
    double64 DiffusionMultiplier( size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double64 CapillaryDiffusionMultiplier(  ) const;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    double64 GravityTerm() const;
    double64 GravityMultiplier_G( ) const;
    double64 GravityMultiplier_dGds( ) const;

    /// G multiplier
    double64 G() const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double64 AdvectionMultiplier() const;

    /// linearized diffusion multiplier for large-timestep calculations
    virtual double64 DiffusionCharacteristic( size_t phase ) const;

    /// linearized fractional flow derivative
    virtual double64 ShockSpeed() const;
    virtual double64 ShockHeight() const;
    void ShockSpeedHeight( double64& speed, double64& height) const;

    /// maximum value of previous derivative
    virtual double64 MaxFractionalFlowDerivative( ) const;

    /// maximum value of dpcdS
    double64 MaxCapillaryPressure( size_t phase = 1U ) const;

    /// always of the wetting phase by convention
    virtual double64  EffectiveSaturation() const;
    virtual double64  SeffToSw() const;
    virtual double64  SeffToSw( double64 seff) const;

    /// fractional flow
    double64 f_Phase( size_t phase ) const;

    /// relative permeabilities
    virtual double64 krn_Phase() const = 0;
    virtual double64 krw_Phase() const = 0;

    /// derivatives of relative permeabilities
    virtual double64 dkrnds_Phase() const;
    virtual double64 dkrwds_Phase() const;

    virtual double64 dlnds( ) const;
    virtual double64 dlwds( ) const;

    /// capillary pressure (limit this to 4e7, the max strength of the rock)
    /// do this by computing seff for which pc=4e7, then use this seff as
    /// a limiting value
    virtual double64 pc_Phase( ) const = 0;

    /// capillary pressure derivatives (treat seff as for previous function)
    virtual double64 dpcds_Phase( ) const = 0;

    /// inverse capillary pressure function
    virtual double64 Sw_Phase( double64 pc_Phase ) const;

    /// derivatives of inverse capillary pressure function
    virtual double64 dsdpc_Phase( double64 pc_Phase ) const;
                                  
    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    virtual double64 dfds( ) const;

    /// derivatives of gravitational flow (advection multipliers)
    virtual double64 dGds( ) const;
    
    /// Numerical derivatives
    double64  mobility_w_at( double64 ) const;
    double64  mobility_n_at( double64 ) const;
    double64  fw_at( double64 ) const;
    double64  fn_at( double64 ) const;
    double64  krw_at( double64 ) const;
    double64  krn_at( double64 ) const;
    double64  dkrwds_at( double64 ) const;
    double64  dkrnds_at( double64 ) const;
    double64  dfds_at( double64 se ) const;
    double64  G_at( double64 se ) const;
    double64  pc_at( double64 se ) const;
    double64  dpcds_at( double64 se) const;
    double64  Sw_at( double64 se) const;
    double64  dsdpc_at( double64 se) const;

    virtual double64 dlwds_numerical(  double64 h = 0.001 ) const;
    virtual double64 dlnds_numerical(  double64 h = 0.001 ) const;
    virtual double64 dkrwds_numerical( double64 h = 0.001 ) const;
    virtual double64 dkrnds_numerical( double64 h = 0.001 ) const;
    virtual double64 dfds_numerical (  double64 h = 0.001 ) const;
    virtual double64 dGds_numerical (  double64 h = 0.000001 ) const;
    virtual double64 dpcds_numerical(  double64 h = 0.00001 ) const;

    double64 spline_value( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;
    double64 spline_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;
    double64 spline_second_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;

    /// default is the wetting phase
    virtual void Out( size_t phase=1U ) const;
    
  protected:

    TwoPhaseModel();

    csmp::Index       perm_key_,    // permeability = scalar element property
                      sat_key_,     // saturation of the wetting phase (nodal property)
                      snr_key_,     // irreducible saturation of non-wetting phase (element property)
                      swr_key_,     // irreducible saturation of wetting phase (element property)
                      ift_key_,     // interfacial tension key
                      mun_key_,     // viscosity of non-wetting phase (nodal property)
                      muw_key_,     // viscosity of wetting phase (nodal property)
                      rhw_key_,     // density of wetting phase (nodal property)
                      rhn_key_,     // density of non-wetting phase (nodal property)
                      lt_key_;      // total mobility
    double64          swr_, snr_,   // irreducible saturations
                      k_,           // permeability
                      ift_,         // interfacial tension
                      mun_, muw_,   // viscosities of wetting and non-wetting phase
                      rhn_, rhw_;   // densities
    TensorVariable<dim> K_;         // tensor permeability

    mutable double64  sat_,         // saturation of the wetting phase
                      seff_;        // actual saturation, effective saturation
    const double64    acc_gravity_, // accelaration of gravity
                      tolerance_,   // cut_off value for expensive calculations
                      MAX_CAPILLARY_PRESSURE_,      // maximum permitted capillary pressure
                      MAX_CAPILLARY_PRESSURE_SLOPE_;// maximum permitted capillary pressure slope
    const bool        interpolate_fluid_properties_;
    const bool        sw_ro_mu_placement_;  // NODE=true ELEMENT=false
    bool              tensor_permeability_; // TENSOR=true SCALAR=false
};

template<size_t dim>
inline double64 csmp::TwoPhaseModel<dim>::Swr() const
  {
    return swr_;
  }

template<size_t dim>
inline double csmp::TwoPhaseModel<dim>::Snr() const
  {
    return snr_;
  }

template<size_t dim>
inline void csmp::TwoPhaseModel<dim>::Swr( double wettingResidual )
  {
     swr_ = wettingResidual;
  }

template<size_t dim>
inline void csmp::TwoPhaseModel<dim>::Snr( double wettingNonResidual )
  {
    snr_ = wettingNonResidual;
  }


template<size_t dim>
inline void TwoPhaseModel<dim>::SaturationWettingPhase( double64 s_wetting )
{ sat_ = s_wetting; }

template<size_t dim>
inline void TwoPhaseModel<dim>::ViscosityNonWettingPhase( double64 visc )
{ mun_ = visc; }

template<size_t dim>
inline void TwoPhaseModel<dim>::ViscosityWettingPhase( double64 visc )
{ muw_ = visc; }

template<size_t dim>
inline void TwoPhaseModel<dim>::DensityNonWettingPhase( double64 dens )
{ rhn_ = dens; }

template<size_t dim>
inline void TwoPhaseModel<dim>::DensityWettingPhase( double64 dens )
{ rhw_ = dens; }


template<size_t dim>
inline double64 TwoPhaseModel<dim>::Permeability() const
{ return k_; }

template<size_t dim>
inline void TwoPhaseModel<dim>::Permeability( double64 permeability )
  { k_ = permeability; }

template<size_t dim>
inline TensorVariable<dim> TwoPhaseModel<dim>::TensorPermeability() const
  { return K_; }

template<size_t dim>
inline void TwoPhaseModel<dim>::TensorPermeability( TensorVariable<dim> permeability )
  { K_ = permeability; }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::ViscosityNonWettingPhase() const
{ return mun_; }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::ViscosityWettingPhase() const
{ return muw_; }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::DensityNonWettingPhase() const
{ return rhn_; }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::DensityWettingPhase() const
{ return rhw_; }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::MaxCapillaryPressure( size_t phase ) const
{ 
   assert( phase == 1U or phase == 2U );
   if ( phase == 2U ) return MAX_CAPILLARY_PRESSURE_;
   return static_cast<double64>(0.); 
}


template<size_t dim>
inline double64 TwoPhaseModel<dim>::Saturation( size_t phase ) const
{ 
   assert( phase == 1U or phase == 2U );
   if ( phase == 2U ) return static_cast<double64>(1.) - sat_;
   return sat_;
}


template<size_t dim>
inline double64 TwoPhaseModel<dim>::MobilityPhase( size_t phase ) const 
 {
    assert( phase == 1U or phase == 2U );
    if ( phase == 1U ) return krw_Phase() / muw_;
    return krn_Phase() / mun_;
 }


template<size_t dim>
inline double64 TwoPhaseModel<dim>::TotalMobilityMultiplier() const 
 {
    return krn_Phase() / mun_ + krw_Phase() / muw_;
 }

template<size_t dim>
double64 TwoPhaseModel<dim>::TotalMobility() const 
 {
    return k_ * TotalMobilityMultiplier();
 }

template<size_t dim>
inline double64 TwoPhaseModel<dim>::ViscosityRatio() const
{
  return mun_ / muw_;
}

/** 
    calculates lambda_t (total mobility) = sum of phase mobilities

    @todo SKM averages tensor properties, but should create tensor mobilities in stead
*/
template<size_t dim>
inline void TwoPhaseModel<dim>::Initialize( const Element<dim>& e )
 {
    swr_ = e.Read( swr_key_ );
    snr_ = e.Read( snr_key_ );

    if( tensor_permeability_){
        e.Read( perm_key_, K_);
        // TODO: Skm: fix these unwanted averages of the tensor k
        k_ = K_.Trace()/static_cast<double64>(dim);
    }else{
        k_ = e.Read( perm_key_ );
        K_.operator=( VectorVariable<dim>(PLAIN, k_ ) );
    }
    
    // if properties are discretized on element
    if (!sw_ro_mu_placement_) {
        sat_ = e.Read( sat_key_ );
        mun_ = e.Read( mun_key_ );
        muw_ = e.Read( muw_key_ );
        rhn_ = e.Read( rhn_key_ );
        rhw_ = e.Read( rhw_key_ );
    }
    
 } 
 
 
 
template<size_t dim>
inline void TwoPhaseModel<dim>::InitializeForBaryCenter( const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        e.N_AtBaryCenter( e.FE()->NRST );
        InterpolateNodeProperties( e );
    }
 }  


template<size_t dim>
inline void TwoPhaseModel<dim>::InitializeForNode( const Element<dim>& e,
                                                   size_t fem_node )
 {
    if (sw_ro_mu_placement_) {
        sat_ = e.N(fem_node)->Read( sat_key_ );
        mun_ = e.N(fem_node)->Read( mun_key_ );
        muw_ = e.N(fem_node)->Read( muw_key_ );
        rhn_ = e.N(fem_node)->Read( rhn_key_ );
        rhw_ = e.N(fem_node)->Read( rhw_key_ );
    }

 }  


template<size_t dim>
inline void TwoPhaseModel<dim>::InitializeForIntegrationPoint( size_t ip, 
                                                               const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        e.N_AtIntegrationPoint( ip, e.FE()->NRST );
        InterpolateNodeProperties( e );
    }
    
 } // end 
 
 
 
 
template<size_t dim>
inline void TwoPhaseModel<dim>::InitializeForFacetIntegrationPoint( size_t facet, size_t ip,
                                                                    const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        (e).N_AtFacetIntegrationPoint( facet, ip );  
        InterpolateNodeProperties( e );
    }

 }  



template<size_t dim>
inline void TwoPhaseModel<dim>::InitializeForSectorIntegrationPoint( size_t sector, size_t ip,
                                                                     const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        (e).N_AtSectorIntegrationPoint( sector, ip );  
        InterpolateNodeProperties( e );
    }

 } 





 

/**
 
The effective saturation is always computed w.r.t. the wetting phase. 
(compare with Helmig 1997, p. 71 and note that phase 1 is the wetting 
phase; note also that Initialize() must be called first to get input 
parameters like residual saturations).  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::EffectiveSaturation() const 
 {
    return seff_ = std::min( std::max( (sat_ - swr_) / (1. - swr_ - snr_), 0. ), 1. );
    
 } // end EffectiveSaturation


/// return wetting-phase saturation based on effective saturation
template<size_t dim>
inline double64 TwoPhaseModel<dim>::SeffToSw() const
 {
    return sat_ = seff_*(1. - swr_ - snr_) + swr_;

 }

/// return wetting-phase saturation based on effective saturation ( useful for numerical calculations )
template<size_t dim>
inline double64 TwoPhaseModel<dim>::SeffToSw( double64 seff ) const
 {
    return sat_ = seff*(1. - swr_ - snr_) + swr_;

 }

/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::f_Phase( size_t phase ) const
 {
    assert( phase == 1U or phase == 2U );
      if ( phase == 1U )
        return (krw_Phase() / muw_) / TotalMobilityMultiplier();

      return   (krn_Phase() / mun_) / TotalMobilityMultiplier();
 }

      
                                  
/**
Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn 
and de Neef (1998). Note that Initialize() must be called first.  
 */
template<size_t dim>
inline double64 TwoPhaseModel<dim>::G() const
 {
    const double64 lambda_w(krw_Phase() / muw_),
             lambda_n(krn_Phase() / mun_);
             
	return (lambda_w * lambda_n) / (lambda_w + lambda_n);
 }
 



/**
Returns the diffusion coefficient for the phase of interest. If not 
overloaeded, the hydraulic conductivity is returned.  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::DiffusionMultiplier( size_t phase ) const
{
     assert( phase == 1U or phase == 2U );
     return k_ / ( (phase==1u) ? muw_ : mun_ );
} 
 


/**
See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
permeability in direction of flow  x  lambda_overbar  x pc-gradient.  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::CapillaryDiffusionMultiplier( ) const
{
   return k_ * G() * dpcds_Phase( );
} 



/**
 
If not overloaded, this returns the derivative of the fractional flow
function at the current saturation of the wetting phase (see Helmig, 1997, 
p. 108, eqn. 3.74, term 2 (first part).  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::AdvectionMultiplier( ) const
 {
    return dfds();
 } 

template<size_t dim>
inline double64 TwoPhaseModel<dim>::GravityTerm() const
{
  // note that the projected gravity acts opposite the y-axis
  const double64 k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);

    // economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double64>(0.);

    // else compute result using G saturation derivative
    return k_g_drho;
}

/** See Sebastian Geiger's thesis (2004), closed form.
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::GravityMultiplier_G( ) const
{
  // note that the projected gravity acts opposite the y-axis
  const double64 k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);
	
	// economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double64>(0.);

	// else compute result using G saturation derivative
	return k_g_drho * G();	
} 


 
/**
 
Computes multiplier for advection gravity coefficient. The divergence 
lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997, 
p. 108, eqn. 3.74, term 2 (second part).  
*/
template<size_t dim>
inline double64 TwoPhaseModel<dim>::GravityMultiplier_dGds( ) const
{
    // SKM flow equations worked out with Adrian
    const double64 k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);
	
	// economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double64>(0.);

	// else compute result using G saturation derivative
    return k_g_drho * dGds(  );
} 

template<size_t dim>
inline csmp::Index TwoPhaseModel<dim>::WettingPhaseSaturationKey() const
{
  return sat_key_;
}

template<size_t dim>
inline csmp::Index TwoPhaseModel<dim>::NonWettingPhaseSaturationKey() const
{
    throw csmp::Exception( ERROR, "TwoPhaseModel<dim>::NonWettingPhaseSaturationKey()",
                           "It appears someone decided this should be the irreducible non-wet. phase saturation. Wrong." );
    return snr_key_;
}

template<size_t dim>
inline csmp::Index TwoPhaseModel<dim>::TotalMobilityKey() const
{
  return lt_key_;
}

} // end namespace csmp

#endif 





















