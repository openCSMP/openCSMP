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

enum { WETTING_PHASE = 1U, NONWETTING_PHASE = 2U };

/** base class for 2-phase flow models

base class for 2-phase flow models (Linear, Brooks Corey, Van Genuchten, Richards)
now you can inherit BrooksCoreyWetting, BrooksCoreyNonWetting... 
 
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
 
    /// read and modify element properties when the saturation functions are history dependent
    virtual void InitializeAndStore( Element<dim>& e );

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
  
    /// fractional flow of water at water saturation sw
    double64 fw_at( double64 sw ) const;
    double64 fn_at( double64 sw ) const;

    /// relative permeabilities
    virtual double64 krn_Phase() const = 0;
    virtual double64 krw_Phase() const = 0;

    /// derivatives of relative permeabilities
    virtual double64 dkrnds_Phase() const;
    virtual double64 dkrwds_Phase() const;

    virtual double64 dlnds() const;
    virtual double64 dlwds() const;

    /// capillary pressure (limit this to 4e7, the max strength of the rock)
    /// do this by computing seff for which pc=4e7, then use this seff as
    /// a limiting value
    virtual double64 pc_Phase() const = 0;

    /// capillary pressure derivatives (treat seff as for previous function)
    virtual double64 dpcds_Phase() const = 0;

    /// inverse capillary pressure function
    virtual double64 Sw_Phase( double64 pc_Phase ) const;

    /// derivatives of inverse capillary pressure function
    virtual double64 dsdpc_Phase( double64 pc_Phase ) const;
                                  
    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    virtual double64 dfds() const;

    /// derivatives of gravitational flow (advection multipliers)
    virtual double64 dGds() const;
    
    /// Numerical derivatives
    // mobility ratio?
    double64  mobility_w_at( double64 ) const;
    double64  mobility_n_at( double64 ) const;
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

} // end namespace csmp

#endif 





















