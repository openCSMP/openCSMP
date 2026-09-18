// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_TWO_PHASE_MODEL_H
#define CSMP_TWO_PHASE_MODEL_H

#include "Node.h"
#include "Element.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "TensorVariable.h"
#include "TensorVariable1.h"
#include "TensorVariable2.h"

namespace csmp {

enum { WETTING_PHASE = 1U, NONWETTING_PHASE = 2U };

/** base class for 2-phase flow models

base class for 2-phase flow models (Linear, Brooks Corey, Van Genuchten, Richards)
now you can inherit BrooksCoreyWetting, BrooksCoreyNonWetting... 
 
 */
template<uint32_t dim>
class TwoPhaseModel {
  public:
    /// for fixed values of fluid viscosities and densities
    TwoPhaseModel( const PropertyDatabase<dim>& database,
                   double viscosity_nw, double viscosity_w,
                   double density_nw, double density_w,
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
                                    uint32_t node );
                                    
    /// finite element integration points
    virtual void InitializeForIntegrationPoint( uint32_t ip,
                                                const Element<dim>& e );
    /// finite volume facets
    virtual void InitializeForFacetIntegrationPoint( uint32_t facet, uint32_t ip,
                                                     const Element<dim>& e );
    /// finite volume sectors
    virtual void InitializeForSectorIntegrationPoint( uint32_t sector, uint32_t ip,
                                                      const Element<dim>& e );
     
    virtual void InterpolateNodeProperties( const Element<dim>& e );

    /// to set properties on the fly for compressible flow modelling
    void SaturationWettingPhase( double s_wetting );
    void ViscosityNonWettingPhase( double visc );
    void ViscosityWettingPhase( double visc );
    void DensityNonWettingPhase( double dens );
    void DensityWettingPhase( double dens );

    /// Return of keys for calculation of two phase mobilities in transport
    csmp::Index WettingPhaseSaturationKey() const;
    csmp::Index NonWettingPhaseSaturationKey() const;
    csmp::Index TotalMobilityKey() const;

    /// accessors
    double     Permeability()  const;
    void       Permeability( double permeability );
    TensorVariable<dim>   TensorPermeability()  const;
    void       TensorPermeability( TensorVariable<dim> permeability );
    double     Swr() const;
    double     Snr() const;
    void       Swr( double swr );
    void       Snr( double snr );    

    /// fluid properties
    double   ViscosityNonWettingPhase() const;
    double   ViscosityWettingPhase() const;
    double   DensityNonWettingPhase() const;
    double   DensityWettingPhase() const;
    double   Saturation( uint32_t phase ) const;
    double   MobilityPhase( uint32_t phase ) const; // kr_i/mu_i
    double   ViscosityRatio() const;

    /// sum of all phase mobilities * k
    double TotalMobility() const;
    
    /// sum of all phase mobilities
    double TotalMobilityMultiplier() const;

    // multiplier for diffusion coefficient in the case of non-linear diffusion
    // uses viscosity(phase) 
    double DiffusionMultiplier( size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double CapillaryDiffusionMultiplier(  ) const;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    double GravityTerm() const;
    double GravityMultiplier_G( ) const;
    double GravityMultiplier_dGds( ) const;

    /// G multiplier
    double G() const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double AdvectionMultiplier() const;

    /// linearized diffusion multiplier for large-timestep calculations
    virtual double DiffusionCharacteristic( size_t phase ) const;

    /// linearized fractional flow derivative
    virtual double ShockSpeed() const;
    virtual double ShockHeight() const;
    void ShockSpeedHeight( double& speed, double& height) const;

    /// maximum value of previous derivative
    virtual double MaxFractionalFlowDerivative( ) const;

    /// maximum value of dpcdS
    double MaxCapillaryPressure( uint32_t phase = 1U ) const;

    /// always of the wetting phase by convention
    virtual double  EffectiveSaturation() const;
    virtual double  SeffToSw() const;
    virtual double  SeffToSw( double seff) const;

    /// fractional flow
    double f_Phase( uint32_t phase ) const;
  
    /// fractional flow of water at water saturation sw
    double fw_at( double sw ) const;
    double fn_at( double sw ) const;

    /// relative permeabilities
    virtual double krn_Phase() const = 0;
    virtual double krw_Phase() const = 0;

    /// derivatives of relative permeabilities
    virtual double dkrnds_Phase() const;
    virtual double dkrwds_Phase() const;

    virtual double dlnds() const;
    virtual double dlwds() const;

    /// capillary pressure (limit this to 4e7, the max strength of the rock)
    /// do this by computing seff for which pc=4e7, then use this seff as
    /// a limiting value
    virtual double pc_Phase() const = 0;

    /// capillary pressure derivatives (treat seff as for previous function)
    virtual double dpcds_Phase() const = 0;

    /// inverse capillary pressure function
    virtual double Sw_Phase( double pc_Phase ) const;

    /// derivatives of inverse capillary pressure function
    virtual double dsdpc_Phase( double pc_Phase ) const;
                                  
    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    virtual double dfds() const;

    /// derivatives of gravitational flow (advection multipliers)
    virtual double dGds() const;
    
    /// Numerical derivatives
    // mobility ratio?
    double  mobility_w_at( double ) const;
    double  mobility_n_at( double ) const;
    double  krw_at( double ) const;
    double  krn_at( double ) const;
    double  dkrwds_at( double ) const;
    double  dkrnds_at( double ) const;
    double  dfds_at( double se ) const;
    double  G_at( double se ) const;
    double  pc_at( double se ) const;
    double  dpcds_at( double se) const;
    double  Sw_at( double se) const;
    double  dsdpc_at( double se) const;

    virtual double dlwds_numerical(  double h = 0.001 ) const;
    virtual double dlnds_numerical(  double h = 0.001 ) const;
    virtual double dkrwds_numerical( double h = 0.001 ) const;
    virtual double dkrnds_numerical( double h = 0.001 ) const;
    virtual double dfds_numerical (  double h = 0.001 ) const;
    virtual double dGds_numerical (  double h = 0.000001 ) const;
    virtual double dpcds_numerical(  double h = 0.00001 ) const;

    /// default is the wetting phase
    virtual void Out( uint32_t phase=1U ) const;
    
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
    double            swr_, snr_,   // irreducible saturations
                      k_,           // permeability
                      ift_,         // interfacial tension
                      mun_, muw_,   // viscosities of wetting and non-wetting phase
                      rhn_, rhw_;   // densities
    TensorVariable<dim> K_;         // tensor permeability

    mutable double    sat_,         // saturation of the wetting phase
                      seff_;        // actual saturation, effective saturation
    const double      acc_gravity_, // accelaration of gravity
                      tolerance_,   // cut_off value for expensive calculations
                      MAX_CAPILLARY_PRESSURE_,      // maximum permitted capillary pressure
                      MAX_CAPILLARY_PRESSURE_SLOPE_;// maximum permitted capillary pressure slope
    const bool        interpolate_fluid_properties_;
    const bool        sw_ro_mu_placement_;  // NODE=true ELEMENT=false
    bool              tensor_permeability_; // TENSOR=true SCALAR=false
};

} // end namespace csmp

#endif 





















