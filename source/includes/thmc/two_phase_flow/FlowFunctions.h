//
//  FlowFunctions.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 7/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLOW_FUNCTIONS_H
#define CSMP_FLOW_FUNCTIONS_H


#include "BrooksCoreySaturationFunctions.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Region;
template<size_t> class Model;
  
///
/**
    @brief 2-phase flow functions

    Captures all relevant input values from the model which are then used to compute the
    relevant constitutive relations at the point of interest in the element.
    
    The variables that are specific to each relperm model and capillary pressure curve are stored
    there. All other variables are stored here.
    
    @section example Example
    ...
    Initialize( model_domain, e );
    InitializeForBaryCenter( e );
    double lt = TotalMobility();
    
    @attention the phases are numbered, phase 1 is the water phase (usually the wetting phase)
    
    @attention for simplicity and efficiency, fluid property values are always computed at the nodes by respective EOS modules.
    Subsequently these values are interpolated to the points of interest.
    
    Code serves as an example of how static polymorphism can be used to implement constitutive relationships for multiphase flow.
*/
template<size_t dim>
class FlowFunctions : public VariableSet2PhaseSlightlyCompressible,              ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
                public BrooksCoreySaturationFunctions<dim,FlowFunctions>,  ///< placeholder for saturation function model
                public Fluid<dim,FlowFunctions> {                          ///< placeholder for fluids module / EOS interface
  public:
    FlowFunctions();
                
    // INITIALISATION (relperm models etc.) - maybe put this into a separate policy
                
    /// reads region and element properties snr, swr, k, mtrl_param, fluid_param and interpolates to target placement
    void Initialize( const Region<dim>&, const Element<dim>& e );
    
    /// interpolates node properties to element barycenter
    void InitializeForBaryCenter( const Element<dim>& e );

    /// finite element integration points
    void InitializeForIntegrationPoint( size_t ip, const Element<dim>& e );
    /// finite volume facets
    void InitializeForFacetIntegrationPoint( size_t facet, size_t ip, const Element<dim>& e );
    /// finite volume sectors
    void InitializeForSectorIntegrationPoint( size_t sector, size_t ip, const Element<dim>& e );
                
    // OUTPUTS
    
    double64 Mobility( size_t phase ) const;
    
    /// sum of all phase mobilities * k
    double64 TotalMobility() const;
    
    /// sum of all phase mobilities
    double64 TotalMobilityMultiplier() const;
    
    /// l1 * l2 / l1 + l2
    double64 MobilityProduct() const;

    // multiplier for diffusion coefficient in the case of non-linear diffusion
    // uses viscosity(phase) 
    double64 DiffusionMultiplier( size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double64 CapillaryDiffusionMultiplier() const;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    double64 GravityTerm() const;
    double64 GravityMultiplier_G() const;
    double64 GravityMultiplier_dGds() const;

    /// G multiplier
    double64 G() const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double64 AdvectionMultiplier() const;

    /// maximum value of previous derivative
    double64 MaxFractionalFlowDerivative() const;

    /// maximum value of dpcdS
    double64 MaxCapillaryPressure( size_t phase = 1U ) const;

    /// always of the wetting phase by convention
    double64  EffectiveSaturation() const;

    /// fractional flow
    double64 f_Phase( size_t phase ) const;
                  
    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds() const;

    /// derivatives of gravitational flow (advection multipliers)
    double64 dGds() const;
    
    // BASIC ACCESSORS - element properties are read from model, others are stored.
    double64 Permeability() const;
    double64 Porosity() const;
    double64 Density( size_t phase ) const;
    double64 Viscosity( size_t phase ) const;
    double64 IrreducibleWaterSaturation() const;
    double64 ResidualSaturation() const;
    /// defaults to water phase
    double64 Saturation( size_t phase=1 ) const;
    
  private:
    csmp::Index pf0_key_, pf1_key_, T_key_,  ///< pressures at previous and current time levels
                s1_key_, s2_key_,            ///< convention that water is the first phase
                rho1_key_, rho2_key,         ///< fluid densities (kg/m3)
                mu1_key_, mu2_key;           ///< dynamic fluid viscosities (Pa.s)
  
    double64  P_, T_, sw_;
    double64  rhow_, rhon_, muw_, mun_;
    
};

} // end csmp

#endif /* CSMP_FLOW_FUNCTIONS_H */
