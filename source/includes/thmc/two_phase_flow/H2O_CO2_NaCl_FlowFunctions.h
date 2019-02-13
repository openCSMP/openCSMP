//
//  H2O_CO2_NaCl_FlowFunctions.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai and Qi Shao on 07/January/2019.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_H2O_CO2_NACl_FLOW_FUNCTIONS_H
#define CSMP_H2O_CO2_NACl_FLOW_FUNCTIONS_H

#include "CSMP_definitions.h"

namespace csmp {

/**
    @brief H2O_CO2_NaCl_FlowFunctions: Saturation functions and fluid properties combined into flow functions for the H2O + CO2+ NaCl system
 
    Use in conjunction with the PhaseStateFinder_H2O_CO2_NaCl for compositional simulation
    in a mass-based transport framework. Plugs into a suitable FlowFunctionsModule
    within the Colleoli class framework for a multiphase-multicomponent transport scheme.

    Generate all relevant input values for transport scheme from SaturationFunctions and Fluid modules,
    using the constitutive relationships specified therein.
    Gets the fluid properties by interpolation from the Element nodes to its barycentre.
    
    Expects that the water saturation, the fluid densities and viscosities to be placed on the element nodes,
    while all material (k, phi, bcp etc.), and rocktype related parameters are placed on the element.
 
    @note all parameters are weighted by the element thickness which is 1 for iso-dimensional elements.
 
    @author Stephan Matthai
    @date 7/1/2019
 
    @discussion Phase state cannot be interpolated between nodes !
      If we read the system state from the element,
      in the case of an element-based phase state calculation, it may be
      entirely different from that computed at the nodes and therefore,
      we may try to fetch densities and viscosities which do not exist ?!

*/
// TODO: use these to generate values at nodes as well (no interpolations required!)
// TODO: does the capillary diffusion multiplier need multiplication with density
// TODO: create tensor permeability versions of these functions
template<size_t dim, template<size_t> class USER>
class H2O_CO2_NaCl_FlowFunctions {
  public:
  
    // ELEMENT BASED CALCULATIONS with '_at' versions for computations involving facet and sector integration point saturations
  
    /// water saturation interpolated to element barycentre; use Read( key_sH2O ) to get nodal value
    double64 Sw( Element<dim>* const e ) const;
      
    /// non-wetting phase saturation interpolated to element barycentre; use Read( key_sCO2 ) to get nodal value
    double64 Snw( Element<dim>* const e ) const;

    /// density * lambda = kri(sw)/mu_i  of the phase i: 0 for water, 1 for the non-wetting phase
    double64 Mobility( Element<dim>* const, size_t phase ) const;
    
    /// density * lambda_i: i=0 for water, 1 for the non-wetting phase (using prescribed sw)
    double64 Mobility_at( Element<dim>* const, size_t phase, double64 sw ) const;

    /// d lambda_i / dsw
    double64 MobilityDerivative( Element<dim>* const, size_t phase ) const;
 
    /// d lambda_i / dsw
    double64 MobilityDerivative_at( Element<dim>* const, size_t phase, double64 sw ) const;
                        
    /// lambda_t: sum of phase-mobility * density products
    double64 TotalMobility( Element<dim>* const ) const;
    
    /// lambda_t: sum of phase-mobility * density products
    double64 TotalMobility_at( Element<dim>* const, double64 sw ) const;
    
    /// lambda overbar: mobility product l_overbar = (li * rhow * lj * rhonw) / (li*rhow + lj*rhonw),  also known as G
    double64 MobilityProduct( Element<dim>* const ) const;

    /// lambda overbar: mobility product l_overbar = (li * rhow * lj * rhonw) / (li*rhow + lj*rhonw),  also known as G
    double64 MobilityProduct_at( Element<dim>* const, double64 sw ) const;

    /// d lambda overbar / dsw also known as dGds
    double64 MobilityProductDerivative( Element<dim>* const, bool  evaluate_numerically=false ) const;
                        
    double64 MobilityProductDerivative_at( Element<dim>* const, double64 sw ) const;

     /// fractional mass flow; 0=water, 1=non-wetting phase
    double64 f(Element<dim>* const, size_t phase ) const;
    
     /// fractional mass flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    double64 f_at(Element<dim>* const, size_t phase, double64 sw ) const;
    
    /// permeability
    double64 Permeability( Element<dim>* const ) const;

    /// derivative of fractional flow w.r.t water saturation (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds( Element<dim>* const, size_t phase ) const;
    
    double64 dfds_at( Element<dim>* const, size_t phase, double64 sw ) const;
    
    /// maximum value of fractional flow / saturation derivative
    double64 MaxFractionalFlowDerivative( Element<dim>* const, size_t phase ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
    double64 AdvectionMultiplier( Element<dim>* const ) const;

    /// outputs the shock wave celerity for mass flow
    double64 ShockSpeed( Element<dim>* const ) const;
    
    /// outputs the characteristic water saturation at the shock (tied to rock properties in viscous dominated flow)
    double64 ShockHeight( Element<dim>* const ) const;
    
    /// calculated shock height and mass flow rate (all-in-one function)
    void ShockSpeedAndHeight( Element<dim>* const, double64& speed, double64& height ) const;
    
    /// speed (m/s) of the shock front in current cell
    double64 ShockFrontVelocity( Element<dim>* const ) const;

    /// k * kri(sw)/mi * rho_i^2 projected onto the dip vector of the current element; writes result to dip vector
    void GravityTerm(Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
    
    void GravityMultiplier_phase( Element<dim>* const, VectorVariable<dim>& dip_vec, size_t phase ) const;
  
    /// gravity multiplier, gmult = k * kri(sw)/mi * rho_i^2 * G (=mobility product); writes result on dip vector
    void GravityMultiplier_G( Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
  
    /// gravity multiplier, gmult = k * kri(sw)/mi * rho_i^2 * saturation derivative of the mobility product dG/ds, writes result on dip vectpor
    void GravityMultiplier_dGds( Element<dim>* const, VectorVariable<dim>& dip_vec ) const;

    ///  mass diffusion coefficient for CO2 in water saturated porous medium; use phase=1
    double64 DiffusionMultiplier( Element<dim>* const, size_t phase ) const;
    
    /// mass diffusion coefficient for the non-linear diffusion of saturation due to the saturation dependent dpc/ds
    double64 CapillaryDiffusionMultiplier( Element<dim>* const ) const;
  
    /// k * dpds(sw) * mobility_product(sw), @attention value can be zero.
    double64 CapillaryDiffusionMultiplier_at( Element<dim>* const, double64 sw ) const;

    // TODO: used when? double64 CapillaryDiffusionMultiplier_Phase( Element<dim>* const, size_t phase ) const;

  
    // NODE-BASED COMPUTATIONS using element parameters, but saturations and fluid properties from the current node

    /// density * lambda = kri(sw)/mu_i  of the phase i: 0 for water, 1 for the non-wetting phase
    double64 Mobility( Element<dim>* const, size_t phase, size_t node ) const;

    /// d lambda_i / dsw
    double64 MobilityDerivative( Element<dim>* const, size_t phase, size_t node ) const;
 
    /// lambda_t: sum of phase-mobility * density products
    double64 TotalMobility( Element<dim>* const, size_t node ) const;
  
    /// lambda overbar: mobility product l_overbar = (li * rhow * lj * rhonw) / (li*rhow + lj*rhonw),  also known as G
    double64 MobilityProduct( Element<dim>* const, size_t node ) const;

    /// d lambda overbar / dsw also known as dGds
    double64 MobilityProductDerivative( Element<dim>* const, size_t node, bool  evaluate_numerically=false ) const;

     /// fractional mass flow; 0=water, 1=non-wetting phase
    double64 f( Element<dim>* const, size_t phase, size_t node ) const;

    /// derivative of fractional flow w.r.t water saturation (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds( Element<dim>* const, size_t phase, size_t node ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
    double64 AdvectionMultiplier( Element<dim>* const, size_t node ) const;

    /// k * kri(sw)/mi * rho_i^2 projected onto the dip vector of the current element; writes result to dip vector
    void     GravityTerm( Element<dim>* const, size_t node, VectorVariable<dim>& dip_vec ) const;

    ///  mass diffusion coefficient for CO2 in water saturated porous medium; use phase=1 and get from Element
  
    /// mass-based coefficient for non-linear diffusion of saturation due to saturation-dependent dpc/ds, @attention value can be zero
    double64 CapillaryDiffusionMultiplier( Element<dim>* const, size_t node ) const;


    // central difference derivatives of saturation and flow functions
  
    double64 dfds_Numerical( Element<dim>* const, size_t phase, double64 delta_s = 0.001 ) const;
      
    double64 dfds_at_Numerical( Element<dim>* const, double64 sw, double64 delta_s = 0.001 ) const;
      
    double64 dGds_Numerical( Element<dim>* const, double64 delta_s = 0.000001 ) const;
      
    double64 dlwds_Numerical( Element<dim>* const, double64 delta_s = 0.001 ) const;
      
    double64 dlnds_Numerical( Element<dim>* const, double64 delta_s = 0.001 ) const;

           
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
      
    // NON-STANDARD INTERFACES
    /// quick look at which phases exist on all nodes, returns true if aqueous and carbonic phases are continuous and there is no salt
    bool ContinuousPhases( Element<dim>* const, bool& aqueous, bool& carbonic, bool& halite ) const;
  
    /// just the saturation of the aqueous phase
    double64 InterpolateAqueousPhaseSaturation( const Element<dim>* const ) const;
  
    /// interpolates saturations of all mobile phases to barycenter and returns their sum which should be 1 - halite_saturation
    double64 InterpolateSaturations( Element<dim>* const, double64& sw, double64& snw ) const;
  
    /// returns aqueous phase properties at element barycentre, if there is a value at at least a single node, else returns false
    bool InterpolateAqueousPhase( Element<dim>* const, double64& sw, double64& rhow, double64& muw ) const;

    /// returns carbonic phase properties at element barycentre, if there is a value at at least a single node, else returns false
    bool InterpolateCarbonicPhase( Element<dim>* const, double64& snw, double64& rhon, double64& mun ) const;
  
    /// returns halite saturation = volume fraction of salt
    double64 InterpolateSystem( Element<dim>* const,
                                double64& sw, double64& rhow, double64& muw,
                                double64& snw, double64& rhon, double64& mun ) const;

    /// calculating the  Inflection point
    double64 InflectionPointSaturation( Element<dim>* const ) const;
    
    /// calculating the Tangent Point
    double64 TangentPointSaturation( Element<dim>* const ) const;
    
    /// calculating  Buckley Leverett function... this has to be zero at shock point
    double64 TangentOfFractionalFlowFunction( Element<dim>* const, double64 S) const;
    
    /// finding root of  Buckley Leverett function... which results is shock point saturation.
    double64 FindRootSecantMethod( Element<dim>* const, double64 S1, double64 S2) const;
  
  private:
    mutable double64 sw_,   ///< saturation aqueous phase
                     sn_,   ///< saturation carbonic phase
                     rhow_, ///< density aqueous phase (kg/m3)
                     rhon_, ///< density carbonic phase
                     muw_,  ///< viscosity aqueous phase (Pa.s)
                     mun_;  ///< viscosity carbonic phase
};

} // end csmp

#endif /* CSMP_H2O_CO2_NACl_FLOW_FUNCTIONS_H */
