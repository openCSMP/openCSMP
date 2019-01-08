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
 
*/
template<size_t dim, template<size_t> class USER>
class H2O_CO2_NaCl_FlowFunctions {
  public:
    /// water saturation interpolated to element barycentre; use Read( key_sH2O ) to get nodal value
    double64 Sw( const Element<dim>* const e ) const;
      
    /// density * lambda = kri(sw)/mu_i  of the phase i: 0 for water, 1 for the non-wetting phase
    double64 Mobility( const Element<dim>* const, size_t phase ) const;
    
    /// density * lambda_i: i=0 for water, 1 for the non-wetting phase (using prescribed sw)
    double64 Mobility_at( const Element<dim>* const, size_t phase, double64 sw ) const;

    /// d lambda_i / dsw
    double64 MobilityDerivative( const Element<dim>* const, size_t phase ) const;
 
    /// d lambda_i / dsw
    double64 MobilityDerivative_at( const Element<dim>* const, size_t phase, double64 sw ) const;
                        
    /// lambda_t: sum of phase-mobility * density products
    double64 TotalMobility(  const Element<dim>* const ) const;
    
    /// lambda_t: sum of phase-mobility * density products
    double64 TotalMobility_at(  const Element<dim>* const, double64 sw ) const;
    
    /// lambda overbar: mobility product l_overbar = (li * rhow * lj * rhonw) / (li*rhow + lj*rhonw),  also known as G
    double64 MobilityProduct(  const Element<dim>* const ) const;

    /// d lambda overbar / dsw also known as dGds
    double64 MobilityProductDerivative( const Element<dim>* const, bool  evaluate_numerically=false ) const;
                        
    double64 MobilityProductDerivative_at( const Element<dim>* const, double64 sw ) const;

     /// fractional mass flow; 0=water, 1=non-wetting phase
    double64 f( const Element<dim>* const, size_t phase ) const;
    
     /// fractional mass flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    double64 f_at( const Element<dim>* const, size_t phase, double64 sw ) const;
    
    /// permeability
    double64 Permeability(  const Element<dim>* const ) const;

    /// derivative of fractional flow w.r.t water saturation (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds(  const Element<dim>* const, size_t phase ) const;
    
    double64 dfds_at( const Element<dim>* const, double64 sw ) const;
    
    /// maximum value of fractional flow / saturation derivative
    double64 MaxFractionalFlowDerivative(  const Element<dim>* const ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
    double64 AdvectionMultiplier(  const Element<dim>* const ) const;

    /// outputs the shock wave celerity for mass flow
    double64 ShockSpeed( const Element<dim>* const ) const;
    
    /// outputs the characteristic water saturation at the shock (tied to rock properties in viscous dominated flow)
    double64 ShockHeight( const Element<dim>* const ) const;
    
    /// calculated shock height and mass flow rate (all-in-one function)
    void ShockSpeedAndHeight( const Element<dim>* const, double64& speed, double64& height ) const;
    
    /// speed (m/s) of the shock front in current cell
    double64 ShockFrontVelocity( const Element<dim>* const ) const;

    /// k * kri(sw)/mi * rho_i^2 projected onto the dip vector of the current element; writes result to dip vector
    void GravityTerm( const Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
  
    /// gravity multiplier, gmult = k * kri(sw)/mi * rho_i^2 * G (=mobility product); writes result on dip vector
    void GravityMultiplier_G( const Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
  
    /// gravity multiplier, gmult = k * kri(sw)/mi * rho_i^2 * saturation derivative of the mobility product dG/ds, writes result on dip vectpor
    void GravityMultiplier_dGds( const Element<dim>* const, VectorVariable<dim>& dip_vec ) const;

    /// average mass diffusion coefficient for CO2 in the aqueous phase
    double64 DiffusionMultiplier(  const Element<dim>* const, size_t phase ) const;
    
    /// mass diffusion coefficient for the non-linear diffusion of saturation due to the saturation dependent dpc/ds
    double64 CapillaryDiffusionMultiplier( const Element<dim>* const ) const;

    /// to model cappillary diffusion of the non-wetting, simply use a negative sign on the multiplier computed with previous function
    //double64 CapillaryDiffusionMultiplier( const Element<dim>* const, size_t phase ) const;

    // central difference derivatives of saturation and flow functions
  
    double64 dfds_Numerical(  const Element<dim>* const, size_t phase, double64 delta_s = 0.001 ) const;
      
    double64 dfds_at_Numerical(  const Element<dim>* const, double64 sw, double64 delta_s = 0.001 ) const;
      
    double64 dGds_Numerical(  const Element<dim>* const, double64 delta_s = 0.000001 ) const;
      
    double64 dlwds_Numerical(  const Element<dim>* const, double64 delta_s = 0.001 ) const;
      
    double64 dlnds_Numerical(  const Element<dim>* const, double64 delta_s = 0.001 ) const;

           
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
      
    // NON-STANDARD INTERFACES

    /// calculating the  Inflection point
    double64 InflectionPointSaturation( const Element<dim>* const ) const ;
    
    /// calculating the Tangent Point
    double64 TangentPointSaturation( const Element<dim>* const ) const ;
    
    /// calculating  Buckley Leverett function... this has to be zero at shock point
    double64 TangentOfFractionalFlowFunction( const Element<dim>* const, double64 S) const ;
    
    /// finding root of  Buckley Leverett function... which results is shock point saturation.
    double64 FindRootSecantMethod( const Element<dim>* const, double64 S1, double64 S2) const ;

};

} // end csmp

#endif /* CSMP_H2O_CO2_NACl_FLOW_FUNCTIONS_H */
