//
//  TwoPhaseFlowFunctions.h
//  CSMP_GitHub
//
//  Created by Mahyar Madadi  on 16/July/2018.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_TWO_PHASE_FLOW_FUNCTIONS_H
#define CSMP_TWO_PHASE_FLOW_FUNCTIONS_H

#include "CSMP_definitions.h"

namespace csmp {

/**
    @brief 2-phase multiphase flow functions for immiscible displacement of slightly compressible fluids.

    Generates all relevant input values for 2-phase immiscible (saturation-based volume) transport scheme.
    Uses the SaturationFunctions and Fluid module and the
    constitutive relationships specified therein.
*/
template<uint32_t dim, template<uint32_t> class USER>
class TwoPhaseFlowFunctions {
  public:
    /// current water saturation initialised inside of the model
    double Sw( Element<dim>* const e ) const;
      
    /// lambda parameter: 0 for water, 1 for the non-wetting phase
    double Mobility( Element<dim>* const, size_t phase ) const;
    
    /// lambda parameter: 0 for water, 1 for the non-wetting phase (using prescribed sw)
    double Mobility_at( Element<dim>* const, size_t phase, double sw ) const;

    /// d lambda_i / dsw
    double MobilityDerivative( Element<dim>* const, size_t phase, bool evaluate_numerically=false ) const;
 
    double MobilityDerivative_at( Element<dim>* const, size_t phase, double sw ) const;
                        
    /// lambda_t: sum of phase mobilities
    double TotalMobility(  Element<dim>* const ) const;
    
    /// lambda_t: sum of phase mobilities (using prescribed sw)
    double TotalMobility_at(  Element<dim>* const, double sw ) const;
    
    /// lambda overbar: l1 * l2 / l1 + l2 = mobility product / total mobility also known as G
    double MobilityProduct(  Element<dim>* const ) const;

    /// d lambda overbar / dsw also known as dGds
    double MobilityProductDerivative( Element<dim>* const, bool  evaluate_numerically=false ) const;
                        
    double MobilityProductDerivative_at( Element<dim>* const, double sw ) const;

     /// fractional flow; 0=water, 1=non-wetting phase
    double f( Element<dim>* const, size_t phase ) const;
    
     /// fractional flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    double f_at( Element<dim>* const, size_t phase, double sw ) const;
    
    /// Permeability
    double Permeability(  Element<dim>* const ) const;

    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double dfds(  Element<dim>* const, size_t phase ) const;
    
    double dfds_at( Element<dim>* const, double sw ) const;
    
    /// maximum value of previous derivative
    double MaxFractionalFlowDerivative(  Element<dim>* const ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double AdvectionMultiplier(  Element<dim>* const ) const;

    /// outputs the shock wave celerity
    double ShockSpeed(  Element<dim>* const ) const;
    
    /// outputs the water saturation at the non-wetting phase shock front
    double ShockHeight( Element<dim>* const ) const;
    
    /// calculated shock height and speed
    void     ShockSpeedAndHeight( Element<dim>* const, double& speed, double& height ) const;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    void GravityMultiplier( Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
    
    void GravityMultiplier_G( Element<dim>* const, VectorVariable<dim>& dip_vec ) const;
    
    void GravityMultiplier_dGds( Element<dim>* const, VectorVariable<dim>& dip_vec ) const;

    /// multiplier for diffusion coefficient in the case of non-linear diffusion uses viscosity(phase)
    double DiffusionMultiplier(  Element<dim>* const, size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double CapillaryDiffusionMultiplier( Element<dim>* const ) const;
    
    double CapillaryDiffusionMultiplier_Phase(  Element<dim>* const, size_t phase ) const;

    /// calculating the Shock velocity
    double ShockFrontVelocity( Element<dim>* const ) const ;
    
    double dfds_Numerical(  Element<dim>* const, size_t phase, double delta_s = 0.001 ) const;
      
    double dfds_at_Numerical(  Element<dim>* const, double sw, double delta_s = 0.001 ) const;
      
    double dGds_Numerical(  Element<dim>* const, double delta_s = 0.000001 ) const;
      
    double dlwds_Numerical(  Element<dim>* const, double delta_s = 0.001 ) const;
      
    double dlnds_Numerical(  Element<dim>* const, double delta_s = 0.001 ) const;

           
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
      
    // NON-STANDARD INTERFACES

    /// calculating the  Inflection point
    double InflectionPointSaturation( Element<dim>* const ) const ;
    
    /// calculating the Tangent Point
    double TangentPointSaturation( Element<dim>* const ) const ;
    
    /// calculating  Buckley Leverett function... this has to be zero at shock point
    double TangentOfFractionalFlowFunction( Element<dim>* const, double S) const ;
    
    /// finding root of  Buckley Leverett function... which results is shock point saturation.
    double FindRootSecantMethod( Element<dim>* const, double S1, double S2) const ;
};

} // end csmp

#endif /* CSMP_TWO_PHASE_FLOW_FUNCTIONS_H */
