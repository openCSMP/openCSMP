//
//  TwoPhaseFlowFunctions.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_TWO_PHASE_FLOW_FUNCTIONS_H
#define CSMP_TWO_PHASE_FLOW_FUNCTIONS_H

#include "CSMP_definitions.h"

namespace csmp {

/**
    @brief 2-phase multiphase flow functions

    Generates all relevant input values for transport scheme from relperm model and fluids module
    using the constitutive relationships specified therein.    
    
    The results are written to the model at the respective variable placements.
*/
template<size_t dim, template<size_t> class USER>
class TwoPhaseFlowFunctions {
  public:
    /// current water saturation initialised inside of the model
    template<class TARGET_PLACEMENT>
    double64 Sw( const TARGET_PLACEMENT& p ) const { return p.Obtain(this->key_sH2O); }
      
    /// lambda parameter: 0 for water, 1 for the non-wetting phase
    template<class TARGET_PLACEMENT>
    double64 Mobility( const TARGET_PLACEMENT&, size_t phase ) const;
    
    /// lambda parameter: 0 for water, 1 for the non-wetting phase (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 Mobility_at( const TARGET_PLACEMENT&, size_t phase, double64 ) const;

    /// d lambda_i / dsw
    template<class TARGET_PLACEMENT>
    double64 MobilityDerivative( const TARGET_PLACEMENT&, size_t phase, bool evaluate_numerically=false ) const;
 
    template<class TARGET_PLACEMENT>
    double64 MobilityDerivative_at( const TARGET_PLACEMENT&, size_t phase, double64 ) const;
                        
    /// lambda_t: sum of phase mobilities
    template<class TARGET_PLACEMENT>
    double64 TotalMobility(  const TARGET_PLACEMENT& ) const;
    
    /// lambda_t: sum of phase mobilities (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 TotalMobility_at(  const TARGET_PLACEMENT&, double64 ) const;
    
    /// lambda overbar: l1 * l2 / l1 + l2 = mobility product / total mobility also known as G
    template<class TARGET_PLACEMENT>
    double64 MobilityProduct(  const TARGET_PLACEMENT& ) const;

    /// d lambda overbar / dsw also known as dGds
    template<class TARGET_PLACEMENT>
    double64 MobilityProductDerivative( const TARGET_PLACEMENT&, bool  evaluate_numerically=false ) const;
                        
    template<class TARGET_PLACEMENT>
    double64 MobilityProductDerivative_at( const TARGET_PLACEMENT&, double64 ) const;

     /// fractional flow; 0=water, 1=non-wetting phase
    template<class TARGET_PLACEMENT>
    double64 f( const TARGET_PLACEMENT&, size_t phase ) const;
    
     /// fractional flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 f_at( const TARGET_PLACEMENT&, size_t phase, double64 ) const;
    
    /// Permeability
    template<class TARGET_PLACEMENT>
    double64 Permeability(  const TARGET_PLACEMENT& ) const;    

    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    template<class TARGET_PLACEMENT>
    double64 dfds(  const TARGET_PLACEMENT&, size_t phase, bool evaluate_numerically=false ) const;
    
    template<class TARGET_PLACEMENT>
    double64 dfds_at( const TARGET_PLACEMENT&, double64 sw ) const;
    
    /// maximum value of previous derivative
    template<class TARGET_PLACEMENT>
    double64 MaxFractionalFlowDerivative(  const TARGET_PLACEMENT& ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
    template<class TARGET_PLACEMENT>
   	double64 AdvectionMultiplier(  const TARGET_PLACEMENT&, bool evaluate_numerically=false ) const;

    /// outputs the saturation of the desired phase at the shock front
    template<class TARGET_PLACEMENT>
    double64 ShockSaturation( const TARGET_PLACEMENT&, size_t phase, bool evaluate_numerically=false ) const;

    /// output Shock Speed
//    template<class TARGET_PLACEMENT>
//    double64 ShockSpeed(  const TARGET_PLACEMENT& ) const;
    
    /// outputs the water saturation at the shock front
//    template<class TARGET_PLACEMENT>
//    double64 ShockHeight( const TARGET_PLACEMENT& ) const;
    
    /// calculated shock height and speed
//    template<class TARGET_PLACEMENT>
//    void     ShockSpeedHeight( const TARGET_PLACEMENT&, double64& speed, double64& height ) const;

    /// calculating the Shock velocity
//    template<class TARGET_PLACEMENT>
//double64 ShockFrontVelocity( const TARGET_PLACEMENT& ) const ;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    template<class TARGET_PLACEMENT>
    double64 GravityTerm(  const TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 GravityMultiplier_G(  const TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 GravityMultiplier_dGds(  const TARGET_PLACEMENT&, bool evaluate_numerically=false ) const;

    /// multiplier for diffusion coefficient in the case of non-linear diffusion uses viscosity(phase)
    template<class TARGET_PLACEMENT>
    double64 DiffusionMultiplier(  const TARGET_PLACEMENT&, size_t phase ) const;
    
    /// driven by capillary pressure gradient
    template<class TARGET_PLACEMENT>
    double64 CapillaryDiffusionMultiplier( const TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 CapillaryDiffusionMultiplier_Phase(  const TARGET_PLACEMENT&, size_t phase ) const;
    
  
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
      
    template<class TARGET_PLACEMENT>
    double64 dfds_Numerical(  const TARGET_PLACEMENT&, size_t phase, double64 h = 0.001 ) const;
      
    template<class TARGET_PLACEMENT>
    double64 dfds_at_Numerical(  const TARGET_PLACEMENT&, double64 sw, double64 h = 0.001 ) const;
      
    template<class TARGET_PLACEMENT>
    double64 dGds_Numerical(  const TARGET_PLACEMENT&, double64 h = 0.000001 ) const;
      
    template<class TARGET_PLACEMENT>
    double64 dlwds_Numerical(  const TARGET_PLACEMENT&, double64 h = 0.001 ) const;
      
    template<class TARGET_PLACEMENT>
    double64 dlnds_Numerical(  const TARGET_PLACEMENT&, double64 h = 0.001 ) const;
  
    template<class TARGET_PLACEMENT>
    double64 InflectionPointSaturation( const TARGET_PLACEMENT& p ) const;
};

} // end csmp

#endif /* CSMP_FLOW_FUNCTIONS_BC_HYSTERETIC_H */
