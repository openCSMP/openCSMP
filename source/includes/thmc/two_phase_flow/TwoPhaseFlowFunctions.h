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
    @brief 2-phase multiphase flow functions

    Generates all relevant input values for transport scheme from relperm model and fluids module
    using the constitutive relationships specified therein.    
    
    The results are written to the model at the respective variable placements.
*/
template<size_t dim, template<size_t> class USER>
class TwoPhaseFlowFunctions {
  public:
    /// current water saturation initialised inside of the model
    double64 Sw( const Element<dim>* const e ) const;
      
    /// lambda parameter: 0 for water, 1 for the non-wetting phase
    double64 Mobility( const Element<dim>* const, size_t phase ) const;
    
    /// lambda parameter: 0 for water, 1 for the non-wetting phase (using prescribed sw)
    double64 Mobility_at( const Element<dim>* const, size_t phase, double64 sw ) const;

    /// d lambda_i / dsw
    double64 MobilityDerivative( const Element<dim>* const, size_t phase, bool evaluate_numerically=false ) const;
 
    double64 MobilityDerivative_at( const Element<dim>* const, size_t phase, double64 sw ) const;
                        
    /// lambda_t: sum of phase mobilities
    double64 TotalMobility(  const Element<dim>* const ) const;
    
    /// lambda_t: sum of phase mobilities (using prescribed sw)
    double64 TotalMobility_at(  const Element<dim>* const, double64 sw ) const;
    
    /// lambda overbar: l1 * l2 / l1 + l2 = mobility product / total mobility also known as G
    double64 MobilityProduct(  const Element<dim>* const ) const;

    /// d lambda overbar / dsw also known as dGds
    double64 MobilityProductDerivative( const Element<dim>* const, bool  evaluate_numerically=false ) const;
                        
    double64 MobilityProductDerivative_at( const Element<dim>* const, double64 sw ) const;

     /// fractional flow; 0=water, 1=non-wetting phase
    double64 f( const Element<dim>* const, size_t phase ) const;
    
     /// fractional flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    double64 f_at( const Element<dim>* const, size_t phase, double64 sw ) const;
    
    /// Permeability
    double64 Permeability(  const Element<dim>* const ) const;

    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds(  const Element<dim>* const, size_t phase ) const;
    
    double64 dfds_at( const Element<dim>* const, double64 sw ) const;
    
    /// maximum value of previous derivative
    double64 MaxFractionalFlowDerivative(  const Element<dim>* const ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double64 AdvectionMultiplier(  const Element<dim>* const ) const;

    /// outputs the shock wave celerity
    double64 ShockSpeed(  const Element<dim>* const ) const;
    
    /// outputs the water saturation at the non-wetting phase shock front
    double64 ShockHeight( const Element<dim>* const ) const;
    
    /// calculated shock height and speed
    void     ShockSpeedAndHeight( const Element<dim>* const, double64& speed, double64& height ) const;
    
    /// multipliers for gravity-driven flow (advection multiplier and source term)
    double64 GravityTerm(  const Element<dim>* const ) const;
    
    double64 GravityMultiplier_G(  const Element<dim>* const ) const;
    
    double64 GravityMultiplier_dGds(  const Element<dim>* const ) const;

    /// multiplier for diffusion coefficient in the case of non-linear diffusion uses viscosity(phase)
    double64 DiffusionMultiplier(  const Element<dim>* const, size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double64 CapillaryDiffusionMultiplier( const Element<dim>* const ) const;
    
    double64 CapillaryDiffusionMultiplier_Phase(  const Element<dim>* const, size_t phase ) const;

    /// calculating the Shock velocity
    double64 ShockFrontVelocity( const Element<dim>* const ) const ;
    
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

#endif /* CSMP_TWO_PHASE_FLOW_FUNCTIONS_H */
