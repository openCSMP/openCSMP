//
//  FlowFunctions.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 7/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLOW_FUNCTIONS_H
#define CSMP_FLOW_FUNCTIONS_H

#include "SaturationFunction.h"
#include "BrooksCoreySaturationFunctions.h"
#include "Variables_TwoPhaseFlow.h"
#include "Fluid.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Region;
template<size_t> class Model;

/// auxiliary functions for spline interpolation
double64 spline_value( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
double64 spline_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
double64 spline_second_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
  
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
    
    TODO: specify through template parameter for what PLACEMENT/ipoint the flow functions shall be initialised
*/
template<size_t dim>
class FlowFunctions : public variables::Variables_TwoPhaseFlow,              ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
                      public SaturationFunction<dim,FlowFunctions>,
                      public BrooksCoreySaturationFunctions<dim,FlowFunctions>,  ///< placeholder for saturation function model
                      public Fluid<dim,FlowFunctions> {                          ///< placeholder for fluids module / EOS interface
  public:
    FlowFunctions();
    ~FlowFunctions();
    
    /// current water saturation initialised inside of the model
    double64 Sw() const;
    
    /// lambda parameter: 0 for water, 1 for the non-wetting phase
    double64 Mobility( double64 sw, size_t phase ) const;

    /// d lambda_i / dsw
    double64 MobilityDerivative( double64 sw, size_t phase, bool evaluate_numerically=true ) const;
 
    /// lambda_t: sum of phase mobilities
    double64 TotalMobility( double64 sw ) const;
    
    /// lambda overbar: l1 * l2 / l1 + l2 = mobility product / total mobility also known as G
    double64 MobilityProduct( double64 sw ) const;

    /// d lambda overbar / dsw also known as dGds
    double64 MobilityProductDerivative( double64 sw, bool evaluate_numerically=true ) const;

     /// fractional flow; 0=water, 1=non-wetting phase
    double64 f( double64 sw, size_t phase ) const;

    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    double64 dfds( double64 sw, size_t phase, bool evaluate_numerically=true ) const;
    
    /// maximum value of previous derivative
    double64 MaxFractionalFlowDerivative() const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
   	double64 AdvectionMultiplier( double64 sw, bool evaluate_numerically=true ) const;

    /// linearized fractional flow derivative
    double64 ShockSpeed() const;
    double64 ShockHeight() const;
    void ShockSpeedHeight( double64& speed, double64& height) const;

    /// multipliers for gravity-driven flow (advection multiplier and source term)
    double64 GravityTerm( double64 sw ) const;
    double64 GravityMultiplier_G( double64 sw ) const;
    double64 GravityMultiplier_dGds( double64 sw, bool evaluate_numerically=true ) const;

    /// multiplier for diffusion coefficient in the case of non-linear diffusion uses viscosity(phase)
    double64 DiffusionMultiplier( double64 sw, size_t phase ) const;
    
    /// driven by capillary pressure gradient
    double64 CapillaryDiffusionMultiplier( double64 sw ) const;
    
  private:
    double64 dfds_Numerical( double64 sw, size_t phase, double64 h = 0.001 ) const;
    double64 dGds_Numerical( double64 sw, double64 h = 0.000001 ) const;
    double64 dlwds_Numerical( double64 sw, double64 h = 0.001 ) const;
    double64 dlnds_Numerical( double64 sw, double64 h = 0.001 ) const;
};

} // end csmp

#endif /* CSMP_FLOW_FUNCTIONS_H */
