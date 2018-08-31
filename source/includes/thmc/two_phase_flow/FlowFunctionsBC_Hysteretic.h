//
//  FlowFunctions.h
//  CSMP_GitHub
//
//  Created by Mahyar Madadi on 16/July/2018.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLOW_FUNCTIONS_BC_HYSTERETIC_H
#define CSMP_FLOW_FUNCTIONS_BC_HYSTERETIC_H

#include "VariableSet_CO2GeoSequestration.h"
#include "BrooksCoreySaturationFunctionswithHysteresis.h"
#include "Fluid.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Region;
template<size_t> class Model;

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
    
    @attention the phases are numbered 0..2, phase 0 is the water phase (usually the wetting phase)
    
    @attention for simplicity and efficiency, fluid property values are always computed at the nodes by respective EOS modules.
    Subsequently these values are interpolated to the points of interest.
    
    Code serves as an example of how static polymorphism can be used to implement constitutive relationships for multiphase flow.
    
    TODO: specify through template parameter for what PLACEMENT/ipoint the flow functions shall be initialised
*/
template<size_t dim>
class FlowFunctionsBC_Hysteretic : public variables::VariableSet_CO2GeoSequestration,    ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctionsWithHysteresis<dim,FlowFunctionsBC_Hysteretic>, ///< placeholder for saturation function model
    public Fluid<dim,FlowFunctionsBC_Hysteretic> {                                       ///< placeholder for fluids module / EOS interface
      
  public:
    
    FlowFunctionsBC_Hysteretic(const PropertyDatabase<dim>& db);
    
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

    /// output Shock Speed
    template<class TARGET_PLACEMENT>
    double64 ShockSpeed(  const TARGET_PLACEMENT& ) const;
    
    /// output Shock Height
    template<class TARGET_PLACEMENT>
    double64 ShockHeight( const TARGET_PLACEMENT& ) const;
    
    /// calculated Shock Height and Speed
    template<class TARGET_PLACEMENT>
    void     ShockSpeedHeight( const TARGET_PLACEMENT&, double64& speed, double64& height) const;
    
    /// ???
    template<class TARGET_PLACEMENT>
    double64 ShockSaturation(  const TARGET_PLACEMENT&, size_t, bool evaluate_numerically=false) const;

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


    /// calculating the  Inflection point
    template<class TARGET_PLACEMENT>
    double64 InflectionPointSaturation( const TARGET_PLACEMENT& ) const ;
    
    /// calculating the Tangent Point
    template<class TARGET_PLACEMENT>
    double64 TangentPointSaturation( const TARGET_PLACEMENT& ) const ;
    
    /// calculating the Shock velocity
    template<class TARGET_PLACEMENT>
    double64 ShockFrontVelocity( const TARGET_PLACEMENT& ) const ;
    
    /// calculating  Buckley Leverett function... this has to be zero at shock point
    template<class TARGET_PLACEMENT>
    double64 BuckleyLeverettFunction( const TARGET_PLACEMENT&, double64 S) const ;
    
    /// finding root of  Buckley Leverett function... which results is shock point saturation.
    template<class TARGET_PLACEMENT>
    double64 FindRootSecantMethod( const TARGET_PLACEMENT&, double64 S1, double64 S2) const ;
 
           
  private:
      
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
};

// ALWAYS JUST USE THIS TYPE RATHER THAN THE COMPLEX TEMPLATE
typedef FlowFunctionsBC_Hysteretic<3U>  ACGSS_FlowFunctions;



} // end csmp

#endif /* CSMP_FLOW_FUNCTIONS_BC_HYSTERETIC_H */
