//
//  FlowFunctionsDraft2.h
//  CSMP_GitHub
//
//  Created by Mahyar Madadi on 16/July/2018.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLOW_FUNCTIONS_DRAFT2_H
#define CSMP_FLOW_FUNCTIONS_DRAFT2_H

#include "BrooksCoreySaturationFunctionsWithHysteresis.h"
#include "VariableSet_CO2GeoSequestration.h"
#include "VariablePlacement.h"
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
class FlowFunctionsDraft2 : public variables::VariableSet_CO2GeoSequestration,              ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctionsWithHysteresis<dim,FlowFunctionsDraft2>,  ///< placeholder for saturation function model
    public Fluid<dim,FlowFunctionsDraft2>         {                          ///< placeholder for fluids module / EOS interface
      
  public:
    
    FlowFunctionsDraft2(const PropertyDatabase<dim>& db);
    
    /// current water saturation initialised inside of the model
    template<class TARGET_PLACEMENT>
    double64 Sw( TARGET_PLACEMENT& p ) const { return p.Obtain(this->key_sH2O); }
      
    /// lambda parameter: 0 for water, 1 for the non-wetting phase
    template<class TARGET_PLACEMENT>
    double64 Mobility( TARGET_PLACEMENT&, size_t phase ) const;
    
    /// lambda parameter: 0 for water, 1 for the non-wetting phase (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 Mobility( TARGET_PLACEMENT&, size_t phase, double64 ) const;     

    /// d lambda_i / dsw
    template<class TARGET_PLACEMENT>
    double64 MobilityDerivative( TARGET_PLACEMENT&, size_t phase ) const;
 
    template<class TARGET_PLACEMENT>
    double64 MobilityDerivative( TARGET_PLACEMENT&, size_t phase, double64 ) const;
                        
    /// lambda_t: sum of phase mobilities
    template<class TARGET_PLACEMENT>
    double64 TotalMobility( TARGET_PLACEMENT& ) const;
    
    /// lambda_t: sum of phase mobilities (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 TotalMobility( TARGET_PLACEMENT&, double64 ) const;     
    
    /// lambda overbar: l1 * l2 / l1 + l2 = mobility product / total mobility also known as G
    template<class TARGET_PLACEMENT>
    double64 MobilityProduct( TARGET_PLACEMENT& ) const;

    /// d lambda overbar / dsw also known as dGds
    template<class TARGET_PLACEMENT>
    double64 MobilityProductDerivative( TARGET_PLACEMENT& ) const;
                        
    template<class TARGET_PLACEMENT>
    double64 MobilityProductDerivative( TARGET_PLACEMENT&, double64 ) const;

     /// fractional flow; 0=water, 1=non-wetting phase
    template<class TARGET_PLACEMENT>
    double64 f( TARGET_PLACEMENT&, size_t phase ) const;
    
     /// fractional flow; 0=water, 1=non-wetting phase  (using prescribed sw)
    template<class TARGET_PLACEMENT>
    double64 f( TARGET_PLACEMENT&, size_t phase, double64 ) const;     
    
    /// Permeability
    template<class TARGET_PLACEMENT>
    double64 Permeability( TARGET_PLACEMENT& ) const;    

    /// derivative of fractional flow (used in advection multiplier); note that fw+fn=1, dfw_dsw=dfn_dsn
    template<class TARGET_PLACEMENT>
    double64 dfds( TARGET_PLACEMENT&, size_t phase ) const;
    
    template<class TARGET_PLACEMENT>
    double64 dfds( TARGET_PLACEMENT&, size_t phase, double64 ) const;
                        
    /// maximum value of previous derivative
    template<class TARGET_PLACEMENT>
    double64 MaxFractionalFlowDerivative( TARGET_PLACEMENT& ) const;

    /// multiplier for advection viscosity coefficient in the case of non-linear advection
    template<class TARGET_PLACEMENT>
   	double64 AdvectionMultiplier( TARGET_PLACEMENT& ) const;

    template<class TARGET_PLACEMENT> double64 ShockSpeed( TARGET_PLACEMENT& ) const;
    template<class TARGET_PLACEMENT> double64 ShockHeight( TARGET_PLACEMENT& ) const;

    /// multipliers for gravity-driven flow (advection multiplier and source term)
    template<class TARGET_PLACEMENT>
    double64 GravityTerm( TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 GravityMultiplier_G( TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 GravityMultiplier_dGds( TARGET_PLACEMENT& ) const;

    /// multiplier for diffusion coefficient in the case of non-linear diffusion uses viscosity(phase)
    template<class TARGET_PLACEMENT>
    double64 DiffusionMultiplier( TARGET_PLACEMENT&, size_t phase ) const;
    
    /// driven by capillary pressure gradient
    template<class TARGET_PLACEMENT>
    double64 CapillaryDiffusionMultiplier( TARGET_PLACEMENT& ) const;
    
    template<class TARGET_PLACEMENT>
    double64 CapillaryDiffusionMultiplier_Phase( TARGET_PLACEMENT&, size_t phase ) const;
};

} // end csmp

#endif /* CSMP_FLOW_FUNCTIONS_DARAFT2_H */
