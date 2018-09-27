//
//  CO2H2O_FunctionsModule1.h
//  CSMP_GitHub
//
//  Created by Mahyar Madadi on 16/July/2018.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_CO2_H2O_FUNCTIONS_MODULE_H
#define CSMP_CO2_H2O_FUNCTIONS_MODULE_H

#include "VariableSet_CO2GeoSequestration.h"
#include "BrooksCoreySaturationFunctions.h"
#include "BrooksCoreySaturationFunctionswithHysteresis.h"
#include "ExperimentalSaturationFunctions.h"
#include "TwoPhaseFlowFunctions.h"
#include "Fluid.h"

namespace csmp {

template<size_t> class PropertyDatabase;

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
class CO2H2O_FunctionsModule0 : public variables::VariableSet_CO2GeoSequestration,     ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctions<dim,CO2H2O_FunctionsModule0>,  ///< saturation function model
    public TwoPhaseFlowFunctions<dim,CO2H2O_FunctionsModule0>,                         ///< mobilities etc.
    public Fluid<dim,CO2H2O_FunctionsModule0> {                                        ///< fluids module / EOS interface
      
  public:
    explicit CO2H2O_FunctionsModule0( const PropertyDatabase<dim>& );
};

// USE THIS TYPE RATHER THAN THE COMPLEX TEMPLATE
typedef CO2H2O_FunctionsModule0<3U>  ACGSS_FlowFunctions;


template<size_t dim>
class CO2H2O_FunctionsModule1 : public variables::VariableSet_CO2GeoSequestration,     ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctionsWithHysteresis<dim,CO2H2O_FunctionsModule1>,  ///< saturation function model
    public TwoPhaseFlowFunctions<dim,CO2H2O_FunctionsModule1>,                         ///< mobilities etc.
    public Fluid<dim,CO2H2O_FunctionsModule1> {                                        ///< fluids module / EOS interface
      
  public:
    explicit CO2H2O_FunctionsModule1( const PropertyDatabase<dim>& );
};

// USE THIS TYPE RATHER THAN THE COMPLEX TEMPLATE
typedef CO2H2O_FunctionsModule1<3U>  ACGSS_HystereticFlowFunctions;



template<size_t dim>
class CO2H2O_FunctionsModule2 : public variables::VariableSet_CO2GeoSequestration,
                                public ExperimentalSaturationFunctions<dim,CO2H2O_FunctionsModule2>,
                                public TwoPhaseFlowFunctions<dim,CO2H2O_FunctionsModule2>,
                                public Fluid<dim,CO2H2O_FunctionsModule2> {
      
  public:
    explicit CO2H2O_FunctionsModule2( const PropertyDatabase<dim>& );
};

typedef CO2H2O_FunctionsModule2<3U>  ACGSS_ExperimentalFlowFunctions;



} // end csmp

#endif /* CSMP_CO2_H2O_FUNCTION_MODULE_H */
