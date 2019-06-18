//
//  FlowFunctionsModule.h
//  CSMP_GitHub
//
//  Created by Mahyar Madadi on 16/July/2018.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLOW_FUNCTIONS_MODULE_H
#define CSMP_FLOW_FUNCTIONS_MODULE_H

#include "VariableSet_CO2GeoSequestration.h"
#include "BrooksCoreySaturationFunctions.h"
#include "BrooksCoreySaturationFunctionswithHysteresis.h"
#include "ExperimentalSaturationFunctions.h"
#include "TwoPhaseFlowFunctions.h"
#include "H2O_CO2_NaCl_FlowFunctions.h"
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
class FlowFunctionsModule1 : public variables::VariableSet_CO2GeoSequestration,               ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
                             public BrooksCoreySaturationFunctions<dim,FlowFunctionsModule1>, ///< saturation function model
                             public TwoPhaseFlowFunctions<dim,FlowFunctionsModule1>,          ///< mobilities etc.
                             public Fluid<dim,FlowFunctionsModule1> {                         ///< fluids module / EOS interface
      
   public:
     FlowFunctionsModule1( const PropertyDatabase<dim>&, double64 acc_gravity );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)
};

// USE THIS TYPE RATHER THAN THE COMPLEX TEMPLATE
typedef FlowFunctionsModule1<3U>  ACGSS_SlightlyCompressible2PhaseFlowFunctions;


template<size_t dim>
class FlowFunctionsModule2 : public variables::VariableSet_CO2GeoSequestration,                             ///< all variables in transport scheme
                             public BrooksCoreySaturationFunctionsWithHysteresis<dim,FlowFunctionsModule2>, ///< saturation function model
                             public TwoPhaseFlowFunctions<dim,FlowFunctionsModule2>,                        ///< mobilities etc.
                             public Fluid<dim,FlowFunctionsModule2> {                                       ///< fluids module / EOS interface
      
   public:
     FlowFunctionsModule2( PropertyDatabase<dim>&, double64 acc_gravity );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)
};

typedef FlowFunctionsModule2<3U>  ACGSS_Hysteretic_SlightlyCompressible2PhaseFlowFunctions;



template<size_t dim>
class FlowFunctionsModule3 : public variables::VariableSet_CO2GeoSequestration,
                             public ExperimentalSaturationFunctions<dim,FlowFunctionsModule3>,
                             public TwoPhaseFlowFunctions<dim,FlowFunctionsModule3>,
                             public Fluid<dim,FlowFunctionsModule3> {
      
   public:
     /// initialises saturation functions from file and set range of rocktype values accordingly in property database
     FlowFunctionsModule3( PropertyDatabase<dim>&, double64 acc_gravity, const char* model_name="ACGSS_simulator" );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)
};

typedef FlowFunctionsModule3<3U>  ACGSS_Experimental_SlightlyCompressible2PhaseFlowFunctions;




template<size_t dim>
class FlowFunctionsModule4 : public variables::VariableSet_CO2GeoSequestration,     ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctions<dim,FlowFunctionsModule4>,  ///< saturation function model
    public H2O_CO2_NaCl_FlowFunctions<dim,FlowFunctionsModule4>,      ///< mobilities etc.
    public Fluid<dim,FlowFunctionsModule4> {                          ///< fluids module / EOS interface
      
   public:
     FlowFunctionsModule4( const PropertyDatabase<dim>&, double64 acc_gravity );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)
};

typedef FlowFunctionsModule4<3U>  ACGSS_Compositional_H2O_CO2_NaCl_FlowFunctions;


template<size_t dim>
class FlowFunctionsModule5 : public variables::VariableSet_CO2GeoSequestration,     ///< all variables in transport scheme (and determining the ones that will be included in the initialisation)
    public BrooksCoreySaturationFunctionsWithHysteresis<dim,FlowFunctionsModule5>,  ///< saturation function model
    public H2O_CO2_NaCl_FlowFunctions<dim,FlowFunctionsModule5>,                    ///< mobilities etc.
    public Fluid<dim,FlowFunctionsModule5> {                                        ///< fluids module / EOS interface
      
   public:
     FlowFunctionsModule5( PropertyDatabase<dim>&, double64 acc_gravity );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)};
};

typedef FlowFunctionsModule5<3U>  ACGSS_Hysteretic_Compositional_H2O_CO2_NaCl_FlowFunctions;



template<size_t dim>
class FlowFunctionsModule6 : public variables::VariableSet_CO2GeoSequestration,
                             public ExperimentalSaturationFunctions<dim,FlowFunctionsModule6>,
                             public H2O_CO2_NaCl_FlowFunctions<dim,FlowFunctionsModule6>,
                             public Fluid<dim,FlowFunctionsModule6> {
      
   public:
     /// initialises saturation functions from file and set range of rocktype values accordingly in property database
     FlowFunctionsModule6( PropertyDatabase<dim>&, double64 acc_gravity, const char* model_name="ACGSS_simulator" );
     double64 acceleration_of_gravity_; ///< this must be read and input from the model (key_g)};
};

typedef FlowFunctionsModule6<3U>  ACGSS_Experimental_Compositional_H2O_CO2_NaCl_FlowFunctions;

} // end csmp

#endif /* CSMP_FLOW_FUNCTION_MODULE_H */
