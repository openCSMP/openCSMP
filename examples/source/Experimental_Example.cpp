//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "Experimental_Example.h"
#include "compareFloats.h"
#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "PDE_IntegratorExperimental.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include <iostream>
#include "InputDataManager.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"

// include any header files that you need here...
#include "SteadyStateDiffusor.h"
#include "CSMP_highLevelUtilities.h"
#include "ComputationalSettings.h"
#include "IAPWS_H2OPropertiesVisitor.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "StressesAndStrains.h"
#include "StressesAndStrains2.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "SinglePhaseVelocityVisitor.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }

// OTHER RELATIONSHIPS

/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.
     
     Systems Modelling and Design model that calculates the effective stress in a dam, 
     using the gravitational loading and the plane stress assumption.
     
     input models:
     - slope_model1
*/
void Experimental_Example::Run()
{
    cout <<"\nHello World and size of uint_fast32_t: "<< sizeof(uint_fast32_t) << endl;
    cout <<"\nauto i{0}: "<< sizeof(uint_fast32_t) << endl;
    vector<double> doubs(1e9,2e-4);
    for ( auto j{0}; j<doubs.size(); j++ ) {
          doubs[j] = 2.3;
          cout <<"\n\tsize of vector loop variable j: "<< sizeof(j) << endl;
          if ( j == 1 ) break;
      }
      
#if 0
    // ODLING 720 x 720 meter
    string         variables_file("CSMP-1phase-variables.txt");
    ANSYS_Model2D  model( "odling720x720", variables_file.c_str(), false, true, true );
    // boolean flags set reading to: 2) default prop.values, 3) group prop.values, 4) essential conditions for box-shaped model
    InputDataManager<2>().ConfigureFromFile( model, "Fluid_Flower",
                                             false,           // region name from parameter range
                                             true,            // default property values
                                             true,            // regional property values
                                             true,            // boundary conditions for box-shaped model
                                             true );          // essential conditions for regions
                                           // default: boundary conditions for arbitrary-shaped model
#endif

    // FLUID FLOWER TESTCASE
//#if 0
    string  variables_file("DES_2phase_variables.txt");
    ANSYS_Model2D  model( "Fluid_Flower", variables_file.c_str(), false, true, true );
    
    // Assignment of material properties, initial conditions, and boundary
    InputDataManager<2>  model_configuration;

    // boolean flags set reading to: 2) default prop.values, 3) group prop.values, 4) essential conditions for box-shaped model
    model_configuration.ConfigureFromFile( model, "Fluid_Flower",
                                           false,           // region name from parameter range
                                           true,            // default property values
                                           true,            // regional property values
                                           true,            // boundary conditions for box-shaped model
                                           true );          // essential conditions for regions
                                           // default: boundary conditions for arbitrary-shaped model
//#endif

} // end Run

} // csmp
