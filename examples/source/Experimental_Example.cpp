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

const size_t DIM(3U);

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
      enum{dim=2};
      const double ym (1000.0), pr(0.3), P0(10.0);
      bool quarterpoint = false;
  
      // Model configuration:
      ANSYS_Model2D  model( "Fluid_Flower", false, true, false );    // Constractor for empty variables
      for ( auto& E : model.Region("FRACTURE").CellVector() ) {
          std::cout << "Nbrs -> " << E->ConnectedNeighbors() << std::endl;
          if (E->ConnectedNeighbors() == 1){
              E->Out();
          }
      }
/*
  // 1. create point property mapper - Fluid_Flower-points is a csv file
  PointPropertyToCellMapper2D  mapper( "Fluid_Flower-points" );
  mapper.Out();

  // 2. Bring mode thickness data in interpolate them across the model
  mapper.MapPointDataToElements( ansys_model, "Model", "thickness" );
  vtk_output.OutputDataToVTK( ansys_model, "Fluid_Flower-thickness", "thickness",  0 );

  // 3. Output the interpolated values into another CSV file
  mapper.MapNodeToPointData( ansys_model, "Model", "thickness" );
  mapper.OutputPointDataToCSV_File( "Fluid_Flower-interpolated_thickness", "thickness" );
*/
} // end Run

} // csmp
