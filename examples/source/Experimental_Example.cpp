// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//

#include "Experimental_Example.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

// TODO: include the declarations that you need here !
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"

#include "ANSYS_Model2D.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "RegionMonitor.h"

#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "ConstantFactor.h"            // interrelation to calculate hydraulic conductivity

#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"  // conductance matrix (LHS)
#include "NumIntegral_NT_op_N_dS.h"    // boundary integral
#include "NumIntegral_NT_lhsop_N_dV.h" // capacitance matrix LHS
#include "NumIntegral_NT_rhsop_N_dV.h" // capacitance matrix (lumped) RHS

#include "VelocityAndVolumeFlux.h"     // post-processing of Darcy velocity

#include "VTU_Interface.h"


using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Simulation of Physical Processes" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }
  

// total systems compressilicity, ct (1/Pa) = (confined aquifer) storativity, S, but in pressure form
double storativity( double phi, double beta_f, double beta_r ) {
   return phi * beta_f + (1.-phi) * beta_r;
}

/**
      Use this Experimental_Example  to experiment with the functionality of csmp
      TODO: write your code into the following Run() member function of 'Experimental_Example'
 */
void Experimental_Example::Run()
{
  constexpr uint32_t dim{2};
  string example_name  = "Experimental_Example";
  string model_name    = "Jura-slope1";
  string variable_file = "CSMP_field_scale_mechanics_variables.txt";
  string config_file = model_name;
  //find the name of this current example source file
  string file_name = GetExampleFileName(__FILE__);
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(example_name, model_name, variable_file, config_file);

  Model<dim>  model( model_name, variable_file );
  printModelDimensions(model);
  model.RegionsOut();
  model.BoundariesOut();


    cout <<"\nmain: That's it..."<< endl;

    // 9.  Output the initial range of the variables
    // ---------------------------------------------
    printRangeOfVariable( model, "fluid pressure" );
    
    cout <<"\n end of program"<< endl;
 
 } // end run


} // csmp
