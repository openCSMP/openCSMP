// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  PressureDependentBlackOilProperties_Example.cpp
//  csmp
//
//  Created by Lukas Mosser on 7/7/14.
//

#include "PressureDependentBlackOilProperties_Example.h"

#include "Triangulator.h"
#include "Model.h"
#include "Region.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

#include "PDE_Integrator.h"
#include "Integral_NT_op_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

#include "VTK_Interface.h"
#include "MatlabInterface.h"
#include "VTU_Interface.h"
#include "CSMP_highLevelUtilities.h"

#include "ConstantFactor.h"
#include "ModelComparator.h"

// example specific include files
#include "UndersaturatedBOCorrelations.h"
#include "BO_LookupVisitor.h"
#include "FluidPropertyLookupTable.h"

using namespace std;

namespace csmp {
  
void PressureDependentBlackOilProperties_Example::Specifications()
{
  SetTitle( "Calculating Pressure dependent properties for Black Oil form Lookup Tables" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "Lukas J. Mosser" );
  AddDescription( "Calculation of pressure dependent black oil properties using Picard iteration." );
  AddDescription( "source in: PressureDependentBlackOilProperties_Example.cpp" );
  AddRequirement( "variables file is: 'tutorial1_variables.txt'");
  AddRequirement( "Triangulator mesher input file: 'tutorial1_input'");
  AddRequirement( "correlation table: 'BlackOilProject.ptb' must also be present in run directory." );
}
 
  
void PressureDependentBlackOilProperties_Example::Run()
  {
    // -------------------------------------
    // 0.0 Set the global variables for CSMP
    // -------------------------------------
    double model_time =  0.0; // time
    
    // -----------------------------------------------------------------------
    // 1.0 Build the CSMP Model from a simple color-coded text input file
    //     that generates a uniform triangular FE mesh. The function
    //     readTextPixelData() reads in the file, generates the FE mesh
    //     (of uniform triangles) and returns it as a VSet
    // -----------------------------------------------------------------------
    auto       vset = readTextPixelData();
    Model<2U>  model( vset, "tutorial1-variables.txt" ); // txt file defines the physical variables to be used in the simulation
    // give the model dimensions
    printModelDimensions( model, true );
      
    
    // -----------------------------------------------------------------------
    // 2.0 Now we apply boundary and initial conditions (this can also be done,
    //     more conveniently, in a configuration file for more realistic runs
    // -----------------------------------------------------------------------
    std::list<string> property_names;
    property_names.push_back("fluid pressure");
    property_names.push_back("viscosity oil");
    property_names.push_back("density oil");
    property_names.push_back("oil viscosity element");
    property_names.push_back("compressibility");
    property_names.push_back("conductivity");
    property_names.push_back("previous fluid pressure");
    property_names.push_back("trial fluid pressure");
    property_names.push_back("nodal compressibility");
    
    
    // -----------------------------------------------------------------------
    // 3.0 PetroleumExperts - PVT table from file
    // -----------------------------------------------------------------------
    /*
    string file_name = "BlackOilProject.ptb";
    std::vector<string> properties;
    properties.push_back("viscosity oil");
    properties.push_back("density oil");
    properties.push_back("compressibility oil");
    //properties.push_back("compressibility oil");
    FluidPropertyLookupTable  table( file_name.c_str() ); // TODO: reading this table does not work under linux ?!
    */
    
    // Oil properties from   - BO_PropertyCalculatorInterface
    // ------------------------------------------------------
    double p = 3.0e7;  // 300 bars
    double T = 120.;   // 120 degrees C
    BO_InputProperties    bo_input_props; // TODO: this struct must be initialised from property table
    BO_PropertyCalculator bo_calculator( bo_input_props );
    bo_calculator.Update( p, T );
  
    // water properties from - ReservoirWaterPropertyCalculatorInterface
    // -----------------------------------------------------------------
    // ReservoirWaterPropertyCalculatorInterface     compressibility, density, viscosity
    ConstantReservoirWaterProperties  water_input_properties( 1.e-09,  1000.0,    0.0001 );
    
    // visitor that computes the property distribution across the model
    // ----------------------------------------------------------------
    // (first the model must be initialised with some p, T distribution)
    BO_LookupVisitor<2>  bovisit( model, bo_calculator, water_input_properties );
    
    // assigning material properties
    model.InputPropertyValue( "porosity",         makeScalar(PLAIN,0.1) );     // always as a fraction
    model.InputPropertyValue( "permeability",     makeScalar(PLAIN,1.0e-15) ); // always in m2 (comment out if heterogeneous k-field is used in input file)
    model.InputPropertyValue( "compressibility",  makeScalar(PLAIN,5.0e-10) ); // for fluid and rock, in Pa-1
    model.InputPropertyValue( "nodal compressibility",  makeScalar(PLAIN,5.0e-10) );
    model.InputPropertyValue( "viscosity oil", makeScalar(PLAIN,1e-03));
    model.InputPropertyValue( "oil viscosity element", makeScalar(PLAIN,1e-03));
    model.InputPropertyValue( "density oil", makeScalar(PLAIN,1e03));
    model.InputPropertyValue( "conductivity", makeScalar(PLAIN, 1.0E-15/1.0E-03));
    // assigning initial conditions
    model.InputPropertyValue( "fluid pressure",      makeScalar(PLAIN,1.0e+07) );
    model.InputPropertyValue( "previous fluid pressure", makeScalar(PLAIN,1.0e+07) );  // always in Pascal
    model.InputPropertyValue( "trial fluid pressure", makeScalar(PLAIN,1.0e+07) );
    model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.0) );      // no sources/sinks (units m3 m-2 s-1)
    // assigning boundary conditions for fluid pressure at the LEFT and RIGHT model boundaries
    // such that a pressure wave travels from left to right through the model
    model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,3.0e+07) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.0e+07) );
    model.InputBoundaryValue( LEFT,  "previous fluid pressure", makeScalar(DIRICH,3.0e+07) );
    model.InputBoundaryValue( RIGHT, "previous fluid pressure", makeScalar(DIRICH,1.0e+07) );
    model.InputBoundaryValue( LEFT,  "trial fluid pressure", makeScalar(DIRICH,3.0e+07) );
    model.InputBoundaryValue( RIGHT, "trial fluid pressure", makeScalar(DIRICH,1.0e+07) );
    
    csmp::Index fluid_pressure_key = model.Database().StorageKey("fluid pressure");
    csmp::Index previous_fluid_pressure_key = model.Database().StorageKey("previous fluid pressure");
    csmp::Index trial_fluid_pressure_key = model.Database().StorageKey("trial fluid pressure");
    csmp::Index conductivity_key = model.Database().StorageKey("conductivity");
    csmp::Index oil_viscosity_element_key = model.Database().StorageKey("oil viscosity element");
    csmp::Index compressibility = model.Database().StorageKey("compressibility");
    csmp::Index permeability_key = model.Database().StorageKey("permeability");
    csmp::Index nodal_compressibility_key = model.Database().StorageKey("compressibility");
    
    //Update Properties at t=0 to Bo Properties at t=0
    model.Accept( bovisit );
    
    Region<2U>& model_domain = model.Region("Model");
    
    double visc_temp = 0.0;
    double permeab_temp = 0.0;
    double cond_temp = 0.0;
    
    model_domain.InterpolateNodeToCellProperty( "viscosity oil", "oil viscosity element" );
    model_domain.InterpolateNodeToCellProperty( "nodal compressibility", "compressibility");
    
    for( auto itel = model_domain.CellsBegin(); itel!=model_domain.CellsEnd(); ++itel )
    {
      visc_temp = (*itel)->Read(oil_viscosity_element_key);
      permeab_temp = (*itel)->Read(permeability_key);
      cond_temp = permeab_temp/visc_temp;
      (*itel)->Store(conductivity_key, makeScalar(PLAIN,cond_temp));
    }
    //End Update Properties
    
    //Collocate the PDE Integrator
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    SAMG_Solver   solver(&settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
    PDE_Integrator<2U,Element>  fluid_pressure(solver);
    
    // LHS stiffness matrix                              operand         basis function    test function
    Integral_dNT_op_dN_dV<2U>  stiffness_matrix( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    
    // LHS mass matrix
    Integral_NT_lhsop_N_dV<2U> mass_matrix_lhs( model.Database(), "compressibility", "fluid pressure", "fluid pressure" );
    
    // RHS mass vector
    Integral_NT_op_N_dV<2U>    mass_matrix_rhs( model.Database(), "compressibility", "fluid pressure" );
    
    // RHS mass vector for integrating source term
    Integral_NT_op_N_dV<2U>    source_term( model.Database(), "fluid volume source", "fluid pressure" );
    
    // mass matrices for dp/dt term must be divided by time increment
    mass_matrix_lhs.MultiplyWithTimeIncrement(true);
    mass_matrix_rhs.MultiplyWithTimeIncrement(true);
    
    // use lumped formulation for all mass matrices (i.e., diagonalise matrices)
    mass_matrix_lhs.LumpedFormulation(true);
    mass_matrix_rhs.LumpedFormulation(true);
    source_term.LumpedFormulation(true);
    
    // evalute source term last
    source_term.AddAccumulateLater();
    
    // now add each FE operation (i.e., PDE Operator) to the FE algorithm
    fluid_pressure.Add( &stiffness_matrix );
    fluid_pressure.Add( &source_term );
    fluid_pressure.Add( &mass_matrix_lhs );
    fluid_pressure.Add( &mass_matrix_rhs );
    
    // -----------------------
    // 5.0 Time Loop Variables
    // -----------------------
    
    // define some constant variables
    const double     hour(3600.0);
    double           time_increment = 2 * hour; // timestep 2 hours
    long	           time;
    
    // set the time increment for the FE algorithm
    fluid_pressure.TimeIncrement( 1.0/time_increment );

    cout << "Output to VTU File" << endl;
    VTU_Interface<2U> vtu(model);
    vtu.OutputDataToVTU("Model",property_names,"Model",model_time);
    
    
    cout << "Before Picard iteration loop" << endl;
    double error    = 0.;
    double l2norm   = 1.0;
    double l2steady = 1000.0;
    double epsilon  = 1E-5;
    double steady_state_criterion = 1.0E+02;
    
    bool converged          = false;
    size_t total_iterations = 0;
    size_t iteration        = 0;
    
    // -----------------------
    // 6.0 Transient loop
    // -----------------------
    while ( l2steady > steady_state_criterion )
    {
      model.CopyReplace("fluid pressure", "previous fluid pressure");
      converged = false;
      
      while(converged == false)
      {
        error = 0.;
        l2norm = 0.;
        
        model.CopyReplace("previous fluid pressure", "fluid pressure");
        model.Apply(fluid_pressure);
        
        //Update the Properties
        model.Accept(bovisit);
        model_domain.InterpolateNodeToCellProperty( "viscosity oil", "oil viscosity element" );
        model_domain.InterpolateNodeToCellProperty( "nodal compressibility", "compressibility");
        
        //Iterate over the elements to update conductivity at each element from oil viscosity and permeability
        for(auto itel = model_domain.CellsBegin(); itel!=model_domain.CellsEnd(); ++itel )
        {
          visc_temp = (*itel)->Read(oil_viscosity_element_key);
          permeab_temp = (*itel)->Read(permeability_key);
          cond_temp = permeab_temp/visc_temp;
          (*itel)->Store(conductivity_key, makeScalar(PLAIN,cond_temp));
        }
        //end Update properties
        
        //Calculate Error between P at t0 +dt with initial properties with P at t0 + dt with properties at P at t0 + dt
        for(auto itnod = model_domain.NodesBegin(); itnod!=model_domain.NodesEnd(); ++itnod )
        {
          error += std::pow(std::fabs((*itnod)->Read(fluid_pressure_key) - (*itnod)->Read(trial_fluid_pressure_key)),2);
        }
        l2norm = std::sqrt( error );
        cout << "l2 norm: " << l2norm << endl;
        //End Error Calculation
        
        
        //Check for Convergence
        if(l2norm > epsilon){ //Compare P(1,1) to P(1,0)
          model.CopyReplace("fluid pressure", "trial fluid pressure");
          cout << "Iteration: "<< iteration << endl;
          ++iteration;
        }
        else{
          cout << "CONVERGED L2Norm: " << l2norm << endl;
          total_iterations += iteration;
          iteration  = 0;
          converged = true;
        }
      }
      vtu.OutputDataToVTU("TransientModel",property_names,"Model",model_time);
      
      //Calculate if steady state criterion met.
      error = 0.;
      l2steady = 1.0E+3;
      for( auto itnod = model_domain.NodesBegin(); itnod!=model_domain.NodesEnd(); ++itnod )
      {
        error += std::pow(std::fabs((*itnod)->Read(fluid_pressure_key) - (*itnod)->Read(previous_fluid_pressure_key)),2);
      }
      l2steady = std::sqrt( error );
      cout << "Steady State Norm: " << l2steady << endl;
      time = static_cast<long>(model_time/hour);
      model_time += time_increment;
    }
    
    cout <<"\nmain: That's it..."<< endl;
  }
  
} // csmp
