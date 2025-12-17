
//
//  RegionMonitor_Example.cpp
//  CSMP_API_library2025
//
//  Created by Stephan Matthai on 10/6/2025
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "RegionMonitor_Example.h"

#include "RegionMonitor.h"

#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"

#include "ANSYS_Model2D.h"
#include "Model.h"
#include "Region.h"

#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "ConstantFactor.h"            // interrelation to calculate hydraulic conductivity

#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"  // conductance matrix (LHS)
#include "NumIntegral_NT_op_N_dV.h"    // fluid volume source
#include "NumIntegral_NT_op_N_dS.h"    // boundary integral
#include "NumIntegral_NT_lhsop_N_dV.h" // capacitance matrix LHS
#include "NumIntegral_NT_op_N_dV.h"    // capacitance matrix (lumped) RHS

#include "VelocityAndVolumeFlux.h"     // post-processing of Darcy velocity

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

#include "VTU_Interface.h"


using namespace std;

namespace csmp {

void RegionMonitor_Example::Specifications()
  {
     SetTitle( "Demonstration of RegionMonitor class" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "SKM" );
     AddDescription( "transient pressure diffusion around a well at a rate that depends of heterogeneous structure");
     AddDescription( "source in: RegionMonitor_Example.cpp" );
     AddRequirement( "any ANSYS model with multiple unique Region objects" );
  }


void RegionMonitor_Example::Run()
 {
  constexpr uint32_t dim{2};
  const string       model_name("Frewens3"); // TODO: choose an ANSYS model available to you!
  ANSYS_Model2D      model( model_name.c_str(), "TransientPressure_Example-variables.txt");


  // 1. Properties from file and hydraulic conductivity from permeability via Interrelation
  // --------------------------------------------------------------------------------------
  const bool region_specifications{false},      // regionname from parameter range
             default_property_values{true},     // default property values
             region_property_values{true},      // regional property values
             box_boundary_conditions{true},     // boundary conditions for box-shaped model
             region_property_conditions{true},  // regional property conditions
             boundary_conditions{false};        // boundary conditions for arbitrary-shaped model
  
  ComputationalSettings settings;
             
  InputDataManager<dim>  model_configuration;
  model_configuration.ConfigureFromFile( model, model.Name(),
                                         region_specifications,       // regionname from parameter range
                                         default_property_values,     // default property values
                                         region_property_values,      // regional property values
                                         box_boundary_conditions,     // boundary conditions for box-shaped model
                                         region_property_conditions,  // regional property conditions
                                         boundary_conditions,         // boundary conditions for arbitrary-shaped model
                                         settings );
  
  const double fluid_viscosity(1.0e-03);
  ConstantFactor<dim,divides>  conductivity( model.Database(),
                                           "conductivity", "permeability",
                                            fluid_viscosity );
  model.Apply( conductivity );
  printRangeOfVariable( model, "conductivity" );



  // 2. Build a transient fluid pressure algorithm using Backward-Euler time-stepping
  // ------------------------------------------------------------------------------------
  // ([C] + dt[K]){p}t+dt = [C]{p}t + dt {Q}t+dt
  PDE_Integrator<2,Element>  transient_pressure;
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver  samg_solver;
  transient_pressure.SetSolver( samg_solver );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
  transient_pressure.SetSolver( linear_solver );
#endif

  NumIntegral_dNT_op_dN_dV<dim> conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                            conductance.MultiplyWithTimeIncrement(true);

  NumIntegral_NT_lhsop_N_dV<dim> capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                             capacitance_lhs.LumpedFormulation(true);

  NumIntegral_NT_op_N_dV<dim> capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                          capacitance_rhs.LumpedFormulation(true);

  NumIntegral_NT_op_N_dV<dim> source( model.Database(), "fluid volume source",  "fluid pressure" );
                                          source.MultiplyWithTimeIncrement(true);
                                          source.AddAccumulateLater();
                                          source.LumpedFormulation(true);

  VelocityAndVolumeFlux<2U>  velocity( model,  "conductivity", "porosity", "fluid pressure", false );

  transient_pressure.Add( &conductance );
  transient_pressure.Add( &capacitance_lhs );
  transient_pressure.Add( &capacitance_rhs );
  transient_pressure.Add( &source );
  transient_pressure.AddPostProcess( &velocity );
  
  
  // TODO: exercise: turn line-element regions called "SHALE_CURVES" into a split boundaries
  // const bool retain_elmts_as_intervening_elements{false};
  // auto new_split_boundaries =  model.CreateSplitBoundaryFrom( "SHALE_CURVES", retain_elmts_as_intervening_elements );
  
  
  // 3.   Variables for transient loop
  // ---------------------------------
  list<string> material_props{"porosity","permeability","storativity"};
  list<string> runtime_variables{"fluid pressure","velocity","volume flux"};
  
  VTU_Interface<dim>  vtu_output( model );
  const size_t zero{ 0 };
  // initial state of the model // vtu_output.OutputDataToVTU( "Model", "hydraulic-conductivity", "conductivity", 0 );
  vtu_output.OutputDataToVTU( (string(model.Name()) + "_mtrl_properties").c_str(), material_props, model.Region( "Model" ), zero );
  
  // looping over unique (non-overlapping) model regions and outputting them individually
  model.RegionsOut();
  model.BoundariesOut();
  cout <<"\n"<<"RegionMonitor_Example: current model regions: "<< endl;
  for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); ++it ) {
    cout <<"\n\t"<< (*it).first;
    vtu_output.OutputDataToVTU( model.Name(), material_props, (*it).first, zero );
  }

  Standard_IO_Handler stdio;
  const double day(86400.);  // 1 day in seconds
  double       model_time{0.}, maxtime{20. * day}, time_increment{5.}, well_pressure{0.};
  size_t       timestep{1}, save_counter{0}, save_frequency{5};
  

  // 4. SETTING UP THE REGION MONITOR to record the integral popery 'volume flux' and range property 'pf'
  // ----------------------------------------------------------------------------------------------------
  string first_integral_property{"volume flux"};
  string first_range_property{"fluid pressure"};
  RegionMonitor<dim> monitor( model, first_integral_property, first_range_property );


  // TODO: exercise: using the ComputationalSettings object, change this loop so that output is created considering config file
  
  // 5. Transient loop: Compute fluid pressure during each time-step and output the results for each time step
  // ----------------------------------------------------------------------------------------------------------
  while ( model_time <= maxtime )
    {
      cout << "\n\nmain: COMPUTING TIMESTEP " << timestep << endl;

      // transient pressure, ramping up time increment step by step
      transient_pressure.TimeIncrement( time_increment );
      model.Apply( transient_pressure );

      // output variables screen
      printRangeOfVariable( model, "fluid pressure" );
      printRangeOfVariable( model, "velocity" );
      printRangeOfVariable( model, "pore velocity" );
      printRangeOfVariable( model, "volume flux" );
      
      // monitoring
      monitor.ScalarPropertyIntegrals( model, model_time );
      monitor.ScalarPropertyRanges( model, model_time );

      // pressure in the well
      well_pressure = model.Region("WELL").Average( "fluid pressure" );
      cout << "\nWell pressure: " << well_pressure << " Pa " << endl;
      // if well pressure is negative, prompt user to continue the simulation
      if ( well_pressure <= 0. )
        if ( (stdio.YesNo("Do you want to continue")) == false ) terminate();


      // output variables file every x steps (defined by user)
      if ( save_counter == save_frequency ) {
           // appending '-monitoring.txt' to model name for the monitoring output file
           monitor.Out( (string(model.Name())+"monitoring.txt").c_str() );
           vtu_output.OutputDataToVTU( (string(model.Name()) + "_simu_output").c_str(),
                                        runtime_variables, model.Region( "Model" ), model_time/day );
           // resetting save counter
           save_counter = 0;
        }

      // Preparing next Time Step
      model_time += time_increment;
      timestep++;
      save_counter++;
      
      // gradually increasing the time increment
      time_increment *= 1.2;

      cout << model_time << "\t" << well_pressure << endl;
      cout << "\nmain: ELAPSED TIME " << model_time / day << " days " << endl;
    }

    // 6.  Output the initial range of the variables
    // ---------------------------------------------
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );
    printRangeOfVariable( model, "pore velocity" );
    printRangeOfVariable( model, "volume flux" );
 
    cout <<"\nmain: That's it..."<< endl;
 
 } // end run


} // csmp
