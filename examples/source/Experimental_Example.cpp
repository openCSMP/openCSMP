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
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(example_name, model_name, variable_file, config_file);

  Model<dim>  model( model_name, variable_file );
  printModelDimensions(model);
  model.RegionsOut();
  model.BoundariesOut();


  // 0. Properties from file and hydraulic conductivity from permeability via Interrelation
  // --------------------------------------------------------------------------------------
  const bool region_specifications{false},      // region name from parameter range
             default_property_values{true},     // default property values
             region_property_values{true},      // regional property values
             box_boundary_conditions{true},     // boundary conditions for box-shaped model
             region_property_conditions{true},  // regional property conditions
             boundary_conditions{true};         // boundary conditions for arbitrary-shaped model
  
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
  
  // computing hydraulic conductivity
  {
  constexpr double fluid_viscosity(1.0e-03);
      ConstantFactor<dim,divides>  conductivity( model.Database(), "conductivity", "permeability", fluid_viscosity );
      model.Apply( conductivity );
      printRangeOfVariable( model, "conductivity" );
  }
  
  // and the total systems compressibility,ct (beta_r * (1-phi) + beta_f * phi) =storativity m3/(m3 * Pa))
  {
      constexpr double beta_w = 4.6e-10; // 4.6 × 10⁻¹⁰ Pa⁻¹ (water)
      
      // Bulk compressibility Kb = 1/K [Pa^-1]
      // Sources: Biot (1941), Detournay & Cheng (1993),
      //          Schanz & Nguyen-Tuan (2011) for Opalinus Clay,
      //          Mavko et al. "Rock Physics Handbook" (2009)

      constexpr double Kb_Dolomite      = 1.4e-11;  // Pa^-1  K ~ 70 GPa
                                                     // Dolomite / dolostone
                                                     // range: 1e-11 to 2e-11

      constexpr double Kb_IronSst       = 5.0e-11;  // Pa^-1  K ~ 20 GPa
                                                     // Iron-cemented sandstone
                                                     // Well cemented, lower porosity
                                                     // range: 2e-11 to 1e-10

      constexpr double Kb_Opalinus_clay = 4.6e-10;  // Pa^-1  K ~ 2.2 GPa
                                                     // Soft clay shale
                                                     // Schanz & Nguyen-Tuan (2011)
                                                     // range: 3e-10 to 8e-10

      constexpr double Kb_RhaetSst      = 1.0e-10;  // Pa^-1  K ~ 10 GPa
                                                     // Rhaetian sandstone
                                                     // Poorly cemented, high porosity
                                                     // range: 5e-11 to 2e-10
                                                     
      // = 0.15 * 4.5e-10 + 4.6e-10 = 5.275e-10 Pa^-1
      model.Region("AMALTHEA_CLAY").InputPropertyValue( "storativity", makeScalar(FIELD_DATA,5.275e-10) );
      model.Region("OPALINUS_CLAY").InputPropertyValue( "storativity", makeScalar(FIELD_DATA,5.275e-10) );
      // = 0.10 * 4.5e-10 + 5.0e-11 = 9.5e-11 Pa^-1
      model.Region("IRON_SANDSTONE").InputPropertyValue( "storativity", makeScalar(FIELD_DATA,9.5e-11) );
      // = 0.20 * 4.5e-10 + 1.0e-10 = 1.9e-10 Pa^-1
      model.Region("RHAET_SANDSTONE").InputPropertyValue( "storativity", makeScalar(FIELD_DATA,1.9e-10) );
      // = 0.05 * 4.5e-10 + 1.4e-11 = 3.65e-11 Pa^-1
      model.Region("REEF_DOLOMITE").InputPropertyValue( "storativity", makeScalar(FIELD_DATA,3.65e-11) );
      printRangeOfVariable( model, "storativity" );
  }

  // 1. Initial conditions (air pressure: varying 960–1040 hPa over the year at sea level.
  // -------------------------------------------------------------------------------------
  constexpr double initial_barometric_pressure{ 104000. }; // a high
  // ignoring the gravity-induced increase of pressure with depth we start with a "fluid pressure"
  // equal to the airpressure to explore how long it takes before atmospheric pressure fluctuations are felt inside the mountain
  model.InputPropertyValue( "fluid pressure", makeScalar(ANY,initial_barometric_pressure) );
  // now we drop the atmospheric pressure, simulating a low, and compute has fast this change diffuses into the model
  model.Boundary("IRREGULAR").InputPropertyValue( "fluid pressure", makeScalar(DIRICH,96000.) );


  // 2. Build a transient fluid pressure algorithm using Backward-Euler time-stepping
  // --------------------------------------------------------------------------------
  // ([S]/dt + [C]){p}t+dt = [S]/dt {p}t + {Q}t+dt
  PDE_Integrator<dim,Element>  transient_pressure;
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver  samg_solver;
  transient_pressure.SetSolver( samg_solver );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
  transient_pressure.SetSolver( linear_solver );
#endif

  NumIntegral_dNT_lhsop_dN_dV<dim> conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
 
  NumIntegral_NT_lhsop_N_dV<dim> capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                                  capacitance_lhs.LumpedFormulation(true); // NB: volume lumping (works only for simplex elements)
                                                  capacitance_lhs.MultiplyWithTimeIncrement(true);

  NumIntegral_NT_rhsop_N_dV<dim> capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                                  capacitance_rhs.LumpedFormulation(true);
                                                  capacitance_rhs.MultiplyWithTimeIncrement(true);

  NumIntegral_NT_rhsop_N_dV<dim> source( model.Database(), "fluid volume source",  "fluid pressure" );
                                         source.AddAccumulateLater();
                                         source.LumpedFormulation(true);

  transient_pressure.Add( &conductance );
  transient_pressure.Add( &capacitance_lhs );
  transient_pressure.Add( &capacitance_rhs );
  transient_pressure.Add( &source );
  
  // 3.   Variables for transient loop
  // ---------------------------------
  list<string> material_props{"porosity","permeability","storativity"};
  list<string> runtime_variables{"fluid pressure"};
  
  VTU_Interface<dim>  vtu_output( model );
  // initial state of the model // vtu_output.OutputDataToVTU( "Model", "hydraulic-conductivity", "conductivity", 0 );
  vtu_output.OutputDataToVTU( (string(model.Name()) + "_mtrl_properties").c_str(), material_props, model.Region( "Model" ), 0 );

  Standard_IO_Handler stdio;
  const double day(86400.);  // 1 day in seconds, max time 1 million years.
  double       model_time{0.}, maxtime{1e6 * 365. * day}, time_increment{5.}, well_pressure{0.};
  size_t       timestep{1}, save_counter{0}, save_frequency{5};

  cout << "\nExperimentalExample: enter after how many steps you would like to save the results (1 = every step): " << endl;
  cin  >> save_frequency;

  // setting up pressure monitoring for all model regions
  string first_integral_property{"volume flux"};
  string first_range_property{"fluid pressure"};
  RegionMonitor<dim> monitor( model, first_integral_property, first_range_property );


  // TODO: change this output loop so that output is created only at those times enlisted in configuration file
  
  // 4. Transient loop: Compute fluid pressure during each time-step and output the results for each time step
  // ----------------------------------------------------------------------------------------------------------
  while ( model_time <= maxtime )
    {
      cout << "\n\nmain: COMPUTING TIMESTEP " << timestep << endl;

      // screening of computed input variables
      printRangeOfVariable( model, "fluid pressure" );

      // transient pressure computation
      transient_pressure.TimeIncrement( 1. / time_increment );
      model.Apply( transient_pressure );

      // screening output variables
      printRangeOfVariable( model, "fluid pressure" );
      
      monitor.ScalarPropertyIntegrals( model, model_time );
      monitor.ScalarPropertyRanges( model, model_time );

      // output variables file every x steps (defined by user)
      if ( save_counter == save_frequency ) {
           vtu_output.OutputDataToVTU( (string(model.Name()) + "_output").c_str(),
                                       runtime_variables, model.Region( "Model" ),
                                       static_cast<long>(model_time/day) );
           save_counter = 0;
        }

      // Preparing next Time Step
      model_time += time_increment;
      timestep++;
      save_counter++;
      
      // gradually increasing the time increment
      time_increment *= 1.2;

      // writing the fluid pressure to the txt file
      cout << model_time << "\t" << well_pressure << endl;
      cout << "\nmain: ELAPSED TIME " << model_time / day << " days " << endl;
    }

    monitor.Out( model.Name() );
    cout <<"\nmain: That's it..."<< endl;

    // 9.  Output the initial range of the variables
    // ---------------------------------------------
    printRangeOfVariable( model, "fluid pressure" );
    
    cout <<"\n end of program"<< endl;
 
 } // end run


} // csmp
