#include "TransientPressure_Example.h"

// CSMP Files
#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "TRIANGLE_Interface.h"
#include "VSetConverter.h"

// visualization of results
#include "JPEG_RegionInterface.h"
#include "VTK_Interface.h"

// fluid pressure algorithm and velocity computation
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "VelocityAndVolumeFlux.h"

// interrelation
#include "ConstantFactor.h"

// utilities
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"
#include "ModelTime.h"

using namespace std;

namespace csmp{

void TransientPressure_Example::Specifications()
{
  SetTitle( "Transient pressure diffusion around a well (2D)" );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "source in: TransientPressure_Example.cpp" );
  AddDescription( "transient calculation constant rate draw down, 2D" );
  AddDescription( "source in: TransientPressure_Example.cpp" );
  AddRequirement( "file set well.1, configuration(example4.txt)" );
}




/** *****************************************************************************************

  // (4) 2D transient calculation constant rate draw down from a pumping well near a highly
  //     permeable fault. Model uses a 'Triangle' generated input mesh (file set well.1)
  //
  //     Use example_viewer.tcl to visualize vtk output
  //
  //     Vary the draw down rate and observe at which rate the fluid pressure becomes negative
  //     and the computation breaks down. Observe further how the presence of a highly-
  //     permeable fault zone changes the pressure diffusion
  // *****************************************************************************************

*/
void TransientPressure_Example::Run()
{
  double64& model_time( ModelTime::Instance().modelTime );
  model_time = 0.;

  // 1.0 Building quadratic triangular FE mesh
  // -----------------------------------------
  TRIANGLE_Interface  mesh_interface;
  VSet<2U>            mesh_container;
  VSetConverter<2U>   mesh_converter;

  char  file_name[200];
  cout<< "\n2D Transient flow around a pumping well (example 4) "<< endl;
  cout <<"\nmain: Enter name of 'Triangle' input file set: ";
  cin >> file_name;
  mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );
  mesh_converter.ConvertLinearToQuadraticTriangles( mesh_container );

  // 2.0 Building the Region named "model"
  // set boolean for isoparametric elements to true
  // -----------------------------------------------
  Model<2U>   model( mesh_container, "example4.txt", true );

  // 3.0 Input the initial Conditions
  // --------------------------------
  model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.0) );
  model.InputPropertyValue( "porosity",            makeScalar(PLAIN,0.2) );
  model.InputPropertyValue( "storativity",         makeScalar(PLAIN,2.0e-9) );
  model.InputPropertyValue( "fluid pressure",      makeScalar(PLAIN,1.0e+7) );

  // 4.0 Forming various unique Groups
  // ---------------------------------
  cout << "\nForming the region 'fault zone' " << endl;
  model.FormRegionFrom( "fault zone", "permeability", 1.0e-19, 1.0e-16, true );
  cout << "\nForming the region 'well' " << endl;
  model.FormRegionFrom( "well", "permeability", 1.0e-06, 1.0e-04, true );

  // 4.1 Changing the storativity and porosity in the region objects
  // ---------------------------------------------------------------
  // low permeability fault
  model.Region("fault zone").InputPropertyValue( "porosity", makeScalar(PLAIN,0.01) );
  model.Region("fault zone").InputPropertyValue( "storativity", makeScalar(PLAIN, 8.0e-10) );

  // well
  model.Region("well").InputPropertyValue( "porosity", makeScalar(PLAIN,1.0) );
  model.Region("well").InputPropertyValue( "storativity", makeScalar(PLAIN, 7.0e-10) );
  model.Region("well").InputPropertyValue( "permeability", makeScalar(PLAIN, 2.0e-11) );

  // 4.2 Prompt the user to change the fault zone to be highly-permeable
  // -------------------------------------------------------------------
  Standard_IO_Handler  stdio;
  if ( stdio.YesNo("Do you want the fault zone to be highly permeable") ) {
      model.Region("fault zone").InputPropertyValue( "porosity", makeScalar(PLAIN,0.3) );
      model.Region("fault zone").InputPropertyValue( "storativity", makeScalar(PLAIN,6.0e-10) );
      model.Region("fault zone").InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-7) );
   }


  // 5.0 Calculating hydraulic conductivity from permeability using Interrelation subclass
  // ------------------------------------------------------------------------------------
  const double64 fluid_viscosity(1.0e-03);
  ConstantFactor<2U,divides>  conductivity( model.Database(),
                                                    "conductivity", "permeability",
                                                     fluid_viscosity );
  model.Apply( conductivity );
  printRangeOfVariable( model, "conductivity" );


  // 6.0 Assign farfield Dirichlet boundary conditions for the fluid pressure
  // ------------------------------------------------------------------------
  model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,1.0e+7) );
  model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.0e+7) );


  // 7.0 Let the user specify a drawdown rate and change the source/sink term in the well area
  // -----------------------------------------------------------------------------------------
  ScalarVariable  pumping_rate;
  cout << "\nEnter the pumping rate (Units: m3 m-2 s-1, try 1.0e-05 to start with. Positve rate: injection, negative rate: extraction) " << endl;
  cin  >> pumping_rate();
  model.Region("well").InputPropertyValue( "fluid volume source", pumping_rate );


  // 8.0 Variables for transient loop
  // --------------------------------
  VTK_Interface<2U>  vtk_output;
  // only write the results of the fault zone to jpg files because they will be too large otherwise
  JPEG_RegionInterface  jpg_output( model, "fault zone" );

  const double64 day(86400.);  // 1 year in seconds
  double64       maxtime(20. * day), time_increment(0.5 * day), well_pressure;
  size_t         timestep(1), save_counter(1), save_frequency;
  cout << "\nEnter after how many steps you would like to save the results (1 = every step) " << endl;
  cin  >> save_frequency;

  // 9.0 Opening a file for writing the fluid pressure in the well to a txt file
  // ---------------------------------------------------------------------------
  char  outfile[200];
  strcpy( outfile, "well-pressure" );
  strcat( outfile, ".txt" );
  ofstream ofs;
  ofs.open( outfile, ios::out|ios::trunc );
  ofs << "Fluid pressure in the well with an extraction rate of " << pumping_rate << " m3 m-2 s-1" <<  endl;
  ofs << "Time [sec]\tPressure [Pa] " << endl;
  ofs << "0\t1.0e+07 " << endl; // initial condition




  // 10.0 Build a transient fluid pressure algorithm using Backward-Euler time-stepping
  // ----------------------------------------------------------------------------------
  // ([C] + dt[K]){p}t+dt = [C]{p}t + dt {Q}t+dt
  PDE_Integrator<2U,Region>  transient_pressure;
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver  samg_solver;
  transient_pressure.SetSolver(&samg_solver);
#else
  CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
  transient_pressure.SetSolver(&linear_solver);
#endif

  NumIntegral_dNT_op_dN_dV<2U,Element<2U> >  conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                         conductance.MultiplyWithTimeIncrement(true);

  NumIntegral_NT_lhsop_N_dV<2U,Element<2U> > capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                         capacitance_lhs.LumpedFormulation(true);

  NumIntegral_NT_op_N_dV<2U,Element<2U> >    capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                         capacitance_rhs.LumpedFormulation(true);

  NumIntegral_NT_op_N_dV<2U,Element<2U> >    source( model.Database(), "fluid volume source",  "fluid pressure" );
                                         source.MultiplyWithTimeIncrement(true);
                                         source.AddAccumulateLater();
                                         source.LumpedFormulation(true);

  VelocityAndVolumeFlux<2U,Element<2U> >    velocity( model,  "conductivity", "porosity", "fluid pressure", false );

  transient_pressure.Add( &conductance );
  transient_pressure.Add( &capacitance_lhs );
  transient_pressure.Add( &capacitance_rhs );
  transient_pressure.Add( &source );
  transient_pressure.AddPostProcess( &velocity );
  transient_pressure.TimeIncrement( time_increment );


  // 11. Transient loop: Compute fluid pressure during each time-step and output the results for each time step
  // ----------------------------------------------------------------------------------------------------------
  while ( model_time <= maxtime )
    {
      cout << "\n\nmain: COMPUTING TIMESTEP " << timestep << endl;

      // transient pressure
      model.Apply( transient_pressure );

      // output variables screen
      printRangeOfVariable( model, "fluid pressure" );
      printRangeOfVariable( model, "velocity" );
      printRangeOfVariable( model, "pore velocity" );
      printRangeOfVariable( model, "volume flux" );

      // pressure in the well
      well_pressure = model.Region("well").Average( "fluid pressure" );
      cout << "\nWell pressure: " << well_pressure << " Pa " << endl;

      // if well pressure is negative, prompt user to continue the simulation
      if ( well_pressure <= 0. )
        if ( (stdio.YesNo("Do you want to continue")) == false ) terminate();


      // output variables file every x steps (defined by user)
      if ( save_counter == save_frequency ) {
          vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure",    timestep );
          vtk_output.OutputDataToVTK( model, "velocity",       "velocity",          timestep);
          vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",       timestep );

          jpg_output.OutputRegionDataToJPG( model, "fluid-pressure", "fluid pressure", timestep );
          jpg_output.OutputRegionDataToJPG( model, "volume-flux",    "volume flux",    timestep );

          save_counter = 0;
        }

      // Preparing next Time Step
      model_time += time_increment;
      timestep++;
      save_counter++;

      // writing the fluid pressure to the txt file
      ofs << model_time << "\t" << well_pressure << endl;

      cout << "\nmain: ELAPSED TIME " << model_time / day << " days " << endl;

    }

  cout <<"\nmain: That's it..."<< endl;
} // Run()

} // csmp

