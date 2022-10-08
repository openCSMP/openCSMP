#include "DESAdvectionDiffusion3D_Example.h"

#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "ComputationalSettings.h"
#include "InputDataManager.h"
#include "ANSYS_Model3D.h"
#include "ModelTopology.h"
#include "VTK_Interface.h"

// integration od PDEs and post-processing
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// interrelations
#include "ConstantFactor.h"

// DES finite-volume transport scheme
#include "DESAdvectionDiffusion.h"

#include "MeshDiagnostics.h"

using namespace std;

namespace csmp {

void DESAdvectionDiffusion3D_Example::Specifications()
  {
    SetTitle( "3D Passive advection of concentration using DES (discrete event simulation)");
    SetDifficulty( 3 );
    SetCategory( "Simulation of Physical Processes" );
    AddAuthor( "Qi Shao" );
    AddDescription( "3D passive tracer advection by both discrete event simulation (DES) and time-driven simulation (TDS)." );
    AddDescription( "Their outputs are written to VTK files and and their efficiency are compared." );
    AddRequirement( "source files: 'DESAdvectionDiffusion3D_Example.cpp' and '*.h'" );
    AddRequirement( "fracs4 (CSMP native binary files), frac4-configuration.txt, DES_variables.txt");
  }


/** 
    3D passive tracer advection via CSMP's DES transport method 
    combining finite elements (for pressure) with finite volumes (for advection of concentration profile)

    Use model 'fracs4' (CSMP native binary files, frac4-configuration.txt, DES_variables.txt) as input file suites.
*/
void DESAdvectionDiffusion3D_Example::Run()
{
 // ------------------------------------------------------------
 // 1. Load CSMP native format model
 // ------------------------------------------------------------
  string model_name;
  cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'fracs4':"<<endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "fracs4";

  //find the name of current example source file
  string file_name = GetExampleFileName(__FILE__);
  string variable_file = "DES_variables.txt";
  string config_file = "fracs4(DES_tracer)";
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);
  //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model<3U>  model3D(model_name, variable_file);
  printModelDimensions( model3D, true );


 // ------------------------------------------------------------
 // 2. checking the mesh quality
 // ------------------------------------------------------------
  MeshDiagnostics<3U>  mesh_check;
  double             vol_min, vol_max;
  mesh_check.ElementVolumeRange( model3D, vol_min, vol_max );
  cout <<"\nmain: element volume range: "<< vol_min <<" to "<< vol_max << endl;
  cout <<"\tlarge element-volume range (>10^3) can pose problems for the linear solver.2\n";


 // ------------------------------------------------------------
 // 3. configuring the model
 // ------------------------------------------------------------
  InputDataManager<3U>  model_configuration;
//  model_configuration.ConfigureFromFile( model3D, model_name.c_str(),
//                                         false, true, true, true, false );

  // the boolean variables determine which blocks in the input file shall be read
    ComputationalSettings  run_settings;
    model_configuration.ConfigureFromFile( model3D,
                                           config_file.c_str(),
                                           false, 
                                           true,   // 2) default prop.values
                                           true,   // 3) group prop.values
                                           true,   // 4) essential box-boundary conditions
                                           true,   // 5) essential flags
                                           true,   // 6) boundary conditions
                                           run_settings );

    model3D.InputPropertyValue( "nodal concentration source",  makeScalar(PLAIN,0.0) );      // no sources/sinks for solute (units kg m-3 s-1)
    
    Standard_IO_Handler  stdio;
    printRangeOfVariable( model3D, stdio, "permeability" );

    // visualizing the input permeability and boundary conditions
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model3D, "permeability", "permeability", 0 );
    vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );


 // -----------------------------------------------------------------------
 // 4. hydraulic conductivity and other interrelations
 // -----------------------------------------------------------------------
  const double  fluid_viscosity(1.0e-03);
  ConstantFactor<3U,divides>  conductivity( model3D.Database(),
                                           "conductivity", "permeability",
                                            fluid_viscosity );
  model3D.Apply( conductivity );
  printRangeOfVariable( model3D, "conductivity" );

  vtk_output.OutputDataToVTK( model3D, "conductivity", "conductivity", 0 );


 // -----------------------------------------------------------------------
 // 5. computing a steady-state fluid pressure distribution in the model
 // -----------------------------------------------------------------------
  SteadyStateDiffusor<3U,Region> steady_state_pressure( model3D,
                                                        "conductivity", "fluid pressure",
                                                        "fluid volume source" );
  // postprocessing of pressure gradients and flow velocities
  VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model3D, "conductivity", "porosity", "fluid pressure" );
  steady_state_pressure.AddPostProcess( &postpro0 );

  // the calculation of fluid pressure
  steady_state_pressure.ComputeSteadyState( model3D.Region("Model") );

  // results: the pore velocity is the Darcy velocity divided by the porosity
  printRangeOfVariable( model3D, stdio, "fluid pressure" );
  printRangeOfVariable( model3D, stdio, "velocity" );
  printRangeOfVariable( model3D, stdio, "pore velocity" );

  vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 1 );
  vtk_output.OutputDataToVTK( model3D, "velocity",       "velocity",       1 );


 // -----------------------------------------------------------------------
 // 6. Demonstration of the generic transport algorithm
 //   (here you can compare different schemes with one another and
 //    overstep CFL to see how this adds numerical diffusion to the solution)
 // -----------------------------------------------------------------------
  printModelDimensions( model3D, true );

  double time_interval;
  cerr <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds,suggested value: 6.0e6)";
  cin >> time_interval;
  
  double Courant_multiplier;
  cerr <<"\nEnter CFL multiplier (suggested value: 0.3 ~ 0.7)";
  cin >> Courant_multiplier;

  double PEP_factor = 1.0;  
  size_t n_threads = 1;
  VTK_Interface<3U>	vtkOut;
  double model_time(0.);  
    
  bool  DES = stdio.YesNo("Do you want to solve the advection equation with DES (y=DES, n=TDS)");  
  
  if (DES)
  {  
  // -----------------------------------------------------------------------
  // 7. Transient loop using DES simulation
  // -----------------------------------------------------------------------
  cerr <<"\n\nmain: Starting DES simulation: "<<endl; 
  clock_t	T_begin= clock();
  DESAdvectionDiffusion<3U>	DES_advector( model3D, "Model", Courant_multiplier, PEP_factor, false);
  
  for ( int i=0; i<20; ++i ) {
      DES_advector.AdvectVariable_DES( model_time+time_interval/20., n_threads);
      model_time += time_interval/20.;
      vtkOut.OutputDataToVTK( model3D, "DES_concentration", "concentration", i+1 );
      vtkOut.OutputDataToVTK( model3D, "DES_update_count", "update count", i+1 );
      vtkOut.OutputDataToVTK( model3D, "DES_variation_rate_count", "rate count", i+1 );
      vtkOut.OutputDataToVTK( model3D, "DES_schedule_count", "schedule count", i+1 );
      vtkOut.OutputDataToVTK( model3D, "DES_synchronization_count", "synchronize count", i+1 );
   }
   cerr <<"\nmain: Finished DES simulation, using time (sec): "<< (clock() - T_begin)/double(CLOCKS_PER_SEC) << endl;
  }

  else
  {
  // -----------------------------------------------------------------------
  // 8. Transient loop using TDS simulation
  // -----------------------------------------------------------------------
  cerr <<"\n\nmain: Starting TDS simulation: "<<endl; 
  clock_t	T_begin= clock();
  DESAdvectionDiffusion<3U>	TDS_advector( model3D, "Model", Courant_multiplier, PEP_factor, false); 

  for ( int i=0; i<20; ++i ) {
      TDS_advector.AdvectVariable_TDS( time_interval/20., n_threads);
      model_time += time_interval/20.;
      vtkOut.OutputDataToVTK( model3D, "TDS_concentration", "concentration", i+1 );
   }
   cerr <<"\nmain: Finished TDS simulation, using time (sec): "<< (clock() - T_begin)/double(CLOCKS_PER_SEC) << endl;   
   }
   

   cerr <<"\nmain: That's it."<< endl;

   fs::current_path("../../example_inputs/");

} // end Run



} // csmp


