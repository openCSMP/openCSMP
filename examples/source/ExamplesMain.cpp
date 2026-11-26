// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Examples.h"
#include "ExampleSuite.h"
#include "Exception.h"
#include "ErrorHandler.h"

// GUI
#ifdef EXAMPLE_SUITE_WITH_GUI
#define COMMAND_LINE_EXAMPLES
#include <QApplication>
#include "ExampleSuiteMainWindow.h"
#endif

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

using namespace csmp;
using namespace std;


namespace {

/// replaces a leading '~' by the user's home directory (the shell does not do this for environment variables)
filesystem::path expandHome( const string& directory )
 {
    if ( !directory.empty() && directory[0] == '~' && ( directory.size() == 1U || directory[1] == '/' ) ) {
         const char* home = getenv( "HOME" );
         if ( home != nullptr ) return filesystem::path( string(home) + directory.substr(1) );
      }
    return filesystem::path( directory );
 }

/**
    Makes the 'example_inputs/' directory of the examples workspace the current working directory.

    The examples read their inputs relative to this directory and write their results to
    'example_outputs/<example name>/' next to it. The workspace is, in this order:

      1. the directory given with the command-line option  --workspace <dir>  (or --workspace=<dir>),
      2. the directory in the environment variable CSMP_EXAMPLES_DIR,
      3. the parent of the current working directory, if the latter is named 'example_inputs'
         (this is the setup in which an IDE has been given this working directory),
      4. the directory chosen when openCSMP was configured with CMake (CSMP_EXAMPLES_WORKSPACE;
         default: ~/openCSMP/examples).

    @return false if the workspace does not contain an 'example_inputs/' directory
*/
bool selectExamplesWorkspace( int argc, char* argv[] )
 {
    string workspace, chosen_by;

    for ( int i=1; i<argc; ++i ) {
         const string argument( argv[i] );
         if ( argument == "--workspace" && i+1 < argc ) {
              workspace = argv[i+1];
              chosen_by = "command-line option --workspace";
            }
         else if ( argument.rfind( "--workspace=", 0 ) == 0 ) {
              workspace = argument.substr( string("--workspace=").size() );
              chosen_by = "command-line option --workspace";
            }
       }
    if ( workspace.empty() ) {
         const char* from_environment = getenv( "CSMP_EXAMPLES_DIR" );
         if ( from_environment != nullptr && string(from_environment).size() > 0U ) {
              workspace = from_environment;
              chosen_by = "environment variable CSMP_EXAMPLES_DIR";
           }
      }
    if ( workspace.empty() && filesystem::current_path().filename() == "example_inputs" ) {
         workspace = filesystem::current_path().parent_path().string();
         chosen_by = "current working directory";
      }
#ifdef CSMP_EXAMPLES_WORKSPACE_DEFAULT
    if ( workspace.empty() ) {
         workspace = CSMP_EXAMPLES_WORKSPACE_DEFAULT;
         chosen_by = "default set when openCSMP was configured";
      }
#endif

    error_code error;
    const filesystem::path workspace_path = filesystem::absolute( expandHome( workspace ), error );
    const filesystem::path inputs         = workspace_path / "example_inputs";

    if ( workspace.empty() || !filesystem::is_directory( inputs, error ) ) {
         cout <<"\nopenCSMP examples: no example inputs found in the examples workspace\n"
              <<"   workspace: "<< ( workspace.empty() ? string("(none)") : workspace_path.string() );
         if ( !chosen_by.empty() ) cout <<"   ("<< chosen_by <<")";
         cout <<"\n   expected:  "<< inputs.string() <<"\n\n"
              <<"Create the workspace from your openCSMP build directory with\n\n"
              <<"   cmake --build . --target setup_examples\n\n"
              <<"or copy the directory 'examples/example_inputs' of the repository to a place outside\n"
              <<"of the repository and start the examples with\n\n"
              <<"   openCSMP-examples --workspace <directory that contains example_inputs>\n"<< endl;
         return false;
      }

    filesystem::current_path( inputs );
    cout <<"\nopenCSMP examples workspace ("<< chosen_by <<")\n"
         <<"   inputs are read from:   "<< inputs.string() <<"\n"
         <<"   results are written to: "<< ( workspace_path / "example_outputs" / "<example name>" ).string() <<"\n"<< endl;
    return true;
 }

} // anonymous namespace


int main( int argc, char* argv[] )
 {
    // the examples run in a workspace outside of the repository, see selectExamplesWorkspace()
    if ( !selectExamplesWorkspace( argc, argv ) ) return 1;

    try {
     // if PETSc is chosen as the default linear solver
#ifdef CSMP_WITH_PETSC_SOLVER
      PetscInitializeNoArguments();
#endif
      // instantiating suite
      ExampleSuite examplesSuite( "CSMP++ EXAMPLES", &cout );

      // CSMP FOR BEGINNERS
            examplesSuite.RegisterExample( new CplusplusInCSMP_Example() );
            examplesSuite.RegisterExample( new BuildingACSMPModel_Example );
            examplesSuite.RegisterExample( new CSMP_PhysicalVariables_Example() );
            examplesSuite.RegisterExample( new SteadyStateHeatConduction_Example() );
            examplesSuite.RegisterExample( new InputDataManager_Example() );
            examplesSuite.RegisterExample( new AtmosphericPressureDiffusion_Example() );
            examplesSuite.RegisterExample( new MeshConvergence_Example() );
      /*
            examplesSuite.RegisterExample( new CSMP_Variables_Example() );
            examplesSuite.RegisterExample( new CSMP_Model_Example() );
            examplesSuite.RegisterExample( new ModelSubdomainTypes_Example() );
      */

      // this is where your examples go (they will be ordered alphabetically when listed)
  
      // C++ EXAMPLES
            examplesSuite.RegisterExample( new CuriouslyRecurringTemplate_Example() );
            examplesSuite.RegisterExample( new PolicyBased_Example() );
            examplesSuite.RegisterExample( new TDDUnitTest_Example() );
  
      // APPLIED MATH / NUMERICAL METHODS
            examplesSuite.RegisterExample( new DenseMatrix_Example() );
            examplesSuite.RegisterExample( new ErrorMetric_Example() );
            examplesSuite.RegisterExample( new LinearSolver_Example() );
            examplesSuite.RegisterExample( new StokesDiscrepancyMeasure_Example() );
            examplesSuite.RegisterExample( new StokesDiscrepancyMeasureQuadratic_Example() );
            examplesSuite.RegisterExample( new TimeSteppingApproaches_Example() );
            examplesSuite.RegisterExample( new PoreRadiusAnd_Pc_Example() );
#if defined(CSMP_WITH_PETSC_SOLVER)
            examplesSuite.RegisterExample( new PETSc_Example() );
#endif

      // SOFTWARE INTERFACES
            // interfaces with other software tools
            examplesSuite.RegisterExample( new Interfaces_Example() ); // example to generate csmp native format models from various interfaces
            examplesSuite.RegisterExample( new ModelANSYS_Example() );
            examplesSuite.RegisterExample( new Triangulator_Example() );
            examplesSuite.RegisterExample( new QuadrilateratorToCSMPbinary_Example() );
            examplesSuite.RegisterExample( new EclipseMeshInterface_Example() );
            examplesSuite.RegisterExample( new Intrepid_Example() );
            examplesSuite.RegisterExample( new RhinoMesh_Example() );
            examplesSuite.RegisterExample( new SKUA_Example() );
            
      // SOFTWARE FUNCTIONALITY
            examplesSuite.RegisterExample( new VariablesBasic_Example() );
            examplesSuite.RegisterExample( new Variables_Example() ); //deprecate (this is essentially the unit test
            examplesSuite.RegisterExample( new TemplatizedIndex_Example() );
            examplesSuite.RegisterExample( new VariableManagement_Example() );
            examplesSuite.RegisterExample( new Averaging_Example() );
            examplesSuite.RegisterExample( new PermeabilityTensor_Example() );
            examplesSuite.RegisterExample( new Visitor_Example() );
            examplesSuite.RegisterExample( new Region_Example() );
            examplesSuite.RegisterExample( new RegionProperties_Example() );
            examplesSuite.RegisterExample( new RegionMonitor_Example() );
            examplesSuite.RegisterExample( new StatisticalAnalyzer_Example() );

      // SIMULATION OF PHYSICAL PROCESSES
            examplesSuite.RegisterExample( new PressureDiffusion_Example() );
            examplesSuite.RegisterExample( new TransientPressure_Example() );
            examplesSuite.RegisterExample( new ParallelPlateFracture_Example() );
            examplesSuite.RegisterExample( new EffectiveStressDilatation2D_Example() ); // TODO: incomplete! - refactor
            examplesSuite.RegisterExample( new StreamFunction_Example() );
            examplesSuite.RegisterExample( new TopographyDrivenFlow_Example() );
            examplesSuite.RegisterExample( new QuadraticPressure_parallelPlatePermeability_Example() );
            examplesSuite.RegisterExample( new PassiveAdvectionOfTracer_Example() );
            examplesSuite.RegisterExample( new DESAdvectionDiffusion2D_Example() );
            examplesSuite.RegisterExample( new DESAdvectionDiffusion3D_Example() );
            examplesSuite.RegisterExample( new DES2PhaseSlightlyCompressibleFlow_Example() );
            examplesSuite.RegisterExample( new DES2PhaseFlowWithSplitBoundary_Example() );
            examplesSuite.RegisterExample( new LinearElasticity_Example() );
            examplesSuite.RegisterExample( new SlopeMechanics_Example() );
            examplesSuite.RegisterExample( new PoroElasticity_Example() );
            examplesSuite.RegisterExample( new TemperatureDensityPressure_Example() );
            examplesSuite.RegisterExample( new ThermalConvectionETHZ_Example() );
            examplesSuite.RegisterExample( new Geothermal_Example() );
            examplesSuite.RegisterExample( new UG4_ProMeshOutput_Example() );
            examplesSuite.RegisterExample( new StaggeredGridStokesSolver_Example() );
            examplesSuite.RegisterExample( new PressureDependentBlackOilProperties_Example() );
            examplesSuite.RegisterExample( new Tractions_Example() );

      // CVFEM
            examplesSuite.RegisterExample( new CVFEM_fault_well_lithium_example() );
            examplesSuite.RegisterExample( new CVFEM_topo_magma_air_example() );
            examplesSuite.RegisterExample( new CVFEM_geothermal_cooling_magma_chamber_example<2U>() );

      // TUTORIALS
            examplesSuite.RegisterExample( new Tutorial1_Example() );
            examplesSuite.RegisterExample( new Tutorial2_Example() );
            examplesSuite.RegisterExample( new Tutorial3_Example() );
            examplesSuite.RegisterExample( new Tutorial4_Example() );
            
      // WRITE YOUR OWN EXAMPLES
            examplesSuite.RegisterExample( new Experimental_Example() );

            // calling the suite either as stdIO or GUI
#ifndef EXAMPLE_SUITE_WITH_GUI
            examplesSuite.Run();
#else
            QApplication app( argc, argv );
            PL::ExampleSuiteMainWindow exampleSuiteMainWindow;
            exampleSuiteMainWindow.show();
            exampleSuiteMainWindow.SetSuite( &examplesSuite );
            return app.exec();
#endif


    } //--- EXCEPTION HANDLING---
    catch( bad_alloc& ba ) {
         cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
      }
    catch( bad_cast& ba ) {
         cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
      }
    catch( bad_exception& ba ) {
         cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
      }
    catch( bad_typeid& ba ) {
         cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
      }
    catch( ios_base::failure& ba ) {
         cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
      }
    // standard logic errors
    catch( domain_error& ba ) {
         cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
      }
    catch( invalid_argument& ba ) {
         cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
      }
    catch( length_error& ba ) {
         cout <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
       }
    catch( out_of_range& ba ) {
         cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
      }
    // runtime errors
    catch( overflow_error& ba ) {
         cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( range_error& ba ) {
         cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( underflow_error& ba ) {
         cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( Exception& ba ) {
         cout <<"\nException: Exception raised: "<< ba.What() << endl;
         cout <<"\nDiagnostics:"<< endl;
         ba.Out();
      }

#ifndef NDEBUG
    cout <<"\nHit return to end program."<< endl;
    getchar();
#endif

#ifdef CSMP_WITH_PETSC_SOLVER
    PetscFinalize();
#endif

    return 0;

} // end main


