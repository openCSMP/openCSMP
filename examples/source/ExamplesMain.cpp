
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

using namespace csmp;
using namespace std;


#ifdef COMMAND_LINE_EXAMPLES
int main( int argc, char* argv[] )
#else
int main()
#endif
 {
    try {
     // if PETSc is chosen as the default linear solver
#ifdef CSMP_WITH_PETSC_SOLVER
      PetscInitializeNoArguments();
#endif
      // instantiating suite
      ExampleSuite examplesSuite( "Open-CSMP++ EXAMPLES", &cout );

      // write your own
      examplesSuite.RegisterExample( new Experimental_Example() );

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

      // TUTORIALS
            examplesSuite.RegisterExample( new Tutorial1_Example() );
            examplesSuite.RegisterExample( new Tutorial2_Example() );
            examplesSuite.RegisterExample( new Tutorial3_Example() );
            examplesSuite.RegisterExample( new Tutorial4_Example() );
            

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


