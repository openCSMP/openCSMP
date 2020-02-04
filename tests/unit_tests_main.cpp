// Disable wingdi.h because it steps on our toes
#define NOGDI

#include <iostream>
#include <string>

#define CATCH_CONFIG_RUNNER
#define RUNNING_UNDER_CATCH // read: http://hiltmon.com/blog/2014/10/26/simple-c-plus-plus-testing-with-catch-in-xcode/ tutorial how to use with XCode

#include "Test.h" // includes catch.hpp
#include "TestSuite.h"
#include "CSMP_definitions.h"
#include "Exception.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Exception.h"
#endif

#include "ScalarVar_Test.h"
#include "VectorVar_Test.h"
#include "VectorVar_Test1.h"
#include "VectorVar_Test2.h"
#include "TensorVar_Test.h"
#include "TensorVar_Test1.h"
#include "TensorVar_Test2.h"
#include "ArrayVariable_Test.h"
#include "Variables_Test.h"
#include "PropertyStorageSpeed_Test.h"

#include "CommandLineParser_Test.h"
#include "GenericSingleton_Test.h"

#include "IsnanIsinf_Test.h"
#include "CubicSpline_Test.h"
#include "DenseMatrix_Test.h"
#include "Matrix_Test.h"
#include "SparseMatrix_Test.h"
#include "CompressedRowMatrix_Test.h"
#include "FibonacciHeap_Test.h"

#include "Point_Test.h"
#include "Node_Test.h"
#include "Element_Test.h"
#include "Face_Test.h"
#include "InterFace_Test.h"
#include "MeshManager_Test.h"

#include "FiniteElement_Test.h"
#include "FiniteElement_Test.h"
#include "IsoparametricQuadraticTetrahedron_Test.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricQuadraticTetrahedron.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricQuadraticTriangle.h"
#include "LinearCuboid_Test.h"

#include "FiniteVolumeStencil_Test.h"
#include "FiniteVolumePolicy_Test.h"
#include "FluxMismatch_Test.h"
#include "FV_Parameter_Test.h"
#include "IsoparametricLinearHexahedron_Test.h"

#include "PropertyDatabase_Test.h"
#include "Index_Test.h"
#include "PropertyData_Test.h"
#include "ColorPalette_Test.h"

#include "InputDataManager_Test.h"
#include "LocalVariableStorage_Test.h"
#include "Parameter_Test.h"

#include "ModelTopology_Test.h"
#include "StatisticalAnalyzer_Test.h"
#include "RegionMonitor_Test.h"
#include "Visitor_TestSuite.h"

#include "ModelSubDomain_Test.h"
#include "BoundaryInterface_Test.h"
#include "Boundary_Test.h"
#include "Region_Test.h"
#include "Box_Test.h"
#include "ANSYS_SplitBoundaryMatch_Test.h"
#include "ANSYS_Model2D_Test.h"
#include "ANSYS_Model3D_Test.h"

#include "BinaryFileInterface_Test.h"
#include "VTU_Interface_Test.h"
#include "FEM_Data_Test.h"
#include "VData_Test.h"
#include "VSet_Test.h"
#include "VSet_TestCase.h"
#include "PropertyHandle_Test.h"
#include "PropertyAtPointVisitor_Test.h"
#include "ModelComparator_Test.h"

#include "MathOperatorLHS_Test.h"
#include "MathOperatorRHS_Test.h"
#include "Operand_Test.h"
#include "PDE_Integrator_Test.h"

#include "FluxMismatch_Test.h"
#include "TwoPhaseModel_TestSuite.h"
#include "ExponentialTransferFunction_Test.h"

#include "CSMP_VariableBenchmarking_Test.h"
#include "GenericFiniteVolumeTransport_Test.h"

// new tests 2017 onwards
#include "LinearCuboid_Test.h"
#include "SplitBoundary_Test.h"
#include "SplitBoundaryInterface_Test.h"

using namespace std;
using namespace csmp;

/**  Unit Test Development (status: 1/3/2017, XCode OS X)

     interdependent1:   fails on volume-flux integral balance of face fluxes
     new_developments:  in progress, not ready yet
     
     TODO: BoxTest - in some cases perimeter nodes do not seem to be recognised correctly
     TODO: ModelSubDomain - add tests to cover complete functionality
 
     @section Additional Code Coverage Required

     - ModelSubDomain_Test - extend to cover extensive functionality
     - Box_Test - add unit-normal tests for prism and hexahedral mesh configurations
     - Model
     - MeshManager
     - FiniteElementManager
     - PDE_Integrator - extend to test assembly of solution matrix for systems
     - PropertyConstraints
     - RegionInterface
     - BoundaryInterface
     - ANSYS_Interface
     - InterFace_Test - create
     - IsNan_Test - refactor as it is built on false premises
     - VSet_Test - add functionality that checks storing and reading of properties
     
     - refactor   VTU_Interface_Test
     - refactor   CompressedRowMatrix - include elimitation etc.
     
     @section Failing Tests
     
     - IndexTracker - problems with variables created at runtime
     - Boundary_Test (TODO: separate BoundaryInterface functionality)
     - SplitBoundary_Test
     - InterFace_Test
     - Box_Test
          - TODO: boundary ChangePropertyStatus( INTERIOR) fails on boundary
 
     - TODO: delete non-unique regions with boundary names before saving the model to disk so that 
             they do not get stored and brought back when the model is rebuild
     
     - AnsysModel3D - when reconstructed from file volumetric elements seem to have surface element neighbors
     
     @section Comments
     - after Boundary construction, the parent regions are moved to non-unique, but are kept, is this what we want?
*/

TEST_CASE("CSMP unit tests", "[CSMP unit tests]")
 {
  const bool verbose(false);

  const bool test_fundamentals(true),
  test_interdependent1(true),
  test_interdependent2(true),
  test_composite(true),
  test_refactoring(true),
  test_new_developments(false);

  long fails_fundamentals(0),
  fails_interdependent1(0),
  fails_interdependent2(0),
  fails_composite(0),
  fails_new_developments(0),
  total_failures(0);
  
  try {
    cout <<"\nunit_test_main: running tests..."<< endl;
    if ( test_fundamentals ) {
      cout <<"\n1. underpinning functionality: running tests..."<< endl;
      TestSuite basic("CSMP-fundamental-unit test suite", &cout );
      
      // Auxiliaries
      basic.addTest( new IsnanIsinf_Test() );
      basic.addTest( new GenericSingleton_Test() );
// RUNS WITH EXCEPTIONS    basic.addTest( new CommandLineParser_Test() );
      
      // Data storage tests
      basic.addTest( new LocalVariableStorage_Test() );
      basic.addTest( new PropertyDatabase_Test());
      basic.addTest( new Index_Test());
      basic.addTest( new Parameter_Test());
      basic.addTest( new PropertyData_Test());
      
      // Model
      basic.addTest( new Node_Test() );
      basic.addTest( new Element_Test(verbose));
      basic.addTest( new Face_Test() );
      basic.addTest( new InterFace_Test() );
      basic.addTest( new ModelSubDomain_Test() );       
      basic.addTest( new Region_Test(false) );          
      basic.addTest( new BoundaryInterface_Test(false) ); 
      //mesh manager    
      basic.addTest(new MeshManager_Test(true));
      basic.addTest(new BoundaryInterface_Test(true));
      basic.addTest(new Boundary_Test());
      basic.addTest(new SplitBoundaryInterface_Test());
      basic.addTest(new SplitBoundary_Test());
      
      // Variable tests
      basic.addTest( new Point_Test());
      basic.addTest( new ScalarVariable_Test());
      basic.addTest( new VectorVariable_Test());
      basic.addTest( new VectorVariable_Test1());
      basic.addTest( new VectorVariable_Test2());
      basic.addTest( new TensorVariable_Test());
      basic.addTest( new TensorVariable_Test1());
      basic.addTest( new TensorVariable_Test2());
      basic.addTest( new ArrayVariable_Test());
      basic.addTest( new Variables_Test("FracBox")); //- requires refactoring of InterFace / SplitBoundary functionality
      
      // utilities tests
      basic.addTest( new Matrix_Test(verbose) );
      basic.addTest( new DenseMatrix_Test(verbose) );
      basic.addTest( new SparseMatrix_Test(verbose) );
      basic.addTest( new CompressedRowMatrix_Test() );
      basic.addTest( new CubicSpline_Test() );
      basic.addTest( new FibonacciHeap_Test() );
      
      basic.addTest( new FiniteVolumeStencil_Test(verbose));
      // also compares speed of mapping facet areas and normals versus computing them
      basic.addTest( new FiniteVolumePolicy_Test());
      basic.addTest( new FV_Parameter_Test());
      
      // interfaces / containers
      basic.addTest( new VData_Test() );
      basic.addTest( new FEM_Data_Test());
      basic.addTest( new PropertyData_Test() );
      basic.addTest( new VSet_Test() );
//      basic.addTest( new VSet_TestCase("UnitSquareFracs_irregular",true /* verbose */) );
      basic.addTest( new ColorPalette_Test() );
      
// FAILS TO RUN	  basic.addTest(new ANSYS_SplitBoundaryMatch_Test(true));

      // Running unit tests and reporting errors
      basic.run();
      fails_fundamentals = basic.report();
      basic.free();
      cerr << "\nunit_tests_main: 1. CSMP fundamentals: Total unit test failures: " << fails_fundamentals << endl;
    }
    
    // FINITE ELEMENTS + MATH OPERATORS
    if ( test_interdependent1 ) {
      cout <<"\n2. partially interdependent functionality: running tests..."<< endl;
      TestSuite interdependent1("CSMP-interdependent1-unit test suite", &cout );
      interdependent1.addTest( new FEM_Data_Test());
      // TODO: add the tests for all elements
      // finite elements
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(4), "IsoparametricLinearTetrahedron4P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(3,3), "IsoparametricLinearTriangle3D3IP.txt", verbose ) ); // 3D case 3 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,3), "IsoparametricLinearTriangle3IP.txt", verbose ) );   // 2D case 3 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,4), "IsoparametricLinearTriangle4IP.txt", verbose ) );   // 2D case 4 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt", verbose ) );  // 3D case 3 integration point
      // non-standard element tests
      interdependent1.addTest( new IsoparametricQuadraticTetrahedron_Test(verbose) ); // FAIL - flux balance on constant velocity projected on sides
      // volume conservation of distorted hexahedra - fails for certain deformation modes, highlighting limitations of this elements
      interdependent1.addTest(new IsoparametricLinearHexahedron_Test(verbose));
      // straight sided analytically integrated elements
      interdependent1.addTest( new LinearCuboid_Test(false) );  // TODO: add automatic test of interpolation function derivative matrix
      // math operators etc.
      interdependent1.addTest( new Operand_Test() );
      interdependent1.addTest( new MathOperatorLHS_Test());
      interdependent1.addTest( new MathOperatorRHS_Test());
      // TODO: add test of assembly of matrix for systems, elimination of boundary conditions etc.
      interdependent1.addTest( new PDE_Integrator_Test() );
      interdependent1.addTest( new LinearCuboid_Test(false) );
      
      // running unit tests and reporting errors
      interdependent1.run();
      fails_interdependent1 = interdependent1.report();
      interdependent1.free();
      cerr << "\nunit_tests_main: 2. CSMP interdependent-functionality1: Total unit test failures: " << fails_interdependent1 << endl;
    }
    
    if ( test_interdependent2 ) {
      cout <<"\n3. Model-related interdependent functionality: running tests..."<< endl;
      TestSuite interdependent2("CSMP-interdependent2-unit test suite", &cout );
      interdependent2.addTest( new ModelTopology_Test() );
      // model
      interdependent2.addTest( new Box_Test() );                  // XCode OK (SKM) but does not test hexahedral or prism element meshes
      interdependent2.addTest( new ANSYS_Model2D_Test() );        // XCode OK (SKM)
      interdependent2.addTest( new InputDataManager_Test());      // XCode OK (SKM)
      interdependent2.addTest( new ANSYS_Model3D_Test() );        // XCode OK (SKM)
      interdependent2.addTest( new PropertyHandle_Test() );       // XCode OK (SKM)
      // interfaces
      interdependent2.addTest( new BinaryFileInterface_Test() );  // XCode OK (SKM)
      interdependent2.addTest( new VTU_Interface_Test() );
      interdependent2.addTest( new StatisticalAnalyzer_Test() );
      // running unit tests and reporting errors
      interdependent2.run();
      fails_interdependent2 = interdependent2.report();
      interdependent2.free();
      cerr << "\nunit_tests_main: 3. CSMP Model-related, interdependent-functionality2: Total unit test failures: " << fails_interdependent2 << endl;
    }
    
    if ( test_composite ) {
      cout <<"\n4. Composite-dependent functionality: running tests..."<< endl;
      TestSuite composite("CSMP-dependent-unit test suite", &cout );
      // misc
      composite.addTest( new RegionMonitor_Test() );
      composite.addTest( new Variables_Test("FracBox") );

      // constitutive relationships TODO: create a separate test section for this
      
      /// Property data search tests
      Visitor_TestSuite visitorTests( composite );
      visitorTests.run();
      composite.addTest( new PropertyAtPointVisitor_Test(verbose) ); // PASS
      
      /// Two phase flow tests
      TwoPhaseModel_TestSuite  twoPhaseModelTests( composite );
      twoPhaseModelTests.run();
      composite.addTest( new ExponentialTransferFunction_Test() );
      
      // running unit tests and reporting errors
      composite.run();
      fails_composite = composite.report();
      composite.free();
      cerr << "\nunit_tests_main: 4. CSMP-dependent-functionality: Total unit test failures: " << fails_composite << endl;
    }
    
    // tests related to code that is currently being refactored
    if ( test_refactoring ) {
      cout <<"\n5. Refactored and new code functionality: running tests..."<< endl;
      TestSuite refactored("CSMP-refactored code unit-test suite", &cout );
      
      // TODO: review and get these tests to run (in this sequence)
      // update composite.addTest( new ModelComparator_Test() ); // crashes on PropertyData
      // basic.addTest( new VariableBenchmarking_Test() ); - needs redesign, tests tensor with random numbers
      // basic.addTest( new PropertyStorageSpeed_Test( &cout )); // needs redesign, calls Eigenvectors on random numbers

      refactored.run();
      long nFail = refactored.report();
      refactored.free();
      cerr << "\nunit_tests_main: 5. CSMP refactored and new functionality: Total unit test failures: " << nFail << endl;
    }
    
    // tests related to code that is currently being refactored
    if ( test_new_developments ) {
      cout <<"\n5. Refactored and new code functionality: running tests..."<< endl;
      TestSuite new_developments("new tests of the CSMP base library", &cout );
      
      // EVERYTHING THAT PERTAINS TO REFACTORED TRANSPORT SCHEME
      // TODO: broken   composite.addTest( new FluxMismatch_Test() );
      
      new_developments.addTest( new GenericFiniteVolumeTransport_Test() );
      
      // running unit tests and reporting errors
      new_developments.run();
      fails_new_developments = new_developments.report();
      new_developments.free();
      cerr << "\nunit_tests_main: 6. New functionality: Total unit test failures: " << fails_new_developments << endl;
    }
    
    cerr << "\nunit_tests_main: Total unit test failures: ";
    total_failures = fails_fundamentals + fails_interdependent1 + fails_interdependent2 + fails_composite + fails_new_developments;
    cerr << total_failures << endl;
    if ( fails_fundamentals > 0 )     cerr <<"\nfundamental functionality tests failed.";
    if ( fails_interdependent1 > 0 )  cerr <<"\ninterpedendent (basic) functionality tests failed.";
    if ( fails_interdependent2 > 0 )  cerr <<"\ninterpedendent (advanced) functionality tests failed.";
    if ( fails_composite > 0 )        cerr <<"\ncomposite functionality tests failed.";
    if ( fails_new_developments > 0 ) cerr <<"\nnew development tests failed.";
    cerr << endl << endl;
    
  } // Exception handling (warnings etc. are caught at a much lower level)
  catch( bad_alloc& ba ) {
    cerr <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( bad_cast& ba ) {
    cerr <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( bad_exception& ba ) {
    cerr <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( bad_typeid& ba ) {
    cerr <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( ios_base::failure& ba ) {
    cerr <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
    system("pause");
  }
  // standard logic errors
  catch( domain_error& ba ) {
    cerr <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( invalid_argument& ba ) {
    cerr <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( length_error& ba ) {
    cerr <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( out_of_range& ba ) {
    cerr <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
    system("pause");
  }
  // runtime errors
  catch( overflow_error& ba ) {
    cerr <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( range_error& ba ) {
    cerr <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( underflow_error& ba ) {
    cerr <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
    system("pause");
  }
  catch( Exception& ba ) {
#ifdef __GNUC__
    // this is a fix for gcc name demangling.
    const std::type_info  &ti = typeid(ba);
    int status;
    char* realname = abi::__cxa_demangle(ti.name(), 0, 0, &status );
    cout<<"\nException: Exception raised by: "<<realname<<endl;
#else
    cout<<"\nException: Exception raised by: "<< typeid(ba).name() << endl;
#endif
    cout <<"\nDiagnostics:"<< endl;
    ba.Out();
  }
#ifdef CSMP_WITH_SAMG_SOLVER
  catch( SAMG_Exception& ba ) {
    cout <<"\nSAMG_Exception: "<< ba.what() << endl;
    system("pause");
  }
#endif
  
  std::exit(total_failures);

} // end TEST_CASE



// -----------------------------------------
// RUNNING THE SUITES OF TESTS THROUGH CATCH
// -----------------------------------------
int main(int argc, char* argv[])
{
  Catch::Session session;

  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
    return returnCode;


  // Run Catch tests.
  session.run();
} // end main

