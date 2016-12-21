#include <iostream>
#include <string>

#include "CSMP_definitions.h"
#include "Exception.h"

#include "Test.h"
#include "TestSuite.h"

// read: http://hiltmon.com/blog/2014/10/26/simple-c-plus-plus-testing-with-catch-in-xcode/  as a tutorial how to use with XCode
#include "catch159.h"

#include "ScalarVar_Test.h"
#include "VectorVar_Test.h"
#include "VectorVar_Test1.h"
#include "VectorVar_Test2.h"
#include "TensorVar_Test.h"
#include "TensorVar_Test1.h"
#include "TensorVar_Test2.h"

#include "CommandLineParser_Test.h"
#include "GenericSingleton_Test.h"

#include "IsnanIsinf_Test.h"
#include "ErrorFunction_Test.h"
#include "CubicSpline_Test.h"
#include "DenseMatrix_Test.h"
#include "Matrix_Test.h"
#include "SparseMatrix_Test.h"
#include "CompressedRowMatrix_Test.h"

#include "Point_Test.h"
#include "Node_Test.h"
#include "Element_Test.h"
#include "Face_Test.h"
//#include "Edge_Test.h"				// jc: not exist

#include "FiniteElement_Test.h"
#include "IsoparametricQuadraticTetrahedron_Test.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricQuadraticTetrahedron.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricQuadraticTriangle.h"

#include "FiniteVolumeStencil_Test.h"
#include "FiniteVolumeTraits_Test.h"
#include "FluxMismatch_Test.h"
#include "FV_Parameter_Test.h"
#include "IsoparametricLinearHexahedron_Test.h"

#include "PropertyDatabase_Test.h"
#include "Index_Test.h"
#include "PropertyData_Test.hpp"

#include "InputDataManager_Test.h"
#include "LocalVariableStorage_Test.h"
#include "Parameter_Test.h"

#include "ModelTopology_Test.h"
#include "StatisticalAnalyzerTest.h"
#include "RegionMonitorTest.h"
#include "Visitor_TestSuite.h"

#include "Boundary_Test.h"
#include "Region_Test.h"
#include "Box_Test.h"
#include "SplitBoundary_Test.h"
#include "ANSYS_Model2D_Test.h"
#include "ANSYS_Model3D_Test.h"

#include "BinaryFileInterface_Test.h"
#include "VTU_Interface_Test.h"
#include "FEM_Data_Test.h"
#include "VData_Test.h"
#include "VSet_Test.h"
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

using namespace std;
using namespace csmp;

/**  SKM 16/12/2016 - status report (XCode Mac)

     fundamentals:      24 real precision issues in FiniteVolumeTraits, VSet test needs to be rewritten for variable storage in PropertyData container
     
     interdependent1:   fails on volume-flux integral balance of face fluxes
     
     interdependent2:   fails when re-reading the 3D Model in PropertyData for a tensor variable placed on the sector integration point
 
     composite:
     
     @section Missing Unit Tests
     
     Triage needed to generate order:
 
     - ModelTopology
     - Model
     - MeshManager
     - FiniteElementManager
     - ModelSubDomain
     - PDE_Integrator
     - MathOperatorLHS
     - MathOperatorRHS
     - PropertyConstraints
     - RegionInterface
     - BoundaryInterface
     - Point
     - ANSYS_Interface
     
     - refactor   Box_Test
     - refactor   VTU_Interface_Test
     - refactored CompressedRowMatrix
     
     
     @section Failing Tests
     
     - FiniteVolumeTraits_Test (TODO: urgent fixes needed)
     - Boundary_Test (TODO: separate BoundaryInterface functionality)
     - SplitBoundary_Test
     - InterFace_Test
     -
*/
int main()
{
  long nFail(0);
  const bool test_fundamentals(false),
             test_interdependent1(false),
             test_interdependent2(true),
             test_composite(false),
             test_refactoring(false);
  try {
        cout <<"\nunit_test_main: running tests..."<< endl;
        if ( test_fundamentals ) {
              cout <<"\n1. underpinning functionality: running tests..."<< endl;
              TestSuite basic("CSMP-fundamental-unit test suite", &cout );
          
              // Auxiliaries
              basic.addTest( new GenericSingleton_Test() );
              basic.addTest( new CommandLineParser_Test() );
              basic.addTest( new IsnanIsinf_Test() );

              // Data storage tests
              basic.addTest( new LocalVariableStorage_Test() );
              basic.addTest( new PropertyDatabase_Test());
              basic.addTest( new Index_Test());
              basic.addTest( new Parameter_Test());

              // Model
              basic.addTest( new Node_Test() );
              basic.addTest( new Element_Test());
              basic.addTest( new Face_Test() );

              // Variable tests
              basic.addTest( new Point_Test());
              basic.addTest( new ScalarVariable_Test());
              basic.addTest( new VectorVariable_Test());
              basic.addTest( new VectorVariable_Test1());
              basic.addTest( new VectorVariable_Test2());
              basic.addTest( new TensorVariable_Test());
              basic.addTest( new TensorVariable_Test1());
              basic.addTest( new TensorVariable_Test2());

              // Math utilities tests
              basic.addTest( new Matrix_Test() );
              basic.addTest( new DenseMatrix_Test() );
              basic.addTest( new SparseMatrix_Test() );
              basic.addTest( new CompressedRowMatrix_Test() );
              basic.addTest( new CubicSpline_Test() );
              basic.addTest( new ErrorFunction_Test() );

              basic.addTest( new FiniteVolumeStencil_Test());
              basic.addTest( new FiniteVolumeTraits_Test()); //also compares the speed between mapping the facet areas and normals and computing them.
              basic.addTest( new FV_Parameter_Test());
          
              // interfaces / containers
// add ModelTopology_Test
              basic.addTest( new VData_Test() );
              basic.addTest( new FEM_Data_Test());
              basic.addTest( new PropertyData_Test() );
          
              // Running unit tests and reporting errors
              basic.run();
              nFail = basic.report();
              basic.free();
              cerr << "\nunit_tests_main: 1. CSMP fundamentals: Total unit test failures: " << nFail << endl;
          }
    
        if ( test_interdependent1 ) {
              cout <<"\n2. partially interdependent functionality: running tests..."<< endl;
              TestSuite interdependent1("CSMP-interdependent1-unit test suite", &cout );
              interdependent1.addTest( new FEM_Data_Test());
              // finite elemenents
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt" ) );
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt" ) );
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(4), "IsoparametricLinearTetrahedron4P.txt" ) );
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(3,3), "IsoparametricLinearTriangle3D3IP.txt" ) ); // 3D case 3 integration points
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,3), "IsoparametricLinearTriangle3IP.txt" ) );   // 2D case 3 integration points
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,4), "IsoparametricLinearTriangle4IP.txt" ) );   // 2D case 4 integration points
              interdependent1.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt" ) );  // 3D case 3 integration point
              // non-standard element tests
              interdependent1.addTest( new IsoparametricQuadraticTetrahedron_Test(false) ); // FAIL - flux balance on constant velocity projected on sides
              // volume conservation of distorted hexahedra - fails for certain deformation modes, highlighting limitations of this elements
              interdependent1.addTest(new IsoparametricLinearHexahedron_Test());
              // math operators etc.
              interdependent1.addTest( new Operand_Test() );
              interdependent1.addTest( new MathOperatorLHS_Test());
              interdependent1.addTest( new MathOperatorRHS_Test());
              interdependent1.addTest( new PDE_Integrator_Test() );
              // running unit tests and reporting errors
              interdependent1.run();
              nFail = interdependent1.report();
              interdependent1.free();
              cerr << "\nunit_tests_main: 2. CSMP interdependent-functionality1: Total unit test failures: " << nFail << endl;
          }

        if ( test_interdependent2 ) {
              cout <<"\n3. Model-related interdependent functionality: running tests..."<< endl;
              TestSuite interdependent2("CSMP-interdependent2-unit test suite", &cout );
              // model
              interdependent2.addTest( new Box_Test() );
              interdependent2.addTest( new ANSYS_Model2D_Test() );
              interdependent2.addTest( new InputDataManager_Test());
              interdependent2.addTest( new ANSYS_Model3D_Test() );           // FAIL - PropertyData (tensor, sector-ip) when model is re-imported from binary file
              interdependent2.addTest( new Region_Test() );                  // SKM OK
              interdependent2.addTest( new ModelTopology_Test() );
              interdependent2.addTest( new PropertyHandle_Test() );          // SKM OK
              // interfaces
              interdependent2.addTest( new VTU_Interface_Test() );
              interdependent2.addTest( new StatisticalAnalyzerTest() );
              // running unit tests and reporting errors
              interdependent2.run();
              nFail = interdependent2.report();
              interdependent2.free();
              cerr << "\nunit_tests_main: 3. CSMP Model-related, interdependent-functionality2: Total unit test failures: " << nFail << endl;
          }
    
        if ( test_composite ) {
              cout <<"\n4. Composite-dependent functionality: running tests..."<< endl;
              TestSuite composite("CSMP-dependent-unit test suite", &cout );
              // misc
              // composite.addTest( new PropertyAtPointVisitor_Test() ); // PASS
              composite.addTest( new BinaryFileInterface_Test() );  // FAIL on assert
              // composite.addTest( new FluxMismatch_Test() );
              composite.addTest( new RegionMonitorTest() );         // FAILS - tolerance issues?
              composite.addTest( new ModelComparator_Test() );      // crashes on PropertyData
              // constitutive relationships
              composite.addTest( new ExponentialTransferFunction_Test() );

              /// Property data search tests
              Visitor_TestSuite visitorTests( composite );
              visitorTests.run();

              /// Two phase flow tests
              TwoPhaseModel_TestSuite twoPhaseModelTests( composite );
              twoPhaseModelTests.run();

             // running unit tests and reporting errors
              composite.run();
              nFail = composite.report();
              composite.free();
              cerr << "\nunit_tests_main: 4. CSMP-dependent-functionality: Total unit test failures: " << nFail << endl;
          }

        // tests related to code that is currently being refactored
        if ( test_refactoring ) {
              cout <<"\n5. Refactored and new code functionality: running tests..."<< endl;
              TestSuite refactored("CSMP-refactored code unit-test suite", &cout );
              refactored.addTest( new Boundary_Test() );
              // composite.addTest( new SplitBoundary_Test() );
             // running unit tests and reporting errors


              // refactored.addTest( new PropertyData_Test() ); // retested: OK - includes vectors, tensors, arrays
              refactored.run();
              nFail = refactored.report();
              refactored.free();
              cerr << "\nunit_tests_main: 5. CSMP refactored and new functionality: Total unit test failures: " << nFail << endl;
          }

        cerr << "\nunit_tests_main: Total unit test failures: " << nFail << endl;

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
    cerr <<"\nException: Exception raised: "<< ba.What() << endl;
    cerr <<"\nDiagnostics:"<< endl;
    ba.Out();
    system("pause");
  }

  // tell operating system that no error occurred (by returning 0 as opposed to
  // some error number)
  return nFail;

} // end main

