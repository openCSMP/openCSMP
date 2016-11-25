#include <iostream>
#include <string>

#include "CSMP_definitions.h"
#include "Exception.h"

#include "Test.h"
#include "TestSuite.h"

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

//#include "PropertyDatabase_Test.h"	// jc: error LNK2001: unresolved external symbol
//#include "Index_Test.h"				// jc: error LNK2001: unresolved external symbol
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

int main()
{
  long nFail(0);
  try {

    cout <<"\nunit_test_main: running tests..."<< endl;
    TestSuite s("CSMP unit test suite (class-level functionality only", &cout );

    s.addTest( new PropertyHandle_Test() );          // SKM OK
    s.addTest( new CompressedRowMatrix_Test() );     // SKM OK
    s.addTest( new Box_Test() );                     // SKM OK
    //s.addTest( new BinaryFileInterface_Test() );
    //s.addTest( new InputDataManager_Test());
    s.addTest( new Region_Test() );                  // SKM OK

//    s.addTest( new Boundary_Test() );
    s.addTest( new ANSYS_Model2D_Test() );
    //s.addTest( new FluxMismatch_Test()); // TODO: Adriana's FV consistency check (needs to be turned into a test)

    /// High level utilities tests
    s.addTest( new GenericSingleton_Test() );
    s.addTest( new CommandLineParser_Test() );
    s.addTest( new IsnanIsinf_Test() );

    /// Data storage tests
    s.addTest( new LocalVariableStorage_Test() );
    //s.addTest( new PropertyDatabase_Test());		// jc: error LNK2001: unresolved external symbol
    //s.addTest( new Index_Test());					// jc: error LNK2001: unresolved external symbol
    s.addTest( new Parameter_Test());

    s.addTest( new VData_Test() );
    s.addTest( new VSet_Test() );

    /// Variable tests
    s.addTest( new ScalarVariable_Test());
    s.addTest( new VectorVariable_Test());
    s.addTest( new VectorVariable_Test1());
    s.addTest( new VectorVariable_Test2());
    s.addTest( new TensorVariable_Test());
    s.addTest( new TensorVariable_Test1());
    s.addTest( new TensorVariable_Test2());

    /// Math utilities tests
    s.addTest( new Point_Test());
    s.addTest( new Matrix_Test() );
    s.addTest( new DenseMatrix_Test() );
    s.addTest( new SparseMatrix_Test() );
    s.addTest( new CubicSpline_Test() );
    s.addTest( new ErrorFunction_Test() );
    s.addTest( new StatisticalAnalyzerTest() );

    /// Low level geometry tests
    s.addTest( new Node_Test() );
    s.addTest( new Element_Test());
    s.addTest( new Face_Test() );
    //s.addTest( new Edge_Test() );					// jc: not exist

    // High level geometry test's
    s.addTest( new ModelTopology_Test() );
    s.addTest( new ANSYS_Model3D_Test() );

    // Monitoring
    s.addTest( new RegionMonitorTest() );

    // Visualization tools
    s.addTest( new VTU_Interface_Test() );

    /// Finite Volume test's
    s.addTest( new FiniteVolumeStencil_Test());
    s.addTest( new FiniteVolumeTraits_Test());//also compares the speed between mapping the facet areas and normals and computing them.
    s.addTest( new FV_Parameter_Test());

    /// Finite Elements test's
    s.addTest( new FEM_Data_Test());
    s.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt" ) );
    s.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(4), "IsoparametricLinearTetrahedron4P.txt" ) );
    s.addTest( new FiniteElement_Test( new IsoparametricQuadraticTetrahedron(), "IsoparametricQuadraticTetrahedron.txt" ) );
    s.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(3,3), "IsoparametricLinearTriangle3D3IP.txt" ) ); // 3D case 3 integration points
    s.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,3), "IsoparametricLinearTriangle3IP.txt" ) );   // 2D case 3 integration points
    s.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,4), "IsoparametricLinearTriangle4IP.txt" ) );   // 2D case 4 integration points
    s.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt" ) );  // 3D case 3 integration point
    
    s.addTest( new IsoparametricQuadraticTetrahedron_Test(false) );
	  
    /// volume conservation of distorted hexahedra - fails for certain deformation modes, highlighting limitations of this elements
    s.addTest(new IsoparametricLinearHexahedron_Test());

    /// Property data search tests
    s.addTest( new PropertyAtPointVisitor_Test() );
    s.addTest( new ModelComparator_Test() );

    /// Visitor tests
    Visitor_TestSuite visitorTests( s );
    visitorTests.run();

    /// Two phase flow tests
    s.addTest( new ExponentialTransferFunction_Test() );

    TwoPhaseModel_TestSuite twoPhaseModelTests( s );
    twoPhaseModelTests.run();

    /// Algorithm tests
    s.addTest( new Operand_Test() );
    s.addTest( new MathOperatorLHS_Test());
    s.addTest( new MathOperatorRHS_Test());
    s.addTest( new PDE_Integrator_Test() );

    // Running unit tests and reporting errors
    s.run();
    nFail = s.report();
    s.free();
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

