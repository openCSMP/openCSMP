// Disable wingdi.h because it steps on our toes
#define NOGDI

#define CATCH_CONFIG_MAIN
#include "catch.hpp"


#include "ScalarVar_Test.h"
#include "VectorVar_Test.h"
#include "VectorVar_Test1.h"
#include "VectorVar_Test2.h"
#include "TensorVar_Test.h"
#include "TensorVar_Test1.h"
#include "TensorVar_Test2.h"
#include "ArrayVariable_Test.h"

#include "CommandLineParser_Test.h"
#include "GenericSingleton_Test.h"

#include "IsnanIsinf_Test.h"
#include "CubicSpline_Test.h"
#include "DenseMatrix_Test.h"
#include "Matrix_Test.h"
#include "SparseMatrix_Test.h"
#include "CompressedRowMatrix_Test.h"
#include "ColorPalette_Test.h"

#include "Point_Test.h"
#include "Node_Test.h"
#include "Element_Test.h"
#include "Face_Test.h"

#include "FiniteElement_Test.h"
#include "IsoparametricQuadraticTetrahedron_Test.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricQuadraticTetrahedron.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricQuadraticTriangle.h"

#include "FiniteVolumeStencil_Test.h"
#include "FiniteVolumePolicy_Test.h"
#include "FiniteElement_Test.h"
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
#include "StatisticalAnalyzer_Test.h"
#include "RegionMonitor_Test.h"
#include "Visitor_TestSuite.h"
#include "Box_Test.h"

#include "ModelSubDomain_Test.h"
#include "BoundaryInterface_Test.h"
#include "Boundary_Test.h"            // sm: needs work
#include "Region_Test.h"
#include "Box_Test.h"
#include "SplitBoundary_Test.h"       // sm: needs work
#include "ANSYS_Model2D_Test.h"
#include "ANSYS_Model3D_Test.h"

#include "BinaryFileInterface_Test.h"
#include "VTU_Interface_Test.h"
#include "FEM_Data_Test.h"
#include "VData_Test.h"
#include "VSet_Test.h"
#include "PropertyHandle_Test.h"
#include "PropertyAtPointVisitor_Test.h"
//#include "ModelComparator_Test.h"     // sm: needs work

#include "MathOperatorLHS_Test.h"
#include "MathOperatorRHS_Test.h"
#include "Operand_Test.h"
#include "PDE_Integrator_Test.h"

//#include "FluxMismatch_Test.h"        // sm: needs work
#include "TwoPhaseModel_TestSuite.h"
#include "ExponentialTransferFunction_Test.h"

using namespace std;
using namespace csmp;

/**  
    SKM revised to include all passing tests 8/3/2017

    @note for code coverage etc., see status report in unit_tests_main.cpp
*/

#define SIMPLE_TEST_SECTION(n)  SECTION(#n) { n##_Test test; test.run(); }
#define SIMPLE_TEST_SECTION_N(n,num)  SECTION(#n) { n##_Test##num test; test.run(); }

TEST_CASE("Auxiliary tests", "[Auxiliaries]") {
    SIMPLE_TEST_SECTION(GenericSingleton)
    SIMPLE_TEST_SECTION(CommandLineParser)
    SIMPLE_TEST_SECTION(IsnanIsinf)
    SIMPLE_TEST_SECTION(ColorPalette)
}

TEST_CASE("Variable database tests", "[Variables]") {
    SIMPLE_TEST_SECTION(LocalVariableStorage)
    SIMPLE_TEST_SECTION(PropertyDatabase)
    SIMPLE_TEST_SECTION(Index)
    SIMPLE_TEST_SECTION(Parameter)
}

TEST_CASE("Variable tests", "[Variable]") {
    SIMPLE_TEST_SECTION(Point);
    SIMPLE_TEST_SECTION(ScalarVariable);
    SIMPLE_TEST_SECTION(VectorVariable);
    SIMPLE_TEST_SECTION_N(VectorVariable,1);
    SIMPLE_TEST_SECTION_N(VectorVariable,2);
    SIMPLE_TEST_SECTION(TensorVariable);
    SIMPLE_TEST_SECTION_N(TensorVariable,1);
    SIMPLE_TEST_SECTION_N(TensorVariable,2);
    SIMPLE_TEST_SECTION(ArrayVariable);
}

TEST_CASE("Math utilities tests", "[MathUtils]") {
    SIMPLE_TEST_SECTION(Matrix)
    SIMPLE_TEST_SECTION(DenseMatrix)
    SIMPLE_TEST_SECTION(SparseMatrix)
    SIMPLE_TEST_SECTION(CompressedRowMatrix)
    SIMPLE_TEST_SECTION(CubicSpline)
}

TEST_CASE("Data containers tests", "[DataContainers]") {
    SIMPLE_TEST_SECTION(VData)
    SIMPLE_TEST_SECTION(VSet)
    SIMPLE_TEST_SECTION(PropertyData)
    SIMPLE_TEST_SECTION(FEM_Data)
    SIMPLE_TEST_SECTION(ModelTopology)
}

TEST_CASE("Model and model functionality tests", "[Model]") {
    SIMPLE_TEST_SECTION(Node)
    SIMPLE_TEST_SECTION(Element)
    SIMPLE_TEST_SECTION(Box)
    SIMPLE_TEST_SECTION(Face)
    SIMPLE_TEST_SECTION(ModelSubDomain)
    SIMPLE_TEST_SECTION(Region)
//    SIMPLE_TEST_SECTION(Boundary)
//    SIMPLE_TEST_SECTION(SplitBoundary)
    SIMPLE_TEST_SECTION(ANSYS_Model2D)
    SIMPLE_TEST_SECTION(ANSYS_Model3D)
    SIMPLE_TEST_SECTION(Box)
    SIMPLE_TEST_SECTION(InputDataManager)
    SIMPLE_TEST_SECTION(PropertyHandle)
 }
// composite.addTest( new ModelComparator_Test() );   crashes on FEM_Data (needs refactoring)

TEST_CASE("Finite element - finite volume integration tests", "[IntegralMethods]") {
    SIMPLE_TEST_SECTION(FiniteVolumePolicy)
    SIMPLE_TEST_SECTION(FiniteVolumeStencil)
    SIMPLE_TEST_SECTION(Operand)
    SIMPLE_TEST_SECTION(MathOperatorLHS)
    SIMPLE_TEST_SECTION(MathOperatorRHS)
    SIMPLE_TEST_SECTION(PDE_Integrator)
}
// composite.addTest( new FluxMismatch_Test() ); // FAILS


TEST_CASE("Interfaces with other software tests", "[Interfaces]") {
    SIMPLE_TEST_SECTION(VTU_Interface)
    SIMPLE_TEST_SECTION(BinaryFileInterface)
 }

TEST_CASE("Analysis of results and integral properties tests", "[Analysis]") {
    SIMPLE_TEST_SECTION(StatisticalAnalyzer)
    SIMPLE_TEST_SECTION(RegionMonitor)
 }


#if 0
    if ( test_interdependent1 ) {
          const bool verbose(false);
          cout <<"\n2. Finite-element functionality: running tests..."<< endl;
          TestSuite interdependent1("Finite element test suite", &cout );
          interdependent1.addTest( new FEM_Data_Test());
          // finite elements
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(4), "IsoparametricLinearTetrahedron4P.txt", verbose ) );
          // 3D case 3 integration points
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(3,3), "IsoparametricLinearTriangle3D3IP.txt", verbose ) );
          // 2D case 3 integration points
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,3), "IsoparametricLinearTriangle3IP.txt", verbose ) );
          // 2D case 4 integration points
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,4), "IsoparametricLinearTriangle4IP.txt", verbose ) );
          // 3D case 3 integration point
          interdependent1.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt", verbose ) );
          // non-standard element tests
          interdependent1.addTest( new IsoparametricQuadraticTetrahedron_Test(verbose) ); // FAIL - flux balance on constant velocity projected on sides
          // volume conservation of distorted hexahedra - fails for certain deformation modes, highlighting limitations of this elements
          interdependent1.addTest(new IsoparametricLinearHexahedron_Test(verbose));
          // running unit tests and reporting errors
          interdependent1.run();
          nFail = interdependent1.report();
          interdependent1.free();
          cerr << "\nunit_tests_main: 2. CSMP interdependent-functionality1: Total unit test failures: " << nFail << endl;
      }

    if ( test_interdependent2 ) {
          const bool verbose(false);
          cout <<"\n3. Model-related interdependent functionality: running tests..."<< endl;
          TestSuite interdependent2("CSMP-interdependent2-unit test suite", &cout );
          // model
          interdependent2.addTest( new ModelSubDomain_Test() );
          interdependent2.addTest( new Region_Test() );
          // interdependent2.addTest( new Boundary_Test() );      // not yet
          // interdependent2.addTest( new SplitBoundary_Test() ); // not yet
          interdependent2.addTest( new ANSYS_Model2D_Test() );
          interdependent2.addTest( new ANSYS_Model3D_Test() );
          interdependent2.addTest( new Box_Test() );
          interdependent2.addTest( new InputDataManager_Test() );
          interdependent2.addTest( new PropertyHandle_Test() );
          // interfaces
          interdependent2.addTest( new VTU_Interface_Test() );
          interdependent2.addTest( new StatisticalAnalyzerTest() );
          // running unit tests and reporting errors
          interdependent2.run();
          nFail = interdependent2.report();
          interdependent2.free();
          cerr << "\nunit_tests_main: 3. CSMP Model-related, interdependent-functionality2: Total unit test failures: " << nFail << endl;
      }

    // visitors
    if ( test_composite ) {
          cout <<"\n4. Composite-dependent functionality: running tests..."<< endl;
          TestSuite composite("CSMP-dependent-unit test suite", &cout );
          // misc

          /// Property data search tests
          Visitor_TestSuite visitorTests( composite );
          visitorTests.run();

          composite.addTest( new PropertyAtPointVisitor_Test() );

          /// Two phase flow tests
          TwoPhaseModel_TestSuite twoPhaseModelTests( composite );
          twoPhaseModelTests.run();
          composite.addTest( new ExponentialTransferFunction_Test() );

         // running unit tests and reporting errors
          composite.run();
          nFail = composite.report();
          composite.free();
          cerr << "\nunit_tests_main: 4. CSMP-dependent-functionality: Total unit test failures: " << nFail << endl;
      }

    cerr << "\nunit_tests_main: Total unit test failures: " << nFail << endl;

#endif

