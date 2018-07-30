// Disable wingdi.h because it steps on our toes
#define NOGDI

#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

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

#include "CommandLineParser_Test.h"
#include "GenericSingleton_Test.h"

#include "IsnanIsinf_Test.h"
#include "CubicSpline_Test.h"
#include "DenseMatrix_Test.h"
#include "Matrix_Test.h"
#include "SparseMatrix_Test.h"
#include "CompressedRowMatrix_Test.h"
#include "ColorPalette_Test.h"
#include "FibonacciHeap_Test.h"

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
#define SIMPLE_TEST_SECTION_NONSTANDARD(n)  SECTION(#n) { n##_Test test(false); test.run(); }
#define SIMPLE_TEST_SECTION_FE(n,p,d) SECTION(#n) {n##_Test test(p, d, false); test.run(); }

TEST_CASE("Auxiliary tests", "[Auxiliaries]") {
    SIMPLE_TEST_SECTION(GenericSingleton)
    SIMPLE_TEST_SECTION(CommandLineParser)
    SIMPLE_TEST_SECTION(IsnanIsinf)
    SIMPLE_TEST_SECTION(ColorPalette)
    SIMPLE_TEST_SECTION(FibonacciHeap)
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
    SIMPLE_TEST_SECTION_NONSTANDARD(SparseMatrix)
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
    SIMPLE_TEST_SECTION_NONSTANDARD(Box) // Unknown exception
    SIMPLE_TEST_SECTION(Face)
    SIMPLE_TEST_SECTION(ModelSubDomain)
    SIMPLE_TEST_SECTION_NONSTANDARD(Region) // Unknown exception
    SIMPLE_TEST_SECTION(Boundary)
    SIMPLE_TEST_SECTION(SplitBoundary)
    SIMPLE_TEST_SECTION(ANSYS_Model2D)
    SIMPLE_TEST_SECTION(ANSYS_Model3D)
    SIMPLE_TEST_SECTION_NONSTANDARD(Box) // Unknown exception
    SIMPLE_TEST_SECTION(InputDataManager) // Unknown exception
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

TEST_CASE("Finite-element functionality", "[FiniteElement]"){
	SIMPLE_TEST_SECTION(FEM_Data)
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricLinearTetrahedron(1),"IsoparametricLinearTetrahedron1P.txt")
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricLinearTetrahedron(4),"IsoparametricLinearTetrahedron4P.txt")
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricLinearTriangle(3,3),"IsoparametricLinearTriangle3D3IP.txt")
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricLinearTriangle(2,3),"IsoparametricLinearTriangle3IP.txt")
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricLinearTriangle(2,4),"IsoparametricLinearTriangle4IP.txt")
	SIMPLE_TEST_SECTION_FE(FiniteElement,new IsoparametricQuadraticTriangle(2),"IsoparametricQuadraticTriangle.txt")
	SIMPLE_TEST_SECTION_NONSTANDARD(IsoparametricQuadraticTetrahedron)
	SIMPLE_TEST_SECTION_NONSTANDARD(IsoparametricLinearHexahedron)
}

TEST_CASE("Model - related interdependent functionality", "[ModelFunctionality]") {
	SIMPLE_TEST_SECTION(ModelSubDomain)
	SIMPLE_TEST_SECTION(Region)
	SIMPLE_TEST_SECTION(Boundary)
	SIMPLE_TEST_SECTION(SplitBoundary)
	SIMPLE_TEST_SECTION(ANSYS_Model2D)
	SIMPLE_TEST_SECTION(ANSYS_Model3D)
	SIMPLE_TEST_SECTION_NONSTANDARD(Box)
	SIMPLE_TEST_SECTION(InputDataManager)
	SIMPLE_TEST_SECTION(PropertyHandle)
	SIMPLE_TEST_SECTION(VTU_Interface)
	SIMPLE_TEST_SECTION(StatisticalAnalyzer)
}

TEST_CASE("Composite-dependent functionality", "[CompositeFunctionality]") {
	SIMPLE_TEST_SECTION_NONSTANDARD(PropertyAtPointVisitor)
	//SIMPLE_TEST_SECTION_NONSTANDARD(ExponentialTransferFunction)
}
