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
#include "MeshManager_Test.h"

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
#include "PropertyData_Test.h"

#include "InputDataManager_Test.h"
#include "LocalVariableStorage_Test.h"
#include "Parameter_Test.h"

#include "ModelTopology_Test.h"
#include "StatisticalAnalyzer_Test.h"
#include "RegionMonitor_Test.h"
#include "Visitor_TestSuite.h"
#include "Box_Test.h"
#include "MeshManager_Test.h"

#include "ModelSubDomain_Test.h"
#include "BoundaryInterface_Test.h"
#include "Boundary_Test.h"
#include "Region_Test.h"
#include "SplitBoundary_Test.h"
#include "ANSYS_Model2D_Test.h"
#include "ANSYS_Model3D_Test.h"

// #include "BinaryFileInterface_Test.h" // alt: not in repository
#include "VTU_Interface_Test.h"
#include "FEM_Data_Test.h"
#include "VData_Test.h"
#include "VSet_Test1.h"
#include "VSet_Test2.h"
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

#define TEST_SECTION(n)				SECTION(#n) { Test* test = new n##_Test(); test->run(); }
#define TEST_SECTION_V(n,v)			SECTION(#n) { Test* test = new n##_Test(v); test->run(); }
#define TEST_SECTION_N(n,num)		SECTION(#n) { Test* test = new n##_Test##num(); test->run(); }
#define TEST_SECTION_AV(n,a,d,v)	SECTION(#n) { Test* test = new n##_Test(a,d,v); test->run(); }

TEST_CASE("Auxiliary tests", "[Auxiliaries]") {	
	TEST_SECTION(IsnanIsinf)
	TEST_SECTION(GenericSingleton)
	TEST_SECTION(CommandLineParser)
}

TEST_CASE("Data storage tests", "[Storages]") {
	TEST_SECTION(LocalVariableStorage)
	TEST_SECTION(PropertyDatabase)
	TEST_SECTION(Index)
	TEST_SECTION(Parameter)
	TEST_SECTION(PropertyData)
}

TEST_CASE("Model tests", "[Models]") {
	TEST_SECTION(Node)
	TEST_SECTION(Element)
	TEST_SECTION(Face)
}

TEST_CASE("Variable tests", "[Variables]") {
	TEST_SECTION(Point);
	TEST_SECTION(ScalarVariable);
	TEST_SECTION(VectorVariable);
	TEST_SECTION_N(VectorVariable, 1);
	TEST_SECTION_N(VectorVariable, 2);
	TEST_SECTION(TensorVariable);
	TEST_SECTION_N(TensorVariable, 1);
	TEST_SECTION_N(TensorVariable, 2);
	TEST_SECTION(ArrayVariable);	
}

TEST_CASE("Utilities tests", "[Utilities]") {
	TEST_SECTION_V(Matrix,false)
	TEST_SECTION_V(DenseMatrix,false)
	TEST_SECTION(SparseMatrix)
  TEST_SECTION(CompressedRowMatrix)
  TEST_SECTION(CubicSpline)
	TEST_SECTION(FibonacciHeap)

	TEST_SECTION_V(FiniteVolumeStencil,false)
	// also compares speed of mapping facet areas and normals versus computing them
	TEST_SECTION(FiniteVolumePolicy)
	TEST_SECTION(FV_Parameter)
}

TEST_CASE("Interfaces / containers tests", "[Interfaces]") {
  TEST_SECTION(VData)
	TEST_SECTION(FEM_Data)
	TEST_SECTION(PropertyData)
 // TEST_SECTION(VSet)
  TEST_SECTION(ColorPalette)
	TEST_SECTION(MeshManager)
}

TEST_CASE("Finite elements and math operators", "[FiniteElements]"){
	TEST_SECTION(FEM_Data)
	TEST_SECTION_AV(FiniteElement,new IsoparametricLinearTetrahedron(1),"IsoparametricLinearTetrahedron1P.txt",false)
	TEST_SECTION_AV(FiniteElement,new IsoparametricLinearTetrahedron(4),"IsoparametricLinearTetrahedron4P.txt",false)
	TEST_SECTION_AV(FiniteElement,new IsoparametricLinearTriangle(3,3),"IsoparametricLinearTriangle3D3IP.txt",false)
	TEST_SECTION_AV(FiniteElement,new IsoparametricLinearTriangle(2,3),"IsoparametricLinearTriangle3IP.txt",false)
	TEST_SECTION_AV(FiniteElement,new IsoparametricLinearTriangle(2,4),"IsoparametricLinearTriangle4IP.txt",false)
	TEST_SECTION_AV(FiniteElement,new IsoparametricQuadraticTriangle(2),"IsoparametricQuadraticTriangle.txt",false)
	TEST_SECTION_V(IsoparametricQuadraticTetrahedron,false)
	TEST_SECTION_V(IsoparametricLinearHexahedron,false)
	TEST_SECTION(Operand)
	TEST_SECTION(MathOperatorLHS)
	TEST_SECTION(MathOperatorRHS)
	//TEST_SECTION(PDE_Integrator) requires a Model as parameter
}

TEST_CASE("Model-related interdependent functionality", "[ModelFunctionality]") {
	TEST_SECTION(ModelTopology)
	TEST_SECTION(Box)
	TEST_SECTION(ModelSubDomain)
	TEST_SECTION(Region)
	TEST_SECTION(BoundaryInterface)
	TEST_SECTION(ANSYS_Model2D)
	TEST_SECTION(InputDataManager)
	TEST_SECTION(ANSYS_Model3D)
	TEST_SECTION(PropertyHandle)
	//TEST_SECTION(BinaryFileInterface)
	TEST_SECTION(VTU_Interface)
	TEST_SECTION(StatisticalAnalyzer)	
  TEST_SECTION(Boundary)
  TEST_SECTION(SplitBoundary)
}

TEST_CASE("Composite-dependent functionality", "[CompositeFunctionality]") {
	TEST_SECTION(RegionMonitor)
	TEST_SECTION(ExponentialTransferFunction)
	TEST_SECTION_V(PropertyAtPointVisitor,false)
	//TEST_SECTION_V(TwoPhaseModel,composite);
}
