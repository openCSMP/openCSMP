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

#include "GenericFiniteVolumeTransport_Test.h"
#include "FiniteVolumeTransportBasics_Test.h"

// new tests from SKM
#include "ANSYS_SplitBoundaryMatch_Test.h"

#include "Vset_TestCase.h"

using namespace std;
using namespace csmp;

TEST_CASE("Development tests", "[Dev]") {

    // ModelSubDomain_Test test;

    // Vset_TestCase test;
    // GenericFiniteVolumeTransport_Test test;
    FiniteVolumeTransportBasics_Test test;
    // FiniteVolumeStencil_Test test;
    // ANSYS_SplitBoundaryMatch_Test  skm_test;
    test.run();

 //   ANSYS_Model2D_Test test;
 //   test.run();
}

