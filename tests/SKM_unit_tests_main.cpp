// Disable wingdi.h because it steps on our toes
#define NOGDI

#include "CSMP_definitions.h"
#include "Exception.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Exception.h"
#endif

#include "Test.h"
#include "TestSuite.h"

// fundamentals
#include "binaryReadWrite_Test.h"
#include "Colony_Test.h"
#include "ColorPalette_Test.h"
#include "ScalarVar_Test.h"
#include "VectorVar_Test.h"
#include "VectorVar_Test1.h"
#include "VectorVar_Test2.h"
#include "TensorVar_Test.h"
#include "TensorVar_Test1.h"
#include "TensorVar_Test2.h"
#include "ArrayVariable_Test.h"
#include "FlaggedArrayVariable_Test.h"
#include "Variables_Test.h"
#include "GenericSingleton_Test.h"
#include "ConvexPolygon_Test.h"
#include "vsetMakers.h"

// math
#include "IsnanIsinf_Test.h"
#include "CubicSpline_Test.h"
#include "DenseMatrix_Test.h"
#include "Matrix_Test.h"
#include "SparseMatrix_Test.h"
#include "CompressedRowMatrix_Test.h"
#include "FibonacciHeap_Test.h"
#include "geometricCalculations_Test.h"
#include "DynamicArray_Test.h"
#include "ColorPalette_Test.h"

// Model with variable storage
#include "Parameter_Test.h"
#include "PropertyDatabase_Test.h"
#include "Index_Test.h"
#include "INDEXandVariables_Test.h"
#include "LocalVariableStorage_Test.h"
#include "Box_Test.h"
#include "Point_Test.h"
#include "Node_Test.h"
#include "NodeManifold_Test.h"
#include "NodeFunctions_Test.h"
#include "Element_Test.h"
#include "ElementPolicyIntegrity_Test.h"
#include "FiniteElementPolicy_Test.h"
#include "FiniteElementManager_Test.h"
#include "Face_Test.h"
#include "InterFace_Test.h"
#include "MeshManager_Test.h"
#include "MeshManagementUtilities_Test.h"
#include "NodeManifoldManager_Test.h"
#include "ModelBasics_Test.h"
#include "parentElementStorageOptions_Test.h"
#include "ModelSubDomain_Test.h"
#include "BoundaryInterface_Test.h"
#include "Boundary_Test.h"
#include "Region_Test.h"
#include "SplitBoundary_Test.h"
#include "SplitBoundaryInterface_Test.h"
#include "SplitBoundaryTest_JK.h"
#include "CopyReplaceVisitor_Test.h"
// model manipulation and property retrieval
#include "IntegrationPointToNodePropertyVisitor_Test.h"
#include "PointPropertyToCellMapper2D_Test.h"
#include "PropertyHandle_Test.h"
#include "PropertyHandle_MathTest.h"
#include "PropertyAtPointVisitor_Test.h"

// finite elements
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricQuadraticTriangle.h"
#include "FiniteElement_Test.h"
#include "FiniteElement_Test2.h" // misc. tests
#include "IsoparametricQuadraticTetrahedron_Test.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPyramid_Test.h"
#include "IsoparametricLinearHexahedron_Test.h"
#include "LinearCuboid_Test.h"

// finite volumes
#include "ExplicitTransport_Test.h"
#include "FiniteVolumeStencil_Test.h"
#include "FiniteVolumePolicy_Test.h"
#include "FV_Parameter_Test.h"
#include "FluxMismatch_Test.h"
#include "GFVT_ParametricSpaceComputation_Test.h"
#include "GenericFiniteVolumeTransport_Test.h" // TODO: from Andrew Bromage: test fails for non-simplex elements

// interfaces to other software
// importing
#include "CommandLineParser_Test.h"
#include "VData_Test.h"
#include "VSet_Test1.h"
#include "VSet_Test2.h"
#include "vsetMakers_Test.h"
#include "PropertyData_Test.h"
#include "FEM_Data_Test.h"
#include "ModelTopology_Test.h"
#include "InputDataManager_Test.h"
#include "SKUA_FiniteElementMeshInterface_Test.h"
#include "ANSYS_Model2D_Test.h"
#include "ANSYS_Model3D_Test.h"
#include "ANSYS_SplitBoundaryMatch_Test.h"
#include "TRIANGLE_Interface_Test.h"
// TODO: test Quadrilaterator
// exporting
#include "UG4_UGX_FileExport_Test.h"
#include "VTU_Interface_Test.h"
#include "ModelComparator_Test.h"

// FE/FV integration
#include "Operand_Test.h"
#include "MathOperatorLHS_Test.h"
#include "MathOperatorRHS_Test.h"
// PDE operators
#include "Integral_dNT_lhsop_dN_dV_Test.h"
#include "Integral_dNT_lhsop_dN_NT_v_dN_dV_Test.h"
#include "Integral_NT_lhsop_N_dV_Test.h"
#include "Integral_dNT_rhsop_dN_dV_Test.h"
#include "Integral_dNT_rhsop_dV_Test.h"

#include "Integral_var_NT_lhsop_N_dV_Test.h"
#include "Integral_var_NT_rhsop_N_dV_Test.h"

#include "NumIntegral_NT_lhsop_N_dV_Test.h"
#include "NumIntegral_NT_rhsop_N_dV_Test.h"
#include "NumIntegral_dNT_lhsop_dN_dV_Test.h"
#include "NumIntegral_dNT_rhsop_dN_dV_Test.h"
#include "NumIntegral_dNT_op_dV_Test.h"
#include "NumIntegral_DNT_op_DN_NT_v_DN_dV_Test.h"
// PDE integrators
#include "PDE_Integrator_Test.h"
#include "PDE_Integrator_Transient_Test.h"
#include "PDE_Integrator_Computation_Test.h"

// constitutive relationships
#include "TwoPhaseModel_TestSuite.h"
#include "ExponentialTransferFunction_Test.h"
#include "TwoPhaseModelWithHysteresis_Test.h"
#include "H2O_CO2_NaCl_FlowFunctions.h"
#include "FlowFunctionsModule.h"

// analysis
#include "PropertyConstraints_Test.h"
#include "StatisticalAnalyzer_Test.h"
#include "RegionMonitor_Test.h"
#include "Visitor_TestSuite.h"

// performance tests
#include "VariableStorageSpeed_Test.h"
#include "FiniteVolumeStencilSpeed_Test.h"
#include "PropertyStorageSpeed_Test.h"
#include "CSMP_VariableBenchmarking_Test.h"
#include "JaggedArray3D_Comparison_Test.h"
#include "AccumulationSpeedProfiling_Test.h"
#include "ExactVersusNumericIntegrationSpeed_Test.h"
#include "FiniteVolumeStencilSpeed_Test.h"

// integration tests
#include "Geothermal_1D_VVCase.h"
#include "DirichletPressureBoxModel_VVCase.h"
#include "GravityInducedFluidPressure_Test.h"
#include "SplitBoundaryPressureDiffusion_Test.h"

using namespace std;
using namespace csmp;

constexpr bool COMPREHENSIVE_TESTING = false;

/**  Unit Test Development
 
@todo  Model illustrate functionality with Example(s)
 
@todo Additional Code Coverage Required: PDE_Integrator: assembly of solution matrix for systems

@todo Base majority of tests on small CSMP native models avoiding all the disk read/write

@todo  IsNan_Test - refactor as it is built on false premises

@todo  refactor   CompressedRowMatrix - include elimitation etc.

*/
int main()
 {
  constexpr bool verbose(false);

  const bool is_comp = (COMPREHENSIVE_TESTING == true) ? true : false;

  const bool test_fundamentals     = is_comp;
  const bool test_interdependent1  = is_comp;
  const bool test_interdependent2  = is_comp;
  const bool test_interfaces       = is_comp;
  const bool test_composite        = is_comp;
  const bool test_refactoring      = !is_comp; // True only when NOT comprehensive
  const bool test_new_developments = false;    // Always false for now
 
  long    fails_fundamentals(0),
          fails_interdependent1(0),
          fails_interdependent2(0),
          fails_interfaces(0),
          fails_composite(0),
          fails_new_developments(0),
          total_failures(0);

  try {
    cout <<"\nunit_test_main: running tests..."<< endl;

    // =========================================================================================================
    //
    //             TESTING REFACTORED CODE
    //
    // =========================================================================================================
    if ( test_refactoring ) {
        cout <<"\n5. Refactored and new code functionality: running tests..."<< endl;
        TestSuite refactored("CSMP-refactored code unit-test suite", &cout );

// TODO:       refactored.addTest( new IsoparametricLinearPyramid_Test() ); // pyramid is correct but test fails
// TODO:       refactored.addTest( new BoundaryInterface_Test() );          // - test alternative BREP flagging approaches
       

// TODO: work in progress:        refactored.addTest( new PropertyConstraints_Test() );
 //       refactored.addTest( new SplitBoundaryTest_JK() ); // tested OK


// TODO: does not run transient problem yet; compare analytic with num integrals
//        refactored.addTest( new PDE_Integrator_Transient_Test() ); // OK (Anne-Laure Tertois) granite_model1
//          refactored.addTest( new PDE_Integrator_Computation_Test() ); // OK
// TODO: test PDE_Integrator with periodic boundary conditions


//        refactored.addTest( new PDE_Integrator_Test() ); // refactored and passed 2/2/26
//        refactored.addTest( new NodeManifold_Test() ); refactored and passed 3/5/26
//        refactored.addTest( new SplitBoundary_Test() ); passed 2/5/2026 (only 2D version tested)
//        refactored.addTest( new ModelSubDomain_Test() );  passed 1/5/2026 (detected need to improve CreatSplitBoundaryBetween() method
//        refactored.addTest( new VData_Test() ); passed: 27/4/2026
//        refactored.addTest( new VSet_Test2() ); passed: 27/4/2026
//       refactored.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt", verbose ) );



// Exact integration
// -----------------
// LHS
//          refactored.addTest( new Integral_dNT_lhsop_dN_dV_Test(true) ); // tested: OK
//          refactored.addTest( new Integral_dNT_lhsop_dN_NT_v_dN_dV_Test(true) ); // advection-dispersion, tested: OK
//          refactored.addTest( new Integral_NT_lhsop_N_dV_Test(true) );     // tested: OK
//          refactored.addTest( new Integral_var_NT_lhsop_N_dV_Test(true) ); // tested: OK
// RHS
//          refactored.addTest( new Integral_dNT_rhsop_dN_dV_Test(true) ); // tested: OK
//          refactored.addTest( new Integral_dNT_rhsop_dV_Test(true) );    // gradient operand, tested: OK
//          refactored.addTest( new Integral_var_NT_rhsop_N_dV_Test(true) ); // tested: OK
          
// Numeric integration
// -------------------
// LHS
//          refactored.addTest( new NumIntegral_dNT_op_dN_NT_v_dN_dV_Test(true) ); // advection-dispersion, tested: OK
//          refactored.addTest( new NumIntegral_dNT_op_dV_Test(true) ); // tested: OK
//          refactored.addTest( new NumIntegral_dNT_lhsop_dN_dV_Test(true) ); // tested: OK
//          refactored.addTest( new NumIntegral_NT_lhsop_N_dV_Test(true) ); // tested: OK
// RHS
//          refactored.addTest( new NumIntegral_dNT_rhsop_dN_dV_Test(true) ); // tested: OK
//          refactored.addTest( new NumIntegral_NT_rhsop_N_dV_Test(true) ); // tested: OK


        // SplitBoundary related testing
        // -----------------------------
        // OK: refactored.addTest( new ModelSubDomain_Test() );
        
        // read and write SplitBoundary to file (SplitBasic22: TestWriteModelToDiskAndReadBackWithInterfaces())
        // OK: refactored.addTest( new ModelBasics_Test() );
        
        // creation of 2D SplitBoundary during simulation and reading and writing from file
        // OK: refactored.addTest( new ANSYS_Model2D_Test() );
        
        // creates 3D model with multiple split boundaries 'ModelDykeAllLayersSplit' writing it to disk and bringing it back and comparing them
        // OK: refactored.addTest( new ANSYS_Model3D_Test() );
        
        // creating a model that was split already in ANSYS, matching up node-matched but disconnected boundaries
        // OK: refactored.addTest( new ANSYS_SplitBoundaryMatch_Test() );
        
        // insert lower-dim fracture into split boundary and test it
        // OK: LFEM refactored.addTest( new SplitBoundary_Test() ); // 2D only
        
        // 2 and 3D testing of creation methods for split boundaries: TODO: revisit correctness and reinstate all component tests
        // OK: refactored.addTest( new SplitBoundaryInterface_Test() );

//        refactored.addTest( new MeshManagementUtilities_Test() ); // TODO: complete this test
        
        /*
            Compares physical space with parametric space computations
            Not using test framework yet but printing everything to std::cerr
        */
//        refactored.addTest( new GFVT_ParametricSpaceComputation_Test() );  // TODO: understand why test is failing for non-simplex elements
//        ex         refactored.addTest( new GenericFiniteVolumeTransport_Test() );


// TODO: no satisfactorily fast access yet; compare with access of a fictious rock-type to determine whether improvement would pay off
// tested: 6/1/26: no speed-up from extra inlining, complications when attempting to remove macros in LocalVariableStorageArithmetic
// added new method to read vecs, tensors and arrays using declarative programming
//        refactored.addTest( new VariableStorageSpeed_Test() );
                
        
   // TODO: uncomment and fix failing tests listed below
   // fail - needs refactoring      interdependent1.addTest( new Integral_var_NT_lhsop_N_dV_Test( verbose ) );
   // fail - needs refactoring      interdependent1.addTest( new Integral_var_NT_rhsop_N_dV_Test( verbose ) );

        // actually running the test
        refactored.run();
        long nFail = refactored.report();
        total_failures += nFail;
        refactored.free();
        cerr << "\nunit_tests_main: 5. CSMP refactored and new functionality: Total unit test failures: " << nFail << endl;

    } // end refactoring



    // =========================================================================================================
    //
    //             FUNDAMENTALS
    //
    // =========================================================================================================
    cout <<"\nunit_test_main: running comprehensive suite of tests in ";
#ifndef NDEBUG
    cout <<"DEBUG mode..."<< endl;
#else
    cout <<"RELEASE mode..."<< endl;
#endif
    auto t0 = chrono::high_resolution_clock::now();

    if ( test_fundamentals ) {
      cout <<"\n"<<"1. underpinning functionality: running tests..."<< endl;
      TestSuite basic("CSMP-fundamental-unit test suite", &cout );
      
      // Auxiliaries
      basic.addTest( new Colony_Test() );
      basic.addTest( new IsnanIsinf_Test() );
      basic.addTest( new GenericSingleton_Test() );
      basic.addTest( new ColorPalette_Test() );
      
      // Data storage tests
      basic.addTest( new LocalVariableStorage_Test() );
      basic.addTest( new PropertyDatabase_Test());
      basic.addTest( new Index_Test());
      basic.addTest( new Parameter_Test());
      basic.addTest( new PropertyData_Test());
      
      // Model
      basic.addTest( new Node_Test() );
      basic.addTest( new NodeManifold_Test() );
      basic.addTest( new Element_Test());
      basic.addTest( new Face_Test() );
      basic.addTest( new InterFace_Test() );
      
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
      basic.addTest( new FlaggedArrayVariable_Test());
      
      // utilities tests
      basic.addTest( new Matrix_Test(verbose) );
      basic.addTest( new DenseMatrix_Test(verbose) );
      basic.addTest( new SparseMatrix_Test() );
      basic.addTest( new CompressedRowMatrix_Test() );
      basic.addTest( new CubicSpline_Test() );
      basic.addTest( new FibonacciHeap_Test() );
      basic.addTest( new DynamicArray_Test() );
      basic.addTest( new geometricCalculations_Test() );
      
      basic.addTest( new FiniteVolumeStencil_Test(verbose));
      // also compares speed of mapping facet areas and normals versus computing them
      basic.addTest( new FiniteVolumePolicy_Test());
      basic.addTest( new FV_Parameter_Test());
      
      // interfaces / containers
      basic.addTest( new FEM_Data_Test());
      basic.addTest( new Parameter_Test() );
      basic.addTest( new PropertyData_Test() );
      basic.addTest( new VData_Test() );
      basic.addTest( new VSet_Test1() );
      basic.addTest( new ColorPalette_Test() );
      
      // Running unit tests and reporting errors
      basic.run();
      fails_fundamentals = basic.report();
      total_failures += fails_fundamentals;
      basic.free();
      cerr << "\nunit_tests_main: 1. CSMP fundamentals: test failures: " << fails_fundamentals << endl;
    }
    


    // =========================================================================================================
    //
    //             INTERDPENDENT1: FINITE ELEMENTS + MATH OPERATORS
    //
    // =========================================================================================================
    if ( test_interdependent1 ) {
      cout <<"\n"<<"2. partially interdependent functionality: running tests..."<< endl;
      TestSuite interdependent1("CSMP-interdependent1-unit test suite", &cout );
      // simplex finite elements
      // -----------------------
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(1), "IsoparametricLinearTetrahedron1P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTetrahedron(4), "IsoparametricLinearTetrahedron4P.txt", verbose ) );
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(3,3), "IsoparametricLinearTriangle3D3IP.txt", verbose ) ); // 3D case 3 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,3), "IsoparametricLinearTriangle3IP.txt", verbose ) );   // 2D case 3 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricLinearTriangle(2,4), "IsoparametricLinearTriangle4IP.txt", verbose ) );   // 2D case 4 integration points
      interdependent1.addTest( new FiniteElement_Test( new IsoparametricQuadraticTriangle(2), "IsoparametricQuadraticTriangle.txt", verbose ) );  // 3D case 3 integration point
      // non-standard element tests
      interdependent1.addTest( new IsoparametricQuadraticTetrahedron_Test(verbose) ); // TODO: no flux balance for constant velocity projected on sides
      // other elements
      // --------------
      interdependent1.addTest( new LinearCuboid_Test(verbose) );
      // volume conservation of distorted hexahedra - fails for certain deformation modes, highlighting limitations of this elements
      interdependent1.addTest(new IsoparametricLinearHexahedron_Test(verbose));
      // interdependent1.addTest(new IsoparametricLinearPyramid_Test(verbose)); // TODO: pyramid produces correct output but test is broken
      // finite element policy
      interdependent1.addTest( new FiniteElementPolicy_Test() ); // needs refactor, no ANSYS model required
      interdependent1.addTest( new ElementPolicyIntegrity_Test() );
      // math operators etc.
      interdependent1.addTest( new Operand_Test() ); // needs refactor, no ANSYS model required
      interdependent1.addTest( new MathOperatorLHS_Test());
      interdependent1.addTest( new MathOperatorRHS_Test());
      // finite element stencil manager
      interdependent1.addTest( new FiniteElementManager_Test() );

// Exact integration
// -----------------
// LHS
      interdependent1.addTest( new Integral_dNT_lhsop_dN_dV_Test(true) ); // tested: OK
      interdependent1.addTest( new Integral_dNT_lhsop_dN_NT_v_dN_dV_Test(true) ); // advection-dispersion, tested: OK
      interdependent1.addTest( new Integral_NT_lhsop_N_dV_Test(true) );     // tested: OK
      interdependent1.addTest( new Integral_var_NT_lhsop_N_dV_Test(true) ); // tested: OK
// RHS
      interdependent1.addTest( new Integral_dNT_rhsop_dN_dV_Test(true) ); // tested: OK
      interdependent1.addTest( new Integral_dNT_rhsop_dV_Test(true) );    // gradient operand, tested: OK
      interdependent1.addTest( new Integral_var_NT_rhsop_N_dV_Test(true) ); // tested: OK
          
// Numeric integration
// -------------------
// LHS
      interdependent1.addTest( new NumIntegral_dNT_op_dN_NT_v_dN_dV_Test(true) ); // advection-dispersion, tested: OK
      interdependent1.addTest( new NumIntegral_dNT_op_dV_Test(true) ); // tested: OK
      interdependent1.addTest( new NumIntegral_dNT_lhsop_dN_dV_Test(true) ); // tested: OK
      interdependent1.addTest( new NumIntegral_NT_lhsop_N_dV_Test(true) ); // tested: OK
// RHS
      interdependent1.addTest( new NumIntegral_dNT_rhsop_dN_dV_Test(true) ); // tested: OK
      interdependent1.addTest( new NumIntegral_NT_rhsop_N_dV_Test(true) ); // tested: OK
      
      // variables but now tested on a full-fledged model
      interdependent1.addTest( new Variables_Test("FracBox") );
      interdependent1.addTest( new VariableStorageSpeed_Test() );

      // running unit tests and reporting errors
      interdependent1.run();
      fails_interdependent1 = interdependent1.report();
      total_failures += fails_interdependent1;
      interdependent1.free();
      cerr << "\nunit_tests_main: 2. CSMP interdependent-functionality1: test failures: " << fails_interdependent1 << endl;
    }



    // =========================================================================================================
    //
    //             INTERDPENDENT2: MODEL, BOUNDARY, SPLIT-BOUNDARY
    //
    // =========================================================================================================
    if ( test_interdependent2 ) {
      // making a test model for the following two tests
      VSet<3U> vset;
      ModelTopology topo = create_FracBox( vset );
      const bool convert_side_surfs_into_boundaries{true};
      Model<3U> model( topo, vset, "CSMP-variables.txt", convert_side_surfs_into_boundaries );
      // creating the test suite
      cout <<"\n"<<"3. Model-related interdependent functionality: running tests..."<< endl;
      TestSuite interdependent2("CSMP-interdependent2-unit test suite", &cout );
      interdependent2.addTest( new INDEXandVariables_Test() );
      interdependent2.addTest( new ModelTopology_Test() ); // TODO: tests only minor functionality
      interdependent2.addTest( new MeshManager_Test() );
      interdependent2.addTest( new ModelBasics_Test() );
      interdependent2.addTest( new VSet_Test2() );
      interdependent2.addTest( new NodeFunctions_Test() );
      // model
      interdependent2.addTest( new Box_Test() );
      interdependent2.addTest( new ModelSubDomain_Test() );
      interdependent2.addTest( new NodeManifoldManager_Test() );
      interdependent2.addTest( new Region_Test() );
      interdependent2.addTest( new BoundaryInterface_Test() );
      interdependent2.addTest( new Boundary_Test() );
      interdependent2.addTest( new SplitBoundaryInterface_Test() );
      interdependent2.addTest( new SplitBoundary_Test() );
      interdependent2.addTest( new ANSYS_SplitBoundaryMatch_Test() );
      interdependent2.addTest( new PropertyHandle_Test() );    
      interdependent2.addTest( new PropertyHandle_MathTest() );
      interdependent2.addTest( new IntegrationPointToNodePropertyVisitor_Test() ); // insufficient accuracy for IsoLinPyra
      interdependent2.addTest( new CopyReplaceVisitor_Test( &model ) );
      // interfaces
      interdependent2.addTest( new InputDataManager_Test());
      interdependent2.addTest( new ANSYS_Model3D_Test() );
      interdependent2.addTest( new ANSYS_Model2D_Test() );
      interdependent2.addTest( new VTU_Interface_Test() );
      interdependent2.addTest( new StatisticalAnalyzer_Test() );
      // running unit tests and reporting errors
      interdependent2.run();
      fails_interdependent2 = interdependent2.report();
      total_failures += fails_interdependent2;
      interdependent2.free();
      cerr << "\nunit_tests_main: 3. CSMP Model-related, interdependent-functionality2: test failures: " << fails_interdependent2 << endl;
    }



    // =========================================================================================================
    //
    //             INTERFACES TO OTHER SOFTWARE: ANSYS, VTK, SKUA
    //
    // =========================================================================================================
    if ( test_interfaces ) {
        cout <<"\n"<<"4. Input & output interfaces of CSMP: running tests..."<< endl;
        TestSuite interfaces("CSMP-refactored code unit-test suite", &cout );

        interfaces.addTest( new TRIANGLE_Interface_Test() ); // FAILED ASSERTION ON COLLOCATED NODES
        interfaces.addTest( new ANSYS_SplitBoundaryMatch_Test() );   // WORKS 11/11/2024
        interfaces.addTest( new ANSYS_Model2D_Test() );              // WORKS 11/11/2024
        interfaces.addTest( new ANSYS_Model3D_Test() );              // WORKS 22/10/2024

        interfaces.addTest( new VTU_Interface_Test() );
        interfaces.addTest( new UG4_UGX_FileExport_Test() );

// TODO: write and create tests for interfaces: Gmsh, GeoModeller, Eclipse, Matlab, Maple, JPEG, Rhino, Tecplot etc. here
// TODO: SKUA NEEDS RECREATION of binary files:  composite.addTest( new SKUA_FiniteElementMeshInterface_Test() );
    
        // actually running the test
        interfaces.run();
        total_failures += fails_interfaces;
        interfaces.free();
        cerr << "\nunit_tests_main: 4. input & output interfaces: Total unit test failures: " << fails_interfaces << endl;
        
    } // end refactoring



    // =========================================================================================================
    //
    //             COMPOSITE FUNCTIONALITY
    //
    // =========================================================================================================
    if ( test_composite ) {
      // creating test models
      VSet<2U> vset2D;
      ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset2D );
      vset2D.RemoveData("element variable"); // not needed here
      bool treat_domains_as_regions{true};
      Model<2U> model2D( topo, vset2D, "CSMP-1phase-variables.txt", treat_domains_as_regions );
      VSet<3U> vset3D;
      topo = create_FracBox( vset3D );
      Model<3U> model3D( topo, vset3D, "CSMP-1phase-variables.txt", treat_domains_as_regions );
      
      cout <<"\n4. Composite-dependent functionality: running tests..."<< endl;
      TestSuite composite("CSMP-dependent-unit test suite", &cout );

      // vistors
      Visitor_TestSuite visitorTests( composite );
      visitorTests.run();
      composite.addTest( new PropertyAtPointVisitor_Test(verbose) ); // PASS
      composite.addTest( new PointPropertyToCellMapper2D_Test() );
      composite.addTest( new PropertyConstraints_Test() );
      composite.addTest( new PropertyStorageSpeed_Test( &cout ) );

      // computations
      composite.addTest( new PDE_Integrator_Test( model2D ) );
      composite.addTest( new PDE_Integrator_Transient_Test() );
      composite.addTest( new PDE_Integrator_Computation_Test() );
      // TODO: test periodic BCs
      // misc
      composite.addTest( new RegionMonitor_Test() );
      composite.addTest( new ModelComparator_Test() );

      // constitutive relationships
      composite.addTest( new ExponentialTransferFunction_Test() );
      composite.addTest( new PropertyStorageSpeed_Test( &clog ) ); // luxury extra

      // constitutive relations
      TwoPhaseModel_TestSuite  twoPhaseModelTests( composite );
      twoPhaseModelTests.run();
      composite.addTest( new TwoPhaseModelwithHysteresis_Test() );
      // running unit tests and reporting errors
      composite.run();
      fails_composite = composite.report();
      total_failures += fails_composite;
      composite.free();
      cerr << "\nunit_tests_main: 4. CSMP-dependent-functionality: test failures: " << fails_composite << endl;
    }
        
    // tests related to code that is currently being refactored
    if ( test_new_developments ) {
      cout <<"\n5. Refactored and new code functionality: running tests..."<< endl;
      TestSuite new_developments("new tests of the CSMP base library", &cout );
      
      new_developments.addTest( new SplitBoundaryPressureDiffusion_Test() );
      
      // running unit tests and reporting errors
      new_developments.run();
      fails_new_developments = new_developments.report();
      total_failures += fails_new_developments;
      new_developments.free();
      cerr << "\nunit_tests_main: 6. New functionality: test failures: " << fails_new_developments << endl;
    }
    
    auto t1 = chrono::high_resolution_clock::now();
	  cout <<"\n\t"<<"Time taken to run the comprensive suite of unit tests: "<< chrono::duration_cast<chrono::seconds>(t1-t0).count() << " seconds." << endl;
   
    cout << "\nunit_tests_main: Total unit test failures: ";
    total_failures = fails_fundamentals + fails_interdependent1 + fails_interdependent2 + fails_interfaces + fails_composite + fails_new_developments;
    cerr << total_failures << endl;
    if ( fails_fundamentals > 0 )     cerr <<"\nfundamental functionality tests failed.";
    if ( fails_interdependent1 > 0 )  cerr <<"\ninterpedendent (basic) functionality tests failed.";
    if ( fails_interdependent2 > 0 )  cerr <<"\ninterpedendent (advanced) functionality tests failed.";
    if ( fails_interfaces > 0 )       cerr <<"\ninput & output interfaces tests failed.";
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

  return 0;

} // end main


// ECMOR24 - speed stuff

// refactored.addTest( new ExactVersusNumericIntegrationSpeed_Test() );
//  refactored.addTest( new AccumulationSpeedProfiling_Test() );
// refactored.addTest( new VariableBenchmarking_Test() );  // FAIL

// refactored.addTest( new PropertyStorageSpeed_Test( &cout ) );
// refactored.addTest( new VariableStorageSpeed_Test() );
// refactored.addTest( new JaggedArray3D_Comparison_Test() );
// refactored.addTest( new FiniteVolumeStencilSpeed_Test() );
