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
#include "CopyReplaceVisitor_Test.h"
// model manipulation and property retrieval
#include "IntegrationPointToNodePropertyVisitor_Test.h"
#include "PointPropertyToCellMapper2D_Test.h"
#include "PropertyHandle_Test.h"
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
// exporting
#include "UG4_UGX_FileExport_Test.h"
#include "VTU_Interface_Test.h"
#include "ModelComparator_Test.h"

// FE/FV integration
#include "PDE_Integrator_Test.h"
#include "Operand_Test.h"
#include "MathOperatorLHS_Test.h"
#include "MathOperatorRHS_Test.h"
#include "Integral_var_NT_lhsop_N_dV_Test.h"
#include "Integral_var_NT_rhsop_N_dV_Test.h"
#include "PDE_Integrator_Test.h"
#include "PDE_Integrator_Transient_Test.h"

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

//  number-of-args  arg-strings
int main( int argc, char* argv[] )
{
    bool verbose = true;
    bool run_all = false;
    std::string target_suite = "";

    // 1. Basic Command Line Parser
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--suite" && i + 1 < argc) {
            target_suite = argv[++i];
        } else if (arg == "--all") {
            run_all = true;
        } else if (arg == "--quiet") {
            verbose = false;
        }
    }

    // Default to comprehensive if no suite is specified
    if (target_suite.empty()) {
        run_all = true;
    }

    // 2. Map flags based on CLI arguments
    const bool test_fundamentals     = run_all || (target_suite == "fundamentals");
    const bool test_interdependent1  = run_all || (target_suite == "interdependent1");
    const bool test_interdependent2  = run_all || (target_suite == "interdependent2");
    const bool test_interfaces       = run_all || (target_suite == "interfaces");
    const bool test_composite        = run_all || (target_suite == "composite");
    const bool test_refactoring      = run_all || (target_suite == "refactored");
    const bool test_new_developments = run_all || (target_suite == "new");

    long  fails_fundamentals(0),
          fails_interdependent1(0),
          fails_interdependent2(0),
          fails_interfaces(0),
          fails_composite(0),
          fails_new_developments(0),
          total_failures(0);

    try {
        if (run_all) cout <<"\nunit_test_main: running ALL tests..."<< endl;

        // =====================================================================
        // REFACTORED CODE
        // =====================================================================
        if ( test_refactoring ) {
        cout <<"\n1. Refactored and new code functionality: running tests..."<< endl;
        TestSuite refactored("CSMP-refactored code unit-test suite", &cout );

        // add your test of the suggested refactoring here
        //refactored.addTest( new PDE_Integrator_Test() );

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
      // finite element policy
      interdependent1.addTest( new FiniteElementPolicy_Test() ); // needs refactor, no ANSYS model required
      interdependent1.addTest( new ElementPolicyIntegrity_Test() );
      // math operators etc.
      interdependent1.addTest( new Operand_Test() ); // needs refactor, no ANSYS model required
      interdependent1.addTest( new MathOperatorLHS_Test());
      interdependent1.addTest( new MathOperatorRHS_Test());
      // finite element stencil manager
      interdependent1.addTest( new FiniteElementManager_Test() );

// fail - TODO: legacy code needs refactoring      interdependent1.addTest( new Integral_var_NT_lhsop_N_dV_Test( verbose ) );
// fail - TODO: legacy code needs refactoring      interdependent1.addTest( new Integral_var_NT_rhsop_N_dV_Test( verbose ) );
      
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
      interdependent2.addTest( new ModelSubDomain_Test() );      // TODO: runs but does not test core functionality
      interdependent2.addTest( new NodeManifoldManager_Test() );
      interdependent2.addTest( new Region_Test() );
      interdependent2.addTest( new BoundaryInterface_Test() );
      interdependent2.addTest( new Boundary_Test() );
      interdependent2.addTest( new SplitBoundaryInterface_Test() );
      interdependent2.addTest( new SplitBoundary_Test() );
      interdependent2.addTest( new ANSYS_SplitBoundaryMatch_Test() );
      interdependent2.addTest( new PropertyHandle_Test() );    // TODO: refactor without ANSYS model
      interdependent2.addTest( new IntegrationPointToNodePropertyVisitor_Test() ); // insufficient accuracy for IsoLinPyra
      interdependent2.addTest( new CopyReplaceVisitor_Test( &model ) );
      // interfaces
      interdependent2.addTest( new InputDataManager_Test());
      interdependent2.addTest( new ANSYS_Model3D_Test() );
      interdependent2.addTest( new ANSYS_Model2D_Test() );
      interdependent2.addTest( new VTU_Interface_Test() ); // TODO: fails for split boundaries, needs refactoring
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
      // TODO: add test of assembly of matrix for systems, elimination of boundary conditions etc.
      composite.addTest( new PDE_Integrator_Test( model2D ) ); // TODO: not comprehensive enough, haha!
      composite.addTest( new PDE_Integrator_Transient_Test() ); // TODO: not comprehensive enough
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
      
       new_developments.addTest( new PropertyStorageSpeed_Test( &cout ) );
       new_developments.addTest( new VariableStorageSpeed_Test() );
       new_developments.addTest( new JaggedArray3D_Comparison_Test() );
       new_developments.addTest( new FiniteVolumeStencilSpeed_Test() );
       new_developments.addTest( new AccumulationSpeedProfiling_Test() );
       new_developments.addTest( new ExactVersusNumericIntegrationSpeed_Test() );

  //     new_developments.addTest( new VariableBenchmarking_Test() );  // FAIL
  //     new_developments.addTest( new SplitBoundaryPressureDiffusion_Test() );

   
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


