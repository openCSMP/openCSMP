#include <iostream>
#include <string>
#include "CSMP_definitions.h"
#include "Region.h"
#include "Exception.h"
#include "Test.h"
#include "TestSuite.h"

// FEM tests
#include "IncompressibleSinglePhaseFlowFEM_VVCase.h"
#include "IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase.h"
#include "IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase.h"

// FEFV Tests
#include "CFL_Calculation_VVCase.h"
#include "InflowOutflowCalculation_VVCase.h"

#include "IncompressibleTwoPhaseFlowFractures_Viscous_VVCase.h"

// Mechanics Tests
#include "LinearElasticityA_VVCase.h"
#include "BoreHole_stability2D_VVCase.h"
#include "BoreHole_stability_VerticalWell3D_VVCase.h"
#include "BoreHole_stability_InclinedWell3D_VVCase.h"

// Geothermal Tests
#include "Geothermal_pseudo1D_VVCase.h"
#include "Geothermal_2D_VVCase.h"
#include "Geothermal_3D_VVCase.h"


using namespace std;
using namespace csmp;

int main(int argc, char **argv )
{
    cout << "Received " << argc << " arguments...\n";
    for (int i=0; i<argc; i++)
        cout << "argument " << i << ": " << argv[i] << endl;
    long nFail(0);
    bool restart(false);
    for (int i=0; i<argc; i++)
    {
        if (string(argv[i])=="restart")
            restart=true;

    }


    try {

        //argv[1] is the Test name (also configuration files), and argv[2] is the project name (i.e. geometry)
        TestSuite s("VV Cases Test Suite", &cout );// Create test suite
		std::string test_name (argv[1]);

		if( test_name.find("_TestSuite") == std::string::npos ){

            //*****************************************************************************************************************
            // Tests list

            //Please add your test to the list imitating the ones shown below.  If you are beginning to implement a test,
			//please add it in the second section for verification cases under construction.
			//Tests in Section 1 are online.

			//*****************************************************************
			// SECTION 1 - HERE WE START WITH THE TESTS THAT ARE ACTUALLY ONLINE IN BUILDBOT
			//*****************************************************************

            //*****************************************************************************************************************
			//**** Finite Element, PDE_Integrator-related tests.
			s.addTest("IncompressibleSinglePhaseFlowFEM_VVCase", new IncompressibleSinglePhaseFlowFEM_VVCase(argv[2])); //By: Julian
			s.addTest("IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase", new IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase(argv[2])) ; //By: Julian
			s.addTest("IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase", new IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase(argv[2])) ; //By: Julian

			//*****************************************************************************************************************
			//**** Two-Phase flow tests. FEFV - In some cases, please scroll to the right to read more details about each test!
			//s.addTest("InflowOutflowCalculation_VVCase", new InflowOutflowCalculation_VVCase<2U>(argv[2])); //By: Philipp & Georg, modified by Julian
			s.addTest("InflowOutflowCalculation_VVCase", new InflowOutflowCalculation_VVCase<3U>(argv[2])); //By: Philipp & Georg, modified by Julian


            s.addTest("IncompressibleTwoPhaseFlowFractures_Viscous_VVCase_2D_BC_imp_1st", new IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<2U>(argv[2], "Brooks Corey Model","1st order")); //By: Christine

			//*****************************************************************
			// SECTION 2 - TESTS BELOW ARE STILL UNDER CONSTRUCTION
			//*****************************************************************

            //*****************************************************************************************************************
            //**** Mechanics Tests
            s.addTest("Linear Elasitcity 2D Beam Displacement", new LinearElasticityA_VVCase(argv[2])); //By: Philipp L.
            s.addTest("BoreHole_stability2D_VVCase", new BoreHole_stability2D_VVCase(argv[2])); //By: Mokhles M.
            s.addTest("BoreHole_stability_VerticalWell3D_VVCase", new BoreHole_stability_VerticalWell3D_VVCase(argv[2])); //By: Mokhles M.
            s.addTest("BoreHole_stability_InclinedWell3D_VVCase", new BoreHole_stability_InclinedWell3D_VVCase(argv[2])); //By: Mokhles M.

            //*****************************************************************************************************************
            //**** Geothermal Tests
            s.addTest("Geothermal_pseudo1D_VVCase", new Geothermal_pseudo1D_VVCase(argv[2])); //By: Alina
            s.addTest("Geothermal_2D_VVCase", new Geothermal_2D_VVCase(argv[2])); //By: Alina
            s.addTest("Geothermal_3D_VVCase", new Geothermal_3D_VVCase(argv[2])); //By: Alina

            //*****************************************************************************************************************
            //**** Two Phase Flow Tests

            s.addTest("CFL_Calculation_VVCase", new CFL_Calculation_VVCase<1U>(argv[2])); //By: Julian
            s.addTest("CFL_Calculation_VVCase", new CFL_Calculation_VVCase<2U>(argv[2])); //By: Julian
            s.addTest("CFL_Calculation_VVCase", new CFL_Calculation_VVCase<3U>(argv[2])); //By: Julian


            //*****************************************************************************************************************
            // end of Tests list


			s.RunSpecificTest(argv[1]);
			s.FreeAllButSpecificTest(argv[1]);
			nFail = s.report();
			if(!s.IsInSuite(argv[1]))
			{
				cout<<s.getName()<<" test : "<<string(argv[1])<<" not found. Please add test to vv_cases_main.cpp"<<endl;
				nFail=1;
			}
			s.free();

	}else{

			// will	be filled with test suites		
			s.run();
			nFail = s.report();
			s.free();
	}



    }
    // Exception handling (warnings etc. are caught at a much lower level)
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
    catch( TestSuiteError& ba ) {
        cout <<"\nException: Exception raised: "<< ba.what() << endl;
    }

    //returning 0 will mean that the test passed
    return nFail;

} // end main


