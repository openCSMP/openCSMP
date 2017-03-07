#include <iostream>
#include <string>
#include "CSMP_definitions.h"
#include "Region.h"
#include "Exception.h"
#include "Test.h"
#include "TestSuite.h"

#include "CsmpIntro1_TestCase.h"

#include "Vset_Test.h"
#include "Variables_Test.h"


using namespace std;

using namespace csmp;

int main(int argc, char **argv)
{
    cout << "Received " << argc << " arguments...\n";
    for (int i=0; i<argc; i++)
        cout << "argument " << i << ": " << argv[i] << endl;
    long nFail(0);

    try {

        //argv[1] is the Test name (also configuration files), and argv[2] is the project name (i.e. geometry)
        TestSuite s("Test Cases Test Suite", &cout );// Create test suite

        //Please add your test to the list imitating the ones shown below.  If you are beginning to implement a test,
        //please add it in the second section for test cases under construction.
        //Tests in Section 1 are online.

        //*****************************************************************
        // SECTION 1 - HERE WE START WITH THE TESTS THAT ARE ACTUALLY ONLINE IN BUILDBOT
        //*****************************************************************
        s.addTest("Variable_Test", new Variables_Test(argv[2]));
//        s.addTest("Vset_Test", new VSet_Test(argv[2]));
        s.addTest("Vset_Test", new VSet_Test());

        //*****************************************************************
        // SECTION 2 - TESTS BELOW ARE STILL UNDER CONSTRUCTION
        //*****************************************************************


        s.RunSpecificTest(argv[1]);
        s.FreeAllButSpecificTest(argv[1]);
        nFail = s.report();
        if(!s.IsInSuite(argv[1]))
        {
            cout<<s.getName()<<" test : "<<string(argv[1])<<" not found. Please add test to vv_cases_main.cpp"<<endl;
            nFail=1;
        }
        s.free();

    }
    // Exception handling (warnings etc. are caught at a much lower level)
    catch( bad_alloc& ba ) {
        cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( bad_cast& ba ) {
        cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( bad_exception& ba ) {
        cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( bad_typeid& ba ) {
        cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( ios_base::failure& ba ) {
        cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
        system("pause");
    }
    // standard logic errors
    catch( domain_error& ba ) {
        cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( invalid_argument& ba ) {
        cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( length_error& ba ) {
        cout <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( out_of_range& ba ) {
        cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
        system("pause");
    }
    // runtime errors
    catch( overflow_error& ba ) {
        cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( range_error& ba ) {
        cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( underflow_error& ba ) {
        cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
        system("pause");
    }
    catch( Exception& ba ) {
        cout <<"\nException: Exception raised: "<< ba.What() << endl;
        cout <<"\nDiagnostics:"<< endl;
        ba.Out();
        system("pause");
    }
    catch( TestSuiteError& ba ) {
        cout <<"\nException: Exception raised: "<< ba.what() << endl;
    }

    //returning 0 will mean that the test passed
    return nFail;

} // end main





