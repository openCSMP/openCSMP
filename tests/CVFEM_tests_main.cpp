// CVFEM_tests_main.cpp
//
// Runs 1D and 2D Benchmarks to test base functionality of CVFEM used for hydrothermal system modelling

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#ifdef CSMP_WITH_PETSC_SOLVER
#include <petscsys.h>
#endif

#include "Test.h"
#include "TestSuite.h"
#include "Exception.h"
#include "SAMG_Exception.h"

#include "CVFEM_1D_VVCase.h"
#include "CVFEM_2D_VVCase.h"

#include <stdexcept>
#include <typeinfo>
#ifdef __GNUC__
#include "cxxabi.h"
#endif

using namespace std;
using namespace csmp;


// Single source of truth for every known 1D/2D test name — used both to build the
// suites below and to validate the user's --test-selection input against typos.
const std::vector<std::string> ALL_TEST_FILES_1D = {
    "test-101",
    "test-101-wg",
    "test-102",
    "test-102-wg",
    "test-103",
    "test-103-wg",
    "test-104",
    "test-104-wg",
};

const std::vector<std::string> ALL_TEST_FILES_2D = {
    "test-201a",
    "test-201b",
    "test-202a",
    "test-202b",
    "test-203a",
    "test-203b",
};


std::string getTestDataPath_1D(const std::string& benchmarkDir, const std::string& filename) {
    return benchmarkDir + "/1D_Benchmarks/" + filename;
}

std::string getTestDataPath_2D(const std::string& benchmarkDir, const std::string& filename) {
    return benchmarkDir + "/2D_Benchmarks/" + filename + "/" + filename;
}

std::string getMeshPath2D(const std::string& benchmarkDir) {
    return benchmarkDir + "/2D_Benchmarks/2D_Benchmarks";
}

// Splits a comma-separated string into trimmed, non-empty tokens, e.g.
// "test-201a, test-202b" -> {"test-201a", "test-202b"}
std::vector<std::string> parseCommaSeparatedList(const std::string& input) {
    std::vector<std::string> tokens;
    std::stringstream ss(input);
    std::string token;
    while (std::getline(ss, token, ',')) {
        const size_t start = token.find_first_not_of(" \t");
        const size_t end   = token.find_last_not_of(" \t");
        if (start != std::string::npos)
            tokens.push_back(token.substr(start, end - start + 1));
    }
    return tokens;
}

// Returns allTests unchanged if runAll is true; otherwise returns only the entries of
// allTests that also appear in requestedTests, preserving allTests' original order.
std::vector<std::string> filterRequestedTests(const std::vector<std::string>& allTests,
                                              const std::vector<std::string>& requestedTests,
                                              bool runAll) {
    if (runAll)
        return allTests;

    std::vector<std::string> filtered;
    for (const std::string& name : allTests) {
        if (std::find(requestedTests.begin(), requestedTests.end(), name) != requestedTests.end())
            filtered.push_back(name);
    }
    return filtered;
}


int main( int argc, char* argv[] )
{
#ifdef CSMP_WITH_PETSC_SOLVER
    PetscInitializeNoArguments();
#endif

    bool verbose = true;

    // 1. Minimal command-line parsing
    //    (kept simple since there's only one suite here; extend if needed)
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--quiet") verbose = false;
    }

    long total_failures = 0;

    //! USER INPUT
    //! Define main path
    std::string benchmarkDir;
    std::cerr << "\nEnter the path to the benchmark data directory "
                 "\n(the folder containing 1D_Benchmarks/ and 2D_Benchmarks/): ";
    while ( !std::getline(std::cin, benchmarkDir) || benchmarkDir.empty() )
    {
        std::cerr << "Invalid input. Please enter a non-empty path: ";
    }
    // benchmarkDir = "/home/jkoepping/Documents/0_modelling/openCSMP_2026_github/openCSMP-v3/tests/data-tests_CVFEM"; // JK


    // define which suite(s) to run:
    uint32_t suite_selection = 0;
    std::cerr << "\nSelect which simulations to run:"
                 "\n  0 = 1D only"
                 "\n  1 = 2D only"
                 "\nEnter 0 or 1: ";
    while ( !(std::cin >> suite_selection) || suite_selection > 1 )
    {
        std::cerr << "Invalid input. Please enter 0, or 1: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    const bool run1D = (suite_selection == 0);
    const bool run2D = (suite_selection == 1);


    // define test_type:
    uint32_t test_type = 0;
    std::cerr << "\nSelect test type:"
                 "\n  0 = QUICK    — only test model initialisation and output after the first iteration"
                 "\n  1 = EXTENDED — run the whole simulation and compare initial, first-iteration, and final state"
                 "\n      NOTE that the EXTENDED version takes a couple of hours but also allows for more detailed comparison with Benchmark data."
                 "\nEnter 0 or 1: ";
    while ( !(std::cin >> test_type) || (test_type != 0 && test_type != 1) )
    {
        std::cerr << "Invalid input. Please enter 0 or 1: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // discard leftover newline from test_type, once, before any getline below

    // define which specific tests to run:
    bool runAllTests_1D = false, runAllTests_2D = false;
    std::vector<std::string> requestedTests_1D, requestedTests_2D;
    if ( run1D )
    {
        std::cerr << "\nEnter the tests to run, comma-separated (test-101, test-101-wg, test-102, test-102-wg, test-103, test-103-wg, test-104, test-104-wg),"
                     "\nor type 'all' to run every test in the selected suite(s): ";

        std::string testSelectionInput;
        std::getline(std::cin, testSelectionInput);
        {
            const size_t start = testSelectionInput.find_first_not_of(" \t");
            const size_t end   = testSelectionInput.find_last_not_of(" \t");
            testSelectionInput = (start != std::string::npos) ? testSelectionInput.substr(start, end - start + 1) : "";
        }
        runAllTests_1D = (testSelectionInput == "all" || testSelectionInput == "All" || testSelectionInput == "ALL");
        requestedTests_1D = runAllTests_1D ? std::vector<std::string>{} : parseCommaSeparatedList(testSelectionInput);
    }
    if ( run2D )
    {
        std::cerr << "\nEnter the tests to run, comma-separated (test-201a, test-201b, test-202a, test-202b, test-203a, test-203b),"
                     "\nor type 'all' to run every test in the selected suite(s): ";

        std::string testSelectionInput;
        std::getline(std::cin, testSelectionInput);
        {
            const size_t start = testSelectionInput.find_first_not_of(" \t");
            const size_t end   = testSelectionInput.find_last_not_of(" \t");
            testSelectionInput = (start != std::string::npos) ? testSelectionInput.substr(start, end - start + 1) : "";
        }
        // std::string testSelectionInput = "test-202a"; // JK
        runAllTests_2D = (testSelectionInput == "all" || testSelectionInput == "All" || testSelectionInput == "ALL");
        requestedTests_2D = runAllTests_2D ? std::vector<std::string>{} : parseCommaSeparatedList(testSelectionInput);
    }

    try {
        if ( run1D )
        {
            const std::vector<std::string> configFiles1D = filterRequestedTests( ALL_TEST_FILES_1D, requestedTests_1D, runAllTests_1D );

            cerr << "\nCVFEM_tests_main: running 1D test suite (" << configFiles1D.size() << " tests)...\n" << endl;

            // 2. Build the suite
            TestSuite suite1D( "New CVFEM test suite", &cerr );

            std::vector<std::string> configPaths1D;
            configPaths1D.reserve( configFiles1D.size() );


            for ( const std::string& filename : configFiles1D )
            {
                configPaths1D.push_back( getTestDataPath_1D(benchmarkDir, filename) );
                suite1D.addTest( filename.c_str(), new CVFEM_1D_VVCase( configPaths1D.back().c_str(), test_type ) );
            }

            // 3. Run, report, clean up
            suite1D.run();
            total_failures = suite1D.report();
            suite1D.free();
        }

        // 2D
        if ( run2D )
        {
            const std::vector<std::string> configFiles2D = filterRequestedTests( ALL_TEST_FILES_2D, requestedTests_2D, runAllTests_2D );

            cerr << "\nCVFEM_tests_main: running 2D test suite (" << configFiles2D.size() << " tests)...\n" << endl;

            TestSuite suite2D( "New CVFEM 2D test suite", &cerr );

            const std::string meshPath2D = getMeshPath2D(benchmarkDir);

            std::vector<std::string> configPaths2D;
            configPaths2D.reserve( configFiles2D.size() );
            for ( const std::string& filename : configFiles2D )
            {
                configPaths2D.push_back( getTestDataPath_2D(benchmarkDir, filename) );
                // argv[0] is conventionally the program name (unused here); argv[1] = mesh, argv[2] = config file
                const char* argv2D[] = { "CVFEM_2D_VVCase", meshPath2D.c_str(), configPaths2D.back().c_str() };
                suite2D.addTest( filename.c_str(),
                                new CVFEM_2D_VVCase( 3, const_cast<char**>(argv2D), test_type ) );
            }

            // 3. Run, report, clean up
            suite2D.run();
            total_failures += suite2D.report();
            suite2D.free();
        }

        if ( run1D && !runAllTests_1D )
        {
            // Warn about any requested 1D test name that matched nothing at all — likely a typo.
            for ( const std::string& requested : requestedTests_1D )
            {
                const bool matches1D = std::find(ALL_TEST_FILES_1D.begin(), ALL_TEST_FILES_1D.end(), requested) != ALL_TEST_FILES_1D.end();
                if ( !matches1D )
                    cerr << "\nWARNING: requested 1D test \"" << requested << "\" does not match any known test — check for typos." << endl;
            }
        }
        if ( run2D && !runAllTests_2D )
        {
            // Warn about any requested 2D test name that matched nothing at all — likely a typo.
            for ( const std::string& requested : requestedTests_2D )
            {
                const bool matches2D = std::find(ALL_TEST_FILES_2D.begin(), ALL_TEST_FILES_2D.end(), requested) != ALL_TEST_FILES_2D.end();
                if ( !matches2D )
                    cerr << "\nWARNING: requested 2D test \"" << requested << "\" does not match any known test — check for typos." << endl;
            }
        }

        cerr << "\nCVFEM_tests_main: Total test failures: "
             << total_failures << endl;
    }


    // 4. Error handling
    catch ( const TestSuiteError& e ) {
        cerr << "\nTestSuiteError: " << e.what() << endl;
    }
    catch ( Exception& ba ) {
#ifdef __GNUC__
        const std::type_info& ti = typeid(ba);
        int status;
        char* realname = abi::__cxa_demangle( ti.name(), 0, 0, &status );
        cerr << "\nException: Exception raised by: " << realname << endl;
#else
        cout << "\nException: Exception raised by: " << typeid(ba).name() << endl;
#endif
        cerr << "\nDiagnostics:" << endl;
        ba.Out();
    }
#ifdef CSMP_WITH_SAMG_SOLVER
    catch ( SAMG_Exception& ba ) {
        cerr << "\nSAMG_Exception: " << ba.what() << endl;
    }
#endif
    catch ( const std::exception& e ) {
        cerr << "\nstd::exception: " << e.what() << endl;
    }

#ifdef CSMP_WITH_PETSC_SOLVER
    PetscFinalize();
#endif

    std::exit( total_failures );
    return 0;

} // end main