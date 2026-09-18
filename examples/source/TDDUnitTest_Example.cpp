// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TDDUnitTest_Example.h"

#include "CSMP_definitions.h"
#include "Exception.h"
#include "compareFloats.h"
#include "TestSuite.h"
#include "Test.h"

using namespace std;

namespace csmp {

// class to test
class Unit_c {
  public:
    int Foo1( int a ) const { return a/a; }
    int Foo2( int a ) const { return a*a; }
  };

  // @note this is a Test not an example
class FooFunction_Test : public Test {
  public:
    ~FooFunction_Test() = default;
    void run() override final { TestFunctions(); }

  private:
    void TestFunctions() {
         Unit_c sut;
         _info("TestFunctions: Testing the member functions of class Unit_c...");
         _equal( sut.Foo1( 5 ), 1, 1.0e-6 );
         _test( approximatelyEqual( sut.Foo1(5), 1 ) );
         _test( sut.Foo2( 5 ) > 5 );
         // provoked failure
         _test( sut.Foo2( 5 ) < 20 );
      }
};

  void TDDUnitTest_Example::Specifications()
  {
    SetTitle( "Test-Driven Development=TDD of unit tests" );
    SetDifficulty( 1 );
    SetCategory( "C++" );
    AddAuthor( "SKM" );
    AddDescription( "source in: TDDUnitTest_Example.cpp" );
    AddDescription( "how to use a unit test for test driven development" );
  }




/**
    Demonstrates the testing macros defined in 'Test.h'
*/
void TDDUnitTest_Example::Run()
  {
     cout <<"\n"<<"Verification of the functions foo1() and foo2(): running tests...";
     
     TestSuite example_suite("CSMP examples unit-test suite (analogous to testing mechanism used in tests/unit/", &cout );

     example_suite.addTest( new FooFunction_Test() );

     // running the test
     example_suite.run();
     long nFail = example_suite.report();
     example_suite.free();
     cout << "\n"<<"TDDUnitTest_Example::Run: Total unit test failures: " << nFail << endl;
  }

} // csmp
