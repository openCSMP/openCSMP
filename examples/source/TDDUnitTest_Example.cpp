#include "TDDUnitTest_Example.h"

// Use CATCH_CONFIG_MAIN if you want Catch to supply main().
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace std;

namespace csmp {

  // class to test
  class Unit_c
  {
  public:
    int foo1( int a ) { return a/a; }
    int foo2( int a ) { return a*a; }
  };

  // this comes first(failing)
  struct UnitTest_c
  {
      // system under test
      Unit_c sut;

      void test1() { CHECK(sut.foo1( 5 ) == 1); }
      void test2() { CHECK(sut.foo1( 5 ) > 5); }
  };

// Basic test
TEST_CASE("Test 1", "") {
    REQUIRE_NOTHROW(UnitTest_c().test1());
}

// Test with a fixture
TEST_CASE_METHOD(UnitTest_c, "Test 2", "[create]") {
    REQUIRE_NOTHROW(test1());
}

// Correctly failing test
TEST_CASE("Test 3", "[!shouldfail]") {
    REQUIRE_NOTHROW(UnitTest_c().test2());
}

void TDDUnitTest_Example::Specifications()
{
  SetTitle( "TDD - (unit) Test-Driven Development of software" );
  SetDifficulty( 1 );
  SetCategory( "C++" );
  AddAuthor( "P. Lang" );
  AddDescription( "source in: TDDUnitTest_Example.cpp" );
  AddDescription( "how to use a unit test for test driven development" );
}

void TDDUnitTest_Example::Run()
{
    // If you don't use CATCH_CONFIG_MAIN, this is the minimal
    // infrastructure required to run a test.

    Catch::Session session;
    session.run();
} // Run()

} // csmp
