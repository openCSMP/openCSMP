#include "TDDUnitTest_Example.h"

#include "TestSuite.h"

using namespace std;

namespace csmp {

  // implementation follos
  class Unit_c
  {
  public:
    int foo1( int a ) { return a/a; }
    int foo2( int a ) { return a*a; }
  };

  // this comes first(failing)
  class  UnitTest_c : public Test
  {
  public:
    virtual void run() { _test( Unit_c().foo1( 5 ) == 1 );
                         _test( Unit_c().foo2( 5 ) > 5 ); }
  };

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
  TestSuite s("TDD SUITE", &cout);

  s.addTest( new UnitTest_c() );

  s.run();
  s.report();

} // Run()

} // csmp
