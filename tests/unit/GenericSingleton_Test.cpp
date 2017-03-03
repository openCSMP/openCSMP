#include "GenericSingleton_Test.h"

using namespace std;

namespace csmp{

  void accessSingletonGloballyWithinFile();
  void accessSingletonGloballyExternalFile();

void GenericSingleton_Test::run()
{
    TestSingleton& single( TestSingleton::Instance() );

    single.i = 4;
    _test( single.i == 4 );

    accessSingletonGloballyWithinFile();
    _test( single.i == 1000 );

    accessSingletonGloballyExternalFile();
    _test( single.i == 500 );

}

void accessSingletonGloballyWithinFile()
{
  TestSingleton& single( TestSingleton::Instance() );
  single.i = 1000;
}

} // csmp
