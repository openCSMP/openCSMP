#include "GenericSingleton_Test.h"

using namespace std;

namespace csmp{

void accessSingletonGloballyExternalFile()
{
  TestSingleton& single( TestSingleton::Instance() );
  single.i = 500;
}

} // csmp
