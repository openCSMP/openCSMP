// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "GenericSingleton_Test.h"

using namespace std;

namespace csmp {

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
