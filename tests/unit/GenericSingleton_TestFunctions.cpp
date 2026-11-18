// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "GenericSingleton_Test.h"

using namespace std;

namespace csmp {

void accessSingletonGloballyExternalFile()
{
  TestSingleton& single( TestSingleton::Instance() );
  single.i = 500;
}

} // csmp
