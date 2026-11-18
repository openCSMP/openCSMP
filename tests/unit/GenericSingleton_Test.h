// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GENERIC_SINGLETON_TEST_H
#define GENERIC_SINGLETON_TEST_H

#include "Test.h"
#include "GenericSingleton.h"

namespace csmp {

class GenericSingleton_Test : public Test
{
public:
  virtual void run();
};


 

class TestSingleton : public GenericSingleton<TestSingleton>
{
  friend class GenericSingleton<TestSingleton>;
 
  private:
    TestSingleton() : i(3) {}

  public:
    int i;
};

void accessSingletonGloballyWithinFile();
void accessSingletonGloballyExternalFile();

} // csmp

#endif // GENERIC_SINGLETON_TEST_H
