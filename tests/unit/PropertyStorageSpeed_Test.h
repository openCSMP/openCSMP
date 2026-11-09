// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PROPERTY_STORAGE_SPEED_TEST_H
#define CSMP_PROPERTY_STORAGE_SPEED_TEST_H

#include "Test.h"
#include "PropertyDatabase.h"
#include "LocalVariableStorage.h"

namespace csmp {

/// speed benchmark for read-write operations
/**
@todo This is actually a test case, plus it contains no test criteria/statements
*/
class  PropertyStorageSpeed_Test : public Test {
  public:
    explicit PropertyStorageSpeed_Test( std::ostream* osptr );
    virtual void run();
    
    const bool verbose_ = false; // should be off to get diagnostic value
};

} // end csmp

#endif
