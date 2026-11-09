// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FINITE_ELEMENT_TEST2_H
#define FINITE_ELEMENT_TEST2_H

#include "Test.h"

namespace csmp {

/// Misc tests relating to non-simplex elements
class FiniteElement_Test2 : public Test
{
  public:
    FiniteElement_Test2() = default;
    
    virtual void run();
    
  private:
    const bool verbose_ = true;
};

}

#endif // FINITE_ELEMENT_TEST2_H
