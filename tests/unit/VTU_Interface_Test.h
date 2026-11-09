// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_VTU_INTERFACE_TEST_H
#define CSMP_VTU_INTERFACE_TEST_H

#include "Test.h"

namespace csmp {

class VTU_Interface_Test : public Test
{
  public:
    VTU_Interface_Test();

    virtual void run();
};

} // csmp

#endif // CSMP_VTU_INTERFACE_TEST_H
