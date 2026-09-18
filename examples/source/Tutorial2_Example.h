// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TUTORIAL2_EXAMPLE_H
#define TUTORIAL2_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// Classic advection - diffusion equation (ADE) solver
class  Tutorial2_Example : public Example
{

public:

    virtual void Run();
    virtual void Specifications();

};

} // csmp

#endif // TUTORIAL2_EXAMPLE_H
