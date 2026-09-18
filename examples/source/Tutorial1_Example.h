// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TUTORIAL1_EXAMPLE_H
#define TUTORIAL1_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// transient pressure diffusion through a porous medium
class  Tutorial1_Example : public Example
{

public:

    virtual void Run();
    virtual void Specifications();

};

} // csmp

#endif // TUTORIAL1_EXAMPLE_H
