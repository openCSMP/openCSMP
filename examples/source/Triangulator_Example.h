// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_TRIANGULATOR_EXAMPLE_H
#define CSMP_TRIANGULATOR_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
   Use Triangulator to make 2D mesh and perform a simple pressure diffusion computation on it.
*/
class  Triangulator_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // CSMP_TRIANGULATOR_EXAMPLE_H
