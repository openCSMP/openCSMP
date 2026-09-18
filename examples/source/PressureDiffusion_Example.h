// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PRESSURE_DIFFUSION_EXAMPLE_H
#define PRESSURE_DIFFUSION_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  PressureDiffusion_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // PFDIFFUSION_EXAMPLE_H
