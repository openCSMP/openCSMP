// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef QUADRATIC_PRESSURE_PARALLEL_PLATE_PERMEABILITY_EXAMPLE_H
#define QUADRATIC_PRESSURE_PARALLEL_PLATE_PERMEABILITY_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// pressure computation using FEM with quadratic basis functions and the parallel-plate model for fracture permeability
class  QuadraticPressure_parallelPlatePermeability_Example : public Example {

  public:
  
    virtual void Run();
    virtual void Specifications();
    
};

} // csmp

#endif // QUADRATIC_PRESSURE_PARALLEL_PLATE_PERMEABILITY_EXAMPLE_H
