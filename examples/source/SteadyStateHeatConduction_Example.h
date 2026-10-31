// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STEADY_STATE_HEAT_CONDUCTION_EXAMPLE_H
#define STEADY_STATE_HEAT_CONDUCTION_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive introduction to solving a PDE with CSMP.

    Solves steady-state heat conduction in a layered geological model
    with Dirichlet (fixed temperature) and Neumann (heat flux) boundary
    conditions. No equation of state or transient effects.
*/
class SteadyStateHeatConduction_Example : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~SteadyStateHeatConduction_Example() override = default;
};

} // csmp

#endif // STEADY_STATE_HEAT_CONDUCTION_EXAMPLE_H

