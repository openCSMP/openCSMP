// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ATMOSPHERIC_PRESSURE_DIFFUSION_EXAMPLE_H
#define ATMOSPHERIC_PRESSURE_DIFFUSION_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive introduction to transient pressure diffusion,
           verified against an analytical solution in 1D before
           applying to a 2D field-scale geological cross-section.

    Demonstrates how atmospheric pressure fluctuations diffuse into
    rock, illustrating the role of hydraulic diffusivity and
    storativity across contrasting geological formations.
*/
class AtmosphericPressureDiffusion_Example : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~AtmosphericPressureDiffusion_Example() override = default;
};

} // csmp

#endif // ATMOSPHERIC_PRESSURE_DIFFUSION_EXAMPLE_H

