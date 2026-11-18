// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PHYSICAL_VARIABLES_EXAMPLE_H
#define CSMP_PHYSICAL_VARIABLES_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive introduction to CSMP physical variable types,
           their flags, placements, interpolation, extrapolation,
           and assignment to DenseMatrix objects.
*/
class CSMP_PhysicalVariables_Example : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~CSMP_PhysicalVariables_Example() override = default;
};

} // csmp

#endif // CSMP_PHYSICAL_VARIABLES_EXAMPLE_H

