// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef BUILDING_A_CSMP_MODEL_EXAMPLE_H
#define BUILDING_A_CSMP_MODEL_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive introduction to building and querying a CSMP Model.

    Loads a mesh from CSMP binary files, queries its structure,
    accesses regions and boundaries, reads and writes variables,
    and introduces PropertyHandle.
*/
class BuildingACSMPModel_Example final : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~BuildingACSMPModel_Example() override = default;
};

} // csmp

#endif // BUILDING_A_CSMP_MODEL_EXAMPLE_H
