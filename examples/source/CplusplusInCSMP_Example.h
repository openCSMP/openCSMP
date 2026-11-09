// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CPLUSPLUS_IN_CSMP_EXAMPLE_H
#define CPLUSPLUS_IN_CSMP_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Introductory example demonstrating the C++ language features
           used throughout CSMP simulation code.

    No mesh or file input is required. All demonstrations use
    in-memory objects only.
*/
class CplusplusInCSMP_Example final : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~CplusplusInCSMP_Example() override = default;
};

} // csmp

#endif // CPLUSPLUS_IN_CSMP_EXAMPLE_H
