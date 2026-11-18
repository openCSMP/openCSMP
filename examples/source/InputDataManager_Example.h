
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INPUT_DATA_MANAGER_EXAMPLE_H
#define INPUT_DATA_MANAGER_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive introduction to model configuration via
           block-structured text files read by InputDataManager.

    Demonstrates all six configuration blocks, the PropertyDatabase
    variable range checking, whitespace rules, and ComputationalSettings.
*/
class InputDataManager_Example : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~InputDataManager_Example() override = default;
};

} // csmp

#endif // INPUT_DATA_MANAGER_EXAMPLE_H
