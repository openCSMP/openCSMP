// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SOLVER_SETTINGS_H
#define SOLVER_SETTINGS_H

#include "CSMP_definitions.h"

// TODO: deprecate since this is redundant and overly prescriptive because not all solvers have settings

namespace csmp {

class SolverSettings {
public:
    SolverSettings();
    virtual ~SolverSettings();
};

} // end namespace csmp

#endif //SOLVERSETTINGS_H
