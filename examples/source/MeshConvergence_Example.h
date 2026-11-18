// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MESH_CONVERGENCE_EXAMPLE_H
#define MESH_CONVERGENCE_EXAMPLE_H

#include "Example.h"

namespace csmp {

/**
    @brief Interactive demonstration of FEM mesh convergence using
           the method of manufactured solutions.

    Solves the Poisson equation with a known exact solution on
    successively refined meshes, measuring the L2 error norm at
    each refinement level for both linear and quadratic triangles.
*/
class MeshConvergence_Example : public Example {
public:
    void Run()            override;
    void Specifications() override;
    ~MeshConvergence_Example() override = default;
};

} // csmp

#endif // MESH_CONVERGENCE_EXAMPLE_H

