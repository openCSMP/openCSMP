// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_LINEAR_PRISM_TEST_H
#define ISOPARAMETRIC_LINEAR_PRISM_TEST_H

#include "Test.h"

namespace csmp {

/**
    Unit test of the 6-noded isoparametric linear prism (wedge) finite element.

    The checks rely on properties that every correct implementation must have, not on
    numbers recorded from this implementation:
     - topology: numbers of nodes, faces and segments, face types and face nodes
     - interpolation functions: value 1 at the own node and 0 at the other nodes, partition of unity,
       derivatives summing to zero
     - numerical integration: weights add up to the volume of the reference element, integration
       points lie inside of it
     - volume of the reference prism and of scaled, translated and sheared prisms
     - derivatives in global coordinates reproduce the gradient of a linear field
     - unit normals to the faces point outward (formerly FiniteElement_Test2)
*/
class IsoparametricLinearPrism_Test : public Test
{
  public:
    IsoparametricLinearPrism_Test() = default;

    void run() override;

  private:
    void Topology_Test();
    void InterpolationFunctions_Test();
    void NumericalIntegration_Test();
    void Volume_Test();
    void GlobalDerivatives_Test();
    void OutwardPointingFaceNormals_Test();
};

}

#endif // ISOPARAMETRIC_LINEAR_PRISM_TEST_H
