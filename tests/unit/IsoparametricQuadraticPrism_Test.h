// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_QUADRATIC_PRISM_TEST_H
#define ISOPARAMETRIC_QUADRATIC_PRISM_TEST_H

#include "Test.h"

namespace csmp {

/**
    Unit test of the 18-noded isoparametric quadratic prism (wedge) finite element:
    6 corner nodes, 9 nodes on the middle of the edges and 3 nodes in the centre of the quadrilateral faces.

    The checks rely on properties that every correct implementation must have, not on
    numbers recorded from this implementation:
     - topology: numbers of nodes, faces and segments, corner nodes of the faces
     - interpolation functions: value 1 at the own node and 0 at the other nodes, partition of unity,
       derivatives summing to zero and agreeing with finite differences, exact interpolation of a
       quadratic field
     - numerical integration: weights add up to the volume of the reference element, integration
       points lie inside of it
     - volume of the reference prism and of scaled, translated and sheared prisms
     - derivatives in global coordinates reproduce the gradient of a linear field

    @see IsoparametricLinearPrism_Test for the corresponding test of the linear element
*/
class IsoparametricQuadraticPrism_Test : public Test
{
  public:
    IsoparametricQuadraticPrism_Test() = default;

    void run() override;

  private:
    void Topology_Test();
    void InterpolationFunctions_Test();
    void NumericalIntegration_Test();
    void Volume_Test();
    void GlobalDerivatives_Test();
};

}

#endif // ISOPARAMETRIC_QUADRATIC_PRISM_TEST_H
