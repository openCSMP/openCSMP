// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_QUADRATIC_PYRAMID_TEST_H
#define ISOPARAMETRIC_QUADRATIC_PYRAMID_TEST_H

#include "Test.h"

namespace csmp {

/**
 @brief Unit test suite for IsoparametricQuadraticPyramid.

 Tests:
   1. TopologyTest                  — face/segment/node connectivity
   2. PartitionOfUnityTest          — sum(N_i) = 1, Kronecker delta
   3. NodeRecoveryTest              — ParametricToPhysical at node coords
   4. VolumeTest                    — reference and distorted volumes
   5. ConsistencyWithLinearPyramid  — quadratic reproduces linear results
   6. ExtrapolationTest             — linear field extrapolated exactly
*/
class IsoparametricQuadraticPyramid_Test : public Test
{
public:
    IsoparametricQuadraticPyramid_Test() = default;
    ~IsoparametricQuadraticPyramid_Test() override = default;

    void run() override;

private:
    bool TopologyTest();
    bool PartitionOfUnityTest();
    bool NodeRecoveryTest();
    bool VolumeTest();
    bool ConsistencyWithLinearPyramidTest();
    bool ExtrapolationTest();
    
    constexpr static bool verbose_ = true;
};

} // namespace csmp

#endif // ISOPARAMETRIC_QUADRATIC_PYRAMID_TEST_H

