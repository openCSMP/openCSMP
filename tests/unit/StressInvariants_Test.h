/*
 *  StressInvariants_Test.h
 *  csmp_core
 *
 *  Created to test StressInvariants for correctness across a range
 *  of stress states including hydrostatic, uniaxial, pure shear,
 *  asymmetric, and rotated tensors.
 */

#ifndef CSMP_STRESS_INVARIANTS_TEST_H
#define CSMP_STRESS_INVARIANTS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

/**
    @brief Unit test for StressInvariants.

    Covers:
    - Hydrostatic stress state (Lode angle undefined — throws)
    - Uniaxial compression
    - Triaxial compression (Lode angle = +pi/6)
    - Triaxial extension  (Lode angle = -pi/6)
    - Pure shear          (Lode angle = 0)
    - 2D plane stress
    - 2D plane strain
    - Asymmetric tensor (off-diagonal terms non-zero)
    - Rotated tensor (principal stresses invariant under rotation)
    - Yield envelope evaluation (Mohr-Coulomb and smooth variant)
    - Sign convention: compression positive throughout
*/
class StressInvariants_Test : public Test {
  public:
    explicit StressInvariants_Test( bool verbose = false );
    ~StressInvariants_Test() = default;

    virtual void run();

  private:
    /// Hydrostatic state: all principal stresses equal, Lode angle undefined
    bool TestHydrostatic();

    /// Uniaxial compression: one non-zero principal stress
    bool TestUniaxialCompression();

    /// Triaxial compression: sigma2 == sigma3, Lode angle == +pi/6
    bool TestTriaxialCompression();

    /// Triaxial extension: sigma1 == sigma2, Lode angle == -pi/6
    bool TestTriaxialExtension();

    /// Pure shear: mean stress zero, Lode angle == 0
    bool TestPureShear();

    /// 2D plane stress: sigma_zz == 0
    bool TestPlaneStress2D();

    /// 2D plane strain: sigma_zz == nu*(sigma_xx + sigma_yy)
    bool TestPlaneStrain2D();

    /// Asymmetric off-diagonal terms: invariants must still be correct
    bool TestAsymmetricTensor();

    /// Rotated tensor: principal stresses invariant under coordinate rotation
    bool TestRotatedTensor();

    /// Mohr-Coulomb yield envelope at known Lode angles
    bool TestMohrCoulombYieldEnvelope();

    /// Smooth Mohr-Coulomb yield envelope
    bool TestSmoothMohrCoulombYieldEnvelope();

    /// Principal stress ordering: sigma1 >= sigma2 >= sigma3
    bool TestPrincipalStressOrdering();

  private:
    bool   verbose_;
    double tolerance_; ///< absolute tolerance for floating point comparisons
};

} // end csmp

#endif /* CSMP_STRESS_INVARIANTS_TEST_H */

