/*
 *  MohrCoulombFailure_Visitor_Test.h
 *  csmp_core
 */

#ifndef CSMP_MOHR_COULOMB_FAILURE_VISITOR_TEST_H
#define CSMP_MOHR_COULOMB_FAILURE_VISITOR_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

/**
    @brief Unit test for MohrCoulombFailure_Visitor.

    Builds a FracBox model, assigns overburden stress that increases
    linearly with depth, and verifies that the visitor correctly
    identifies failure at the model sides where deviatoric stress
    is highest.

    @par Test setup

    The FracBox model occupies a unit cube. Overburden stress is
    assigned at each element barycentre as:

        σ_zz = ρ·g·z   (vertical, compression positive)
        σ_xx = σ_yy = K₀·σ_zz   (horizontal, with K₀ = 0.5)

    where z is the depth below the model top. This produces a stress
    state that is purely isotropic at the model centre and increasingly
    deviatoric toward the sides, where the horizontal confinement is
    reduced relative to the vertical load.

    Material properties assigned:
    - friction angle: 30 degrees
    - cohesion: 1 MPa
    - tensile strength: 0.1 MPa
    - Biot alpha: 1.0
    - fluid pressure: 0 (drained)

    @par Expected outcomes

    - Elements near the top of the model (low overburden) should show
      tensile failure where σ_zz is small and σ_xx = σ_yy = K₀·σ_zz
      is insufficient to prevent tension.
    - Elements at depth with high deviatoric stress should show shear
      failure where F_mc ≥ 0.
    - Elements in the interior under moderate isotropic stress should
      remain stable (F_mc < 0).
*/
class MohrCoulombFailure_Visitor_Test : public Test {
  public:
    explicit MohrCoulombFailure_Visitor_Test( bool verbose = false );
    ~MohrCoulombFailure_Visitor_Test() = default;

    virtual void run();

  private:
    /// Assigns overburden stress to all elements from barycentre depth
    bool TestOverburdenStressAssignment();

    /// Visitor runs without throwing on a correctly initialised model
    bool TestVisitorRuns();

    /// Elements near the top show tensile failure (low confinement)
    bool TestTensileFailureAtTop();

    /// Elements at depth with high deviatoric stress show shear failure
    bool TestShearFailureAtDepth();

    /// Stable elements in the interior have F_mc < 0
    bool TestStableInterior();

    /// Failure variable is stored at the correct placement
    bool TestOutputPlacement();

  private:
    bool verbose_;
};

} // end csmp

#endif /* CSMP_MOHR_COULOMB_FAILURE_VISITOR_TEST_H */

