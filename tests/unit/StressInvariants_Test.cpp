/*
 *  StressInvariants_Test.cpp
 *  csmp_core
 */

#include "StressInvariants_Test.h"
#include "StressInvariants.h"
#include "TensorVariable.h"
#include "CSMP_physical_constants.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Construction
// ----------------------------------------------------------------------------

StressInvariants_Test::StressInvariants_Test( bool verbose )
    : verbose_  ( verbose  ),
      tolerance_( 1.0e-10  )
{
}


// ----------------------------------------------------------------------------
//  run
// ----------------------------------------------------------------------------

void StressInvariants_Test::run()
{
    _test( TestHydrostatic()                  );
    _test( TestUniaxialCompression()          );
    _test( TestTriaxialCompression()          );
    _test( TestTriaxialExtension()            );
    _test( TestPureShear()                    );
    _test( TestPlaneStress2D()                );
    _test( TestPlaneStrain2D()                );
    _test( TestAsymmetricTensor()             );
    _test( TestRotatedTensor()                );
    _test( TestMohrCoulombYieldEnvelope()     );
    _test( TestSmoothMohrCoulombYieldEnvelope());
    _test( TestPrincipalStressOrdering()      );
}


// ----------------------------------------------------------------------------
//  Helpers
// ----------------------------------------------------------------------------

namespace {

/// builds a 3D tensor from six independent components (compression positive)
TensorVariable<3U> makeTensor3( double sxx, double syy, double szz,
                                double sxy, double sxz, double syz )
{
    TensorVariable<3U> ts;
    ts(0,0) = sxx;  ts(1,1) = syy;  ts(2,2) = szz;
    ts(0,1) = sxy;  ts(1,0) = sxy;
    ts(0,2) = sxz;  ts(2,0) = sxz;
    ts(1,2) = syz;  ts(2,1) = syz;
    return ts;
}

/// builds a 2D tensor from three independent components (compression positive)
TensorVariable<2U> makeTensor2( double sxx, double syy, double sxy )
{
    TensorVariable<2U> ts;
    ts(0,0) = sxx;  ts(1,1) = syy;
    ts(0,1) = sxy;  ts(1,0) = sxy;
    return ts;
}

/// returns true if |a - b| <= tol
bool nearlyEqual( double a, double b, double tol )
{
    return std::abs( a - b ) <= tol;
}

} // anonymous namespace


// ----------------------------------------------------------------------------
//  TestHydrostatic
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestHydrostatic()
{
    const double p = 1.0e6;
    const auto ts = makeTensor3( p, p, p, 0., 0., 0. );

    // constructor must not throw for hydrostatic state
    bool threw = false;
    StressInvariants* si_ptr = nullptr;
    try {
        si_ptr = new StressInvariants( ts );
    }
    catch ( ... ) {
        threw = true;
    }
    _test( threw == false );

    if ( si_ptr != nullptr ) {
        // mean stress must equal p
        _test( nearlyEqual( si_ptr->MeanStress(), p, tolerance_ ) );

        // deviatoric stress must be zero
        _test( nearlyEqual( si_ptr->DeviatoricStress(), 0., tolerance_ ) );

        // Lode angle returns sentinel value 0 for hydrostatic state
        _test( nearlyEqual( si_ptr->LodeAngle( ts ), 0., tolerance_ ) );

        // all principal stresses must equal p
        _test( nearlyEqual( si_ptr->MaximumPrincipalStress1(),      p, tolerance_ ) );
        _test( nearlyEqual( si_ptr->IntermediatePrincipalStress2(), p, tolerance_ ) );
        _test( nearlyEqual( si_ptr->LeastPrincipalStress3(),        p, tolerance_ ) );

        if ( verbose_ )
            cout << "TestHydrostatic:"
                 << " p="     << si_ptr->MeanStress()
                 << " q="     << si_ptr->DeviatoricStress()
                 << " theta=" << si_ptr->LodeAngle( ts ) << "\n";

        delete si_ptr;
    }

    // perturbed state: mean stress still close to p
    const double eps = 1.0e-3;
    const auto ts_perturbed = makeTensor3( p + eps, p, p - eps, 0., 0., 0. );
    StressInvariants si_perturbed( ts_perturbed );
    _test( nearlyEqual( si_perturbed.MeanStress(), p, 1.0e-6 ) );

    return true;
}


// ----------------------------------------------------------------------------
//  TestUniaxialCompression
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestUniaxialCompression()
{
    const double q = 3.0e6;
    const auto ts = makeTensor3( q, 0., 0., 0., 0., 0. );
    StressInvariants si( ts );

    // use a scale-relative tolerance for values computed via trig functions
    const double tol = tolerance_ * q;

    _test( nearlyEqual( si.MeanStress(), q / 3., tol ) );
    _test( nearlyEqual( si.MaximumPrincipalStress1(),      q,  tol ) );
    _test( nearlyEqual( si.IntermediatePrincipalStress2(), 0., tol ) );
    _test( nearlyEqual( si.LeastPrincipalStress3(),        0., tol ) );

    if ( verbose_ )
        cout << "TestUniaxialCompression:"
             << " sigma1=" << si.MaximumPrincipalStress1()
             << " sigma2=" << si.IntermediatePrincipalStress2()
             << " sigma3=" << si.LeastPrincipalStress3()
             << " p="      << si.MeanStress() << "\n";

    return true;
}


bool StressInvariants_Test::TestTriaxialCompression()
{
    const double s1 = 5.0e6;
    const double s3 = 2.0e6;
    const auto ts = makeTensor3( s1, s3, s3, 0., 0., 0. );
    StressInvariants si( ts );

    const double tol = tolerance_ * s1;

    _test( nearlyEqual( si.LodeAngle( ts ), -CSMP_PI / 6., tolerance_ ) );
    _test( nearlyEqual( si.MaximumPrincipalStress1(),      s1, tol ) );
    _test( nearlyEqual( si.IntermediatePrincipalStress2(), s3, tol ) );
    _test( nearlyEqual( si.LeastPrincipalStress3(),        s3, tol ) );

    if ( verbose_ )
        cout << "TestTriaxialCompression:"
             << " theta=" << si.LodeAngle( ts )
             << " sigma1=" << si.MaximumPrincipalStress1()
             << " sigma2=" << si.IntermediatePrincipalStress2()
             << " sigma3=" << si.LeastPrincipalStress3() << "\n";

    return true;
}


bool StressInvariants_Test::TestTriaxialExtension()
{
    const double s1 = 5.0e6;
    const double s3 = 2.0e6;
    const auto ts = makeTensor3( s1, s1, s3, 0., 0., 0. );
    StressInvariants si( ts );

    const double tol = tolerance_ * s1;

    _test( nearlyEqual( si.LodeAngle( ts ), +CSMP_PI / 6., tolerance_ ) );
    _test( nearlyEqual( si.MaximumPrincipalStress1(),      s1, tol ) );
    _test( nearlyEqual( si.IntermediatePrincipalStress2(), s1, tol ) );
    _test( nearlyEqual( si.LeastPrincipalStress3(),        s3, tol ) );

    if ( verbose_ )
        cout << "TestTriaxialExtension:"
             << " theta=" << si.LodeAngle( ts )
             << " sigma1=" << si.MaximumPrincipalStress1()
             << " sigma2=" << si.IntermediatePrincipalStress2()
             << " sigma3=" << si.LeastPrincipalStress3() << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestPureShear
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestPureShear()
{
    // Pure shear: sigma_xx = -sigma_yy = tau, sigma_zz = 0
    // Mean stress = 0, Lode angle = 0
    const double tau = 2.0e6;
    const auto ts = makeTensor3( tau, -tau, 0., 0., 0., 0. );
    StressInvariants si( ts );

    _test( nearlyEqual( si.MeanStress(),   0., tolerance_ ) );
    _test( nearlyEqual( si.LodeAngle( ts ), 0., tolerance_ ) );

    if ( verbose_ )
        cout << "TestPureShear:"
             << " p="     << si.MeanStress()
             << " theta=" << si.LodeAngle( ts ) << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestPlaneStress2D
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestPlaneStress2D()
{
    // Plane stress: sigma_zz = 0
    // Compare 2D result with equivalent 3D result (szz = 0)
    const double sxx = 4.0e6;
    const double syy = 2.0e6;
    const double sxy = 1.0e6;

    const auto ts2 = makeTensor2( sxx, syy, sxy );
    const auto ts3 = makeTensor3( sxx, syy, 0., sxy, 0., 0. );

    StressInvariants si2( ts2, PlaneAssumption::PLANE_STRESS );
    StressInvariants si3( ts3 );

    // mean stress and deviatoric stress must match
    _test( nearlyEqual( si2.MeanStress(),       si3.MeanStress(),       tolerance_ ) );
    _test( nearlyEqual( si2.DeviatoricStress(), si3.DeviatoricStress(), tolerance_ ) );

    // principal stresses must match
    _test( nearlyEqual( si2.MaximumPrincipalStress1(),
                        si3.MaximumPrincipalStress1(),      tolerance_ ) );
    _test( nearlyEqual( si2.IntermediatePrincipalStress2(),
                        si3.IntermediatePrincipalStress2(), tolerance_ ) );
    _test( nearlyEqual( si2.LeastPrincipalStress3(),
                        si3.LeastPrincipalStress3(),        tolerance_ ) );

    if ( verbose_ )
        cout << "TestPlaneStress2D:"
             << " p2D="  << si2.MeanStress()
             << " p3D="  << si3.MeanStress()
             << " q2D="  << si2.DeviatoricStress()
             << " q3D="  << si3.DeviatoricStress() << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestPlaneStrain2D
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestPlaneStrain2D()
{
    // Plane strain: sigma_zz = nu * (sigma_xx + sigma_yy)
    // Compare 2D plane strain result with equivalent 3D result
    const double sxx = 4.0e6;
    const double syy = 2.0e6;
    const double sxy = 1.0e6;
    const double nu  = 0.3;
    const double szz = nu * ( sxx + syy );

    const auto ts2 = makeTensor2( sxx, syy, sxy );
    const auto ts3 = makeTensor3( sxx, syy, szz, sxy, 0., 0. );

    StressInvariants si2( ts2, PlaneAssumption::PLANE_STRAIN, nu );
    StressInvariants si3( ts3 );

    _test( nearlyEqual( si2.MeanStress(),       si3.MeanStress(),       tolerance_ ) );
    _test( nearlyEqual( si2.DeviatoricStress(), si3.DeviatoricStress(), tolerance_ ) );

    _test( nearlyEqual( si2.MaximumPrincipalStress1(),
                        si3.MaximumPrincipalStress1(),      tolerance_ ) );
    _test( nearlyEqual( si2.IntermediatePrincipalStress2(),
                        si3.IntermediatePrincipalStress2(), tolerance_ ) );
    _test( nearlyEqual( si2.LeastPrincipalStress3(),
                        si3.LeastPrincipalStress3(),        tolerance_ ) );

    // plane strain mean stress must be larger than plane stress mean stress
    // because sigma_zz > 0 adds to I1
    StressInvariants si2_ps( ts2, PlaneAssumption::PLANE_STRESS );
    _test( si2.MeanStress() > si2_ps.MeanStress() );

    if ( verbose_ )
        cout << "TestPlaneStrain2D:"
             << " szz="         << szz
             << " p_strain="    << si2.MeanStress()
             << " p_stress="    << si2_ps.MeanStress() << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestAsymmetricTensor
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestAsymmetricTensor()
{
    // A tensor with significant off-diagonal shear terms.
    // The invariants must still satisfy known relationships:
    //   sigma1 + sigma2 + sigma3 = sigma_xx + sigma_yy + sigma_zz = I1
    //   3 * MeanStress = I1
    const double sxx = 5.0e6;
    const double syy = 3.0e6;
    const double szz = 1.0e6;
    const double sxy = 2.0e6;
    const double sxz = 1.5e6;
    const double syz = 0.5e6;

    const auto ts = makeTensor3( sxx, syy, szz, sxy, sxz, syz );
    StressInvariants si( ts );

    // sum of principal stresses must equal I1 = sxx + syy + szz
    const double I1 = sxx + syy + szz;
    const double sum_principal = si.MaximumPrincipalStress1()
                               + si.IntermediatePrincipalStress2()
                               + si.LeastPrincipalStress3();

    _test( nearlyEqual( sum_principal, I1, 1.0e-6 ) );

    // mean stress = I1/3
    _test( nearlyEqual( si.MeanStress(), I1 / 3., tolerance_ ) );

    // Lode angle must be in [-pi/6, pi/6]
    _test( si.LodeAngle( ts ) >= -CSMP_PI / 6. - tolerance_ );
    _test( si.LodeAngle( ts ) <=  CSMP_PI / 6. + tolerance_ );

    if ( verbose_ )
        cout << "TestAsymmetricTensor:"
             << " I1="            << I1
             << " sum_principal=" << sum_principal
             << " theta="         << si.LodeAngle( ts ) << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestRotatedTensor
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestRotatedTensor()
{
    // Rotate a diagonal stress tensor by 45 degrees about the z-axis.
    // The principal stresses and invariants must be unchanged.
    //
    // Original: sigma1 = 5 MPa, sigma2 = 3 MPa, sigma3 = 1 MPa (diagonal)
    // After 45-degree rotation about z:
    //   sigma_xx' = (sigma1 + sigma2)/2 = 4 MPa
    //   sigma_yy' = (sigma1 + sigma2)/2 = 4 MPa
    //   sigma_xy' = (sigma1 - sigma2)/2 = 1 MPa
    //   sigma_zz' = sigma3 = 1 MPa (unchanged)
    const double s1 = 5.0e6;
    const double s2 = 3.0e6;
    const double s3 = 1.0e6;

    // original diagonal tensor
    const auto ts_diag = makeTensor3( s1, s2, s3, 0., 0., 0. );
    StressInvariants si_diag( ts_diag );

    // rotated tensor (45 degrees about z)
    const double sxx_rot = ( s1 + s2 ) / 2.;
    const double syy_rot = ( s1 + s2 ) / 2.;
    const double szz_rot = s3;
    const double sxy_rot = ( s1 - s2 ) / 2.;

    const auto ts_rot = makeTensor3( sxx_rot, syy_rot, szz_rot,
                                     sxy_rot, 0., 0. );
    StressInvariants si_rot( ts_rot );

    // invariants must be identical
    _test( nearlyEqual( si_diag.MeanStress(),
                        si_rot.MeanStress(),       tolerance_ ) );
    _test( nearlyEqual( si_diag.DeviatoricStress(),
                        si_rot.DeviatoricStress(), tolerance_ ) );
    _test( nearlyEqual( si_diag.LodeAngle( ts_diag ),
                        si_rot.LodeAngle( ts_rot ), tolerance_ ) );

    // principal stresses must be identical
    _test( nearlyEqual( si_diag.MaximumPrincipalStress1(),
                        si_rot.MaximumPrincipalStress1(),      tolerance_ ) );
    _test( nearlyEqual( si_diag.IntermediatePrincipalStress2(),
                        si_rot.IntermediatePrincipalStress2(), tolerance_ ) );
    _test( nearlyEqual( si_diag.LeastPrincipalStress3(),
                        si_rot.LeastPrincipalStress3(),        tolerance_ ) );

    if ( verbose_ )
        cout << "TestRotatedTensor:"
             << " p_diag="  << si_diag.MeanStress()
             << " p_rot="   << si_rot.MeanStress()
             << " q_diag="  << si_diag.DeviatoricStress()
             << " q_rot="   << si_rot.DeviatoricStress() << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestMohrCoulombYieldEnvelope
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestMohrCoulombYieldEnvelope()
{
    // Mohr-Coulomb yield envelope:
    //
    //   K(θ) = cos θ - sin φ · sin θ / √3
    //
    // Sign convention (compression positive, geotechnical):
    //   θ = -π/6 at triaxial compression (σ₁ > σ₂ = σ₃)
    //   θ = +π/6 at triaxial extension   (σ₁ = σ₂ > σ₃)
    //
    // At triaxial compression (θ = -π/6):
    //   K = cos(-π/6) - sin φ · sin(-π/6) / √3
    //     = √3/2 + sin φ / (2√3)   > 1 for φ > 0
    //
    // At triaxial extension (θ = +π/6):
    //   K = cos(π/6) - sin φ · sin(π/6) / √3
    //     = √3/2 - sin φ / (2√3)   < 1 for φ > 0
    //
    // Therefore K_compression > K_extension for φ > 0:
    // the Mohr-Coulomb surface is stronger in compression than extension.

    const double phi     = 30.;
    const double sin_phi = std::sin( phi * CSMP_PI / 180. );
    const double sqrt3   = std::sqrt( 3. );

    // --- triaxial compression: σ₁ > σ₂ = σ₃, θ = -π/6 ---
    const double s1_c = 5.0e6;
    const double s3_c = 2.0e6;
    const auto   ts_c = makeTensor3( s1_c, s3_c, s3_c, 0., 0., 0. );
    StressInvariants si_c( ts_c );

    // verify Lode angle is -π/6
    _test( nearlyEqual( si_c.LodeAngle( ts_c ),
                        -CSMP_PI / 6.,
                        tolerance_ ) );

    const double K_compression_expected = sqrt3 / 2.
                                        + sin_phi / ( 2. * sqrt3 );
    const double K_compression = si_c.MohrCoulombYieldEnvelope( phi );

    _test( nearlyEqual( K_compression, K_compression_expected, tolerance_ ) );

    if ( verbose_ )
        cout << "TestMohrCoulombYieldEnvelope:"
             << " K_compression="     << K_compression
             << " (expected "         << K_compression_expected << ")"
             << " theta_compression=" << si_c.LodeAngle( ts_c )
             << " (expected "         << -CSMP_PI / 6. << ")\n";

    // --- triaxial extension: σ₁ = σ₂ > σ₃, θ = +π/6 ---
    const double s1_e = 5.0e6;
    const double s3_e = 2.0e6;
    const auto   ts_e = makeTensor3( s1_e, s1_e, s3_e, 0., 0., 0. );
    StressInvariants si_e( ts_e );

    // verify Lode angle is +π/6
    _test( nearlyEqual( si_e.LodeAngle( ts_e ),
                        +CSMP_PI / 6.,
                        tolerance_ ) );

    const double K_extension_expected = sqrt3 / 2.
                                      - sin_phi / ( 2. * sqrt3 );
    const double K_extension = si_e.MohrCoulombYieldEnvelope( phi );

    _test( nearlyEqual( K_extension, K_extension_expected, tolerance_ ) );

    if ( verbose_ )
        cout << "TestMohrCoulombYieldEnvelope:"
             << " K_extension="     << K_extension
             << " (expected "       << K_extension_expected << ")"
             << " theta_extension=" << si_e.LodeAngle( ts_e )
             << " (expected "       << +CSMP_PI / 6. << ")\n";

    // --- K_compression > K_extension for φ > 0 ---
    _test( K_compression > K_extension );

    // --- K_compression > 1 for φ > 0 ---
    // the hexagonal Mohr-Coulomb surface extends beyond the unit circle
    // at triaxial compression
    _test( K_compression > 1. );

    // --- K_extension < 1 for φ > 0 ---
    _test( K_extension < 1. );

    // --- both must be positive ---
    _test( K_compression > 0. );
    _test( K_extension   > 0. );

    // --- at φ = 0 (frictionless): K = cos θ regardless of stress state ---
    // cos(-π/6) = cos(+π/6) = √3/2 since cosine is even
    _test( nearlyEqual( si_c.MohrCoulombYieldEnvelope( 0. ),
                        std::cos( CSMP_PI / 6. ),
                        tolerance_ ) );
    _test( nearlyEqual( si_e.MohrCoulombYieldEnvelope( 0. ),
                        std::cos( CSMP_PI / 6. ),
                        tolerance_ ) );
    // both equal √3/2
    _test( nearlyEqual( si_c.MohrCoulombYieldEnvelope( 0. ),
                        si_e.MohrCoulombYieldEnvelope( 0. ),
                        tolerance_ ) );

    // --- pure shear: θ = 0, K = cos(0) - sin φ · sin(0) / √3 = 1 ---
    const double tau   = 2.0e6;
    const auto   ts_ps = makeTensor3( tau, -tau, 0., 0., 0., 0. );
    StressInvariants si_ps( ts_ps );

    _test( nearlyEqual( si_ps.LodeAngle( ts_ps ), 0., tolerance_ ) );
    _test( nearlyEqual( si_ps.MohrCoulombYieldEnvelope( phi ),
                        1.,
                        tolerance_ ) );

    if ( verbose_ )
        cout << "TestMohrCoulombYieldEnvelope:"
             << " K_pure_shear=" << si_ps.MohrCoulombYieldEnvelope( phi )
             << " (expected 1.0)\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestSmoothMohrCoulombYieldEnvelope
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestSmoothMohrCoulombYieldEnvelope()
{
    // Smooth Mohr-Coulomb:
    //   K_param = (3 - sin φ) / (3 + sin φ)
    //   G(θ) = 2·K_param / [(1+K_param) - sin(3θ)·(1-K_param)]
    //   K(θ) = 1/G(θ)
    //
    // At triaxial compression θ = -π/6, sin(3θ) = sin(-π/2) = -1:
    //   G = 2·K_param / [(1+K_param) + (1-K_param)] = 2·K_param / 2 = K_param
    //   K(θ) = 1/K_param = (3 + sin φ) / (3 - sin φ)  > 1 for φ > 0
    //
    // At triaxial extension θ = +π/6, sin(3θ) = sin(+π/2) = +1:
    //   G = 2·K_param / [(1+K_param) - (1-K_param)] = 2·K_param / 2·K_param = 1
    //   K(θ) = 1/G = 1
    //
    // Therefore K_smooth_compression > K_smooth_extension for φ > 0.

    const double phi     = 30.;
    const double sin_phi = std::sin( phi * CSMP_PI / 180. );
    const double K_param = ( 3. - sin_phi ) / ( 3. + sin_phi );

    // triaxial compression: θ = -π/6
    const auto ts_c = makeTensor3( 5.0e6, 2.0e6, 2.0e6, 0., 0., 0. );
    StressInvariants si_c( ts_c );

    _test( nearlyEqual( si_c.LodeAngle( ts_c ), -CSMP_PI / 6., tolerance_ ) );

    const double K_smooth_compression          = si_c.SmoothMohrCoulombYieldEnvelope( phi );
    const double K_smooth_compression_expected = 1. / K_param;  // > 1

    _test( nearlyEqual( K_smooth_compression, K_smooth_compression_expected, tolerance_ ) );

    if ( verbose_ )
        cout << "TestSmoothMohrCoulombYieldEnvelope:"
             << " K_compression=" << K_smooth_compression
             << " (expected "     << K_smooth_compression_expected << ")\n";

    // triaxial extension: θ = +π/6
    const auto ts_e = makeTensor3( 5.0e6, 5.0e6, 2.0e6, 0., 0., 0. );
    StressInvariants si_e( ts_e );

    _test( nearlyEqual( si_e.LodeAngle( ts_e ), +CSMP_PI / 6., tolerance_ ) );

    const double K_smooth_extension          = si_e.SmoothMohrCoulombYieldEnvelope( phi );
    const double K_smooth_extension_expected = 1.;  // exactly 1 at extension

    _test( nearlyEqual( K_smooth_extension, K_smooth_extension_expected, tolerance_ ) );

    if ( verbose_ )
        cout << "TestSmoothMohrCoulombYieldEnvelope:"
             << " K_extension=" << K_smooth_extension
             << " (expected "   << K_smooth_extension_expected << ")\n";

    // K_compression > K_extension for φ > 0
    _test( K_smooth_compression > K_smooth_extension );

    // at φ = 0: K_param = 1, G = 2/(2 - 0) = 1, K(θ) = 1 everywhere
    _test( nearlyEqual( si_c.SmoothMohrCoulombYieldEnvelope( 0. ), 1., tolerance_ ) );
    _test( nearlyEqual( si_e.SmoothMohrCoulombYieldEnvelope( 0. ), 1., tolerance_ ) );

    // pure shear: θ = 0, sin(3θ) = 0
    //   G = 2·K_param / (1+K_param) => K(θ) = (1+K_param)/(2·K_param)
    const double tau   = 2.0e6;
    const auto   ts_ps = makeTensor3( tau, -tau, 0., 0., 0., 0. );
    StressInvariants si_ps( ts_ps );

    _test( nearlyEqual( si_ps.LodeAngle( ts_ps ), 0., tolerance_ ) );

    const double K_smooth_shear_expected = ( 1. + K_param ) / ( 2. * K_param );
    _test( nearlyEqual( si_ps.SmoothMohrCoulombYieldEnvelope( phi ),
                        K_smooth_shear_expected,
                        tolerance_ ) );

    if ( verbose_ )
        cout << "TestSmoothMohrCoulombYieldEnvelope:"
             << " K_pure_shear=" << si_ps.SmoothMohrCoulombYieldEnvelope( phi )
             << " (expected "    << K_smooth_shear_expected << ")\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestPrincipalStressOrdering
// ----------------------------------------------------------------------------

bool StressInvariants_Test::TestPrincipalStressOrdering()
{
    // For any stress state, sigma1 >= sigma2 >= sigma3 must hold
    // (compression positive, so sigma1 is the most compressive)

    // general asymmetric state
    {
        const auto ts = makeTensor3( 5.0e6, 3.0e6, 1.0e6,
                                     2.0e6, 1.5e6, 0.5e6 );
        StressInvariants si( ts );
        _test( si.MaximumPrincipalStress1()
            >= si.IntermediatePrincipalStress2() - tolerance_ );
        _test( si.IntermediatePrincipalStress2()
            >= si.LeastPrincipalStress3() - tolerance_ );
    }

    // state with tensile (negative) least principal stress
    {
        const auto ts = makeTensor3( 3.0e6, 1.0e6, -1.0e6, 0., 0., 0. );
        StressInvariants si( ts );
        _test( si.MaximumPrincipalStress1()
            >= si.IntermediatePrincipalStress2() - tolerance_ );
        _test( si.IntermediatePrincipalStress2()
            >= si.LeastPrincipalStress3() - tolerance_ );
        // least principal stress must be negative (tensile)
        _test( si.LeastPrincipalStress3() < 0. );
    }

    // uniaxial tension (negative compression)
    {
        const auto ts = makeTensor3( -2.0e6, 0., 0., 0., 0., 0. );
        StressInvariants si( ts );
        _test( si.MaximumPrincipalStress1()
            >= si.IntermediatePrincipalStress2() - tolerance_ );
        _test( si.IntermediatePrincipalStress2()
            >= si.LeastPrincipalStress3() - tolerance_ );
    }

    if ( verbose_ )
        cout << "TestPrincipalStressOrdering: all ordering checks passed\n";

    return true;
}

} // end csmp

