/*
 *  StressInvariants.cpp
 *
 *  Created by Stephan Matthai on 10/25/14.
 *  Revised and corrected.
 */

#include <cmath>
#include <limits>
#include "StressInvariants.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

// ============================================================================
//  Internal helpers — file scope only
// ============================================================================

namespace {

struct PrincipalStresses {
    double sigma1, sigma2, sigma3;
};


/**
    @brief Out-of-plane normal stress for 2D problems.

    Returns ν(σ_xx + σ_yy) for plane strain, 0 for plane stress.

    @param ts         2D stress tensor.
    @param condition  Plane condition.
    @param nu         Poisson's ratio (used only for plane strain).
    @return           σ_zz.
*/
double computeSzz( const TensorVariable<2U>& ts,
                   PlaneAssumption           condition,
                   double                    nu ) noexcept
{
    return ( condition == PlaneAssumption::PLANE_STRAIN )
           ? nu * ( ts(0,0) + ts(1,1) )
           : 0.;
}


/**
    @brief Mean stress p = I₁/3.

    Positive in compression (geotechnical convention).

    @param sxx  σ_xx
    @param syy  σ_yy
    @param szz  σ_zz
    @return     p = (σ_xx + σ_yy + σ_zz) / 3
*/
double computeP( double sxx, double syy, double szz ) noexcept
{
    return ( sxx + syy + szz ) / 3.;
}


/**
    @brief Second deviatoric invariant J₂ = ½ s:s.

    Where s = σ - p·I is the deviatoric stress tensor.

        J₂ = ½(s_xx² + s_yy² + s_zz²) + σ_xy² + σ_xz² + σ_yz²

    @param sxx  σ_xx
    @param syy  σ_yy
    @param szz  σ_zz
    @param sxy  σ_xy
    @param sxz  σ_xz
    @param syz  σ_yz
    @return     J₂ ≥ 0
*/
double computeJ2( double sxx, double syy, double szz,
                  double sxy, double sxz, double syz ) noexcept
{
    const double p  = computeP( sxx, syy, szz );
    const double sx = sxx - p;
    const double sy = syy - p;
    const double sz = szz - p;
    return 0.5 * ( sx*sx + sy*sy + sz*sz )
         + sxy*sxy + sxz*sxz + syz*syz;
}


/**
    @brief Third deviatoric invariant J₃ = det(s).

    Where s = σ - p·I is the deviatoric stress tensor.

        J₃ = s_x·s_y·s_z
           - s_x·σ_yz² - s_y·σ_xz² - s_z·σ_xy²
           + 2·σ_xy·σ_xz·σ_yz

    @param sxx  σ_xx
    @param syy  σ_yy
    @param szz  σ_zz
    @param sxy  σ_xy
    @param sxz  σ_xz
    @param syz  σ_yz
    @return     J₃
*/
double computeJ3( double sxx, double syy, double szz,
                  double sxy, double sxz, double syz ) noexcept
{
    const double p  = computeP( sxx, syy, szz );
    const double sx = sxx - p;
    const double sy = syy - p;
    const double sz = szz - p;

    double J3 = sx * sy * sz;
    J3 -= sx * syz * syz;
    J3 -= sy * sxz * sxz;
    J3 -= sz * sxy * sxy;
    J3 += 2. * sxy * sxz * syz;
    return J3;
}


/**
    @brief Lode angle θ from J₂ and J₃ [radians].

    Defined as:

        θ = (1/3) · arcsin( -3√3·J₃ / (2·J₂^(3/2)) )

    Range: [-π/6, π/6].

    Sign convention (compression positive, geotechnical):
    - θ = -π/6 at triaxial compression (σ₁ > σ₂ = σ₃)
    - θ =  0   at pure shear
    - θ = +π/6 at triaxial extension  (σ₁ = σ₂ > σ₃)

    Returns 0 for hydrostatic states (J₂ = 0) where θ is undefined.

    @param J2  Second deviatoric invariant (must be ≥ 0).
    @param J3  Third deviatoric invariant.
    @return    Lode angle in radians.
*/
double computeLodeAngle( double J2, double J3 ) noexcept
{
    if ( J2 < numeric_limits<double>::epsilon() )
        return 0.;  // hydrostatic state — θ undefined, return neutral value

    // clamp argument of arcsin to [-1, 1] to guard against floating point
    // noise pushing it marginally outside the valid range
    const double arg = max( -1., min( 1.,
        ( -3. * sqrt( 3. ) * J3 )
        / ( 2. * pow( J2, 1.5 ) ) ) );

    return ( 1. / 3. ) * asin( arg );
}


/**
    @brief Recovers one principal stress from p, √J₂, and θ.

    Formula (Smith & Griffiths, eq. 6.4):

        σᵢ = p + (2/√3)·√J₂·sin(θ + angle_offset)

    @param p           Mean stress.
    @param sqrt_J2     √J₂.
    @param theta        Lode angle [radians].
    @param sqrt3        √3 (passed to avoid recomputation).
    @param angle_offset Offset in radians: +2π/3 for σ₁, 0 for σ₂, -2π/3 for σ₃.
    @return             Principal stress.
*/
double principalStress( double p,
                        double sqrt_J2,
                        double theta,
                        double sqrt3,
                        double angle_offset ) noexcept
{
    return p + ( 2. / sqrt3 ) * sqrt_J2 * sin( theta + angle_offset );
}


PrincipalStresses computePrincipalStresses( double p,
                                            double sqrt_J2,
                                            double theta,
                                            double sqrt3 ) noexcept
{
    double raw[3] = {
        principalStress( p, sqrt_J2, theta, sqrt3,  2.*CSMP_PI/3. ),
        principalStress( p, sqrt_J2, theta, sqrt3,  0.            ),
        principalStress( p, sqrt_J2, theta, sqrt3, -2.*CSMP_PI/3. )
    };
    std::sort( raw, raw + 3, std::greater<double>() );
    return { raw[0], raw[1], raw[2] };
}

} // anonymous namespace


// ============================================================================
//  Constructors
// ============================================================================

/**
    @note Member initialisers execute in declaration order.
    The order p_, J2_, sqrt_J2_, theta_, q_, sigma1_, sigma2_, sigma3_
    must be preserved in the header declaration.
*/
StressInvariants::StressInvariants( const TensorVariable<2U>& ts,
                                    PlaneAssumption           condition,
                                    double                    nu )
 : sqrt3_            ( sqrt( 3. )      ),
   degrees_to_radians_( CSMP_PI / 180. ),
   p_    ( computeP ( ts(0,0), ts(1,1),
                      computeSzz( ts, condition, nu ) ) ),
   J2_   ( computeJ2( ts(0,0), ts(1,1),
                      computeSzz( ts, condition, nu ),
                      ts(0,1), 0., 0. ) ),
   sqrt_J2_( sqrt( J2_ ) ),
   theta_( computeLodeAngle( J2_,
               computeJ3( ts(0,0), ts(1,1),
                          computeSzz( ts, condition, nu ),
                          ts(0,1), 0., 0. ) ) ),
   q_    ( sqrt( 3. * J2_ ) ),
   sigma1_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma1; }() ),
   sigma2_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma2; }() ),
   sigma3_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma3; }() )
{
}


StressInvariants::StressInvariants( const TensorVariable<3U>& ts )
 : sqrt3_            ( sqrt( 3. )      ),
   degrees_to_radians_( CSMP_PI / 180. ),
   p_    ( computeP ( ts(0,0), ts(1,1), ts(2,2) ) ),
   J2_   ( computeJ2( ts(0,0), ts(1,1), ts(2,2),
                      ts(0,1), ts(0,2), ts(1,2) ) ),
   sqrt_J2_( sqrt( J2_ ) ),
   theta_( computeLodeAngle( J2_,
               computeJ3( ts(0,0), ts(1,1), ts(2,2),
                          ts(0,1), ts(0,2), ts(1,2) ) ) ),
   q_    ( sqrt( 3. * J2_ ) ),
   sigma1_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma1; }() ),
   sigma2_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma2; }() ),
   sigma3_( [&]{ return computePrincipalStresses(p_,sqrt_J2_,theta_,sqrt3_).sigma3; }() )
{
}


// ============================================================================
//  LodeAngle public overloads
// ============================================================================

double StressInvariants::LodeAngle( const TensorVariable<2U>& ts,
                                    double                    sigma_zz ) const noexcept
{
    const double J3 = computeJ3( ts(0,0), ts(1,1), sigma_zz,
                                 ts(0,1), 0., 0. );
    return computeLodeAngle( J2_, J3 );
}


double StressInvariants::LodeAngle( const TensorVariable<3U>& ts ) const noexcept
{
    const double J3 = computeJ3( ts(0,0), ts(1,1), ts(2,2),
                                 ts(0,1), ts(0,2), ts(1,2) );
    return computeLodeAngle( J2_, J3 );
}


// ============================================================================
//  Yield envelope functions
// ============================================================================

/**
    Hexagonal Mohr-Coulomb yield envelope:

        K(θ) = cos θ - sin φ · sin θ / √3

    Zienkiewicz Vol. 2, eq. 4.5.1, p. 233.
*/
double StressInvariants::MohrCoulombYieldEnvelope( double friction_angle_degrees ) const noexcept
{
    const double sin_phi = sin( friction_angle_degrees * degrees_to_radians_ );
    return cos( theta_ ) - ( sin_phi * sin( theta_ ) ) / sqrt3_;
}


/**
    Smooth Mohr-Coulomb yield envelope:

        K_param = (3 - sin φ) / (3 + sin φ)
        G(θ)    = 2·K_param / [(1 + K_param) - sin(3θ)·(1 - K_param)]
        K(θ)    = 1 / G(θ)

    Zienkiewicz Vol. 2, Ch. 4.5.11.
*/
double StressInvariants::SmoothMohrCoulombYieldEnvelope( double friction_angle_degrees ) const noexcept
{
    const double sin_phi = sin( friction_angle_degrees * degrees_to_radians_ );
    const double K_param = ( 3. - sin_phi ) / ( 3. + sin_phi );
    const double G_theta = ( 2. * K_param )
                         / ( ( 1. + K_param )
                           - sin( 3. * theta_ ) * ( 1. - K_param ) );
    return 1. / G_theta;
}

} // end csmp

