/*
 *  StressInvariants.h
 *
 *  Created by Stephan Matthai on 10/25/14.
 *  Revised and corrected.
 */

#ifndef STRESS_INVARIANTS_H
#define STRESS_INVARIANTS_H

#include "TensorVariable.h"

namespace csmp {

/**
    @brief Plane condition for 2D stress analysis.

    Disambiguates the two common 2D assumptions when computing
    stress invariants from a 2D stress tensor.
*/
enum class PlaneAssumption {
    PLANE_STRESS,  ///< out-of-plane normal stress is zero: σ_zz = 0
    PLANE_STRAIN   ///< out-of-plane normal strain is zero: σ_zz = ν(σ_xx + σ_yy)
};


/**
    @brief Standard decompositions of the Cauchy stress tensor into
    engineering invariants following Smith & Griffiths and Zienkiewicz Vol. 2.

    @par Sign convention
    Compressive stresses are **positive**, following geotechnical engineering
    convention. Tensile stresses are negative.

    @par Invariants computed
    Given the symmetric Cauchy stress tensor σ, the following invariants
    are extracted and stored at construction time:

    - p = I₁/3 — mean (isotropic) stress
    - J₂ = ½ s:s — second deviatoric invariant, where s = σ - p·I
    - J₃ = det(s) — third deviatoric invariant
    - θ ∈ [-π/6, π/6] — Lode angle (radians)
    - q = √(3·J₂) — von Mises deviatoric stress

    Principal stresses are recovered as:

        σᵢ = p + (2/√3)·√J₂·sin(θ + 2πi/3)   for i = 1, 0, -1

    giving σ₁ ≥ σ₂ ≥ σ₃ (most to least compressive).

    @par Lode angle sign convention
    Following the geotechnical convention (compression positive):
    - θ = -π/6 at triaxial compression (σ₁ > σ₂ = σ₃)
    - θ =  0   at pure shear
    - θ = +π/6 at triaxial extension  (σ₁ = σ₂ > σ₃)

    @par Efficiency
    All invariants and principal stresses are computed once at construction
    and stored. All accessor methods are inlined and return stored values
    at zero computational cost, making this class suitable for use in
    tight loops over integration points.

    @par References
    - Smith & Griffiths, Programming the Finite Element Method, 4th ed., Ch. 6
    - Zienkiewicz & Taylor, The Finite Element Method, Vol. 2, Ch. 4

    @author S.K. Matthai
*/
class StressInvariants {
  public:

    /**
        @brief Constructs invariants from a 2D stress tensor.

        @param ts         Symmetric 2D Cauchy stress tensor. Components used:
                          σ_xx = ts(0,0), σ_yy = ts(1,1), σ_xy = ts(0,1).
        @param condition  Whether the 2D model is plane stress or plane strain.
                          Under plane strain, σ_zz = ν(σ_xx + σ_yy).
                          Under plane stress, σ_zz = 0.
        @param nu         Poisson's ratio. Required only for plane strain;
                          ignored for plane stress.

        @pre Compressive stresses are positive (geotechnical convention).
    */
    explicit StressInvariants( const TensorVariable<2U>& ts,
                               PlaneAssumption condition = PlaneAssumption::PLANE_STRESS,
                               double nu = 0. );

    /**
        @brief Constructs invariants from a 3D stress tensor.

        @param ts  Symmetric 3D Cauchy stress tensor. Components used:
                   σ_xx = ts(0,0), σ_yy = ts(1,1), σ_zz = ts(2,2),
                   σ_xy = ts(0,1), σ_xz = ts(0,2), σ_yz = ts(1,2).

        @pre Compressive stresses are positive (geotechnical convention).
    */
    explicit StressInvariants( const TensorVariable<3U>& ts );

    /**
        @brief Lode angle θ for a 2D stress state [radians].

        Recomputes the Lode angle from the supplied tensor and the stored
        √J₂. Returns 0 for hydrostatic states (J₂ = 0).

        @param ts        2D stress tensor.
        @param sigma_zz  Out-of-plane normal stress. Zero for plane stress;
                         ν(σ_xx + σ_yy) for plane strain.
        @return Lode angle in radians, range [-π/6, π/6].
    */
    double LodeAngle( const TensorVariable<2U>& ts,
                      double                    sigma_zz ) const noexcept;

    /**
        @brief Lode angle θ for a 3D stress state [radians].

        Recomputes the Lode angle from the supplied tensor and the stored
        √J₂. Returns 0 for hydrostatic states (J₂ = 0).

        @param ts  3D stress tensor.
        @return Lode angle in radians, range [-π/6, π/6].
    */
    double LodeAngle( const TensorVariable<3U>& ts ) const noexcept;

    /**
        @brief Mean stress p = I₁/3 [Pa].

        Positive in compression (geotechnical convention).
        Precomputed at construction — zero cost to call.
    */
    inline double MeanStress() const noexcept;

    /**
        @brief Von Mises deviatoric stress q = √(3·J₂) [Pa].

        Always non-negative. Precomputed at construction — zero cost to call.
    */
    inline double DeviatoricStress() const noexcept;

    /**
        @brief Maximum principal stress σ₁ [Pa].

        Most compressive stress under geotechnical sign convention.
        Satisfies σ₁ ≥ σ₂ ≥ σ₃.
        Precomputed at construction — zero cost to call.
    */
    inline double MaximumPrincipalStress1() const noexcept;

    /**
        @brief Intermediate principal stress σ₂ [Pa].

        Satisfies σ₁ ≥ σ₂ ≥ σ₃.
        Precomputed at construction — zero cost to call.
    */
    inline double IntermediatePrincipalStress2() const noexcept;

    /**
        @brief Least principal stress σ₃ [Pa].

        Least compressive (or most tensile) stress.
        Satisfies σ₁ ≥ σ₂ ≥ σ₃.
        Precomputed at construction — zero cost to call.
    */
    inline double LeastPrincipalStress3() const noexcept;

    /**
        @brief Hexagonal Mohr-Coulomb yield envelope K(θ).

        Angular factor of the Mohr-Coulomb yield surface on the
        deviatoric plane:

            K(θ) = cos θ - sin φ · sin θ / √3

        See Zienkiewicz Vol. 2, Ch. 4.5.1, p. 233.

        @param friction_angle_degrees  Friction angle φ in degrees.
        @return K(θ), dimensionless.
    */
    double MohrCoulombYieldEnvelope( double friction_angle_degrees ) const noexcept;

    /**
        @brief Smooth (rounded) Mohr-Coulomb yield envelope K(θ).

        Avoids the corners of the hexagonal Mohr-Coulomb surface:

            K_param = (3 - sin φ) / (3 + sin φ)
            G(θ)    = 2·K_param / [(1 + K_param) - sin(3θ)·(1 - K_param)]
            K(θ)    = 1 / G(θ)

        See Zienkiewicz Vol. 2, Ch. 4.5.11.

        @param friction_angle_degrees  Friction angle φ in degrees.
        @return K(θ), dimensionless.
    */
    double SmoothMohrCoulombYieldEnvelope( double friction_angle_degrees ) const noexcept;

  private:

    /// @name Precomputed constants
    /// @{
    const double sqrt3_;               ///< √3
    const double degrees_to_radians_;  ///< π/180
    /// @}

    /// @name Stress invariants — computed once in constructor
    /// @{
    const double p_;       ///< mean stress I₁/3
    const double J2_;      ///< second deviatoric invariant ½ s:s
    const double sqrt_J2_; ///< √J₂, precomputed for principal stress recovery
    const double theta_;   ///< Lode angle [radians]
    const double q_;       ///< von Mises stress √(3·J₂)
    /// @}

    /// @name Derived quantities — precomputed once, returned at zero cost
    /// @{
    const double sigma1_; ///< maximum principal stress
    const double sigma2_; ///< intermediate principal stress
    const double sigma3_; ///< least principal stress
    /// @}
};


// ----------------------------------------------------------------------------
//  Inlines
// ----------------------------------------------------------------------------

/// Mean stress p = I₁/3. Positive in compression (geotechnical convention).
inline double StressInvariants::MeanStress() const noexcept
{
    return p_;
}

/// Von Mises deviatoric stress q = √(3·J₂). Always non-negative.
inline double StressInvariants::DeviatoricStress() const noexcept
{
    return q_;
}

/// Maximum principal stress σ₁. Most compressive under geotechnical convention.
inline double StressInvariants::MaximumPrincipalStress1() const noexcept
{
    return sigma1_;
}

/// Intermediate principal stress σ₂.
inline double StressInvariants::IntermediatePrincipalStress2() const noexcept
{
    return sigma2_;
}

/// Least principal stress σ₃. Least compressive or most tensile.
inline double StressInvariants::LeastPrincipalStress3() const noexcept
{
    return sigma3_;
}

} // end csmp

#endif /* defined(STRESS_INVARIANTS_H) */
