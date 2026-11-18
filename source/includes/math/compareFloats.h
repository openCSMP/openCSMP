// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_COMPARE_FLOATS_H
#define CSMP_COMPARE_FLOATS_H

#include <limits>
#include <cmath>
#include <algorithm>

namespace csmp
{

/**
 @file compareFloats.h "math/compareFloats.h"

 @brief Robust floating-point comparison utilities for double precision values.

 @author Thomas Driesner, ETH Zürich (original)
 @author openCSMP project (revised)

 @section motivation Motivation
 Direct equality comparison of floating-point values is unreliable due to
 rounding errors accumulated during computation. This header provides a
 selection of inline comparison functions that account for numerical precision,
 following the formulations described in:

   Knuth, D.E., "The Art of Computer Programming", Vol. 2, §4.2.2

 and discussed at:
   https://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison

 @section conventions Conventions

 All relative-tolerance functions scale by the larger of the two magnitudes
 unless otherwise noted. This means:

 - `approximatelyEqual`   — looser:  passes if |a-b| ≤ max(|a|,|b|) · ε
 - `essentiallyEqual`     — stricter: passes if |a-b| ≤ min(|a|,|b|) · ε

 Both include an absolute-tolerance floor of `std::numeric_limits<double>::min()`
 to handle comparisons involving zero or near-zero values correctly.

 @note `!definitelyGreaterThan(a,b)` does NOT imply `definitelyLessThan(a,b)`.
 Values within the tolerance band satisfy neither. Use `approximatelyEqual`
 or `approximatelyEqualAbsRel` to test the middle case explicitly.

 @section history History
 - Oct. 2012: ported to CSMP++ by Thomas Driesner
 - 2026:      revised — zero handling, noexcept, Doxygen, std::fabs
*/


// ============================================================================
//  Machine-epsilon versions
// ============================================================================

/**
 @brief Loose relative equality test scaled by the larger magnitude.

 Returns true if @p a and @p b agree to within one ULP of the larger value:

 @f[
   |a - b| \leq \max(|a|,|b|) \cdot \varepsilon_{\text{machine}}
 @f]

 An absolute floor of `std::numeric_limits<double>::min()` is applied so that
 comparisons involving zero or subnormal values do not spuriously fail.

 @param a  First operand.
 @param b  Second operand.
 @return   true if the values are approximately equal.
*/
inline bool approximatelyEqual( double a, double b ) noexcept
{
    const double diff  = std::fabs( a - b );
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return diff <= std::max( scale * std::numeric_limits<double>::epsilon(),
                             std::numeric_limits<double>::min() );
}


/**
 @brief Strict relative equality test scaled by the smaller magnitude.

 Returns true if @p a and @p b agree to within one ULP of the smaller value:

 @f[
   |a - b| \leq \min(|a|,|b|) \cdot \varepsilon_{\text{machine}}
 @f]

 This is stricter than `approximatelyEqual` — two values must agree to the
 precision of the less significant one. An absolute floor of
 `std::numeric_limits<double>::min()` is applied for near-zero safety.

 @param a  First operand.
 @param b  Second operand.
 @return   true if the values are essentially equal.
*/
inline bool essentiallyEqual( double a, double b ) noexcept
{
    const double diff  = std::fabs( a - b );
    const double scale = std::min( std::fabs(a), std::fabs(b) );
    return diff <= std::max( scale * std::numeric_limits<double>::epsilon(),
                             std::numeric_limits<double>::min() );
}


/**
 @brief Returns true if @p a is strictly greater than @p b by more than
 one ULP of the larger magnitude.

 @f[
   (a - b) > \max(|a|,|b|) \cdot \varepsilon_{\text{machine}}
 @f]

 @note `!definitelyGreaterThan(a,b)` does NOT imply `definitelyLessThan(a,b)`.
 Values within the tolerance band satisfy neither predicate.

 @param a  First operand.
 @param b  Second operand.
 @return   true if a is definitely greater than b.
*/
inline bool definitelyGreaterThan( double a, double b ) noexcept
{
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return (a - b) > std::max( scale * std::numeric_limits<double>::epsilon(),
                                std::numeric_limits<double>::min() );
}


/**
 @brief Returns true if @p a is strictly less than @p b by more than
 one ULP of the larger magnitude.

 @f[
   (b - a) > \max(|a|,|b|) \cdot \varepsilon_{\text{machine}}
 @f]

 @note `!definitelyLessThan(a,b)` does NOT imply `definitelyGreaterThan(a,b)`.
 Values within the tolerance band satisfy neither predicate.

 @param a  First operand.
 @param b  Second operand.
 @return   true if a is definitely less than b.
*/
inline bool definitelyLessThan( double a, double b ) noexcept
{
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return (b - a) > std::max( scale * std::numeric_limits<double>::epsilon(),
                                std::numeric_limits<double>::min() );
}


// ============================================================================
//  User-defined tolerance versions
// ============================================================================

/**
 @brief Loose relative equality test with a user-supplied relative tolerance.

 Returns true if:

 @f[
   |a - b| \leq \max(|a|,|b|) \cdot \varepsilon
 @f]

 An absolute floor of `std::numeric_limits<double>::min()` is applied for
 near-zero safety.

 @param a        First operand.
 @param b        Second operand.
 @param epsilon  Relative tolerance (e.g. 1e-9).
 @return         true if the values are approximately equal within @p epsilon.
*/
inline bool approximatelyEqual( double a, double b, double epsilon ) noexcept
{
    const double diff  = std::fabs( a - b );
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return diff <= std::max( scale * epsilon,
                             std::numeric_limits<double>::min() );
}


/**
 @brief Strict relative equality test with a user-supplied relative tolerance.

 Returns true if:

 @f[
   |a - b| \leq \min(|a|,|b|) \cdot \varepsilon
 @f]

 An absolute floor of `std::numeric_limits<double>::min()` is applied for
 near-zero safety.

 @param a        First operand.
 @param b        Second operand.
 @param epsilon  Relative tolerance (e.g. 1e-9).
 @return         true if the values are essentially equal within @p epsilon.
*/
inline bool essentiallyEqual( double a, double b, double epsilon ) noexcept
{
    const double diff  = std::fabs( a - b );
    const double scale = std::min( std::fabs(a), std::fabs(b) );
    return diff <= std::max( scale * epsilon,
                             std::numeric_limits<double>::min() );
}


/**
 @brief Returns true if @p a is strictly greater than @p b by more than
 @p epsilon times the larger magnitude.

 @f[
   (a - b) > \max(|a|,|b|) \cdot \varepsilon
 @f]

 @note `!definitelyGreaterThan(a,b,eps)` does NOT imply
 `definitelyLessThan(a,b,eps)`.

 @param a        First operand.
 @param b        Second operand.
 @param epsilon  Relative tolerance.
 @return         true if a is definitely greater than b.
*/
inline bool definitelyGreaterThan( double a, double b, double epsilon ) noexcept
{
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return (a - b) > std::max( scale * epsilon,
                                std::numeric_limits<double>::min() );
}


/**
 @brief Returns true if @p a is strictly less than @p b by more than
 @p epsilon times the larger magnitude.

 @f[
   (b - a) > \max(|a|,|b|) \cdot \varepsilon
 @f]

 @note `!definitelyLessThan(a,b,eps)` does NOT imply
 `definitelyGreaterThan(a,b,eps)`.

 @param a        First operand.
 @param b        Second operand.
 @param epsilon  Relative tolerance.
 @return         true if a is definitely less than b.
*/
inline bool definitelyLessThan( double a, double b, double epsilon ) noexcept
{
    const double scale = std::max( std::fabs(a), std::fabs(b) );
    return (b - a) > std::max( scale * epsilon,
                                std::numeric_limits<double>::min() );
}


// ============================================================================
//  Combined absolute + relative tolerance version
// ============================================================================

/**
 @brief Equality test combining both a relative and an absolute tolerance.

 Returns true if:

 @f[
   |a - b| \leq \max\!\bigl(\,\text{relTol} \cdot \max(|a|,|b|),\; \text{absTol}\,\bigr)
 @f]

 This is the most robust function for general use because it handles all
 three problematic cases:

 - **Large values** — the relative term dominates, scaling with magnitude.
 - **Near-zero values** — the absolute term provides a floor, preventing
   spurious failures when both values are close to zero.
 - **Mixed scale** — whichever tolerance is larger governs, so neither
   term can force a false negative.

 Recommended defaults:
 - `relTol = 1e-9`  for typical FEM residuals
 - `absTol = 1e-12` for near-zero physical quantities (pressures, strains)

 @param a       First operand.
 @param b       Second operand.
 @param relTol  Relative tolerance (default 1e-9).
 @param absTol  Absolute tolerance floor (default 1e-12).
 @return        true if the values are equal within the combined tolerance.
*/
inline bool approximatelyEqualAbsRel( double a,
                                       double b,
                                       double relTol = 1e-9,
                                       double absTol = 1e-12 ) noexcept
{
    return std::fabs( a - b )
           <= std::max( relTol * std::max( std::fabs(a), std::fabs(b) ),
                        absTol );
}

} // namespace csmp

#endif // CSMP_COMPARE_FLOATS_H
