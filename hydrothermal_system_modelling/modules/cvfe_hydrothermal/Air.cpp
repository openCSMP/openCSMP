// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <iostream>
#include <cmath>
#include "Air.h"
#include "CSMP_mathUtilities.h"

using namespace std;
namespace csmp
{

// ============================================================================
// AIR THERMODYNAMIC PROPERTIES
// ============================================================================
//
// This class implements dry-air thermodynamic properties for use in
// unsaturated porous-media flow simulations (vadose zone, geothermal, etc.).
//
// IDEAL GAS ASSUMPTION
// --------------------
// All functions treat air as an ideal and dry gas.  This is valid when:
//   - Pressure is well below the critical pressure of air (~37.7 bar).
//   - Temperature is well above the critical temperature of air (~-140°C).
//   For geothermal / vadose-zone applications (p < 10 MPa, T > 0°C) the
//   ideal-gas approximation introduces errors below 0.1% for density,
//   compressibility, and enthalpy, and is fully sufficient.
//
// ENTHALPY REFERENCE
// ------------------
// h = 0 at T = 0°C (273.15 K), consistent with the IAPWS-IF97 water EOS
// which sets h = 0 at the triple point (0.01°C, 611.657 Pa).  The 0.01°C
// offset gives a negligible ~0.04 J/kg discrepancy and is ignored.
// This consistency is essential for the energy budget when air and water
// enthalpies are combined in H_current_.
//
// DRY AIR COMPOSITION ASSUMED
// ----------------------------
// M_AIR = 0.02897 kg/mol corresponds to dry air (78% N2, 21% O2, 1% Ar).
// Humidity would reduce M slightly (M_H2O = 0.018 < M_AIR) and is neglected.
// For a humid-air extension, M_eff = M_AIR*(1-x_v) + M_H2O*x_v where x_v
// is the molar fraction of water vapor — but this is already handled by the
// H2O-NaCl Equilibrator for the vapor phase.
//
// ============================================================================

static constexpr double M_AIR  = 0.02897;    // kg/mol, dry air
static constexpr double T_ZERO = 273.15;      // K, 0°C
static constexpr double R_GAS  = 8.31446;     // J/mol/K, universal gas constant

Air::Air()
    :
    csmp_error( ErrorHandler::Instance() )
{
    cout << "\nAir constructor:\n";
    cout << " Creating a Air with";
    cout << " p-t dependent properties\n";
}


Air::~Air()
{}

double Air::HeatCapacity( double t )
{
    // ── Isobaric heat capacity cp [J/kg/K] ───────────────────────────────────
    //
    // Source: NIST-based polynomial fit for dry air
    //
    // Form:  cp_molar(T) = a0 + a1*T + a2*T² + a3*T³   [J/mol/K]
    //   a0 = 28.11,  a1 = 1.967e-3,  a2 = 4.802e-6,  a3 = -1.966e-9
    //
    // Validity: 0–800°C.  Error < 0.5% in this range.
    //
    // IDEAL GAS: cp is pressure-independent.  For a real gas there is a small
    //   pressure correction (∂cp/∂p)_T = -T*(∂²v/∂T²)_p which is negligible
    //   below ~10 MPa for air.

    double T_K       = t + T_ZERO;
    double cp_molar  = 28.11 + 1.967e-3*T_K + 4.802e-6*T_K*T_K - 1.966e-9*T_K*T_K*T_K;
    return cp_molar / M_AIR;   // J/kg/K
}

double Air::Enthalpy( double t )
{
    // ── Specific enthalpy h [J/kg] ────────────────────────────────────────────
    //
    // h(T) = integral[T_ref → T] cp(T') dT'
    //
    // Obtained by analytically integrating the cp polynomial term by term:
    //   integral[a0 + a1*T + a2*T² + a3*T³] dT
    //     = a0*T + a1/2*T² + a2/3*T³ + a3/4*T⁴  + C
    //
    // REFERENCE: h = 0 at T_ref = 273.15 K (0°C).
    //   The subtraction cp_integral(T_K) - cp_integral(T_ZERO) implements the
    //   definite integral from T_ref to T, cancelling the integration constant C.
    //   This matches the IAPWS-IF97 water reference (triple point at 0.01°C)
    //   to within ~0.04 J/kg — negligible.
    //
    // PRESSURE INDEPENDENCE: for an ideal gas, (∂h/∂p)_T = 0 exactly.
    //   No pressure argument is needed or provided.  For a real gas, an
    //   enthalpy departure function h_dep(T,p) must be added:
    //     h_dep = integral[0 → p] [v - T*(∂v/∂T)_p] dp
    //   For air below 10 MPa this correction is < 200 J/kg (< 0.01% of total h
    //   at geothermal temperatures) and is safely neglected.
    //
    // Validity: same as HeatCapacity, 0–800°C.

    double T_K = t + T_ZERO;
    auto cp_integral = [](double T)
    {
        return (28.11*T + 0.5*1.967e-3*T*T + (1./3.)*4.802e-6*T*T*T - 0.25*1.966e-9*T*T*T*T) / M_AIR;
    };
    return cp_integral(T_K) - cp_integral(T_ZERO);
}

double Air::Density( double t, double p )
{
    // ── Density rho [kg/m³] ───────────────────────────────────────────────────
    //
    // Ideal gas law:  rho = p * M / (R * T)
    //
    // Exact for an ideal gas; valid for air when:
    //   p << p_crit ≈ 3.77 MPa  and  T >> T_crit ≈ 133 K (-140°C).
    //
    // Validity / error vs. real air (from NIST data):
    //   p =  0.1 MPa (1 bar),  0–500°C:  error < 0.01%
    //   p =  1.0 MPa (10 bar), 0–500°C:  error < 0.1%
    //   p = 10.0 MPa (100 bar),  20°C:   error ≈ 3%   ← real-gas corrections needed

    return p * M_AIR / (R_GAS * (t + T_ZERO));
}

double Air::Viscosity( double t )
{
    // ── Dynamic viscosity mu [Pa·s] ───────────────────────────────────────────
    //
    // Sutherland's law:
    //   mu(T) = mu0 * (T/T0)^1.5 * (T0 + C) / (T + C)
    //
    //   mu0 = 1.716e-5 Pa·s   at T0 = 273.15 K (0°C)
    //   C   = 110.4 K          Sutherland constant for air
    //
    // Validity: 0–1500°C.  Error < 2% in this range vs. NIST data.
    //   Below 0°C extrapolates reasonably to -50°C (< 1% error).
    //
    // PRESSURE INDEPENDENCE: for an ideal gas, viscosity is pressure-independent
    //   up to ~10 MPa.

    const double mu0 = 1.716e-5;   // Pa·s at T0 = 273.15 K
    const double C   = 110.4;      // Sutherland constant for air [K]
    double T_K = t + T_ZERO;
    return mu0 * std::pow(T_K/T_ZERO, 1.5) * (T_ZERO + C) / (T_K + C);
}

double Air::Compressibility( double p )
{
    // ── Isothermal compressibility beta [1/Pa] ────────────────────────────────
    //
    // For an ideal gas:  beta = -(1/V)*(dV/dp)_T = 1/p   (exact)
    //
    // Derivation: V = nRT/p → dV/dp = -nRT/p² = -V/p → beta = 1/p.
    //
    // This is also consistent with:  beta = 1/(rho * (dp/drho)_T) = M/(R*T*rho) = 1/p
    //
    // Validity: same as Density — good below ~10 MPa; underestimates compressibility
    //   slightly at higher pressures where Z < 1 (real gas is more compressible
    //   than ideal in the sub-critical region).

    return 1.0 / p;   // [1/Pa]
}

} // namespace csmp
