// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <cmath>
#include <iostream>
#include <sstream>

#include "H2ONaClThermalEquilibrator.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{

// =============================================================================
// CONSTRUCTOR
//
// All external quantities are passed by CONST REFERENCE.  The equilibrator
// therefore always sees the caller's current values without any copy overhead.
// Internal "_eq" copies are taken at the start of each Equilibrate() call so
// that the bisection loop works on a consistent snapshot even if the external
// references change mid-iteration (which cannot happen in the normal visitor
// flow, but is a useful safety invariant).
//
// The "_ini" copies are diagnostic snapshots stored for error reporting only.
//
// NOTE on tmin / tmax initialisation:
//   These bounds must match the validity range of the underlying lookup tables.
//   They are hard-wired here (5 – 1000 °C) and must not be changed without
//   also extending the lookup tables.  icrit_max is derived from them so that
//   the bisection is guaranteed to reach double precision.
// =============================================================================

H2ONaClThermalEquilibrator::H2ONaClThermalEquilibrator(
    const double& external_mass_rock,   // [kg]      rock mass in control volume
    const double& external_cp_rock,     // [J/kg/K]  rock isobaric heat capacity
    const double& external_rho_rock,    // [kg/m³]   rock density

    const double& external_mass_air,    // [kg]      air mass in control volume
    const double& external_cp_air,      // [J/kg/K]  air isobaric heat capacity

    const double& external_porosity,    // [–]       nodal porosity
    const double& external_mass_fluid,  // [kg]      H2O-NaCl mass (no air, no rock)
    const double& external_wt_current,  // [wt%]     NaCl salinity
    const double& external_t_previous,  // [°C]      temperature at previous timestep
    const double& external_p_current,   // [Pa]      current pressure
    const double& external_H_current,   // [J]       total enthalpy (rock+fluid+air)
    const double& external_H_previous,  // [J]       total enthalpy at previous timestep
    const bool&   fixed_external_t,     // if true: skip bisection, use t_fix directly
    const double& t_fix,                // [°C]      fixed temperature (DIRICH nodes)
    const double& external_t_diff,      // [°C]      post-diffusion temperature (fallback)
    const bool&   verbose)
    :
    // ── External references (never copied; always live) ───────────────────────
    mass_rock    (external_mass_rock),
    cp_rock      (external_cp_rock),
    rho_rock     (external_rho_rock),
    mass_air     (external_mass_air),
    cp_air       (external_cp_air),
    phi          (external_porosity),
    mass_fluid   (external_mass_fluid),
    wt_current   (external_wt_current),
    t_previous   (external_t_previous),
    p_current    (external_p_current),
    H_current    (external_H_current),
    H_previous   (external_H_previous),
    t_fix        (t_fix),
    t_diff       (external_t_diff),
    fixed_t      (fixed_external_t),
    verbose      (verbose),

    // ── State flags ───────────────────────────────────────────────────────────
    equilibrated (false),
    fatal        (false),

    // ── Algorithm options ─────────────────────────────────────────────────────
    convergence_speed_up                  (false),
    with_rock_liquidus_solidus            (false),
    tl                                    (-1.0),
    ts                                    (-1.0),
    enable_air_phase                      (false),
    attempt_to_survive_fluid_properties_error(false),

    // ── Convergence parameters ────────────────────────────────────────────────
    convergence_criterion (1.0e-4),  // fractional enthalpy tolerance
    minimum_dt            (20.0),    // [°C] minimum bisection bracket half-width

    // ── Bisection state (initialised to safe dummy values) ────────────────────
    t_eq          (0.0),
    x_current_eq  (1.0),
    h_fluid_eq    (0.0),
    h_fluid_test  (0.0),
    h_air_eq      (0.0),
    h_air_test    (0.0),
    tmin          (5.0),    // lower bound: must match lookup-table validity range
    tmax          (1000.0),  // upper bound: must match lookup-table validity range
    resid         (0.0),
    Hmax          (0.0),
    Hmin          (0.0),

    // ── Working copies of external quantities (snapshotted each Equilibrate call)
    // Initialised slightly offset from the references so that the first
    // call always detects a change and triggers a proper update.
    mass_rock_eq   (mass_rock  - 1.0e-10),
    cp_rock_eq     (cp_rock    - 1.0e-10),
    mass_air_eq    (mass_air   - 1.0e-10),
    cp_air_eq      (cp_air     - 1.0e-10),
    mass_fluid_eq  (mass_fluid - 1.0e-10),
    wt_current_eq  (wt_current - 1.0e-10),
    t_previous_eq  (t_previous - 1.0e-10),
    p_current_eq   (p_current  - 1.0e-10),
    H_current_eq   (H_current  - 1.0e-10),
    H_previous_eq  (H_previous - 1.0e-10),

    // ── Solver temporaries ────────────────────────────────────────────────────
    H_test        (0.0),
    tdummy        (0.0),
    hdummy        (0.0),

    // ── Diagnostic snapshots (for error reporting only) ───────────────────────
    t_previous_ini (0.0),
    mass_rock_ini  (0.0),  cp_rock_ini  (0.0),  rho_rock_ini (0.0),
    mass_air_ini   (0.0),  cp_air_ini   (0.0),
    porosity_ini   (0.0),
    mass_fluid_ini (0.0),  wt_current_ini(0.0), x_current_ini(0.0),
    p_current_ini  (0.0),  H_current_ini (0.0), H_previous_ini(0.0),

    // ── Iteration counters ────────────────────────────────────────────────────
    icrit     (0),
    icrit_max (0),

    // ── Sub-objects ───────────────────────────────────────────────────────────
    // fluid: H2O-NaCl SOWAT lookup.  t and h are fed via tdummy/hdummy so that
    //        ComputeTotalEnthalpyAtTemperature() controls what SOWAT sees.
    //        p and x are the live references; they must never change during bisection.
    fluid(tdummy, p_current_eq, x_current_eq, hdummy, cp_rock, rho_rock, phi, verbose),

    // water_equilibrator: handles the pure-H2O (x=0) limit
    water_equilibrator(mass_rock, cp_rock, rho_rock,
                       mass_air,  cp_air,
                       phi, mass_fluid, t_previous,
                       p_current, H_current, H_previous,
                       fixed_t, t_fix, t_diff),

    csmp_error( ErrorHandler::Instance() )
{
    // Compute icrit_max: the number of bisection halvings needed to shrink
    // the initial bracket (tmax-tmin) down to machine epsilon.
    //   (tmax - tmin) * 0.5^icrit_max  <=  epsilon
    //   icrit_max  =  ceil( log(epsilon/(tmax-tmin)) / log(0.5) )
    // Extra counts: +1 for rounding, +1 safety, +10 to cover the weighted
    // (secant-style) iterations used for the first 10 steps.
    double n_halvings  = log(numeric_limits<double>::epsilon() / (tmax - tmin))
                        / log(0.5);
    icrit_max = static_cast<int>(n_halvings) + 1 + 1 + 10;

    equilibrated = false;
}


H2ONaClThermalEquilibrator::~H2ONaClThermalEquilibrator()
{}


// =============================================================================
// PUBLIC ACCESSORS
// Each triggers a full equilibration before returning the requested properties.
// =============================================================================
Fluidproperties H2ONaClThermalEquilibrator::Liquid() { Equilibrate(); return liquidprops; }
Fluidproperties H2ONaClThermalEquilibrator::Vapor()  { Equilibrate(); return vaporprops;  }
Fluidproperties H2ONaClThermalEquilibrator::Bulk()   { Equilibrate(); return bulkprops;   }
Fluidproperties H2ONaClThermalEquilibrator::Salt()   { Equilibrate(); return saltprops;   }

double H2ONaClThermalEquilibrator::Resid()        { return resid;       }
int    H2ONaClThermalEquilibrator::n_iterations() { return icrit;       }
bool   H2ONaClThermalEquilibrator::Equilibrated() { return equilibrated;}
bool   H2ONaClThermalEquilibrator::Fatal()        { return fatal;       }


// =============================================================================
// ThreePhaseProperties
//
// Primary entry point called by NaClH2OPropertiesVisitorPHX.
// Runs Equilibrate() and writes the four phase-property structs back to the
// caller's variables by reference.  Returns true if convergence was achieved.
// =============================================================================
bool H2ONaClThermalEquilibrator::ThreePhaseProperties(Fluidproperties& bulk_external,
                                                      Fluidproperties& liquid_external,
                                                      Fluidproperties& vapor_external,
                                                      Fluidproperties& salt_external)
{
    Equilibrate();
    bulk_external   = bulkprops;
    liquid_external = liquidprops;
    vapor_external  = vaporprops;
    salt_external   = saltprops;
    return equilibrated;
}


// =============================================================================
// Equilibrate
//
// Core isobaric thermal equilibration.  Finds the temperature T* such that
//
//   H_total(T*) = H_current
//
// where  H_total(T) = mass_rock  * h_rock(T)
//                   + mass_air   * h_air(T)       [if enable_air_phase]
//                   + mass_fluid * h_fluid(T,p,x)
//
// and h_fluid is returned by the H2O-NaCl SOWAT lookup (H2ONaClFluidProperties).
//
// ALGORITHM:
//   The bisection bracket [tmin, tmax] is established by FindInitialValues()
//   (or InitialValues() when convergence_speed_up is true).  The loop then
//   alternates between a weighted-secant step (first 10 iterations, faster
//   convergence in smooth regions) and pure bisection (thereafter, guaranteed
//   convergence).
//
// AIR HANDLING:
//   Air energy is stripped inside ComputeTotalEnthalpyAtTemperature() before
//   passing the remaining enthalpy to SOWAT.  The visitor must NOT subtract
//   air energy before calling this function — doing so would double-strip.
//
// PURE-WATER LIMIT (x = 0):
//   Delegated to water_equilibrator (H2OThermalEquilibrator) which has a
//   dedicated pure-H2O implementation.
//
// FIXED-TEMPERATURE NODES:
//   If fixed_t == true, no bisection is performed.  Fluid properties are
//   evaluated directly at t_fix and returned immediately.
// =============================================================================
void H2ONaClThermalEquilibrator::Equilibrate()
{
    if(verbose) cerr << "***** H2ONaClThermalEquilibrator::Equilibrate() *****\n";

    equilibrated = false;
    fatal        = false;

    // ── Snapshot external references into working copies ─────────────────────
    mass_rock_eq   = mass_rock;     cp_rock_eq   = cp_rock;
    mass_air_eq    = mass_air;      cp_air_eq    = cp_air;
    mass_fluid_eq  = mass_fluid;
    wt_current_eq  = wt_current;
    x_current_eq   = Weight2XNaCl(wt_current);
    p_current_eq   = p_current;
    H_current_eq   = H_current;
    H_previous_eq  = H_previous;

    // ── CASE 1: Pure water (x = 0) ────────────────────────────────────────────
    // Delegate to the dedicated H2O equilibrator.
    //if( essentiallyEqual(x_current_eq, 0.0, numeric_limits<double>::epsilon()) )

    // Benoit test change to 2*epsilon switch, consistent with H2ONaClFluidProperties
    if(definitelyLessThan( x_current_eq, 1.0*numeric_limits<double>::epsilon() ) )
    {
        equilibrated = water_equilibrator.TwoPhaseProperties(bulkprops, liquidprops, vaporprops);
        saltprops.InitToZero();
        FluidPropertiesErrorCheck();
        return;
    }

    // ── CASE 2: Fixed-temperature node (thermal DIRICH boundary) ─────────────
    // No bisection needed.  Compute fluid properties at t_fix directly by
    // setting hdummy = specific H2O-NaCl enthalpy implied by H_current at t_fix.
    if(fixed_t)
    {
        tdummy     = t_fix;
        h_fluid_eq = H_current_eq
                     - (with_rock_liquidus_solidus ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                                                   : mass_rock_eq * rock.Enthalpy(tdummy));
        if(enable_air_phase) h_fluid_eq -= mass_air_eq * air.Enthalpy(tdummy);
        h_fluid_eq /= mass_fluid_eq;
        hdummy      = h_fluid_eq;

        bulkprops   = fluid.BulkProperties();
        liquidprops = fluid.LiquidProperties();
        vaporprops  = fluid.VaporProperties();
        saltprops   = fluid.SaltProperties();
        equilibrated = true;
        return;
    }

    // ── CASE 3: Normal bisection ───────────────────────────────────────────────
    // Snapshot initial state for diagnostic error reporting
    t_previous_eq   = t_previous;
    mass_rock_ini   = mass_rock_eq;   cp_rock_ini  = cp_rock_eq;  rho_rock_ini = rho_rock;
    mass_air_ini    = mass_air_eq;    cp_air_ini   = cp_air_eq;
    mass_fluid_ini  = mass_fluid_eq;  wt_current_ini = wt_current_eq;
    x_current_ini   = x_current_eq;  t_previous_ini = t_previous_eq;
    p_current_ini   = p_current_eq;
    H_current_ini   = H_current_eq;  H_previous_ini = H_previous_eq;

    tmin  = 5.0;
    tmax  = 1000.0;
    resid = 1.0;
    icrit = 0;

    // Establish the initial bracket [tmin, tmax] and evaluate H at both ends
    if(!convergence_speed_up)
    {
        FindInitialValues();
        Hmin = ComputeTotalEnthalpyAtTemperature(tmin);  ErrorCheckHmin(icrit);
        Hmax = ComputeTotalEnthalpyAtTemperature(tmax);  ErrorCheckHmax(icrit);
    }
    else
    {
        InitialValues();
        if(equilibrated) return;   // InitialValues may converge immediately
    }

    // ── Bisection loop ────────────────────────────────────────────────────────
    while(fabs(resid) > 1.0e-10)
    {
        // STEP A: Choose next temperature estimate.
        // For the first 10 iterations use a weighted secant step (faster in
        // smooth regions); afterwards fall back to pure bisection (robust).
        if(icrit < 10
            && !essentiallyEqual(tmin, tmax, numeric_limits<double>::epsilon())
            && !essentiallyEqual(Hmin, Hmax, numeric_limits<double>::epsilon()))
            t_eq = tmin + (H_current_eq - Hmin) / (Hmax - Hmin) * (tmax - tmin);
        else
            t_eq = 0.5 * (tmin + tmax);

        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid  = H_current_eq - H_test;

        // STEP B: Convergence check.
        // |resid| < 1 J is the criterion.  Note: comparison is between fabs values
        // because halite enthalpy can be formally negative, so a plain resid < 1
        // would fire incorrectly on the halite side of the phase diagram.
        if( definitelyLessThan(fabs(resid), 1.0, numeric_limits<double>::epsilon()) )
        {
            // Optional: attempt a small temperature nudge if the converged state
            // has unphysical properties (e.g. negative viscosity near the critical
            // point).  The nudge is at most ±1 °C in 0.05 °C increments.
            if(attempt_to_survive_fluid_properties_error)
            {
                properties_check_fail = false;
                temp_correction       = 0.;
                AssignAllPropertiesViaReport();

                while(properties_check_fail && fabs(temp_correction) < 1.)
                {
                    double sign = (resid >= 0.) ? +1. : -1.;
                    t_eq            += sign * 0.05;
                    temp_correction += sign * 0.05;
                    H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
                    cerr << "temp_correction applied: " << temp_correction
                         << " (resid " << (resid >= 0. ? "+" : "-") << ")\n";
                    AssignAllPropertiesViaReport();
                    cerr << "Properties check " << (properties_check_fail ? "FAIL" : "PASS") << "\n";
                }

                if(properties_check_fail)
                    cerr << "H2ONaClThermalEquilibrator: converged but fluid properties still unphysical\n";
            }

            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }

        // STEP C: Degenerate bracket (tmin == tmax) — apply loose convergence.
        // If the H-range of the collapsed bracket is tiny relative to the fluid
        // enthalpy scale, accept the result as converged (with equilibrated=false
        // to signal the approximate nature of the result to the caller).
        if( essentiallyEqual(tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            double bracket_H_range = (Hmax - Hmin) / mass_fluid_ini;
            if(fabs(bracket_H_range) < 0.1 * fabs(h_fluid_test))
            {
                AssignAllPropertiesViaReport();
                equilibrated = false;
                cerr <<endl<< "Loose convergence: tmin == tmax\n";
            }
            else
            {
                // fatal = true;
                // ErrorConditionsToScreen();
                // csmp_error.Note(FATAL_ERROR,
                //                 "H2ONaClThermalEquilibrator::Equilibrate()",
                //                 "bracket collapsed (tmin==tmax) but H-range too large to accept");
                AssignAllPropertiesViaReport();
                equilibrated = false;
                cerr<<endl<<"Bracket collapsed (tmin==tmax) and H-range should be too large to accept!";

            }
            return;
        }

        // STEP D: Narrow the bracket.
        // The two branches differ only in whether Hmin/Hmax are re-evaluated
        // from scratch (original) or reused from the last ComputeTotal call
        // (speed-up).  The speed-up avoids one extra SOWAT call per iteration.
        if(!convergence_speed_up)
        {
            if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
            {
                tmin = t_eq;
                Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
                ErrorCheckHmin(icrit);
            }
            else
            {
                tmax = t_eq;
                Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
                ErrorCheckHmax(icrit);
            }
        }
        else
        {
            if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
            {
                tmin = t_eq;   Hmin = H_test;   ErrorCheckHmin(icrit);
            }
            else
            {
                tmax = t_eq;   Hmax = H_test;   ErrorCheckHmax(icrit);
            }
        }

        // STEP E: Bracket inversion guard — should never happen, but fatal if it does
        if( definitelyGreaterThan(tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            csmp_error.Note(FATAL_ERROR,
                            "H2ONaClThermalEquilibrator::Equilibrate()",
                            "tmin > tmax after bracket update!");
            return;
        }

        // STEP F: Maximum-iteration guard
        ++icrit;
        if(icrit > icrit_max)
        {
            fluid.InitializeToBogus();
            equilibrated = false;
            fatal        = true;
            return;
        }
    }
}


// =============================================================================
// AssignAllPropertiesViaReport
//
// Called at convergence.  BulkProperties() triggers a final SOWAT evaluation
// at the current (tdummy, p, x, hdummy); the other three phases are extracted
// via the cheaper Report* functions which reuse the same evaluation.
// =============================================================================
void H2ONaClThermalEquilibrator::AssignAllPropertiesViaReport()
{
    bulkprops   = fluid.BulkProperties();          // triggers full SOWAT update
    liquidprops = fluid.ReportLiquidProperties();  // free: reuse last evaluation
    vaporprops  = fluid.ReportVaporProperties();
    saltprops   = fluid.ReportSaltProperties();
    FluidPropertiesErrorCheck();
}


// =============================================================================
// ComputeTotalEnthalpyAtTemperature
//
// Evaluates  H_total(t) = H_rock(t) + H_air(t) + H_fluid(t)
//
// for use inside the bisection loop.  Sets tdummy = t so that the fluid object
// (which holds a const-reference to tdummy) evaluates at the trial temperature.
//
// The specific H2O-NaCl enthalpy is obtained by subtracting the rock and air
// contributions from H_current_eq, then dividing by mass_fluid_eq.  This
// stripped value is assigned to hdummy so that fluid.BulkEnthalpy() receives
// the correct H2O-NaCl-only enthalpy.
//
// AIR ENERGY STRIPPING:
//   mass_air_eq * air.Enthalpy(t) is subtracted before passing to SOWAT.
//   This is the counterpart to the volume/mass shaving done by the visitor.
//   Stripping here and shaving there together make SOWAT completely unaware
//   of the air phase, which is the design intent.
// =============================================================================
double H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t )
{
    tdummy = t;

    // Stripped H2O-NaCl enthalpy [J/kg] passed to SOWAT via hdummy
    hdummy  = H_current_eq;
    hdummy -= with_rock_liquidus_solidus ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                                         : mass_rock_eq * rock.Enthalpy(tdummy);
    if(enable_air_phase) hdummy -= mass_air_eq * air.Enthalpy(tdummy);
    hdummy /= mass_fluid_eq;

    h_fluid_test = fluid.BulkEnthalpy();
    if(enable_air_phase) h_air_test = air.Enthalpy(tdummy);

    if(fluid.Fatal())
    {
        PrintStatusToCerr();
        csmp_error.Note(FATAL_ERROR,
                        "H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature()",
                        "fatal error from H2ONaClFluidProperties");
    }

    // Reconstruct total enthalpy from all constituents
    double H_total = with_rock_liquidus_solidus
                         ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                         : mass_rock_eq * rock.Enthalpy(tdummy);
    H_total += mass_fluid_eq * h_fluid_test;
    if(enable_air_phase) H_total += mass_air_eq * h_air_test;

    return H_total;
}


// =============================================================================
// FindInitialValues  (used when convergence_speed_up == false)
//
// Estimates an initial bracket [tmin, tmax] around the expected solution using
// the rock heat capacity as a proxy for the total system heat capacity.  This
// underestimates the true heat capacity (fluid and air are ignored), so the
// temperature excursion is overestimated — the bracket is guaranteed to contain
// the true solution but may be wider than necessary.  The bisection loop
// tightens it from there.
// =============================================================================
void H2ONaClThermalEquilibrator::FindInitialValues()
{
    H_test = ComputeTotalEnthalpyAtTemperature(t_previous_ini);
    double dH = H_current_eq - H_test;

    // Estimate temperature change using rock heat capacity as a lower bound
    // for the system's total heat capacity (gives an upper bound for dT).
    double dT_estimate = fabs(dH) / (rock.MinimumHeatCapacity() * mass_rock_eq);
    if(mass_rock_eq <= 0.) dT_estimate = 1.;

    if(dH > 0.)
    {
        // Solution is hotter than t_previous
        tmax = std::min(1000., t_previous_ini + dT_estimate + minimum_dt);
        tmin = std::max(5.0,  t_previous_ini                - minimum_dt);
    }
    else if(dH < 0.)
    {
        // Solution is cooler than t_previous
        tmin = std::max(5.0,  t_previous_eq  - dT_estimate - minimum_dt);
        tmax = std::min(1000., t_previous_ini               + minimum_dt);
    }
    else
    {
        // H unchanged: bracket symmetrically around t_previous
        tmin = std::max(5.0,  t_previous_ini - minimum_dt);
        tmax = std::min(1000., t_previous_ini + minimum_dt);
    }

    if(tmin > tmax)
    {
        ErrorConditionsToScreen();
        csmp_error.Note(FATAL_ERROR,
                        "H2ONaClThermalEquilibrator::FindInitialValues()",
                        "tmin > tmax after initial bracket estimation");
    }
}


// =============================================================================
// InitialValues  (used when convergence_speed_up == true)
//
// A smarter bracket initialisation introduced by PW (May 2016).  Checks
// whether the VLH invariant point or the previous temperature already satisfies
// convergence; if not, builds a tighter bracket than FindInitialValues() by
// using the rock heat capacity evaluated at the trial temperatures rather than
// its global minimum.
// =============================================================================
void H2ONaClThermalEquilibrator::InitialValues()
{
    // ── Check VLH invariant point(s) ─────────────────────────────────────────
    auto tryConverge = [&](double t_try) -> bool
    {
        t_eq   = t_try;
        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid  = H_current_eq - H_test;
        if( definitelyLessThan(fabs(resid), convergence_criterion * bulkprops.h,
                               numeric_limits<double>::epsilon()) )
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return true;
        }
        return false;
    };

    if( essentiallyEqual(p_current_eq, fluid.VLH_Pmax()) )
    {
        if( tryConverge(fluid.VLH_Tmax()) ) return;
    }
    else if( definitelyLessThan(p_current_eq, fluid.VLH_Pmax()) )
    {
        if( tryConverge(fluid.VLH_T_low (p_current_eq)) ) return;
        if( tryConverge(fluid.VLH_T_high(p_current_eq)) ) return;
    }

    // ── Check previous temperature ────────────────────────────────────────────
    if( tryConverge(t_previous_eq) ) return;

    // ── Build bracket from the t_previous_eq evaluation ──────────────────────
    // dT is estimated using the rock heat capacity evaluated at two end-points,
    // taking the minimum to give the largest (safest) dT estimate.
    auto minRockCp = [&](double ta, double tb) -> double
    {
        if(with_rock_liquidus_solidus)
            return min(rock.HeatCapacity(ta, tl, ts), rock.HeatCapacity(tb, tl, ts));
        else
            return min(rock.HeatCapacity(ta), rock.HeatCapacity(tb));
    };

    if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
    {
        // t_previous is below solution → tmin = t_previous, expand tmax upward
        tmin = t_eq;   Hmin = H_test;
        double dT = fabs(H_current_eq - H_test) / (mass_rock_eq * minRockCp(t_eq, tmax));
        if(mass_rock_eq <= 0.) dT = 1.;
        tmax = std::min(1000., tmin + dT);
        Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
        while( definitelyLessThan(Hmax, H_current_eq, numeric_limits<double>::epsilon()) )
        {
            tmax = std::min(1000., tmax + dT);
            Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
        }
    }
    else
    {
        // t_previous is above solution → tmax = t_previous, expand tmin downward
        tmax = t_eq;   Hmax = H_test;
        double dT = fabs(H_test - H_current_eq) / (mass_rock_eq * minRockCp(t_eq, tmin));
        if(mass_rock_eq <= 0.) dT = 1.;
        tmin = std::max(5., tmax - dT);
        Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
        while( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
        {
            tmin = std::max(5., tmin - dT);
            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
        }
    }

    ErrorCheckHmin(icrit);
    ErrorCheckHmax(icrit);

    if(tmin > tmax)
    {
        ErrorConditionsToScreen();
        csmp_error.Note(FATAL_ERROR,
                        "H2ONaClThermalEquilibrator::InitialValues()",
                        "tmin > tmax after bracket construction");
    }
}


// =============================================================================
// ErrorCheckHmin / ErrorCheckHmax
//
// Verify that H_current_eq lies within [Hmin, Hmax].  A violation means the
// bracket does not contain the solution, which is a fatal initialisation error
// (typically caused by incorrect domain initialisation).
// =============================================================================
void H2ONaClThermalEquilibrator::ErrorCheckHmin(const int& i)
{
    if( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
    {
        ErrorConditionsToScreen();
        liquidprops = fluid.ReportLiquidProperties();
        vaporprops  = fluid.ReportVaporProperties();
        saltprops   = fluid.ReportSaltProperties();
        PrintStatusToCerr();

        ostringstream msg;
        msg << "iteration #" << i << ": H_current_eq < Hmin at tmin = " << tmin
            << " C.  Check domain initialisation.";
        csmp_error.Note(FATAL_ERROR, "H2ONaClThermalEquilibrator::Equilibrate()", msg.str());
    }
}

void H2ONaClThermalEquilibrator::ErrorCheckHmax(const int& i)
{
    if( definitelyGreaterThan(H_current_eq, Hmax, numeric_limits<double>::epsilon()) )
    {
        ErrorConditionsToScreen();
        liquidprops = fluid.ReportLiquidProperties();
        vaporprops  = fluid.ReportVaporProperties();
        saltprops   = fluid.ReportSaltProperties();
        PrintStatusToCerr();

        ostringstream msg;
        msg << "iteration #" << i << ": H_current_eq > Hmax at tmax = " << tmax
            << " C.  Check domain initialisation.";
        csmp_error.Note(FATAL_ERROR, "H2ONaClThermalEquilibrator::Equilibrate()", msg.str());
    }
}


// =============================================================================
// FluidPropertiesErrorCheck
//
// Checks that all returned Fluidproperties structs have non-negative physical
// quantities.  Called after AssignAllPropertiesViaReport().
// Sets properties_check_fail = true if any check fires; this signals the
// attempt_to_survive loop in Equilibrate() to try a temperature nudge.
// =============================================================================
void H2ONaClThermalEquilibrator::FluidPropertiesErrorCheck()
{
    // Helper lambda to reduce repetition
    auto neg = [](double v){ return definitelyLessThan(v, 0.0, numeric_limits<double>::epsilon()); };

    properties_check_fail =
        neg(bulkprops.t)   || neg(bulkprops.p)   || neg(bulkprops.x)   ||
        neg(bulkprops.wt)  || neg(bulkprops.rho) || neg(bulkprops.beta)||
        neg(bulkprops.s)   || neg(bulkprops.mf)  ||

        neg(liquidprops.t) || neg(liquidprops.p) || neg(liquidprops.x) ||
        neg(liquidprops.wt)|| neg(liquidprops.rho)||neg(liquidprops.h) ||
        neg(liquidprops.cp)|| neg(liquidprops.beta)||neg(liquidprops.s)||
        neg(liquidprops.mf)|| neg(liquidprops.mu)||

        neg(vaporprops.t)  || neg(vaporprops.p)  || neg(vaporprops.x) ||
        neg(vaporprops.wt) || neg(vaporprops.rho)|| neg(vaporprops.h) ||
        neg(vaporprops.cp) || neg(vaporprops.s)  || neg(vaporprops.mf)||
        neg(vaporprops.mu);

    if(properties_check_fail)
    {
        cerr << "H2ONaClThermalEquilibrator::FluidPropertiesErrorCheck() found an error\n";
        PrintStatusToCerr();
    }
}


// =============================================================================
// Report* functions
//
// Low-cost accessors used by the visitor to obtain phase properties at a
// specific (t, p, x, h) without triggering a new equilibration.  These are
// safe to call immediately after BulkProperties() / LiquidProperties() etc.
// have already been evaluated at the same conditions.
// =============================================================================
Fluidproperties H2ONaClThermalEquilibrator::ReportLiquidProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;    t_eq = t;    p_current_eq = p;    x_current_eq = x;
    t_previous_eq = -100.;   // sentinel: prevent Equilibrate() from thinking it succeeded
    return fluid.LiquidProperties();
}

Fluidproperties H2ONaClThermalEquilibrator::ReportVaporProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;    t_eq = t;    p_current_eq = p;    x_current_eq = x;
    t_previous_eq = -100.;
    return fluid.VaporProperties();
}

Fluidproperties H2ONaClThermalEquilibrator::ReportBulkProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;    t_eq = t;    p_current_eq = p;    x_current_eq = x;
    t_previous_eq = -100.;
    return fluid.BulkProperties();
}

Fluidproperties H2ONaClThermalEquilibrator::ReportSaltProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;    t_eq = t;    p_current_eq = p;    x_current_eq = x;
    t_previous_eq = -100.;
    return fluid.SaltProperties();
}


// =============================================================================
// PrintStatusToCerr / ErrorConditionsToScreen
//
// Diagnostic dumps called on error.  PrintStatusToCerr outputs the full
// Fluidproperties state; ErrorConditionsToScreen outputs the bisection state.
// =============================================================================
void H2ONaClThermalEquilibrator::PrintStatusToCerr()
{
    cerr << "tdummy           = " << tdummy          << "\n"
         << "p_current_eq     = " << p_current_eq    << "\n"
         << "x_current_eq     = " << x_current_eq    << "\n"
         << "t_previous_eq    = " << t_previous_eq   << "\n"
         << "tmin             = " << tmin             << "\n"
         << "tmax             = " << tmax             << "\n"
         << "Hmin             = " << Hmin             << "\n"
         << "Hmax             = " << Hmax             << "\n"
         << "H_current_eq     = " << H_current_eq    << "\n"
         << "H_previous_eq    = " << H_previous_eq   << "\n\n";

    auto print = [](const char* label, const Fluidproperties& f)
    {
        cerr << label << ".t      = " << f.t     << "\n"
             << label << ".p      = " << f.p     << "\n"
             << label << ".x      = " << f.x     << "\n"
             << label << ".wt     = " << f.wt    << "\n"
             << label << ".rho    = " << f.rho   << "\n"
             << label << ".h      = " << f.h     << "\n"
             << label << ".cp     = " << f.cp    << "\n"
             << label << ".beta   = " << f.beta  << "\n"
             << label << ".s      = " << f.s     << "\n"
             << label << ".mf     = " << f.mf    << "\n"
             << label << ".mu     = " << f.mu    << "\n"
             << label << ".state  = " << f.state << "\n\n";
    };

    print("bulk",   bulkprops);
    print("liq",    liquidprops);
    print("vap",    vaporprops);
    print("salt",   saltprops);
}

void H2ONaClThermalEquilibrator::ErrorConditionsToScreen()
{
    cerr.setf(ios::scientific);
    cerr.precision(20);
    cerr << "--- H2ONaClThermalEquilibrator::ErrorConditionsToScreen() ---\n"
         << "tmin          = " << tmin               << "\n"
         << "tmax          = " << tmax               << "\n"
         << "Hmin          = " << Hmin               << "\n"
         << "Hmax          = " << Hmax               << "\n"
         << "H_current_eq  = " << H_current_eq       << "\n"
         << "Hmin-H_curr   = " << Hmin-H_current_eq  << "\n"
         << "Hmax-H_curr   = " << Hmax-H_current_eq  << "\n"
         << "t_eq          = " << t_eq               << "\n"
         << "t_eq-tmin     = " << t_eq-tmin          << "\n"
         << "tmax-t_eq     = " << tmax-t_eq          << "\n"
         << "tdummy        = " << tdummy             << "\n"
         << "hdummy        = " << hdummy             << "\n"
         << "p_current_eq  = " << p_current_eq       << "\n"
         << "wt_current_eq = " << wt_current_eq      << "\n"
         << "x_current_eq  = " << x_current_eq       << "\n"
         << "bulk.state    = " << bulkprops.state    << "\n\n"
         << "--- Initial conditions ---\n"
         << "m_rock      = " << mass_rock_ini   << "\n"
         << "cp_rock     = " << cp_rock_ini     << "\n"
         << "m_fluid     = " << mass_fluid_ini  << "\n"
         << "wt          = " << wt_current_ini  << "\n"
         << "x           = " << x_current_ini   << "\n"
         << "t_previous  = " << t_previous_ini  << "\n"
         << "p           = " << p_current_ini   << "\n"
         << "H_current   = " << H_current_ini   << "\n"
         << "H_previous  = " << H_previous_ini  << "\n";
}


// =============================================================================
// Setter functions — propagate options to sub-objects as needed
// =============================================================================

void H2ONaClThermalEquilibrator::SetConvergenceSpeedUpTo(bool boost)
{
    convergence_speed_up = boost;
    water_equilibrator.SetConvergenceSpeedUpTo(boost);
}

void H2ONaClThermalEquilibrator::WithRockLiquidusSolidus(bool b)
{
    with_rock_liquidus_solidus = b;
    fluid.WithRockLiquidusSolidus(b);
    water_equilibrator.WithRockLiquidusSolidus(b);
}

void H2ONaClThermalEquilibrator::SetRockLiquidusSolidusTemperatures(double tl_, double ts_)
{
    tl = tl_;   ts = ts_;
    fluid.SetRockLiquidusSolidusTemperatures(tl_, ts_);
    water_equilibrator.SetRockLiquidusSolidusTemperatures(tl_, ts_);
}

void H2ONaClThermalEquilibrator::AttemptToSurviveFluidPropertiesError(bool b)
{
    attempt_to_survive_fluid_properties_error = b;
    water_equilibrator.AttemptToSurviveFluidPropertiesError(b);
}

void H2ONaClThermalEquilibrator::SetRockHeatCapacity(double mini_cp)
{
    rock.              SetRockHeatCapacity(mini_cp);
    fluid.             SetRockHeatCapacity(mini_cp);
    water_equilibrator.SetRockHeatCapacity(mini_cp);
}

void H2ONaClThermalEquilibrator::SetRockCrystallizationCurve(
    double nu, double sigma1, double latent_heat, double b_coeff, std::string crystallization_curve_)
{
    rock.              SetRockCrystallizationCurve(nu, sigma1, latent_heat, b_coeff, crystallization_curve_);
    fluid.             SetRockCrystallizationCurve(nu, sigma1, latent_heat, b_coeff, crystallization_curve_);
    water_equilibrator.SetRockCrystallizationCurve(nu, sigma1, latent_heat, b_coeff, crystallization_curve_);
}

double H2ONaClThermalEquilibrator::SaturationPressureFromT(double t_C)
{
    return fluid.SaturationPressureFromT(t_C);  // delegates to H2OFluidProperties
}

void H2ONaClThermalEquilibrator::WithAirPhase(bool b)
{
    enable_air_phase = b;
    water_equilibrator.WithAirPhase(b);
    // H2ONaClFluidProperties does not need to know about air:
    // air is hidden from SOWAT by energy-stripping in ComputeTotalEnthalpyAtTemperature
    // and by volume/mass shaving in the visitor before calling ThreePhaseProperties.
}

} // namespace csmp
