// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "H2OThermalEquilibrator.h"
#include "compareFloats.h"
#include <limits>
#include <sstream>

namespace csmp
{

// =============================================================================
// CONSTRUCTOR
//
// All external quantities are passed by CONST REFERENCE — the equilibrator
// always sees the caller's current values without copying overhead.
// Internal "_eq" copies are snapshotted at the start of each Equilibrate()
// call so the bisection loop works on a consistent state.
//
// NOTE on tmin / tmax:
//   These must bracket the validity range of the underlying pure-water lookup
//   tables (currently 5 – 1000 °C).  icrit_max is derived from them so that
//   the bisection is guaranteed to reach double precision.
//
// NOTE on air:
//   This equilibrator handles air via the same energy-stripping mechanism as
//   H2ONaClThermalEquilibrator: air enthalpy is subtracted from H_current_eq
//   before passing the residual to the H2O fluid object, then added back when
//   reconstructing H_total.  Enable with WithAirPhase(true).
// =============================================================================

H2OThermalEquilibrator::H2OThermalEquilibrator(
    const double& external_mass_rock,   // [kg]      rock mass in control volume
    const double& external_cp_rock,     // [J/kg/K]  rock isobaric heat capacity
    const double& external_rho_rock,    // [kg/m³]   rock density

    const double& external_mass_air,    // [kg]      air mass in control volume
    const double& external_cp_air,      // [J/kg/K]  air isobaric heat capacity

    const double& external_porosity,    // [–]       nodal porosity
    const double& external_mass_fluid,  // [kg]      pure-H2O mass
    const double& external_t_previous,  // [°C]      temperature at previous timestep
    const double& external_p_current,   // [Pa]      current pressure
    const double& external_H_current,   // [J]       total enthalpy (rock + fluid + air)
    const double& external_H_previous,  // [J]       total enthalpy at previous timestep
    const bool&   fixed_external_t,     // if true: skip bisection, use t_fix directly
    const double& t_fix,                // [°C]      fixed temperature (DIRICH nodes)
    const double& external_t_diff)      // [°C]      post-diffusion temperature (fallback)
    :
    // ── External references (always live; never copied) ───────────────────────
    mass_rock   (external_mass_rock),
    cp_rock     (external_cp_rock),
    rho_rock    (external_rho_rock),
    mass_air    (external_mass_air),
    cp_air      (external_cp_air),
    phi         (external_porosity),
    mass_fluid  (external_mass_fluid),
    t_previous  (external_t_previous),
    p_current   (external_p_current),
    H_current   (external_H_current),
    H_previous  (external_H_previous),
    t_fix       (t_fix),
    t_diff      (external_t_diff),
    fixed_t     (fixed_external_t),

    // ── State flags ───────────────────────────────────────────────────────────
    equilibrated (false),
    fatal        (false),

    // ── Algorithm options ─────────────────────────────────────────────────────
    convergence_speed_up                     (false),
    with_rock_liquidus_solidus               (false),
    tl                                       (-1.0),
    ts                                       (-1.0),
    enable_air_phase                         (false),
    attempt_to_survive_fluid_properties_error(false),
    properties_check_fail                    (false),
    temp_correction                          (0.0),

    // ── Convergence parameters ────────────────────────────────────────────────
    convergence_criterion (1.0e-4),   // fractional enthalpy tolerance

    // ── Bisection state ───────────────────────────────────────────────────────
    t_eq      (0.0),
    h_fluid_eq(0.0),
    tmin      (5.0),     // lower bound: must match lookup-table validity range
    tmax      (1000.0),  // upper bound: must match lookup-table validity range
    resid     (0.0),
    Hmax      (0.0),
    Hmin      (0.0),
    icrit     (0),
    icrit_max (0),

    // ── Working copies of external quantities (snapshotted each Equilibrate call)
    // Initialised slightly offset so the first call always detects a change.
    mass_rock_eq  (mass_rock  - 1.0e-10),
    cp_rock_eq    (cp_rock    - 1.0e-10),
    mass_air_eq   (mass_air   - 1.0e-10),
    cp_air_eq     (cp_air     - 1.0e-10),
    mass_fluid_eq (mass_fluid - 1.0e-10),
    t_previous_eq (t_previous - 1.0e-10),
    p_current_eq  (p_current  - 1.0e-10),
    H_current_eq  (H_current  - 1.0e-10),
    H_previous_eq (H_previous - 1.0e-10),

    // ── Solver temporaries ────────────────────────────────────────────────────
    H_test  (0.0),
    tdummy  (0.0),
    hdummy  (0.0),

    // ── T-dependent rock heat capacity (currently unused; reserved) ───────────
    cpr_min    (0.0),
    cpr_max    (0.0),
    t_min      (0.0),
    t_max      (0.0),
    cpr_t_dep  (0.0),
    t_dependent_cpr(false),

    // ── Sub-objects ───────────────────────────────────────────────────────────
    // fluid: pure-H2O lookup.  t and h are fed via tdummy/hdummy so that
    //        ComputeTotalEnthalpyAtTemperature() controls what the fluid object sees.
    fluid(tdummy, p_current_eq, hdummy, cp_rock, rho_rock, phi),

    csmp_error( ErrorHandler::Instance() )
{
    // Compute icrit_max: halvings needed to shrink the initial bracket
    // (tmax - tmin) down to machine epsilon.
    //   (tmax - tmin) * 0.5^icrit_max  <=  epsilon
    //   icrit_max = ceil( log(epsilon / (tmax-tmin)) / log(0.5) )
    // +1 for rounding, +1 safety, +10 for the weighted (secant-style) steps.
    double n_halvings = log(numeric_limits<double>::epsilon() / (tmax - tmin))
                        / log(0.5);
    icrit_max = static_cast<int>(n_halvings) + 1 + 1 + 10;

    equilibrated = false;
}


H2OThermalEquilibrator::~H2OThermalEquilibrator()
{}


// =============================================================================
// PUBLIC ACCESSORS
// Each triggers a full equilibration before returning the requested properties.
// =============================================================================
Fluidproperties H2OThermalEquilibrator::Liquid() { Equilibrate(); return liquidprops; }
Fluidproperties H2OThermalEquilibrator::Vapor()  { Equilibrate(); return vaporprops;  }
Fluidproperties H2OThermalEquilibrator::Bulk()   { Equilibrate(); return bulkprops;   }

int  H2OThermalEquilibrator::n_iterations() { return icrit;        }
bool H2OThermalEquilibrator::Equilibrated() { return equilibrated; }
bool H2OThermalEquilibrator::Fatal()        { return fatal;        }


// =============================================================================
// TwoPhaseProperties
//
// Primary entry point called by H2ONaClThermalEquilibrator when salinity is
// zero (pure-water limit).  Runs Equilibrate() and writes the three phase-
// property structs back to the caller by reference.
// Returns true if convergence was achieved.
// =============================================================================
bool H2OThermalEquilibrator::TwoPhaseProperties(Fluidproperties& bulk_external,
                                                Fluidproperties& liquid_external,
                                                Fluidproperties& vapor_external)
{
    Equilibrate();
    bulk_external   = bulkprops;
    liquid_external = liquidprops;
    vapor_external  = vaporprops;
    return equilibrated;
}


// =============================================================================
// Equilibrate
//
// Core isobaric thermal equilibration for pure water.  Finds the temperature
// T* such that:
//
//   H_total(T*) = H_current
//
// where  H_total(T) = mass_rock  * h_rock(T)
//                   + mass_air   * h_air(T)       [if enable_air_phase]
//                   + mass_fluid * h_H2O(T, p)
//
// AIR HANDLING:
//   Identical to H2ONaClThermalEquilibrator: air enthalpy is stripped inside
//   ComputeTotalEnthalpyAtTemperature() before the residual is passed to the
//   pure-H2O fluid object, and added back when reconstructing H_total.
//   The visitor must NOT subtract air energy before calling this function.
//
// ALGORITHM:
//   The bisection bracket is established by InitialValues() (when
//   convergence_speed_up is true) or by direct evaluation at the lookup-table
//   limits (original path).  InitialValues() is phase-aware: it uses the
//   boiling point as a natural bracket boundary below the critical pressure,
//   giving a much tighter initial bracket in two-phase situations.
//
// FIXED-TEMPERATURE NODES:
//   If fixed_t == true, no bisection is performed.  Fluid properties are
//   evaluated directly at t_fix via the (T, h, p) path and returned immediately.
// =============================================================================
void H2OThermalEquilibrator::Equilibrate()
{
    equilibrated = false;
    fatal        = false;

    // ── CASE 1: Fixed-temperature node (thermal DIRICH boundary) ─────────────
    if(fixed_t)
    {
        mass_rock_eq  = mass_rock;    cp_rock_eq  = cp_rock;
        mass_air_eq   = mass_air;     cp_air_eq   = cp_air;
        mass_fluid_eq = mass_fluid;
        t_previous_eq = t_fix;
        p_current_eq  = p_current;
        H_current_eq  = H_current;
        H_previous_eq = H_previous;
        tdummy        = t_fix;

        // Derive the H2O-specific enthalpy by stripping rock (and air) from H_current_eq
        h_fluid_eq  = H_current_eq;
        h_fluid_eq -= with_rock_liquidus_solidus ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                                                 : mass_rock_eq * rock.Enthalpy(tdummy);
        if(enable_air_phase) h_fluid_eq -= mass_air_eq * air.Enthalpy(tdummy);
        h_fluid_eq /= mass_fluid_eq;
        hdummy = h_fluid_eq;

        bulkprops   = fluid.BulkPropertiesFromTHP();
        liquidprops = fluid.LiquidPropertiesFromTHP();
        vaporprops  = fluid.VaporPropertiesFromTHP();
        equilibrated = true;
        return;
    }

    // ── CASE 2: Normal bisection ───────────────────────────────────────────────
    // Snapshot external references into working copies
    mass_rock_eq  = mass_rock;    cp_rock_eq  = cp_rock;
    mass_air_eq   = mass_air;     cp_air_eq   = cp_air;
    mass_fluid_eq = mass_fluid;
    t_previous_eq = t_previous;
    p_current_eq  = p_current;
    H_current_eq  = H_current;
    H_previous_eq = H_previous;

    tmin  = 5.0;
    tmax  = 1000.0;
    resid = 1.0;
    icrit = 0;

    // Establish the initial bracket [tmin, tmax] and evaluate H at both ends
    if(!convergence_speed_up)
    {
        Hmin = ComputeTotalEnthalpyAtTemperature(tmin);  ErrorCheckHmin(icrit);
        Hmax = ComputeTotalEnthalpyAtTemperature(tmax);  ErrorCheckHmax(icrit);
    }
    else
    {
        InitialValues();
        // InitialValues() may converge immediately (sets equilibrated=true and returns)
        if( definitelyLessThan(fabs(resid), convergence_criterion * bulkprops.h,
                               numeric_limits<double>::epsilon()) )
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }
    }

    // ── Bisection loop ────────────────────────────────────────────────────────
    while(fabs(resid) > 1.0e-10)
    {
        // STEP A: Choose next temperature estimate.
        // First 10 iterations: weighted secant step (faster in smooth regions).
        // Thereafter: pure bisection (guaranteed convergence).
        if(icrit < 10)
            t_eq = tmin + (H_current_eq - Hmin) / (Hmax - Hmin) * (tmax - tmin);
        else
            t_eq = 0.5 * (tmin + tmax);

        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid  = H_current_eq - H_test;

        // STEP B: Convergence check.
        // Criterion: |resid| < convergence_criterion * h_fluid  (fractional enthalpy).
        if( definitelyLessThan(fabs(resid), convergence_criterion * bulkprops.h,
                               numeric_limits<double>::epsilon()) )
        {
            // Optional: nudge temperature if converged state has unphysical
            // properties (e.g. negative viscosity near the critical point).
            // Maximum nudge: ±1 °C in 0.05 °C steps.
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
                    cerr << "H2OThermalEquilibrator: converged but fluid properties still unphysical\n";
            }

            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }

        // STEP C: Degenerate bracket (tmin == tmax) — accept as converged.
        // When the bracket collapses to numerical precision without meeting the
        // strict criterion, the current solution is the best achievable.
        if( essentiallyEqual(tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }

        // STEP D: Narrow the bracket.
        // Original path: re-evaluates both Hmin and Hmax every iteration
        // (one extra ComputeTotal call, but avoids any staleness risk).
        // Speed-up path: reuses H_test for the updated bound (saves one call).
        if(!convergence_speed_up)
        {
            if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
                tmin = t_eq;
            else
                tmax = t_eq;

            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);  ErrorCheckHmin(icrit);
            Hmax = ComputeTotalEnthalpyAtTemperature(tmax);  ErrorCheckHmax(icrit);
        }
        else
        {
            if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
            {
                tmin = t_eq;   Hmin = H_test;
            }
            else
            {
                tmax = t_eq;   Hmax = H_test;
            }
            ErrorCheckHmin(icrit);
            ErrorCheckHmax(icrit);
        }

        // STEP E: Bracket inversion guard
        if( definitelyGreaterThan(tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            csmp_error.Note(FATAL_ERROR,
                            "H2OThermalEquilibrator::Equilibrate()",
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
// Called at convergence.  Reads back all three phases from the fluid object
// at the current (tdummy, p, hdummy) state established by the last call to
// ComputeTotalEnthalpyAtTemperature().
// =============================================================================
void H2OThermalEquilibrator::AssignAllPropertiesViaReport()
{
    bulkprops   = fluid.ReportBulkProperties();
    liquidprops = fluid.ReportLiquidProperties();
    vaporprops  = fluid.ReportVaporProperties();
    FluidPropertiesErrorCheck();
}


// =============================================================================
// ComputeTotalEnthalpyAtTemperature
//
// Evaluates  H_total(t) = H_rock(t) + H_air(t) + H_H2O(t)
//
// Sets tdummy = t so the fluid object (which holds a const-ref to tdummy)
// evaluates at the trial temperature.  The H2O-specific enthalpy is derived
// by stripping rock and air contributions from H_current_eq, then dividing
// by mass_fluid_eq; the result is placed in hdummy for BulkPropertiesFromTHP().
//
// AIR ENERGY STRIPPING:
//   mass_air_eq * air.Enthalpy(t) is subtracted before passing to the H2O
//   fluid object, then added back in the return value.  When enable_air_phase
//   is false and mass_air = 0, the air terms vanish automatically.
// =============================================================================
double H2OThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t )
{
    tdummy = t;

    // Stripped H2O enthalpy [J/kg] fed to the fluid object via hdummy
    hdummy  = H_current_eq;
    hdummy -= with_rock_liquidus_solidus ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                                         : mass_rock_eq * rock.Enthalpy(tdummy);
    if(enable_air_phase) hdummy -= mass_air_eq * air.Enthalpy(tdummy);
    hdummy /= mass_fluid_eq;

    bulkprops = fluid.BulkPropertiesFromTHP();  // evaluates at (tdummy, p, hdummy)

    // Reconstruct full H_total from all constituents
    double H_total = with_rock_liquidus_solidus
                         ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts)
                         : mass_rock_eq * rock.Enthalpy(tdummy);
    H_total += mass_fluid_eq * bulkprops.h;
    if(enable_air_phase) H_total += mass_air_eq * air.Enthalpy(tdummy);

    return H_total;
}


// =============================================================================
// CheckForTminTmaxOutOfRange
//
// Guard called after bracket initialisation.  If either bound falls outside
// the lookup-table validity range, the fluid object is poisoned and fatal is
// set so the visitor can handle the failure gracefully.
// =============================================================================
void H2OThermalEquilibrator::CheckForTminTmaxOutOfRange()
{
    if( tmin < 5.0 || tmin > 1000.0 || tmax < 5.0 || tmax > 1000.0 )
    {
        cerr << "H2OThermalEquilibrator: tmin/tmax out of range "
             << "(tmin = " << tmin << ", tmax = " << tmax << ")\n";
        fluid.InitializeToBogus();
        equilibrated = false;
        fatal        = true;
    }
}


// =============================================================================
// ErrorCheckHmin / ErrorCheckHmax
//
// Verify that H_current_eq lies within [Hmin, Hmax].  A violation means the
// bracket does not contain the solution — fatal initialisation error.
// =============================================================================
void H2OThermalEquilibrator::ErrorCheckHmin(const int& i)
{
    if( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
    {
        cerr << "tmin = " << tmin << "\ntmax = " << tmax
             << "\nHmin = " << Hmin << "\nHmax = " << Hmax
             << "\nH_current_eq = " << H_current_eq << "\n";
        ostringstream msg;
        msg << "iteration #" << i << ": H_current_eq < Hmin at tmin = " << tmin
            << " C.  Check domain initialisation.";
        csmp_error.Note(FATAL_ERROR, "H2OThermalEquilibrator::Equilibrate()", msg.str());
    }
}

void H2OThermalEquilibrator::ErrorCheckHmax(const int& i)
{
    if( definitelyGreaterThan(H_current_eq, Hmax, numeric_limits<double>::epsilon()) )
    {
        // Compute useful derived quantities for diagnosis
        const double H_rock_at_tmax  = with_rock_liquidus_solidus
                                          ? mass_rock_eq * rock.Enthalpy(tmax, tl, ts)
                                          : mass_rock_eq * rock.Enthalpy(tmax);
        const double H_air_at_tmax   = enable_air_phase
                                         ? mass_air_eq * air.Enthalpy(tmax) : 0.;
        const double H_fluid_residual = H_current_eq - H_rock_at_tmax - H_air_at_tmax;
        const double h_fluid_at_tmax  = (mass_fluid_eq > 0.)
                                           ? H_fluid_residual / mass_fluid_eq : 0.;

        const double H_rock_at_tprev = with_rock_liquidus_solidus
                                           ? mass_rock_eq * rock.Enthalpy(t_previous_eq, tl, ts)
                                           : mass_rock_eq * rock.Enthalpy(t_previous_eq);
        const double H_air_at_tprev  = enable_air_phase
                                          ? mass_air_eq * air.Enthalpy(t_previous_eq) : 0.;
        const double H_fluid_at_tprev = H_current_eq - H_rock_at_tprev - H_air_at_tprev;
        const double h_fluid_at_tprev = (mass_fluid_eq > 0.)
                                            ? H_fluid_at_tprev / mass_fluid_eq : 0.;

        cerr << "\nH2OThermalEquilibrator::ErrorCheckHmax:"
             << "\n  -- Bracket --"
             << "\n  tmin              = " << tmin
             << "\n  tmax              = " << tmax
             << "\n  t_previous_eq     = " << t_previous_eq
             << "\n  p_current_eq      = " << p_current_eq
             << "\n  -- Enthalpy --"
             << "\n  H_current_eq      = " << H_current_eq
             << "\n  Hmin              = " << Hmin
             << "\n  Hmax              = " << Hmax
             << "\n  H_current_eq-Hmax = " << (H_current_eq - Hmax)
             << "\n  -- At tmax (" << tmax << "C) --"
             << "\n  H_rock            = " << H_rock_at_tmax
             << "\n  H_air             = " << H_air_at_tmax
             << "\n  H_fluid_residual  = " << H_fluid_residual
             << "\n  h_fluid (=hdummy) = " << h_fluid_at_tmax
             << "\n  hdummy (last set) = " << hdummy
             << "\n  -- At t_previous (" << t_previous_eq << "C) --"
             << "\n  H_rock            = " << H_rock_at_tprev
             << "\n  H_air             = " << H_air_at_tprev
             << "\n  H_fluid_residual  = " << H_fluid_at_tprev
             << "\n  h_fluid           = " << h_fluid_at_tprev
             << "\n  -- Masses --"
             << "\n  mass_fluid_eq     = " << mass_fluid_eq
             << "\n  mass_air_eq       = " << mass_air_eq
             << "\n  mass_rock_eq      = " << mass_rock_eq
             << "\n  -- Flags --"
             << "\n  enable_air_phase  = " << enable_air_phase
             << "\n  fixed_t           = " << fixed_t
             << "\n  convergence_su    = " << convergence_speed_up
             << endl;

        ostringstream msg;
        msg << "iteration #" << i << ": H_current_eq > Hmax at tmax = " << tmax
            << " C.  Check domain initialisation.";
        csmp_error.Note(FATAL_ERROR, "H2OThermalEquilibrator::Equilibrate()", msg.str());
    }
}


// =============================================================================
// FluidPropertiesErrorCheck
//
// Checks that all returned Fluidproperties structs have non-negative physical
// quantities.  Sets properties_check_fail = true on failure, which triggers
// the temperature-nudge recovery loop in Equilibrate().
// =============================================================================
void H2OThermalEquilibrator::FluidPropertiesErrorCheck()
{
    auto neg = [](double v){ return definitelyLessThan(v, 0.0, numeric_limits<double>::epsilon()); };

    properties_check_fail =
        neg(bulkprops.t)   || neg(bulkprops.p)  || neg(bulkprops.x) ||
        neg(bulkprops.wt)  || neg(bulkprops.rho)|| neg(bulkprops.cp)||
        neg(bulkprops.beta)|| neg(bulkprops.s)  || neg(bulkprops.mf)||
        neg(bulkprops.mu)  ||

        neg(liquidprops.t) || neg(liquidprops.p) || neg(liquidprops.x)  ||
        neg(liquidprops.wt)|| neg(liquidprops.rho)||neg(liquidprops.h)  ||
        neg(liquidprops.cp)|| neg(liquidprops.beta)||neg(liquidprops.s) ||
        neg(liquidprops.mf)|| neg(liquidprops.mu) ||

        neg(vaporprops.t)  || neg(vaporprops.p)  || neg(vaporprops.x)  ||
        neg(vaporprops.wt) || neg(vaporprops.rho)|| neg(vaporprops.h)  ||
        neg(vaporprops.cp) || neg(vaporprops.beta)||neg(vaporprops.s)  ||
        neg(vaporprops.mf) || neg(vaporprops.mu);

    if(properties_check_fail)
    {
        auto p = [](const char* lbl, const Fluidproperties& f)
        {
            cerr << lbl << ".t = "    << f.t    << "  .p = "   << f.p
                 << "  .rho = "       << f.rho  << "  .h = "   << f.h
                 << "  .cp = "        << f.cp   << "  .beta = "<< f.beta
                 << "  .s = "         << f.s    << "  .mf = "  << f.mf
                 << "  .mu = "        << f.mu   << "\n";
        };
        p("bulk",  bulkprops);
        p("liq",   liquidprops);
        p("vap",   vaporprops);
    }
}


// =============================================================================
// Report* functions
//
// Low-cost accessors: return phase properties at the given (t, p, h) without
// triggering a new equilibration.  The sentinel t_previous_eq = -100 prevents
// any stale convergence assumption if Equilibrate() were called afterwards.
// =============================================================================
Fluidproperties H2OThermalEquilibrator::ReportLiquidProperties(const double& t, const double& p, const double& h)
{
    h_fluid_eq = h;    t_eq = t;    p_current_eq = p;    t_previous_eq = -100.;
    return fluid.LiquidPropertiesFromTHP();
}

Fluidproperties H2OThermalEquilibrator::ReportVaporProperties(const double& t, const double& p, const double& h)
{
    h_fluid_eq = h;    t_eq = t;    p_current_eq = p;    t_previous_eq = -100.;
    return fluid.VaporPropertiesFromTHP();
}

Fluidproperties H2OThermalEquilibrator::ReportBulkProperties(const double& t, const double& p, const double& h)
{
    h_fluid_eq = h;    t_eq = t;    p_current_eq = p;    t_previous_eq = -100.;
    return fluid.BulkPropertiesFromTHP();
}


// =============================================================================
// InitialValues  (used when convergence_speed_up == true)
//
// Constructs a tight initial bracket using prior knowledge of the phase diagram:
//
//   SUPERCRITICAL (p > p_crit):
//     No phase boundary exists.  Start from t_previous and build the bracket
//     by estimating dT from the rock heat capacity (conservative lower bound
//     for total system cp, giving an upper bound for bracket width).
//
//   SUBCRITICAL (p <= p_crit):
//     The boiling point t_boil is a natural bracket boundary.  First check
//     whether the solution lies at t_boil itself.  Then determine which side
//     of the phase boundary t_previous lies on (liquid or vapor) and build the
//     bracket accordingly, clipping to t_boil to avoid straddling it.
//
// In all cases, "to be on the safe side" expansion loops guarantee that
// H_current_eq truly lies within the final [Hmin, Hmax].
// =============================================================================
void H2OThermalEquilibrator::InitialValues()
{
    // Helper: compute dT estimate using the minimum rock heat capacity over
    // [ta, tb].  Using the minimum gives the largest (safest) dT estimate,
    // ensuring the bracket is wide enough.
    auto minRockCp = [&](double ta, double tb) -> double
    {
        if(with_rock_liquidus_solidus)
            return min(rock.HeatCapacity(ta, tl, ts), rock.HeatCapacity(tb, tl, ts));
        else
            return min(rock.HeatCapacity(ta), rock.HeatCapacity(tb));
    };

    // Helper: try t_try as the solution; return true and set equilibrated if it converges
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

    // Helper: compute H at a two-phase boundary point using the saturated
    // liquid or vapor enthalpy.  Used to pin Hmin or Hmax at t_boil.
    auto HatBoilVapor = [&](double t_b) -> double
    {
        tdummy = t_b;
        hdummy = fluid.VaporEnthalpyVL(p_current_eq);
        return with_rock_liquidus_solidus
                   ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts) + mass_fluid_eq * fluid.BulkPropertiesFromTHP().h
                   : mass_rock_eq * rock.Enthalpy(tdummy)         + mass_fluid_eq * fluid.BulkPropertiesFromTHP().h;
    };
    auto HatBoilLiquid = [&](double t_b) -> double
    {
        tdummy = t_b;
        hdummy = fluid.LiquidEnthalpyVL(p_current_eq);
        return with_rock_liquidus_solidus
                   ? mass_rock_eq * rock.Enthalpy(tdummy, tl, ts) + mass_fluid_eq * fluid.BulkPropertiesFromTHP().h
                   : mass_rock_eq * rock.Enthalpy(tdummy)         + mass_fluid_eq * fluid.BulkPropertiesFromTHP().h;
    };

    // ── SUPERCRITICAL branch ──────────────────────────────────────────────────
    if( fluid.AboveCriticalPressureH2O(p_current_eq) )
    {
        if( tryConverge(t_previous_eq) ) return;

        double dT = fabs(H_current_eq - H_test)
                    / (mass_rock_eq * minRockCp(t_eq, H_test < H_current_eq ? tmax : tmin));
        if(mass_rock_eq <= 0.) dT = 1.;

        if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
        {
            // Solution is hotter than t_previous
            tmin = t_eq;   Hmin = H_test;
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
            // Solution is cooler than t_previous
            tmax = t_eq;   Hmax = H_test;
            tmin = std::max(5., tmax - dT);
            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            while( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
            {
                tmin = std::max(5., tmin - dT);
                Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            }
        }

        // ── Safety net (supercritical) ────────────────────────────────────────
        // The branch-specific while loops above expand the bracket using a dT
        // estimate from rock heat capacity. This estimate can be inaccurate near
        // phase boundaries or when mass_rock_eq is small, leaving H_current_eq
        // outside [Hmin, Hmax]. Expand in 1°C steps as a last resort before
        // ErrorCheck fires. Without this, the crash happens inside InitialValues()
        // before control returns to Equilibrate().
        while (definitelyGreaterThan(H_current_eq, Hmax,
                                     numeric_limits<double>::epsilon()) && tmax < 990.)
        {
            tmax = std::min(990., tmax + 1.0);
            Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
        }
        while (definitelyLessThan(H_current_eq, Hmin,
                                  numeric_limits<double>::epsilon()) && tmin > 5.)
        {
            tmin = std::max(5., tmin - 1.0);
            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
        }

        ErrorCheckHmin(icrit);
        ErrorCheckHmax(icrit);
        return;
    }

    // ── SUBCRITICAL branch ────────────────────────────────────────────────────
    // The boiling point t_boil is a known phase-boundary temperature.
    // Check it as a candidate solution first; if not, build the bracket on
    // whichever side of t_boil the current state lies.
    double t_boil = fluid.BoilingPointAtPressure(p_current_eq);

    if( tryConverge(t_boil) ) return;

    if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
    {
        // ── VAPOR field: solution is at T > t_boil ────────────────────────────
        if( tryConverge(t_previous_eq) ) return;

        double dT = fabs(H_current_eq - H_test)
                    / (mass_rock_eq * minRockCp(t_eq, tmax));
        if(mass_rock_eq <= 0.) dT = 1.;

        if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
        {
            // t_previous is also below solution → tmin = max(t_previous, t_boil)
            if(t_eq > t_boil)
            {
                tmin = t_eq;   Hmin = H_test;
            }
            else
            {
                tmin = t_boil;   Hmin = HatBoilVapor(t_boil);
            }
            tmax = std::min(1000., t_eq + dT);
            Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
            while( definitelyLessThan(Hmax, H_current_eq, numeric_limits<double>::epsilon()) )
            {
                tmax = std::min(1000., tmax + dT);
                Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
            }
        }
        else
        {
            // t_previous is above solution → tmax = t_previous
            tmax = t_eq;   Hmax = H_test;
            if((tmax - dT) > t_boil)
            {
                tmin = std::max(0., tmax - dT);
                Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            }
            else
            {
                tmin = t_boil;   Hmin = HatBoilVapor(t_boil);
            }
            while( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
            {
                tmin = std::max(5., tmin - dT);
                Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            }
        }
    }
    else
    {
        // ── LIQUID field: solution is at T < t_boil ───────────────────────────
        if( tryConverge(t_previous_eq) ) return;

        double dT = fabs(H_current_eq - H_test)
                    / (mass_rock_eq * minRockCp(t_eq, tmin));
        if(mass_rock_eq <= 0.) dT = 1.;

        if( definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()) )
        {
            // t_previous is below solution → tmin = t_previous
            tmin = t_eq;   Hmin = H_test;
            if((t_eq + dT) < t_boil)
            {
                tmax = std::min(1000., t_eq + dT);
                Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
            }
            else
            {
                tmax = t_boil;   Hmax = HatBoilLiquid(t_boil);
            }
            while( definitelyLessThan(Hmax, H_current_eq, numeric_limits<double>::epsilon()) )
            {
                tmax = std::min(1000., tmax + dT);
                Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
            }
        }
        else
        {
            // t_previous is above solution → tmax = min(t_previous, t_boil)
            if(t_eq < t_boil)
            {
                tmax = t_eq;   Hmax = H_test;
            }
            else
            {
                tmax = t_boil;   Hmax = HatBoilLiquid(t_boil);
            }
            tmin = std::max(5., t_eq - dT);
            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            while( definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
            {
                tmin = std::max(5., tmin - dT);
                Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            }
        }
    }

    // ── Safety net (subcritical) ──────────────────────────────────────────────
    // Edge case: t_previous sits exactly at a phase boundary (e.g. t_boil at
    // 1 atm = 100.001°C). The branch-specific while loops above may not fully
    // expand the bracket if:
    //   - tryConverge(t_boil) set H_test ≈ H_current_eq but didn't converge
    //   - the dT estimate was too small to cross the latent heat discontinuity
    //   - tmin/tmax were pinned at t_boil via HatBoilLiquid/HatBoilVapor
    // Result: H_current_eq falls just outside [Hmin, Hmax] → ErrorCheck crash.
    // Fix: expand in 1°C steps as last resort before ErrorCheck fires.
    // This is safe because ComputeTotalEnthalpyAtTemperature is monotonic away
    // from phase boundaries, and 1°C steps will cross any boundary within
    // a few iterations.
    while (definitelyGreaterThan(H_current_eq, Hmax,
                                 numeric_limits<double>::epsilon()) && tmax < 990.)
    {
        tmax = std::min(990., tmax + 1.0);
        Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
    }
    while (definitelyLessThan(H_current_eq, Hmin,
                              numeric_limits<double>::epsilon()) && tmin > 5.)
    {
        tmin = std::max(5., tmin - 1.0);
        Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
    }

    ErrorCheckHmin(icrit);
    ErrorCheckHmax(icrit);
}


// =============================================================================
// Setter functions
// =============================================================================

void H2OThermalEquilibrator::SetConvergenceSpeedUpTo(bool boost)
{
    convergence_speed_up = boost;
}

void H2OThermalEquilibrator::WithRockLiquidusSolidus(bool b)
{
    with_rock_liquidus_solidus = b;
}

void H2OThermalEquilibrator::SetRockLiquidusSolidusTemperatures(double tl_, double ts_)
{
    tl = tl_;   ts = ts_;
}

void H2OThermalEquilibrator::AttemptToSurviveFluidPropertiesError(bool b)
{
    attempt_to_survive_fluid_properties_error = b;
}

void H2OThermalEquilibrator::SetRockHeatCapacity(double mini_cp)
{
    rock.SetRockHeatCapacity(mini_cp);
}

void H2OThermalEquilibrator::SetRockCrystallizationCurve(
    double nu, double sigma1, double latent_heat, double b_coeff, std::string crystallization_curve_)
{
    rock.SetRockCrystallizationCurve(nu, sigma1, latent_heat, b_coeff, crystallization_curve_);
}

void H2OThermalEquilibrator::WithAirPhase(bool b)
{
    enable_air_phase = b;
    // No sub-objects to propagate to: the pure-H2O fluid object has no
    // knowledge of air.  Air is handled entirely in ComputeTotalEnthalpyAtTemperature().
}

} // namespace csmp
