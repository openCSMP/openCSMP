// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_PHX_Scheme.h"
#include <chrono>
#include <iostream>

using namespace std;

namespace csmp
{
/** custom constructor */
template <uint32_t dim>
CVFEM_PHX_Scheme<dim>::
    CVFEM_PHX_Scheme(
        Model<dim> &model_ref,
        bool with_gravity,
        bool thermodynamic_density_in_gravity_term,
        const std::vector<std::string> wells,
        bool with_zinc,
        bool with_lithium,
        SolverKind solver_kind,
        SolverKind well_solver_kind,
        const std::map<std::string, WellConfiguration> &well_configs)


    :
    model(model_ref),
    p_ref(model.Database()),

    names(
        thermodynamic_density_in_gravity_term,
        with_zinc,
        with_lithium),

    thermodynamic_density_in_gravity_term_(thermodynamic_density_in_gravity_term),
    verbose(false),
    open_top(false),
    brick_wall_limiter(true),
    well_coupling_rejected(false),
    with_well_calculations(!wells.empty()),
    with_lithium_(with_lithium),

    with_air_(false),                       // set in ActivateAirPhase()

    smooth_timestep_decrease(false),
    counter_smooth_timestep_decrease(0),

    with_split_region_(false), // set in ActivateSplitRegionCoupling()
    open_space_(false),        // set in ActivateSplitRegionCoupling()
    fault_perm_anisotropy_(1.),// set in ActivateSplitRegionCoupling()

    output_CVFEM_timing(false),

    pres_grad(model, names.pres_grad_variables[0].c_str(),
              names.pres_grad_variables[1].c_str(),
              names.pres_grad_variables[2].c_str(),
              names.pres_grad_variables[3].c_str(),
              names.pres_grad_variables[4].c_str()),

    fv_transport_vapor(model, names.lhs_vapor_vector, names.rhs_vapor_vector,
                       names.fv_transport_vapor_variables[0].c_str(),
                       names.fv_transport_vapor_variables[1].c_str(),
                       names.fv_transport_vapor_variables[2].c_str()),

    fv_transport_liquid(model, names.lhs_liquid_vector,
                        names.rhs_liquid_vector,
                        names.fv_transport_liquid_variables[0].c_str(),
                        names.fv_transport_liquid_variables[1].c_str(),
                        names.fv_transport_liquid_variables[2].c_str()),

    fv_transport_air(model, names.lhs_air_vector,
                     names.rhs_air_vector,
                     names.fv_transport_air_variables[0].c_str(),
                     names.fv_transport_air_variables[1].c_str(),
                     names.fv_transport_air_variables[2].c_str()),// dormant until ActivateAirPhase()

    upwind_control(model,
                   fv_transport_liquid,
                   fv_transport_vapor,
                   names.upwind_control_variables[0].c_str(),
                   names.upwind_control_variables[1].c_str(),
                   names.densities,
                   names.relperm_visc,
                   names.saturations,
                   names.time_control,
                   names.velocities,
                   names.pore_velocities),// no air until ActivateAirPhase()

    transport(model,
              upwind_control,
              fv_transport_vapor,
              fv_transport_liquid,
              fv_transport_air),// no air until ActivateAirPhase()

    P_bundle(solver_kind), T_bundle(solver_kind),
    P_FE(P_bundle.Get()),  T_FE(T_bundle.Get()),

    capacitance_lhs(p_ref, names.capacitance_lhs_variables[0].c_str(),
                    names.capacitance_lhs_variables[1].c_str(),
                    names.capacitance_lhs_variables[2].c_str(),
                    "thickness"),

    conductance(p_ref, names.conductance_variables[0].c_str(),
                names.conductance_variables[1].c_str(),
                names.conductance_variables[2].c_str(),
                "thickness"),

    capacitance_rhs(p_ref, names.capacitance_rhs_variables[0].c_str(),
                    names.capacitance_rhs_variables[1].c_str(),
                    "thickness"),

    heat_bottom(p_ref, names.heat_bottom_variables[0].c_str(),
                names.heat_bottom_variables[1].c_str()),

    capacitance_lhs_p(p_ref, names.capacitance_lhs_p_variables[0].c_str(),
                      names.capacitance_lhs_p_variables[1].c_str(),
                      names.capacitance_lhs_p_variables[2].c_str(),
                      "thickness"),

    conductance_p_upwind_liquid(
        p_ref, upwind_control, fv_transport_liquid,
        names.conductance_p_upwind_liquid_variables[0].c_str(),
        names.conductance_p_upwind_liquid_variables[1].c_str(),
        names.conductance_p_upwind_liquid_variables[2].c_str(),
        names.conductance_p_upwind_liquid_variables[3].c_str(),
        names.conductance_p_upwind_liquid_variables[4].c_str(),
        "thickness"),

    conductance_p_upwind_vapor(
        p_ref, upwind_control, fv_transport_vapor,
        names.conductance_p_upwind_vapor_variables[0].c_str(),
        names.conductance_p_upwind_vapor_variables[1].c_str(),
        names.conductance_p_upwind_vapor_variables[2].c_str(),
        names.conductance_p_upwind_vapor_variables[3].c_str(),
        names.conductance_p_upwind_vapor_variables[4].c_str(),
        "thickness"),

    conductance_p_upwind_air(
        p_ref, upwind_control, fv_transport_air,
        names.conductance_p_upwind_air_variables[0].c_str(),
        names.conductance_p_upwind_air_variables[1].c_str(),
        names.conductance_p_upwind_air_variables[2].c_str(),
        names.conductance_p_upwind_air_variables[3].c_str(),
        names.conductance_p_upwind_air_variables[4].c_str(),
        "thickness"),// dormant until ActivateAirPhase()

    capacitance_rhs_p(p_ref, names.capacitance_rhs_p_variables[0].c_str(),
                      names.capacitance_rhs_p_variables[1].c_str(),
                      "thickness"),

    grav_liq(p_ref, upwind_control, fv_transport_liquid,
             names.grav_liq_variables[0].c_str(),
             names.grav_liq_variables[1].c_str(),
             names.grav_liq_variables[2].c_str(),
             names.grav_liq_variables[3].c_str(),
             "thickness"),

    grav_vap(p_ref, upwind_control, fv_transport_vapor,
             names.grav_vap_variables[0].c_str(),
             names.grav_vap_variables[1].c_str(),
             names.grav_vap_variables[2].c_str(),
             names.grav_vap_variables[3].c_str(),
             "thickness"),

    grav_air(p_ref, upwind_control, fv_transport_air,
             names.grav_air_variables[0].c_str(),
             names.grav_air_variables[1].c_str(),
             names.grav_air_variables[2].c_str(),
             names.grav_air_variables[3].c_str(),
             "thickness"),// dormant until ActivateAirPhase()

    source_p(p_ref, names.source_p_variables[0].c_str(),
             names.source_p_variables[1].c_str()),

    source_p2(p_ref, names.source_p2_variables[0].c_str(),
              names.source_p2_variables[1].c_str()),

    // Split boundary stuff
    iface_transferLHS_p(model,"permeability","fluid pressure","fluid pressure"),

    iface_grav(p_ref, "split region gravity mass source", "fluid pressure"),//not used anymore

    iface_transferLHS_t(
        model,
        "thermal conductivity",
        "temperature",
        "temperature"),

    equilibrator_properties(model),

    minimum_pressure (1.e4),//Benoit: lets go down to triple point water pressure, it will be fun

    pressure_limiter_transport(model, "fluid pressure", 0., 1000.e6),
    pressure_limiter_fluid(model, "fluid pressure", minimum_pressure, 490.e6),

    well_models(InstantiateWellModels(model, wells, well_solver_kind, well_configs)),

    largest_timestep(60. * 60. * 24. * 365.), well_source(0.),
    total_tracer_mass_injected(0.), model_time_(0.)

{
    // Eigen is a direct solve and is not meant for reservoir-sized systems.
    RequireReservoirCapable(solver_kind);

    //! suppress output from PDE Integrators
    P_FE.Verbose(false);
    T_FE.Verbose(false);

    // Capture run-start timestamp for the dt-cut log filename (down to the second
    // so back-to-back runs get distinct files). Computed at construction = run
    // start; the file itself is only opened if/when a cut actually occurs.
    {
        std::time_t now = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", std::localtime(&now));
        dt_cut_log_stamp_ = stamp;
    }

    equilibrator_properties.SetThermodynamicDensityInGravityTermTo(thermodynamic_density_in_gravity_term);
    equilibrator_properties.SetThermodynamicDensityInCompressibilityTermTo(true/*thermodynamic_density_in_gravity_term*/);//Making it consistent with GravityTerm

#ifdef CSMP_WITH_SAMG_SOLVER
    // Only meaningful when SAMG is the selected backend; the calls below tune the
    // SAMG settings object, which PETSc and Eigen never read.
    if (solver_kind == SolverKind::SAMG) {

        // ============================================================================
        // SAMG SOLVER SETTINGS — PRESSURE AND TEMPERATURE EQUATIONS
        //
        // SAMG is an algebraic multigrid solver. The settings below control its
        // coarsening strategy, Krylov accelerator, and relaxation scheme.
        // The right choices depend strongly on the physical problem.
        //
        // ── CONTEXT FOR THIS SOLVER ──────────────────────────────────────────────────
        // The pressure equation involves a compressibility matrix CT that can span
        // several orders of magnitude across the domain when air is present:
        //
        //   beta_air  = 1/p  ≈ 1e-5  Pa⁻¹   (ideal gas, highly compressible)
        //   beta_H2O  ≈ 4e-10 Pa⁻¹           (liquid water, nearly incompressible)
        //
        // The ratio beta_air / beta_H2O ~ 25,000 produces a poorly conditioned matrix
        // at air invasion fronts. Standard Gauss-Seidel smoothing degrades badly under
        // such contrast; ILU(0) and stronger Krylov accelerators are preferred.
        //
        // ── nred: COARSENING AGGRESSIVENESS ──────────────────────────────────────────
        //   0   Standard coarsening. Robust, produces more levels and more work per
        //       cycle, but handles anisotropy and heterogeneity well.
        //   1   Most aggressive. Fewer coarse-grid nodes, cheaper per cycle but may
        //       degrade convergence for strongly anisotropic problems (e.g. layered
        //       permeability, elongated elements). Recommended for isotropic problems.
        //   2-3 Intermediate. Good default for mildly anisotropic problems.
        //   4   Least aggressive of the aggressive options.
        //   5-6 Cluster-based coarsening (piecewise constant or multi-pass interp.).
        //       Useful for unstructured meshes with irregular connectivity.
        //   → For geothermal problems with layered stratigraphy: prefer 0 or 2.
        //   → For isotropic/homogeneous problems: 1 is efficient.
        //
        // ── ncgrad: KRYLOV ACCELERATOR ───────────────────────────────────────────────
        // The multigrid cycle is used as a preconditioner for a Krylov method.
        //   0   No Krylov acceleration (multigrid alone). Fastest per iteration but
        //       may not converge for ill-conditioned systems.
        //   1   Preconditioned CG (conjugate gradient). Optimal for symmetric positive
        //       definite (SPD) systems. Requires the matrix to be exactly SPD; any
        //       asymmetry (e.g. upwind advection, compressibility contrasts) can
        //       cause stagnation or divergence.
        //   2   Preconditioned BI-CGSTAB. Handles mildly non-symmetric systems well.
        //       Recommended when air is present, as the compressibility contrast
        //       can break the SPD assumption of CG.
        //   3   Preconditioned GMRES. Most robust for strongly non-symmetric systems.
        //       Higher memory cost (stores the full Krylov basis). Use as a last
        //       resort if BI-CGSTAB still struggles.
        //   → Without air: ncgrad=1 (CG) is efficient and sufficient.
        //   → With air:    ncgrad=2 (BI-CGSTAB) is the recommended default.
        //   → If convergence warnings persist: try ncgrad=3 (GMRES).
        //
        // ── nxtyp: SMOOTHER / RELAXATION SCHEME ──────────────────────────────────────
        // Controls the pre- and post-smoothing step on each multigrid level.
        //   0   Gauss-Seidel. Cheap per sweep, works well for well-conditioned systems
        //       with smooth coefficient variation. Degrades for large off-diagonal
        //       entries or strong material contrasts (e.g. air vs water compressibility).
        //   1   ILU(0). Incomplete LU factorisation with no fill-in. Substantially
        //       better than Gauss-Seidel for heterogeneous or ill-conditioned problems.
        //       Increases memory usage by roughly 2× the matrix storage.
        //       Strongly recommended when air saturations are significant.
        //   2   ILUT. ILU with threshold-based fill-in. More accurate than ILU(0)
        //       but higher memory cost. Rarely needed unless ILU(0) still stagnates.
        //   3   Special box relaxation. Problem-specific; not generally applicable.
        //   5   Gauss-Seidel blockwise. Groups unknowns by element; can help for
        //       coupled problems but adds setup cost.
        //   → Without air: nxtyp=0 (Gauss-Seidel) is the efficient default.
        //   → With air:    nxtyp=1 (ILU(0)) is recommended to handle the
        //                  compressibility contrast at air invasion fronts.
        //
        // ── eps: CONVERGENCE TOLERANCE ───────────────────────────────────────────────
        //   0.0   Disables the relative residual stopping criterion.
        //   >0.0  Stop when ||r|| ≤ eps × ||r0||   (relative criterion).
        //   <0.0  Stop when ||r|| ≤ |eps|           (absolute criterion).
        //   → eps = -1e-12 (current) is an absolute criterion. This is appropriate
        //     for pressure but can be difficult to achieve if the RHS varies by many
        //     orders of magnitude across the domain. If convergence warnings persist,
        //     consider switching to a relative criterion (e.g. eps = 1e-10).
        //
        // ── iswit: SETUP REUSE STRATEGY ──────────────────────────────────────────────
        // Benoit: not sure if doc below is truly relevant to us
        // Controls whether the multigrid hierarchy is rebuilt on each call.
        //   1   Assume matrix is identical to previous call. No setup at all.
        //   2   Reuse coarse grids, interpolation, and Galerkin operators entirely.
        //   3   Reuse coarse grids and interpolation; recompute Galerkin operators.
        //       Best option when the matrix changes slowly (e.g. small timesteps).
        //   4   Full SAMG run but keep memory allocated. Rebuilds everything but
        //       avoids repeated allocation/deallocation overhead.
        //   5   Full SAMG run; release all memory on return. Safest but most expensive.
        //   → For slowly evolving problems (small timesteps, steady state):
        //       iswit=3 is efficient and accurate.
        //   → For rapidly evolving problems (large timesteps, phase changes,
        //       moving air fronts): iswit=4 or 5 to avoid stale coarse grids.
        //   → Current setting (5) is the safest default; switch to 3 or 4 once
        //       the solver is confirmed stable.
        //
        // ── RECOMMENDED SETTINGS SUMMARY ─────────────────────────────────────────────
        //
        //                     Without air     With air (invasion front)
        //   nred              0 or 2          0 or 2
        //   ncgrad            1 (CG)          2 (BI-CGSTAB)
        //   nxtyp             0 (GS)          1 (ILU)
        //   eps               -1e-12          -1e-12 (or 1e-10 relative)
        //   iswit             4               4 (or 3 for small timesteps)
        //
        // ── DIAGNOSING CONVERGENCE WARNINGS ──────────────────────────────────────────
        // SAMG reports two residuals at each solve:
        //   quasi residual  = residual in the preconditioned (transformed) system
        //   real  residual  = residual in the original physical system
        // If these diverge significantly (ratio > ~100), the preconditioner is poorly
        // matched to the physics — typically caused by the compressibility contrast
        // between air and water nodes. The recommended fix is nxtyp=1 (ILU) + ncgrad=2
        // (BI-CGSTAB). If warnings persist, consider capping beta_air to a maximum
        // multiple of Bulk.beta (e.g. 100×) to reduce the condition number at the
        // cost of a slightly slower pressure response at air-bearing nodes.
        // ============================================================================

        // ── Parameter legend (SAMG manual) ──────────────────────────────────────
        //   nred    coarsening: 0 standard | 1-4 aggressive (1 most, 4 least)
        //                       5-6 cluster coarsening
        //   ncgrad  Krylov accelerator: 0 default | 1 CG | 2 BI-CGSTAB | 3 GMRES
        //   nxtyp   smoother: 0 Gauss-Seidel | 1 ILU(0) | 2 ILUT
        //                     3 box relaxation | 5 Gauss-Seidel blockwise
        //   eps     stopping criterion: >0 relative to first residual | <0 absolute
        //   iswit   re-use of decompositions across calls: 5 full run (memory freed)
        //                     4 as 5, memory kept | 3 partial setup | 2-1 no setup
        //   idmp / iout1 / iout2   output verbosity; -1 silences SAMG
        //
        // Both systems currently run the settings recommended above for the
        // air-bearing case: aggressive coarsening, BI-CGSTAB, ILU(0).
        // ────────────────────────────────────────────────────────────────────────

        // Temperature system
        SAMG_Settings& T_samg = T_bundle.SAMGSettings();
        T_samg.Set_idmp (-1);
        T_samg.Set_iout1(-1);
        T_samg.Set_iout2(-1);
        T_samg.Set_nred  (2);       // aggressive coarsening
        T_samg.Set_ncgrad(2);       // BI-CGSTAB
        T_samg.Set_nxtyp (1);       // ILU(0)
        T_samg.Set_eps(1.e-12);

        // Pressure system
        SAMG_Settings& P_samg = P_bundle.SAMGSettings();
        P_samg.Set_idmp (-1);
        P_samg.Set_iout1(-1);
        P_samg.Set_iout2(-1);
        P_samg.Set_nred  (2);       // aggressive coarsening
        P_samg.Set_ncgrad(2);       // BI-CGSTAB
        P_samg.Set_nxtyp (1);       // ILU(0)
        P_samg.Set_eps(1.e-12);
        P_samg.Set_iswit(5);        // complete run, memory released each call

        // NOTE (Benoit DD/MM/YYYY): this block previously set nred/ncgrad/nxtyp
        // twice — first to 1/1/0, then a "TEST" block overwrote them with 2/2/1.
        // Only the second set ever took effect, so the first has been removed
        // along with ~20 lines of commented-out "old settings benoit Nov 2022".
        // To go back to Gauss-Seidel + CG, change the three values above.


    } // if (solver_kind == SolverKind::SAMG)
#endif

#if defined(CSMP_WITH_PETSC_SOLVER)
    if (solver_kind == SolverKind::PETSc) {

        // ========================================================================
        // PETSc SOLVER SETTINGS — PRESSURE AND TEMPERATURE EQUATIONS
        //
        // PETSc pairs a Krylov method (KSP) with a preconditioner (PC). The same
        // reasoning as for SAMG applies: the pressure matrix is poorly conditioned
        // where compressibility contrasts are large (air invasion fronts), so the
        // preconditioner matters more than the Krylov method.
        //
        // ── KSPType ─────────────────────────────────────────────────────────────
        //   "gmres"    robust for non-symmetric systems, higher memory (default)
        //   "bcgs"     BI-CGSTAB — the equivalent of SAMG's ncgrad=2, cheaper
        //   "cg"       only valid if the matrix is exactly SPD; upwind advection
        //              and compressibility contrast generally break that
        //   "preonly"  apply the preconditioner once, no Krylov iteration — use
        //              with a direct PC ("lu") for a direct solve
        //
        // ── PCType ──────────────────────────────────────────────────────────────
        //   "ilu"      incomplete LU (default). Cheap, but a weak preconditioner
        //              for a 10k-node diffusion system: expect many iterations.
        //   "hypre"    BoomerAMG — algebraic multigrid, the closest equivalent to
        //              SAMG and the recommended choice here. Requires a PETSc
        //              build with hypre (check PETSC_HAVE_HYPRE in petscconf.h).
        //   "gamg"     PETSc's native AMG. Always available, usually slower to set
        //              up than hypre but a reasonable fallback.
        //   "lu"       direct factorisation; reference answers only, does not
        //              scale to reservoir sizes.
        //
        // ── Tolerances ──────────────────────────────────────────────────────────
        // SAMG above uses eps = 1e-12 relative to the first residual. The PETSc
        // default rtol is 1e-10, i.e. looser — matched below so the two backends
        // are comparable. atol is left at its default (1e-50, effectively off).
        //
        // NOTE: these are starting points, not tuned values. Enable
        // SetPrintConvergedReason / SetMonitorTrueResidual below when changing
        // them, and compare the resulting fields against a SAMG run.
        // ========================================================================

        PETSc_Settings& T_petsc = T_bundle.PETScSettings();
        T_petsc.SetKSPType("gmres");            // try "bcgs" to mirror SAMG ncgrad=2
        T_petsc.SetPCType("ilu");               // try "hypre" or "gamg" for AMG
        T_petsc.SetRelativeTolerance(1.e-12);   // match SAMG eps
        T_petsc.SetMaximumIterations(1000);
        T_petsc.SetPrintConvergedReason(false);
        T_petsc.SetMonitorTrueResidual(false);

        PETSc_Settings& P_petsc = P_bundle.PETScSettings();
        P_petsc.SetKSPType("gmres");
        P_petsc.SetPCType("ilu");
        P_petsc.SetRelativeTolerance(1.e-12);
        P_petsc.SetMaximumIterations(1000);
        P_petsc.SetPrintConvergedReason(false);
        P_petsc.SetMonitorTrueResidual(false);

    } // if (solver_kind == SolverKind::PETSc)
#endif

    // Eigen exposes no settings: EigenSolver hardcodes a direct SparseLU
    // factorisation. It is rejected for the reservoir by
    // RequireReservoirCapable() above, so nothing to configure here.

    // create PropertyHandles for full reset variables
    for (auto i{0U}; i < names.full_reset_variables.size(); i++)
        reset_properties.push_back(new PropertyHandle<dim>(
            model, ("reset " + names.full_reset_variables[i]).c_str(), SCALAR,
            NODE));

    cfl_dt = current_dt = control_dt = old_dt = largest_timestep;
    timestep_increment_factor = 1.;

    upwind_control.Gravity(with_gravity);

    if (with_gravity) transport.WithGravityComponentAllPhases();// no air until ActivateAirPhase()
    // Temperature diffusion
    capacitance_lhs.LumpedFormulation(true);
    capacitance_rhs.LumpedFormulation(true);

    conductance.MultiplyWithTimeIncrement(true);

    heat_bottom.AddAccumulateLater();
    heat_bottom.MultiplyWithTimeIncrement(true);

    T_FE.Add(&capacitance_lhs);
    T_FE.Add(&capacitance_rhs);
    T_FE.Add(&conductance);
    T_FE.Add(&heat_bottom);

    // Pressure diffusion
    capacitance_lhs_p.LumpedFormulation(true);
    capacitance_rhs_p.LumpedFormulation(true);

    source_p.AddAccumulateLater();
    source_p2.AddAccumulateLater();

    conductance_p_upwind_liquid.MultiplyWithTimeIncrement(true);
    conductance_p_upwind_vapor.MultiplyWithTimeIncrement(true);

    grav_liq.MultiplyWithTimeIncrement(true);
    grav_liq.AddAccumulateLater();
    grav_vap.MultiplyWithTimeIncrement(true);
    grav_vap.AddAccumulateLater();

    if (with_gravity)
    {
        P_FE.Add(&grav_liq);
        P_FE.Add(&grav_vap);
    }

    P_FE.Add(&conductance_p_upwind_liquid);
    P_FE.Add(&conductance_p_upwind_vapor);

    P_FE.Add(&capacitance_lhs_p);
    P_FE.Add(&capacitance_rhs_p);
    P_FE.Add(&source_p);
    P_FE.Add(&source_p2);

    cerr << "\n\n=== CVFEM_PHX_Scheme<" << dim << "> ===";
    cerr << "\n  with_lithium:           " << (with_lithium_         ? "ON"  : "OFF");
    cerr << "\n  with_gravity:           " << (with_gravity          ? "ON"  : "OFF");
    cerr << "\n  with_well_calculations: " << (with_well_calculations? "ON"  : "OFF");
    cerr << "\n  --- linear solvers ---";
    cerr << "\n  pressure system    (P_FE): " << P_bundle.Name();
    cerr << "\n  temperature system (T_FE): " << T_bundle.Name();
    cerr << "\n  well Newton Jacobian:      " << SolverKindName(well_solver_kind)
         << (with_well_calculations ? "" : "   (wells off)");
    cerr << "\n  matrix type (P_FE/T_FE):   SparseMatrix   (pinned: the"
            " split-boundary operator only overrides the SparseMatrix"
            " AssignToGlobal)";
    cerr << "\n  --- PDE operators added to P_FE ---";
    cerr << "\n  conductance_p_upwind_liquid/vapor: always ON";
    cerr << "\n  grav_liq/vap:                      " << (with_gravity ? "ON": "OFF");
    cerr << "\n  iface_transferLHS_p (split BC):    " << (with_split_region_ ? "ON" : "OFF");
    cerr << "\n==============================================\n";

} // end constructor

/** custom destructor */
template <uint32_t dim>
CVFEM_PHX_Scheme <
    dim >::~CVFEM_PHX_Scheme()
{
    for (auto i{0U}; i < names.full_reset_variables.size(); i++)
        delete reset_properties[i];
}

template <uint32_t dim>
std::vector<std::unique_ptr<WellModelPrototype<dim>>>
CVFEM_PHX_Scheme<dim>::InstantiateWellModels(
    Model<dim> &model, const std::vector<std::string> &wells, SolverKind well_solver_kind,
    const std::map<std::string, WellConfiguration> &well_configs)
{
    // Reserve space for well models
    std::vector<std::unique_ptr<WellModelPrototype<dim>>> well_models;
    well_models.reserve(wells.size());

    // Instantiate each well model and add to well_models
    for (const auto &well_name : wells)
    {
        std::cerr << "Instantiating a well model with name: " << well_name
                  << std::endl;
        well_models.emplace_back(
            // Each well gets its own configuration; wells with no entry fall back
            // to the defaults, which is only sensible for a well that needs no
            // completions — so report it rather than let it fail later in
            // WellConfiguration::Validate().
            std::make_unique<WellModelPrototype<dim>>(
                model, well_name, well_solver_kind,
                well_configs.count(well_name) ? well_configs.at(well_name)
                                              : WellConfiguration()));
        if (!well_configs.count(well_name))
            cerr << endl << "WARNING: no WellConfiguration given for well '"
                 << well_name << "'; using defaults, which have no completions.";
    }

    // Initialize injection rates with -1 for each well
    injection_rates = std::vector<double>(wells.size(), -1.0);

    return well_models;
}

/** modifying maximum size of time step */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::SetLargestTimeStep(
    double timestep)
{
    largest_timestep = cfl_dt = current_dt = control_dt = old_dt = timestep;
    upwind_control.SetLargestTimeStep(largest_timestep);
    transport.SetLargestTimeStep(largest_timestep);
} // end SetLargestTimeStep

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ResetLargestTimeStep(
    double timestep)
{
    largest_timestep = timestep;
    upwind_control.SetLargestTimeStep(largest_timestep);
    transport.SetLargestTimeStep(largest_timestep);
} // end ResetLargestTimeStep

/** adjusting timestep */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ChangeTimeStepTo(
    double timestep)
{
    cfl_dt = current_dt = control_dt = old_dt = timestep;
} // end SetLargestTimeStep

/** initialize fluid properties from current PTX conditions */
template <uint32_t dim>
void CVFEM_PHX_Scheme <dim >::InitialFluidPropertiesFromPTX()
{
    equilibrator_properties.InitialPropertiesFromPTX(model);

} // end InitialFluidPropertiesFromPT

/** preparation before transient calculations, calculating pressure gradient and
        updwind nodes from current status */
template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::PrepareTransientCalculations()
{
    model.Accept(pres_grad);
    transport.UpdateProjection();

    // Seed D0 (Upwinder) with a plain direction pass...
    model.Accept(upwind_control);

    // ...then seed D1 (UpwinderTransport) so the first step's transport does not
    // read an all-zero map. Both are built from the same initial velocity field.
    upwind_control.WithVelocity();   // also selects the transport map
    upwind_control.ZeroVelocityOnFlip(false);   // D0 is all-zero here; every facet
    // would read as "reversed"
    model.Accept(upwind_control);

    upwind_control.Reset();
} // end PrepareTransientCalculations

/**
 * Apply — advance the coupled system by one (adaptive) timestep; returns the dt used.
 *
 * Phase order within a step:
 *   1. AdvanceTransientVariables + well advance   (snapshot state into "previous *")
 *   2. Select dt: forced smooth decrease, or gentle growth from old_dt
 *   3. while(do_loop): AdvectionDiffusionLoops (pressure + transport, with reactive
 *      CFL / control_dt / split-boundary cuts) -> temperature solve -> well coupling
 *      -> pressure limiter -> equilibration.
 *   4. Bookkeeping: advance the smooth-decrease counter, old_dt = current_dt
 */
// ═════════════════════════════════════════════════════════════════════════════
//  UPWIND FREEZING WITHIN A TIMESTEP
// ═════════════════════════════════════════════════════════════════════════════
// TWO MODES, selected by upwind_freeze_within_step below.
//
//   false — LEGACY. Both upwind passes re-run on every solve, and the velocity
//           pass runs with `flipping` on (the UPWIND_FLIPPED hysteresis). This
//           is the historical behaviour, kept verbatim and selectable.
//
//   true  — SPLIT MAPS. TWO donor maps, because the pressure equation and the
//           transport equation need different things:
//             D0 = Upwinder          written by the DIRECTION pass, ONCE per
//                                    timestep, and reused by every later solve
//                                    of that step (dt-cut retries, split-boundary
//                                    retries). Read by the pressure operators.
//             D1 = UpwinderTransport written by the VELOCITY pass on EVERY solve
//                                    from the velocities that solve produced.
//                                    Read by transport.
//           The velocity pass also refreshes facet/nodal velocities and the CFL
//           criterion in both modes.
//
// WHY SPLIT THE TWO PASSES THAT WAY — it follows from where each one sits
// relative to the solve:
//
//   * The DIRECTION pass runs BEFORE the solve. After a ResetVariables the model
//     properties are rewound to t^n, but STENCIL_DATA still holds the facet
//     velocities written by UpdateProjection during the REJECTED solve. Deciding
//     donors there means deciding them from velocities that do not match the
//     state they are applied to — stale upwinding. Freezing is what avoids it.
//
//   * The VELOCITY pass runs AFTER UpdateProjection, i.e. on velocities this
//     iteration's solve just produced at this iteration's dt. It is never stale,
//     so there is nothing to gain by freezing it — and freezing it cost
//     correctness: cfl_dt is read from the "courant *" properties that only this
//     pass writes, so a retry compared its new, smaller dt against the REJECTED
//     solve's cfl_dt, passed trivially, and the accepted step's own velocities
//     were never CFL-checked at all.
//
// NOTE on the velocity selection while frozen. DetermineUpwindNodes still picks
// `velocity` from the CURRENT one-sided velocities; only the matrix writes are
// suppressed. That is deliberate. Driving the selection from the frozen donor
// instead would, on a facet whose frozen side has since been gated to zero,
// yield velocity = 0 and contribute NO CFL candidate — silently under-reporting
// the stability limit. Recomputing the active side keeps CFL conservative. The
// price is that on a flipped facet the CFL denominator uses the current donor's
// saturation while assembly uses the frozen one; those facets are near-zero
// velocity by construction, and erring toward a tighter CFL is the safe side.
//
// WHY. Re-evaluating upwind inside a step makes the step NON-REPRODUCIBLE, and
// therefore makes any iteration over it non-convergent. Measured with a
// determinism harness (a loop that ran the pressure solve N times against a
// "test fluid pressure" copy and reported max|dp| per repeat; removed with the
// Picard code, easy to reinstate): repeating
//     solve -> pres_grad -> UpdateProjection -> T solve -> equilibration -> reset
// is BIT-IDENTICAL over 20 repeats when the upwind evaluation is hoisted OUT of
// the loop, and oscillates by +-30 kPa with period 2 when it is left inside.
//
// THE TWO CHANNELS — both INTERNAL C++ STATE, not model properties, so no reset
// list, CopyReplace or snapshot can restore them. (This is why an exhaustive
// rework of transient_variables, and separately adding KgradP and recomputing
// the pressure gradient right after ResetVariables, both failed to remove the
// oscillation: neither touches the state that actually differs.)
//
//   1. STENCIL_DATA facet velocities. UpwindControlVisitor takes its velocity
//      from fv_transport_*.GetFacetNormalVelocity(), i.e. from STENCIL_DATA
//      inside NodeCenteredFiniteVolumeTransport — which UpdateProjection()
//      rewrites every time it runs. So upwind depends on UpdateProjection's
//      output, and a facet whose velocity sits near zero can flip donor between
//      two otherwise identical solves.
//
//   2. UPWIND_FLIPPED hysteresis. The Upwinder matrices carry a one-pass memory
//      of directions that recently changed (the UPWIND_FLIPPED marker in
//      UpwindControlVisitor). They are zeroed only at construction, so
//      consecutive evaluations differ BY DESIGN. Reset() clears two booleans and
//      does not touch this.
//
// A donor flip is DISCRETE: the change in the matrix is finite however small the
// pressure change that caused it. No relaxation can smooth a switch. This is why
// the abandoned sequential-implicit (Picard) experiment could not be rescued by
// damping: with the donor map re-evaluated inside the loop there is no fixed
// point to converge to, and it showed as a residual floor nothing could get
// below, independent of step size or damping factor.
//
// PHYSICALLY this mirrors what transport already does: velocities are frozen for
// the duration of a stage; here the donor choice is frozen for the duration of a
// step. The facets at risk of flipping are by construction those whose velocity
// is near zero — precisely where the donor choice matters least.
//
// Channel 2 in SPLIT-MAP mode: the velocity pass still runs the flipping branch,
// but it now writes D1, never D0. So UPWIND_FLIPPED can appear in the transport
// map and NEVER in the pressure map — the assembled pressure matrix is unaffected
// by the marker either way. flipping therefore stays a live, testable option on
// this path: false gives transport a clean D1 (donor consistent with the flux
// rates); true keeps the legacy veto, which refuses to choose on a reversed facet
// and zeroes its velocity contribution. Measured earlier, with the velocity pass
// writing D0 (flipping disabled vs enabled, same mesh and window): 1320 vs 1418
// dt cuts, median dt within a few percent — the facets it can mark are those
// whose velocity crosses zero between the two passes, so they carry almost no flux.
//
// NOT A COST OF THE SPLIT: D0 is decided before the solve, so the pressure
// matrix upwinds on the previous step's converged field. LEGACY DID THE SAME —
// its direction pass was also pre-solve. The split leaves the pressure side
// unchanged and IMPROVES the transport side: legacy's velocity pass never wrote
// a new donor (it only confirmed or vetoed with UPWIND_FLIPPED), so transport
// also read a lagged map; D1 is genuinely current.
// The residual lag shrinks under retry: each dt cut moves P_N toward P_{N-1},
// so the correct map migrates toward the frozen D0 — the error shrinks, D0
// itself does not change.
//
// Set upwind_freeze_within_step = false to restore LEGACY behaviour (and with it
// the non-reproducibility described above).
static constexpr bool upwind_freeze_within_step = true;

// ── Flip CHECK on the transport donor map (legacy behaviour, default OFF) ────
// When a facet's donor REVERSES between the direction pass's pre-solve choice
// (D0) and the velocity pass's post-solve one, this refuses the change: it writes
// UPWIND_FLIPPED into the transport map instead of a donor, and zeroes that
// facet's velocity, so it drops out of the nodal velocity sum and contributes no
// CFL candidate.
//
// OFF (default): the reversed donor is simply written and the facet keeps its
//   full velocity. Transport gets a clean D1 consistent with its flux rates, and
//   the CFL keeps the facet's contribution — which errs TIGHTER, the safe side.
// ON: the legacy veto. Neutralises reversing facets rather than choosing a side.
//
// Either way the PRESSURE matrix is unaffected: this pass writes D1 only.
static constexpr bool upwind_zero_velocity_on_flip = false;

template <uint32_t dim>
double CVFEM_PHX_Scheme<dim>::Apply()
{

    // New timestep: let the DIRECTION pass run once. Every later solve of this
    // step reuses its donor map (see the UPWIND FREEZING block above). The
    // velocity pass is not latched — it re-runs on every solve.
    upwind_done_this_step_ = false;

    AdvanceTransientVariables();
    if (with_well_calculations)
    {
        for (auto &well : well_models)
        {
            well->Advance_well_variables();
            well->NR_Advance_well_variables();
        }
    }



    // ═════════════════════════════════════════════════════════════════════════════
    //  DETERMINISM TEST  —  temporary; paste into Apply() immediately AFTER
    //  AdvanceTransientVariables(), before the main AdvectionDiffusionLoops() call.
    //
    //  WHAT IT DOES. Repeats the pressure solve N times from the SAME t^n state and
    //  reports how far each repeat lands from the first.
    //
    //  EXPECTED RESULT — this is the whole point, and it needs no code changes to
    //  switch between the two cases, because the mode toggle does the work:
    //
    //    upwind_freeze_within_step = true   the direction pass is latched, so it runs
    //                                       on repeat 0 only and every later repeat
    //                                       reuses that donor map
    //                                       -> max|dp| == 0, BIT-IDENTICAL, all repeats
    //
    //    upwind_freeze_within_step = false  do_upwind is unconditionally true, so the
    //                                       direction pass re-evaluates every repeat
    //                                       from STENCIL_DATA that the PREVIOUS
    //                                       repeat's UpdateProjection rewrote
    //                                       -> max|dp| ~ 3e4 Pa, alternating sign,
    //                                          period 2
    //
    //  NOTE. Do NOT clear upwind_done_this_step_ inside the loop. Leaving it alone is
    //  exactly what makes the test measure the thing we care about: whether holding
    //  the upwind decision fixed is what buys reproducibility.
    //
    //  The harness rewinds the model to t^n when it finishes, so the real step that
    //  follows is unaffected.
    // ═════════════════════════════════════════════════════════════════════════════
    {
        constexpr bool determinism_test    = false;   // flip to true to run
        constexpr int  determinism_repeats = 20;

        if (determinism_test)
        {
            const double dt_entry = current_dt;

            // Snapshot the t^n pressure. NOTE: this harness is now the ONLY user
            // of "stash pressure" and "test fluid pressure" — the well coupling
            // used to stash/restore around its implicit re-solve loop and no
            // longer does. Keep both properties declared for this reason.
            model.CopyReplace("fluid pressure", "stash pressure");

            csmp::Index p_key_  (model.Database().StorageKey("fluid pressure"));
            csmp::Index ref_key_(model.Database().StorageKey("test fluid pressure"));

            // const Region<dim>& flow_ = model.Region("FLOW REGION"); // useful in case flow is only considered in specific regions (e.g., not in an ocean)
            const Region<dim>& flow_ = model.Region("Model");
            typename vector<Node<dim>*>::const_iterator nit;

            cerr << "\n[determinism] START  repeats=" << determinism_repeats
                 << "  freeze=" << (upwind_freeze_within_step ? "ON" : "OFF")
                 << "  dt=" << dt_entry;

            for (int rep = 0; rep < determinism_repeats; ++rep)
            {
                // Rewind STATE only. The upwind latch is deliberately left alone.
                ResetVariables();
                model.CopyReplace("stash pressure", "fluid pressure");
                current_dt = dt_entry;

                AdvectionDiffusionLoops();

                if (rep == 0)
                {
                    model.CopyReplace("fluid pressure", "test fluid pressure");
                    cerr << "\n[determinism] repeat  0 = reference";
                }
                else
                {
                    double max_abs = 0., signed_at_max = 0.;
                    size_t worst = 0;

                    for (nit = flow_.NodesBegin(); nit != flow_.NodesEnd(); ++nit)
                    {
                        const double d = (*nit)->Read(p_key_) - (*nit)->Read(ref_key_);
                        if (fabs(d) > max_abs)
                        {
                            max_abs       = fabs(d);
                            signed_at_max = d;
                            worst         = (*nit)->Idx();
                        }
                    }

                    cerr << "\n[determinism] repeat " << rep
                         << "  max|dp| = "  << max_abs
                         << "  (signed "    << signed_at_max
                         << ", node "       << worst << ")"
                         << (max_abs == 0. ? "   BIT-IDENTICAL" : "");
                }
            }

            // Leave the model at t^n so the real step below is untouched.
            ResetVariables();
            model.CopyReplace("stash pressure", "fluid pressure");
            current_dt             = dt_entry;
            upwind_done_this_step_ = false;

            cerr << "\n[determinism] DONE — proceeding with the real step\n";
        }
    }


    // ── Select this step's dt ─────────────────────────────────────────
    // In "smooth decrease" mode (set after a reactive cut) force a steady 0.5%/step
    // reduction for a while. Otherwise ramp gently back up from last step's dt,
    // capping the per-step growth factor (tighter once dt > 1800 s) and never
    // exceeding largest_timestep.
    if (smooth_timestep_decrease)
        current_dt *= 0.995;
    else
    {
        timestep_increment_factor += 0.00025;
        if (timestep_increment_factor > 1.005 && current_dt > 1800)
            timestep_increment_factor = 1.005;
        if (timestep_increment_factor > 1.01)
            timestep_increment_factor = 1.01;
        current_dt = std::min(old_dt * timestep_increment_factor, largest_timestep);
    }


    bool do_loop(true);
    while (do_loop)
    {
        if (output_CVFEM_timing)
            cerr << endl << "SCHEME START with dt = " << current_dt << endl;
        cerr << "Timestep increment factor = " << timestep_increment_factor;
        cerr << endl << "Smooth dt decrease: ";
        if (smooth_timestep_decrease) cerr << "true "; else cerr << "false ";
        cerr << "\nSmooth counter: " << counter_smooth_timestep_decrease;
        // cerr << endl << "////////////////////////" << endl << endl;

        AdvectionDiffusionLoops();

        if (output_CVFEM_timing)
            cerr << ", TEMPERATURE SOLVE";

        // heat conduction
        auto start = std::chrono::high_resolution_clock::now();

        T_FE.TimeIncrement(current_dt);
        T_FE.IntegrateOver(model, model.Region("Model"));

        auto stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << " (" << duration.count() << ")";
        }

        //! Option to inject mass / energy / tracer
        // PointInjection();


        //! Point injection used in 2D Benchmarks
        //! input = 2D Benchmark tests as per Weis et al., 2014
        //! 202a, 202b, 203a, 203b
        if (with_point_injection_)
            PointInjectionBenchmarks(benchmark_);


        // Well -> reservoir exchange. Leaves a source term for the next pressure solve.
        if (with_well_calculations)
            WellReservoirExplicitCoupling();

        // If the exchange rejected the step (well solve failed, or a source term
        // was impossibly large) it has already cut dt; skip ahead and retry.
        // Otherwise carry on to fluid-rock equilibration.
        if (!well_coupling_rejected)
        {
            // check for pressure below 1atm
            model.MinMaxOf("fluid pressure", min_value, max_value);
            if ((min_value < minimum_pressure || max_value > 490.e6) &&
                brick_wall_limiter)
            {

                // write to Logfile
                std::fstream file("PressureLimiter.log", std::ios::in | std::ios::out); // open file in read and write mode
                file.seekp(0, std::ios::end); // move the write pointer to the end of the file to append new entries
                file << "Pressure limiter fluid";
                file.close();
                cerr << "\nApplying pressure limiter fluid" << endl;
                model.Accept(pressure_limiter_fluid);
            }

            auto start = std::chrono::high_resolution_clock::now();
            FluidRockEquilibration();
            auto stop = std::chrono::high_resolution_clock::now();
            if (output_CVFEM_timing)
            {
                cerr << ", EQUILIBRATION";
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                cerr << " (" << duration.count() << ")";
            }
            do_loop = false;
        }
    }

    // JK26: needs to be conditional ... only if tracer is used
    //    model.CopyReplace("tracer content liquid", "tracer content fluid"); // TEST

    if(model.ContainsRegion("OCEAN"))
    {
        printRangeOfVariable(model,"OCEAN","temperature");
        printRangeOfVariable(model,"OCEAN","salinity");
    }

    if (smooth_timestep_decrease)
        counter_smooth_timestep_decrease++;
    if (counter_smooth_timestep_decrease > 20)
    {
        counter_smooth_timestep_decrease = 0;
        smooth_timestep_decrease = false;
    }

    old_dt = current_dt;
    // PAUSE();
    return current_dt;
} // end Apply

/** book keeping of variables for transient calculations */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::AdvanceTransientVariables()
{
    model.CopyReplace(
        "previous fluid pressure",
        "previous previous fluid pressure"); // needs to be done separately and in
    // the right order!
    for (unsigned i = 0; i < names.transient_variables.size(); i++)
        model.CopyReplace(names.transient_variables[i].c_str(),
                          ("previous " + names.transient_variables[i]).c_str());
    for (unsigned i = 0; i < names.full_reset_variables.size(); i++)
        model.CopyReplace(names.full_reset_variables[i].c_str(),
                          ("reset " + names.full_reset_variables[i]).c_str());
} // end AdvanceTransientVariables

/** reset variables for transient pressure calculations */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ResetVariables()
{
    for (unsigned i = 0; i < names.transient_variables.size(); i++)
        model.CopyReplace(("previous " + names.transient_variables[i]).c_str(),
                          names.transient_variables[i].c_str());
    //    model.CopyReplace("previous previous fluid pressure","previous fluid
    //pressure");//needs to be done separately and in the right order!
    //    model.CopyReplace("previous previous total enthalpy","previous total
    //enthalpy");//needs to be done separately and in the right order!
} // end Reset

// ─────────────────────────────────────────────────────────────────────────────
// ResetAndCutTimestep — shared "a step failed: rewind and retry smaller" routine.
//
// Call sites (all now share this routine; they used to duplicate it with small
// variations):
//   - Apply                   : expected-dp guard — the only one AFTER equilibration
//   - AdvectionDiffusionLoops : mass-based control_dt exceeded (carries transport's
//                               own reason string when it tagged the return)
//   - AdvectionDiffusionLoops : split-boundary leakage requested a cut
//   - PressureLoop            : CFL (Courant) limit exceeded
//   - ApplyWellReservoirExchange : a well solve failed, or a source term was
//                                  impossibly large
//
// Common work: reset the wells, restore the transient variables to the start of
// the step, rescale the explicit pressure source nQ so it stays consistent with
// the smaller dt, assign the new dt and restart the growth ramp.
//
// Per-site differences are passed in:
//   new_dt          reduced timestep to retry with (caller computes it from
//                   control_dt / cfl_dt / current_dt and the site's safety factor)
//   reason          short label for the log line
//   smooth_decrease value assigned to smooth_timestep_decrease
//   reset_counter   if true, counter_smooth_timestep_decrease = 0
//   scale_nfvs      rescale nQ by new_dt/current_dt. NOTE: the rescaling itself
//                   ALWAYS happens — every one of the five call sites takes the
//                   default (true), so nQ is kept consistent with the reduced dt
//                   on every cut. What is dead is only the OFF path: the one site
//                   that once passed false (a split-boundary cut inside
//                   PressureLoop) no longer exists. The parameter survives only
//                   because the method is public API.
//   pause           block on PAUSE() before resetting (debug hook). Also never
//                   overridden by any current caller; kept for interactive use.
//
// NOTE: current_dt is read for the nfvs ratio and then overwritten with new_dt,
//       so callers must pass new_dt computed against the *current* current_dt.
// ─────────────────────────────────────────────────────────────────────────────
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ResetAndCutTimestep(
    double new_dt, const char* reason,
    bool smooth_decrease, bool reset_counter,
    bool scale_nfvs, bool pause)
{
    cerr << endl << "//////////////// RESETTING (" << reason << ") ////////////////" << endl;
    if (pause) { PAUSE(); PAUSE(); }

    // Log the cut while current_dt still holds the PRE-cut value.
    // LogTimestepCut(reason, new_dt);

    if (with_well_calculations)
        for (auto &well : well_models)
            well->Reset();

    ResetVariables();

    if (scale_nfvs)
        nfvs_scaling(new_dt / old_dt);   // keep nQ consistent with the reduced dt

    current_dt                = new_dt;
    timestep_increment_factor = 1.;          // restart the gentle growth ramp
    smooth_timestep_decrease  = smooth_decrease;
    if (reset_counter)
        counter_smooth_timestep_decrease = 0;

    cerr << endl << "Current_dt = " << current_dt;
}

// ── LogTimestepCut ───────────────────────────────────────────────────────────
// Append one CSV row per dt cut. Columns:
//   sim_time_s    : model_time_ reached when the cut happened (sim seconds)
//   sim_time_yr   : same in years (3.15576e7 s/yr) for quick reading
//   reason        : the cut reason string (control_dt exceeded / split boundary
//                   (advection) / CFL / ...) — same string passed to the reset
//   dt_before_s   : current_dt at the moment of the cut (the dt that failed)
//   dt_after_s    : new_dt we are retrying with
//   ratio         : dt_after / dt_before (how hard the cut was)
// Opened lazily on the first cut so a clean run writes no file. Flushed each row
// so a crash still leaves a complete log.
template<uint32_t dim>
void CVFEM_PHX_Scheme<dim>::LogTimestepCut(
    const char* reason, double new_dt)
{
    if (!dt_cut_log_.is_open())
    {
        // Filename: dt_cuts_<label>_<startup-stamp>.csv. <label> from SetCutLogLabel
        // (e.g. output_name); <startup-stamp> captured at construction (run start).
        // If no label was set, omit it (dt_cuts_<stamp>.csv).
        std::string fname = "dt_cuts_";
        if (!dt_cut_log_label_.empty()) fname += dt_cut_log_label_ + "_";
        fname += dt_cut_log_stamp_ + ".csv";

        dt_cut_log_.open(fname.c_str(), std::ios::out | std::ios::trunc);
        if (dt_cut_log_.is_open())
            dt_cut_log_ << "sim_time_s,sim_time_yr,reason,dt_before_s,dt_after_s,ratio\n";
    }
    if (!dt_cut_log_.is_open()) return;   // open failed; skip silently

    const double yr  = model_time_ / 3.15576e7;
    const double rat = (current_dt > 0.) ? new_dt / current_dt : 0.;

    dt_cut_log_ << model_time_ << ',' << yr << ','
                << '"' << reason << '"' << ','
                << current_dt << ',' << new_dt << ',' << rat << '\n';
    dt_cut_log_.flush();
}


/**
 * AdvectionDiffusionLoops — outer transport loop.
 *
 * Each pass: solve pressure (PressureLoop, which itself retries on CFL), then
 * advect mass with that velocity field. AdvectMassConserved returns control_dt,
 * the largest dt the mass-conservative advection would have accepted. If the dt
 * actually used exceeds control_dt the step is invalid → rewind and retry smaller.
 * The split-region leakage path can independently request a cut.
 *
 * Exits only when neither the mass-based (control_dt) nor the split-boundary
 * criterion forces a cut.
 */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::AdvectionDiffusionLoops()
{
    bool do_loop(true);
    bool cut_dt_splitB;

    while (do_loop)
    {
        cut_dt_splitB = false;

        PressureLoop();

        auto start = std::chrono::high_resolution_clock::now();

        // // First-step pressure copy (split boundaries): must precede whoever runs
        // // the leakage — which may now be transport itself (leakage-in-the-loop),
        // // so it runs before AdvectMassConserved rather than inside the scheme's
        // // leakage block below.
        // if ( with_split_region_ && model_time_ < current_dt )
        // {
        //     model.CopyReplace("fluid pressure", "previous fluid pressure");
        //     cerr << endl << "Copying fluid pressure into previous fluid pressure once at initial step!";
        // }

        control_dt = transport.AdvectMassConserved(current_dt);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        if (output_CVFEM_timing) cerr << ", TRANSPORT" << " (" << duration.count() << ")" << endl;

        // Split_region_leakage here do the full calculation of potential gradient (new since pressure has just been solved), etc...
        // Then it performs advection (ComputeGravity_PressureSourceTerm(perform_advection = true)
        // scheme_level mode only: when transport owns the leakage (end_of_stage
        // or per_substep), it has already run inside AdvectMassConserved and
        // this block is skipped; cut_dt_splitB then stays false (the budget
        // clamp never requests cuts).
        if (with_split_region_ && !transport.HandlesLeakage()
            && !(current_dt > control_dt))
        {
            start = std::chrono::high_resolution_clock::now();
            Split_region_leakage->SetTimeStep(current_dt);
            cut_dt_splitB = Split_region_leakage->ComputeGravity_PressureSourceTerm();
        }

        stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << endl << "Leakage: advection, SB LEAKAGE" << " (" << duration.count() << ")" << endl;
        }


        // Noteworthy-but-not-cut stage events -> CSV as a non-reset row
        // (ratio 1), so carryover spikes etc. are greppable next to cuts.
        // if (!(current_dt > control_dt) && !transport.StageEvents().empty())
        // LogTimestepCut(transport.StageEvents().c_str(), current_dt);

        if (current_dt > control_dt)
        {
            // Mass criterion failed: retry at 0.95 * control_dt.
            // If transport tagged this return as a rejection of its own (dry-out
            // redo, or a bounded exit such as the sub-step cap), carry its reason
            // string into the dt-cut CSV so the log distinguishes those from an
            // ordinary over-drain cut; otherwise fall back to a generic label.
            const std::string cut_reason =
                transport.TransportCutReason().empty() ? std::string("control_dt exceeded")
                                                       : transport.TransportCutReason();
            cerr << "!!!current_dt > control_dt!!! current_dt: " << current_dt
                 << ", control_dt: " << control_dt;
            ResetAndCutTimestep(control_dt * 0.95, cut_reason.c_str(),
                                /*smooth_decrease*/true, /*reset_counter*/true);
        }
        else if (cut_dt_splitB)
        {
            // Split-boundary leakage requested a cut: retry at 0.95 * current_dt.
            // Append the leakage's own reason (which mode/node/fraction triggered)
            // so the dt-cut log records WHY the SB cut happened, not just that it did.
            const std::string sb_reason =
                std::string("split boundary (advection): ") +
                Split_region_leakage->GetCutReason();
            ResetAndCutTimestep(current_dt * 0.95, sb_reason.c_str(),
                                /*smooth_decrease*/true, /*reset_counter*/true);
        }
        else
        {
            do_loop = false;   // both criteria satisfied — step accepted
        }
    }
} // end AdvectionDiffusionLoops

/**
 * PressureLoop — inner pressure-diffusion loop.
 *
 * Each pass: rebuild upwind directions, solve the implicit pressure equation,
 * recompute K·gradP and the resulting velocities, then check the CFL (Courant)
 * limit. If current_dt exceeds the CFL-allowed cfl_dt the velocities are too
 * large for a stable explicit transport step → rewind and retry at 0.99·cfl_dt.
 * Exits once current_dt <= cfl_dt.
 */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::PressureLoop()
{
    bool do_loop(true);

    while (do_loop)
    {
        auto start = std::chrono::high_resolution_clock::now();
        const bool do_upwind = !upwind_freeze_within_step || !upwind_done_this_step_;
        if (do_upwind)
        {
            upwind_control.Reset();
            if(model.ContainsRegion("FLOW REGION"))
                model.Region("FLOW REGION").Accept(upwind_control);
            else
                model.Region("Model").Accept(upwind_control);
            upwind_done_this_step_ = true;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << endl << "UPWIND" << (do_upwind ? "" : "(frozen)") << " (" << duration.count() << ")" << endl;
        }

        // Calculate pressure in the reservoir
        start = std::chrono::high_resolution_clock::now();
        P_FE.TimeIncrement(current_dt);
        iface_transferLHS_p.UpdateTimeIncrement(current_dt);

        P_FE.IntegrateOver(model, model.Region("Model"));

        stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << ", PRESSURE SOLVE (" << duration.count() << ")";
        }

        // check for pressure below zero
        model.MinMaxOf("fluid pressure", min_value, max_value);

        if ((min_value < 0. || max_value > 1000.e6) && brick_wall_limiter)
        {
            // write to Logfile
            std::fstream file("PressureLimiter.log", std::ios::in | std::ios::out); // open file in read and write mode
            file.seekp(0, std::ios::end); // move the write pointer to the end of the file to append new entries
            file << "Pressure limiter transport";
            file.close();

            cerr << "\nApplying pressure limiter transport" << endl;
            model.Accept(pressure_limiter_transport);
        }

        // calculate K * grad P
        start = std::chrono::high_resolution_clock::now();

        if(model.ContainsRegion("FLOW REGION"))
            model.Region("FLOW REGION").Accept(pres_grad);
        else
            model.Region("Model").Accept(pres_grad);

        stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << ", PRESSURE GRAD (" << duration.count() << ")";
        }

        // project pressure gradient onto facets
        start = std::chrono::high_resolution_clock::now();
        transport.UpdateProjection();

        if (output_CVFEM_timing) cerr << ", UPDATE PROJ";
        stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << " (" << duration.count() << ")";
        }

        // calculate velocities for Finite Volume calculations
        // or CVFEM_Visitors and determine CFL criterion
        // The VELOCITY pass runs on EVERY iteration, in both modes. It sits AFTER
        // UpdateProjection, so the facet velocities it reads were just produced by
        // THIS iteration's solve at THIS iteration's dt — they are never stale, and
        // freezing them was what let a dt-cut retry test its new dt against the
        // REJECTED solve's cfl_dt.
        //
        // This pass ALWAYS writes the transport map D1 (WithVelocity selects it),
        // in BOTH modes, on every solve including retries. D1 is never frozen:
        // transport's donor map always matches the velocity field whose flux rates
        // DetermineFacetFluxRatesOnce assembles.
        //
        // upwind_freeze_within_step controls ONLY the pressure map D0, by gating
        // whether the DIRECTION pass re-runs above:
        //   true  -> D0 fixed for the step  => every retry assembles the SAME
        //            pressure matrix at a smaller dt, and the solve is reproducible
        //   false -> D0 re-decided each iteration from STENCIL_DATA left by the
        //            REJECTED solve => the historical non-determinism.
        start = std::chrono::high_resolution_clock::now();
        upwind_control.WithVelocity();   // also selects the transport map (D1)
        upwind_control.ZeroVelocityOnFlip(upwind_zero_velocity_on_flip);
        upwind_control.ResetCflDiagnostics();
        if(model.ContainsRegion("FLOW REGION"))
            model.Region("FLOW REGION").Accept(upwind_control);
        else
            model.Region("Model").Accept(upwind_control);

        stop = std::chrono::high_resolution_clock::now();
        if (output_CVFEM_timing)
        {
            cerr << ", UPWIND-VEL" << (upwind_freeze_within_step ? "(donors frozen)" : "");
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cerr << " (" << duration.count() << ")";
        }

        model.ExtrapolateCellToNodeProperty("pore velocity", "nodal pore velocity");

        model.MinMaxOf("courant liquid", cfl_min_l, cfl_max_l);
        model.MinMaxOf("courant vapor", cfl_min_v, cfl_max_v);
        if( with_air_) model.MinMaxOf("courant air", cfl_min_a, cfl_max_a);

        cfl_dt = std::min(cfl_min_l, cfl_min_v);
        if( with_air_) cfl_dt = std::min(cfl_dt, cfl_min_a);

        if (output_CVFEM_timing)
        {
            const char* phase_name = upwind_control.WorstCflPhase() == 0 ? "liquid"
                                     : upwind_control.WorstCflPhase() == 1 ? "vapor"
                                                                           : "air";
            if ( upwind_control.WorstCfl() < largest_timestep )
                cerr << "\n[CFL] bind: elem=" << upwind_control.WorstCflElement()
                     << " phase=" << phase_name
                     << " sat=" << upwind_control.WorstCflSat()
                     << " v_mob=" << upwind_control.WorstCflVel()     // mobility-weighted flux, kr/mu * KgradP
                     << " v_pore=" << upwind_control.WorstCflPoreVel()
                     << " v_raw_Darcy="<< upwind_control.WorstCflVel() / upwind_control.WorstCflRelperm()
                     << " relperm_visc=" << upwind_control.WorstCflRelperm()
                     << " cfl_dt=" << upwind_control.WorstCfl() << endl;
            else
                cerr << "\n[CFL] no binding facet this pass: cfl larger than largest timestep" << endl;
        }

        // We either exit this loop because it has succeeded or retry with smaller
        // timestep
        if (current_dt > cfl_dt)
        {
            // CFL limit exceeded: retry at 0.99 * cfl_dt.
            cerr << endl << "Timestep exceeds reservoir CFL";
            ResetAndCutTimestep(cfl_dt * 0.99, "CFL exceeded",
                                /*smooth_decrease*/true, /*reset_counter*/true);
        }

        else
            do_loop = false;
    }
} // end PressureLoop

template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::WellReservoirExplicitCoupling()
{
    // WELL COUPLING IS EXPLICIT (single pass)
    well_source = 0.;

    for (auto &well : well_models)
    {
        well->SetTimeStep(current_dt);
    }

    well_coupling_rejected = false;

    // Top injection is applied every step for as long as a well is in that mode —
    // it delivers rate * dt each step, so it has to be re-applied. The MODE is
    // what selects it, not the value of the rate: a well leaves top injection by
    // being put under another control mode, never by being passed a rate of zero
    // or less. (Benoit DD/MM/YYYY)
    for (size_t i = 0; i < well_models.size(); ++i)
        if (well_models[i]->Control() == WellControl::TopInjection)
            well_models[i]->Top_Injection(injection_rates[i]);

    // ── PHASE 1 — solve the wells, compute their reservoir source terms ──────
    for (auto &well : well_models)
    {
        well_coupling_rejected = well->Apply(0);   // 0: pass index; there is only one pass
        if (well_coupling_rejected) break;         // well solve failed -> reject the step
    }

    if (!well_coupling_rejected)
    {
        for (auto &well : well_models)
        {
            ComputeWellSourceTerms(/*final_pass*/false, well->getName());
            if (well_coupling_rejected) break;     // source term impossibly large -> reject
        }
    }

    // ── PHASE 2 — accept: store the final source terms ───────────────────────
    if (!well_coupling_rejected)
    {
        for (auto &well : well_models)
        {
            ComputeWellSourceTerms(/*final_pass*/true, well->getName());
        }

        for (auto &well : well_models)
        {
            well->Write_results();
        }
    }

    if (well_coupling_rejected)
    {
        // Rejected: drop the well source so it isn't reused in the retried
        // pressure solve, then rewind and retry at 0.99 * current_dt.

        ResetAndCutTimestep(current_dt * 0.99, "well coupling failed",
                            /*smooth_decrease*/true, /*reset_counter*/true);
    }
}

/** thermal quilibration between fluid and rock */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::FluidRockEquilibration()
{
    equilibrator_properties.SetTimeIncrement(current_dt);
    model.Accept(equilibrator_properties);

} // end FluidRockEquilibration


template<uint32_t dim>
void CVFEM_PHX_Scheme<dim>::OpenBoundaries(
    double wt_top)
{
    open_top = true;
    equilibrator_properties.OpenBoundaries(wt_top);
}

/** switch brick wall limiter for fluid pressure on or off */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::SetBrickWallLimiterTo(
    bool limit)
{
    brick_wall_limiter = limit;
} // end SetBrickWallLimiterTo

/** add further variables for FV calculations */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::AddAdvectionVariable(
    const char *balanced_variable, const char *new_lhs_liquid,
    const char *new_rhs_liquid, const char *new_lhs_vapor,
    const char *new_rhs_vapor)
{
    names.AddAdvectionVariable(balanced_variable, new_lhs_liquid, new_rhs_liquid,
                               new_lhs_vapor, new_rhs_vapor);
    transport.AddAdvectionVariable(new_lhs_liquid, new_rhs_liquid, new_lhs_vapor,
                                   new_rhs_vapor);
} // end AddAdvectionVariable

/** modifying cfl criterion */
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::Adjust_CFL_Criterion(
    double scale_factor, bool take_pore_velocity)
{
    upwind_control.Adjust_CFL_Criterion(scale_factor, take_pore_velocity);
} // end Adjust_CFL_Criterion

/** access function to transient fluxes */
template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::GetFacetFluxFromInsideNodeToOutsideNode(Element<dim> &e,
                                                  unsigned int facet_idx,
                                                  double &flux_liquid,
                                                  double &flux_vapor,
                                                  double &flux_air)

{
    flux_liquid = fv_transport_liquid.GetFacetFlux(e, facet_idx, 0U);
    flux_vapor = fv_transport_vapor.GetFacetFlux(e, facet_idx, 0U);
    if(with_air_) flux_air = fv_transport_air.GetFacetFlux(e, facet_idx, 0U);
}

// PW May 2016
template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::SetEquilibratorConvergenceSpeedUpTo(bool boost)
{
    equilibrator_properties.SetConvergenceSpeedUpTo(boost);
}

// PW Sept 2018
template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::SetAvoidPressureOscillationsAtVLHTo(bool avoid_p_VLH)
{
    equilibrator_properties.SetAvoidPressureOscillationsAtVLHTo(avoid_p_VLH);
}

// BenoitLC add
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::WithRockLiquidusSolidus(bool with_rock_liquidus_solidus)
{
    equilibrator_properties.WithRockLiquidusSolidus(with_rock_liquidus_solidus);
}

// BenoitLC add
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::
    AddFluidContributionToHeatCapacity(
        bool add_fluid_contribution_to_heat_capacity)
{
    equilibrator_properties.AddFluidContributionToHeatCapacity(
        add_fluid_contribution_to_heat_capacity);
}

// BenoitLC add
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::
    AttemptToSurviveFluidPropertiesError(
        bool attempt_to_survive_fluid_properties_error_)
{
    equilibrator_properties.AttemptToSurviveFluidPropertiesError(
        attempt_to_survive_fluid_properties_error_);
}

// BenoitLC add
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::SetRockHeatCapacity(
    double mini_cp)
{
    equilibrator_properties.SetRockHeatCapacity(mini_cp);
}

// BenoitLC add
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::SetRockCrystallizationCurve(
    double nu_coefficient, double sigma1_coefficient,
    double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_)
{
    // In CVFEM_PHX_Scheme, this has to be used with
    // "WithRockLiquidusSolidus" set to "true" to take effect!
    equilibrator_properties.SetRockCrystallizationCurve(
        nu_coefficient, sigma1_coefficient, latent_heat_of_fusion, b_coefficient,
        crystallization_curve_);
}

////BB 2020
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::Initialize_well(
    const std::string &well_name)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            well->Initialize_Well();

            cerr << endl << "Initialized well " << well_name << endl;
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::Initialize_well-> ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::Set_well_top_pressure(const std::string &well_name,
                                const double &pressure)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            // cerr <<endl<< "Set top pressure for well " << well_name <<endl;
            well->Set_well_top_pressure(pressure);
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::Set_well_top_pressure. ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::Set_well_top_temperature(const std::string &well_name,
                                   const double &temperature)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            // cerr <<endl<< "Set top pressure for well " << well_name <<endl;
            well->Set_well_top_temperature(temperature);
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::Set_well_top_pressure. ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::Set_well_water_table(
    const std::string &well_name, const bool &toggle_water_table_displacement,
    const double &initial_water_table_depth)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            cerr << endl << "Set water table for well " << well_name << endl;
            well->Set_water_table(toggle_water_table_displacement,
                                  initial_water_table_depth);
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::Set_well_water_table. ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::Set_target_rate(
    const std::string &well_name, const bool &trigger,
    const double &ext_target_rate)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            cerr << endl<<"Set target production rate for well " << well_name << endl;
            well->Set_target_rate(trigger, ext_target_rate);
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::OutputWellResults. ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::Set_injection_rate(
    const std::string &well_name, const double &ext_injection_rate)
{
    auto it = std::find_if(
        well_models.begin(), well_models.end(),
        [&well_name](const std::unique_ptr<WellModelPrototype<dim>> &well_model)
        {
            return well_model->getName() == well_name;
        });

    if (it != well_models.end())
    {
        // Calculate the index based on the iterator position
        size_t index = std::distance(well_models.begin(), it);

        // Check if the index is within the size of injection_rates
        if (index < injection_rates.size())
        {
            if (ext_injection_rate <= 0.)
                throw csmp::Exception(ERROR,
                                      "CVFEM_PHX_Scheme::Set_injection_rate",
                                      well_name,
                                      "the top-injection rate must be positive. To stop top "
                                      "injection, put the well under another control mode with "
                                      "Set_well_top_pressure() or Set_target_rate().");

            // Store the rate and put the well into top-injection mode. It is
            // re-applied every step for as long as the well stays in that mode.
            injection_rates[index] = ext_injection_rate;
            (*it)->Top_Injection(ext_injection_rate);
            std::cerr << "Updated injection rate for well '" << well_name
                      << "' to: " << ext_injection_rate << std::endl;
        }
        else {
            std::cerr << "Injection rate index out of range for well: " << well_name
                      << std::endl;
            // Handle the error, e.g., throw an exception or use appropriate error
            // handling
        }
    }
    else {
        // Well with the specified name not found
        string message = "Well with the name " + well_name + " does not exist.";
        throw csmp::Exception(
            FATAL_ERROR,
            "CVFEM_PHX_Scheme<dim>::Set_injection_rate. ",
            message);
    }
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::OutputWellResults(
    const std::string &well_name, const std::string &filename)
{
    for (auto &well : well_models)
    {
        if (well->getName() == well_name)
        {
            well->WriteSolutionToFile(filename);
            return;
        }
    }

    // Well with the specified name not found
    string message = "Well with the name " + well_name + " does not exist.";
    throw csmp::Exception(
        FATAL_ERROR,
        "CVFEM_PHX_Scheme<dim>::OutputWellResults. ",
        message);
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::IncludeWellCalculations(
    const bool &ext_with_well)
{
    with_well_calculations = ext_with_well;
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::PointInjection()
{
    ScalarVariable ID;
    csmp::Index n_id_key(model.Database().StorageKey("node ID"));

    csmp::Index mtp_key_(model.Database().StorageKey("previous fluid density"));
    csmp::Index mt_key_(model.Database().StorageKey("fluid density"));
    csmp::Index pore_volume_key_(model.Database().StorageKey("pore volume"));
    csmp::Index bulk_volume_key_(model.Database().StorageKey("bulk volume"));
    csmp::Index sh_key_(model.Database().StorageKey("saturation halite"));
    csmp::Index ml_key_(model.Database().StorageKey("fluid mass liquid"));
    csmp::Index mv_key_(model.Database().StorageKey("fluid mass vapor"));
    csmp::Index bulk_fluid_density_key_(
        model.Database().StorageKey("bulk fluid density"));
    csmp::Index hCl_key_(model.Database().StorageKey("enthalpy content liquid"));
    csmp::Index hCv_key_(model.Database().StorageKey("enthalpy content vapor"));
    csmp::Index msp_key_(model.Database().StorageKey("previous mass salt"));
    csmp::Index lithium_content_liquid_key_(
        model.Database().StorageKey("lithium content liquid"));
    csmp::Index lithium_content_fluid_key_(
        model.Database().StorageKey("lithium content fluid"));
    csmp::Index tracer_content_liquid_key_(
        model.Database().StorageKey("tracer content liquid"));
    csmp::Index tracer_content_fluid_key_(
        model.Database().StorageKey("tracer content fluid"));
    csmp::Index mml_key_(model.Database().StorageKey("liquid mass mobility"));
    csmp::Index trml_key_(model.Database().StorageKey("liquid tracer mobility"));

    ScalarVariable mtp(ANY, 0.);
    //ScalarVariable mf(ANY, 0.);
    ScalarVariable ml(ANY, 0.);
    ScalarVariable mv(ANY, 0.);
    ScalarVariable bulk_fluid_density(ANY, 0.);
    ScalarVariable msp(ANY, 0.);
    ScalarVariable pore_volume(ANY, 0.);
    ScalarVariable sh(ANY, 0.);
    ScalarVariable mass_transfer_rate(ANY, 0.);
    ScalarVariable tracer_mass_transfer_rate(ANY, 0.);
    ScalarVariable energy_transfer_rate_l(ANY, 0.);
    ScalarVariable energy_transfer_rate_v(ANY, 0.);
    ScalarVariable salt_mass_transfer_rate(ANY, 0.);
    ScalarVariable fluid_mass(ANY, 0.);
    ScalarVariable fluid_mass_liquid(ANY, 0.);
    ScalarVariable enthalpy_content_liquid(ANY, 0.);
    ScalarVariable enthalpy_content_vapor(ANY, 0.);

    ScalarVariable lithium_content_liquid(ANY, 0.);
    ScalarVariable lithium_content_fluid(ANY, 0.);
    ScalarVariable tracer_content_liquid(ANY, 0.);
    ScalarVariable tracer_content_fluid(ANY, 0.);


    const Region<dim> &ref = model.Region("Model");
    typename vector<Node<dim> *>::const_iterator nit;

    double max_mass_tracer =
        -1; // COULD SET TRACER CONTROLS (timings + amounts) EXTERNALLY..
    // Negative value means unlimited tracer!!
    double tracer_start_time = 0.;

    cerr << endl << "POINTINJECTION";
    for (nit = ref.NodesBegin(); nit != ref.NodesEnd(); ++nit)
    {
        (*nit)->Read(n_id_key, ID);

        // if (ID() == 24176 or ID() == 24177 or ID() == 24178 or ID() == 24179
        //    or ID() == 24404 or ID() == 30035 or ID() == 30034 or ID() == 30033 or ID() == 30032){

        //1227 for no split test, 1279 for split // Benchmarks: 6067 split, 6991 no split// CENTER  split: 22408 no split: 23215

        if( (*nit)->AtBoundary()==BOTTOM)
        {
            //if (ID() == 25558 /*15924*/){

            //cerr << endl << "At node: " << ID();
            double mass_exchange_term(0.), energy_exchange_term_l(0.),
                energy_exchange_term_v(0.), salt_mass_exchange_term(0.),
                tracer_exchange_term(0.);

            double current_tot_mass_tracer_injected(0.),
                current_tot_mass_injected(0.),
                old_total_mass_tracer_injected =
                total_tracer_mass_injected; // total mass of tracer injected so
            // far

            (*nit)->Read(pore_volume_key_, pore_volume);
            (*nit)->Read(hCl_key_, enthalpy_content_liquid);
            (*nit)->Read(hCv_key_, enthalpy_content_vapor);
            (*nit)->Read(sh_key_, sh);

            // Now we deal with advective terms (potentially including lithium), which
            // occurs only at completions

            (*nit)->Read(mtp_key_, mtp);
            //(*nit)->Read(mf_key_, mf);
            (*nit)->Read(ml_key_, ml);
            (*nit)->Read(mv_key_, mv);
            (*nit)->Read(msp_key_, msp);

            mass_transfer_rate() = -0.000001*pore_volume();//-0.00005;
            //mass_transfer_rate() = -0.0001;
            //mass_transfer_rate() = -0.00001;

            energy_transfer_rate_l() = mass_transfer_rate() * 3500000;
            energy_transfer_rate_v() = 0.;
            salt_mass_transfer_rate() = 0.;

            mass_exchange_term = -mass_transfer_rate() * current_dt /
                                 (pore_volume() /* * (1.-sh()) */); // kg.m-3
            energy_exchange_term_l = -energy_transfer_rate_l() * current_dt /
                                     (pore_volume() /* * (1.-sh()) */); // J.m-3
            energy_exchange_term_v = -energy_transfer_rate_v() * current_dt /
                                     (pore_volume() /* * (1.-sh()) */); // J.m-3
            salt_mass_exchange_term = -salt_mass_transfer_rate() * current_dt /
                                      (pore_volume() /* * (1.-sh()) */); // kg.m-3

            // Production case fora tracer

            current_tot_mass_injected -= mass_transfer_rate() * current_dt; // kg

            // If we did not yet inject "max_mass_tracer" kg of tracer and it is time
            // to inject, we inject tracer

            if ((max_mass_tracer < 0. ||
                 old_total_mass_tracer_injected < max_mass_tracer) &&
                model_time_ >= tracer_start_time)
            {
                // Set tracer concentration here
                // current_tot_mass_tracer_injected -=
                //     mass_transfer_rate() * 0.2 *
                //     current_dt; // Inject 200g tracer per kg fluid
                current_tot_mass_tracer_injected -= mass_transfer_rate() * 1. *
                                                    current_dt;    // Inject 1kg tracer per kg fluid

            } // end control tracer

            // Check total tracer mass injected since start and limit to
            // "max_mass_tracer" kg Calculate total tracer mass injected in this step
            total_tracer_mass_injected += current_tot_mass_tracer_injected;

            if (max_mass_tracer > 0 && total_tracer_mass_injected > max_mass_tracer)
            {
                current_tot_mass_tracer_injected =
                    max_mass_tracer - old_total_mass_tracer_injected;
                total_tracer_mass_injected = max_mass_tracer;
            } // end limit tracer

            // We update transport variables with source terms
            mtp() += mass_exchange_term;
            enthalpy_content_liquid() += energy_exchange_term_l;
            enthalpy_content_vapor() += energy_exchange_term_v;
            msp() += salt_mass_exchange_term;

            // We store the results of advection

            (*nit)->Store(mtp_key_, mtp); // Trying to update liq and vap instead What
            // if reservoir reset due to cfl...
            (*nit)->Store(msp_key_, msp);

            // We store the results of heat conduction, not just at completions!
            (*nit)->Store(hCl_key_, enthalpy_content_liquid);
            (*nit)->Store(hCv_key_, enthalpy_content_vapor);

            // cerr << endl
            //      << "Current mass injected: " << current_tot_mass_injected
            //      << ", tracer: " << current_tot_mass_tracer_injected;

            // We finish dealing with tracer injection. Also, there is no lithium
            // injection case

            // Second loop to inject tracer, distributed based on fractional injection
            // For visualization purposes we output tracer_mass_transfer_rate for
            // production but not injection

            (*nit)->Read(pore_volume_key_, pore_volume);
            (*nit)->Read(tracer_content_liquid_key_, tracer_content_liquid);
            //(*nit)->Read(tracer_content_fluid_key_, tracer_content_fluid);

            double fraction_of_total_injection = 0.;
            if (current_tot_mass_injected > 0.)
                fraction_of_total_injection =
                    -mass_transfer_rate() * current_dt / current_tot_mass_injected;

            tracer_exchange_term = current_tot_mass_tracer_injected *
                                   fraction_of_total_injection /
                                   pore_volume(); // kg.m-3

            // cerr<<endl<<"Tracer mass injected = "<<current_tot_mass_tracer_injected
            // * fraction_of_total_injection<<
            //       ", fraction of total fluid mass injection =
            //       "<<fraction_of_total_injection;

            tracer_content_liquid() += tracer_exchange_term; // kg.m-3
            // tracer_content_fluid()  = tracer_content_liquid(); // kg.m-3
            // We consider that the tracer is only present in the liquid water phase

            // We store results of tracer injection
            (*nit)->Store(tracer_content_liquid_key_, tracer_content_liquid);
            (*nit)->Store(tracer_content_fluid_key_, tracer_content_liquid);

            // end well node loop 2
        }
    }

    // Calculation of total tracer mass in model for info and tracer mobility
    // calculation
    double totalmasstracer(0.); // for monitoring total mass tracer in model
    for (nit = ref.NodesBegin(); nit != ref.NodesEnd(); ++nit)
    {
        (*nit)->Read(ml_key_, ml);
        (*nit)->Read(mtp_key_, mtp);
        (*nit)->Read(pore_volume_key_, pore_volume);
        (*nit)->Read(tracer_content_liquid_key_, tracer_content_liquid);

        totalmasstracer += tracer_content_liquid() * pore_volume();

    } // end model node loop
    cerr << endl << "************Total mass tracer = " << totalmasstracer;
}


template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ActivatePointInjectionBenchmarks(const char* val) {
    std::string s(val);

    if (s == "test-202a" || s == "test-202b" ||
        s == "test-203a" || s == "test-203b")
    {
        with_point_injection_ = true;
        benchmark_ = s;
    }

    std::cerr<<"PointInjection activated: ";
    if (with_point_injection_)  std::cerr<<"true \nBenchmark: "<< s << std::endl;
    else  std::cerr<<"false"<<std::endl;
}


template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::PointInjectionBenchmarks(std::string test) {

    ScalarVariable ID, x;
    csmp::Index mtp_key_(model.Database().StorageKey("previous fluid density"));
    csmp::Index pore_volume_key_(model.Database().StorageKey("pore volume"));
    csmp::Index hCl_key_(model.Database().StorageKey("enthalpy content liquid"));
    csmp::Index hCv_key_(model.Database().StorageKey("enthalpy content vapor"));
    csmp::Index msp_key_(model.Database().StorageKey("previous mass salt"));

    ScalarVariable mtp(ANY, 0.);
    ScalarVariable msp(ANY, 0.);
    ScalarVariable pore_volume(ANY, 0.);
    ScalarVariable mass_transfer_rate(ANY, 0.);
    ScalarVariable energy_transfer_rate_l(ANY, 0.);
    ScalarVariable energy_transfer_rate_v(ANY, 0.);
    ScalarVariable enthalpy_content_liquid(ANY, 0.);
    ScalarVariable enthalpy_content_vapor(ANY, 0.);
    ScalarVariable salt_mass_transfer_rate(ANY, 0.);

    const Region<dim>   &rref(model.Region("HOT_BOTTOM"));
    typename vector<Node<dim> *>::const_iterator nit;

    cerr << endl << "POINTINJECTION";
    const typename vector<Node< dim>*>::const_iterator NodesEnd( rref.NodesEnd() );
    for ( typename vector<Node<dim>*>::const_iterator nit( rref.NodesBegin() ); nit != NodesEnd; ++nit )
    {
        x() = (*nit)->x();
        if (x() == 0.){

            cerr << endl << "At x-coord: " << x();
            double mass_exchange_term(0.),
                energy_exchange_term_l(0.),
                energy_exchange_term_v(0.),
                salt_mass_exchange_term(0.);

            (*nit)->Read(pore_volume_key_, pore_volume);
            (*nit)->Read(hCl_key_, enthalpy_content_liquid);
            (*nit)->Read(hCv_key_, enthalpy_content_vapor);

            // Now we deal with advective terms (potentially including lithium), which
            // occurs only at completions
            (*nit)->Read(mtp_key_, mtp);
            (*nit)->Read(msp_key_, msp);

            //! Benchmark settings
            if (test == "test-202a")
            {
                cerr << "\nBenchmark: " << test << endl;
                mass_transfer_rate() = -0.04; // kg s-1 ... needs to be negative
                salt_mass_transfer_rate() = 0.0;

                energy_transfer_rate_l()  = mass_transfer_rate() * 500000.; // fluid enthalpy J kg-1
                energy_transfer_rate_v()  = 0.;
            }
            else if (test == "test-202b")
            {
                cerr << "\nBenchmark: " << test << endl;
                mass_transfer_rate() = -0.04; // kg s-1 ... needs to be negative
                salt_mass_transfer_rate() = 0.0;

                energy_transfer_rate_l()  = mass_transfer_rate() * 1500000.; // fluid enthalpy J kg-1
                energy_transfer_rate_v()  = 0.;
            }
            else if (test == "test-203a")
            {
                cerr << "\nBenchmark: " << test << endl;
                mass_transfer_rate() = -0.04; // kg s-1 ... needs to be negative
                salt_mass_transfer_rate() = -0.004;

                energy_transfer_rate_l()  = mass_transfer_rate() * 1363800.; // fluid enthalpy J kg-1
                energy_transfer_rate_v()  = 0.;
            }
            else if (test == "test-203b")
            {
                cerr << "\nBenchmark: " << test << endl;
                mass_transfer_rate() = -0.04; // kg s-1 ... needs to be negative
                salt_mass_transfer_rate() = -0.0004;

                energy_transfer_rate_l()  = mass_transfer_rate() * 3000000.; // fluid enthalpy J kg-1
                energy_transfer_rate_v()  = 0.;
            }

            mass_exchange_term = -mass_transfer_rate() * current_dt /
                                 (pore_volume() ); // kg.m-3
            energy_exchange_term_l = -energy_transfer_rate_l() * current_dt /
                                     (pore_volume() ); // J.m-3
            energy_exchange_term_v = -energy_transfer_rate_v() * current_dt /
                                     (pore_volume() ); // J.m-3
            salt_mass_exchange_term = -salt_mass_transfer_rate() * current_dt /
                                      (pore_volume() ); // kg.m-3


            // We update transport variables with source terms
            mtp() += mass_exchange_term;
            enthalpy_content_liquid() += energy_exchange_term_l;
            enthalpy_content_vapor() += energy_exchange_term_v;
            msp() += salt_mass_exchange_term;

            // We store the results of advection
            (*nit)->Store(mtp_key_, mtp); // Trying to update liq and vap instead What
            (*nit)->Store(msp_key_, msp);

            // We store the results of heat conduction, not just at completions!
            (*nit)->Store(hCl_key_, enthalpy_content_liquid);
            (*nit)->Store(hCv_key_, enthalpy_content_vapor);
        }
    }
}
template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::nfvs_scaling(double scaling_factor)
{
    csmp::Index nQ_key_(model.Database().StorageKey("nodal fluid volume source"));
    ScalarVariable nQ(ANY, 0.);

    const Region<dim> &ref = model.Region("Model");
    typename vector<Node<dim> *>::const_iterator nit;

    //test
    scaling_factor = 1;

    for (nit = ref.NodesBegin(); nit != ref.NodesEnd(); ++nit)
    {
        (*nit)->Read(nQ_key_, nQ);

        nQ()*=scaling_factor;

        (*nit)->Store(nQ_key_, nQ);
    }

    cerr<<endl<<"Scaling nQ by "<<scaling_factor;
}


template <uint32_t dim>
void CVFEM_PHX_Scheme <
    dim >::ComputeWellSourceTerms(bool final_pass,
                                 const std::string well_name)
{
    csmp::Index mtp_key_(model.Database().StorageKey("previous fluid density"));
    csmp::Index tot_mass_exchanged_key_(
        model.Database().StorageKey("total mass exchanged"));
    csmp::Index mt_key_(model.Database().StorageKey("fluid density"));
    csmp::Index pore_volume_key_(model.Database().StorageKey("pore volume"));
    csmp::Index bulk_volume_key_(model.Database().StorageKey("bulk volume"));
    csmp::Index sh_key_(model.Database().StorageKey("saturation halite"));
    csmp::Index well_mass_transfer_rate_key_(
        model.Database().StorageKey("well mass transfer rate"));
    csmp::Index well_tracer_mass_transfer_rate_key_(
        model.Database().StorageKey("well tracer mass transfer rate"));
    csmp::Index well_energy_transfer_rate_l_key_(
        model.Database().StorageKey("well energy transfer rate liquid"));
    csmp::Index well_energy_transfer_rate_v_key_(
        model.Database().StorageKey("well energy transfer rate vapor"));
    csmp::Index ml_key_(model.Database().StorageKey("fluid mass liquid"));
    csmp::Index mv_key_(model.Database().StorageKey("fluid mass vapor"));
    csmp::Index bulk_fluid_density_key_(
        model.Database().StorageKey("bulk fluid density"));
    csmp::Index hCl_key_(model.Database().StorageKey("enthalpy content liquid"));
    csmp::Index hCv_key_(model.Database().StorageKey("enthalpy content vapor"));
    csmp::Index well_salt_mass_transfer_rate_key_(
        model.Database().StorageKey("well salt mass transfer rate"));
    csmp::Index msp_key_(model.Database().StorageKey("previous mass salt"));
    csmp::Index Qr_key_(
        model.Database().StorageKey("well reservoir radial heat"));
    csmp::Index completion_key_(model.Database().StorageKey("well completion"));
    csmp::Index lithium_content_liquid_key_(
        model.Database().StorageKey("lithium content liquid"));
    csmp::Index lithium_content_fluid_key_(
        model.Database().StorageKey("lithium content fluid"));
    csmp::Index tracer_content_liquid_key_(
        model.Database().StorageKey("tracer content liquid"));
    csmp::Index tracer_content_fluid_key_(
        model.Database().StorageKey("tracer content fluid"));
    csmp::Index mml_key_(model.Database().StorageKey("liquid mass mobility"));
    // (Benoit 19/08/2026) added for the produced-liquid fraction, see production branch
    csmp::Index mmv_key_(model.Database().StorageKey("vapor mass mobility"));
    csmp::Index trml_key_(model.Database().StorageKey("liquid tracer mobility"));

    ScalarVariable mtp(ANY, 0.);
    ScalarVariable tot_mass_exchanged(ANY, 0.);
    ScalarVariable ml(ANY, 0.);
    ScalarVariable mv(ANY, 0.);
    ScalarVariable bulk_fluid_density(ANY, 0.);
    ScalarVariable msp(ANY, 0.);
    ScalarVariable pore_volume(ANY, 0.);
    ScalarVariable sh(ANY, 0.);
    ScalarVariable mass_transfer_rate(ANY, 0.);
    ScalarVariable tracer_mass_transfer_rate(ANY, 0.);
    ScalarVariable energy_transfer_rate_l(ANY, 0.);
    ScalarVariable energy_transfer_rate_v(ANY, 0.);
    ScalarVariable salt_mass_transfer_rate(ANY, 0.);
    ScalarVariable fluid_mass(ANY, 0.);
    ScalarVariable fluid_mass_liquid(ANY, 0.);
    ScalarVariable enthalpy_content_liquid(ANY, 0.);
    ScalarVariable enthalpy_content_vapor(ANY, 0.);
    ScalarVariable Qr(ANY, 0.);
    ScalarVariable completion(ANY, 0.);

    ScalarVariable lithium_content_liquid(ANY, 0.);
    ScalarVariable lithium_content_fluid(ANY, 0.);
    ScalarVariable tracer_content_liquid(ANY, 0.);
    ScalarVariable tracer_content_fluid(ANY, 0.);
    ScalarVariable mml(ANY, 0.);
    ScalarVariable mmv(ANY, 0.);   // (Benoit 19/08/2026)
    ScalarVariable trml(ANY, 0.);

    const Region<dim> &well_region_ref_ = model.Region(well_name);
    const Region<dim> &ref = model.Region("Model");
    typename vector<Node<dim> *>::const_iterator nit;

    double max_mass_tracer =
        200; // COULD SET TRACER CONTROLS (timings + amounts) EXTERNALLY..
    //max_mass_tracer =-1;
    // Negative value means unlimited tracer!!
    double tracer_start_time = 2 * 24 * 3600;
    tracer_start_time = 0;

    // If implicit, we only iteratively converge the pressure, not the
    // temperature! If implicit, if pressure near well(s) did not converge: we
    // update the mass source due to well-reservoir exchange we repeat pressure
    // solve with this new source (if coupling not set to explicit)

    // If pressure near well(s) converged we compute and store final source terms
    // (same for explicit or implicit)

    // But before anything else, we check if source terms are too large. This
    // would lead to a reset
    for (nit = well_region_ref_.NodesBegin(); nit != well_region_ref_.NodesEnd();
         ++nit)
    {
        (*nit)->Read(completion_key_, completion);
        if (completion() > 0)
        { // we limit computations to completions
            (*nit)->Read(pore_volume_key_, pore_volume);
            (*nit)->Read(sh_key_, sh);
            (*nit)->Read(mtp_key_, mtp);
            (*nit)->Read(well_mass_transfer_rate_key_, mass_transfer_rate);

            well_source = mass_transfer_rate() * current_dt;


            // We limit the source terms so as to avoid, for example, removing more
            // mass than is present in a control volume This also gives smoother
            // results
            double abs_fluid_source = abs(well_source);
            double max_source = 0.75 * mtp() * pore_volume() * (1. - sh());

            if (abs_fluid_source > max_source)
            {
                // If source terms are quite large we first try to smoothly reduce next
                // timesteps without reset

                cerr << endl << "SOURCE TERM FROM WELL " << well_name << " A BIT TOO LARGE";
                cerr << endl
                     << "Mass transfered: " << well_source
                     << ", vs Max: " << max_source;

                if (abs_fluid_source > mtp() * pore_volume() * (1. - sh()))
                {
                    // If source terms are definitelly too large we will reset and cut
                    // timestep
                    cerr << endl << "SOURCE TERM FROM WELL " << well_name << " DEFINITELY TOO LARGE";

                    well_coupling_rejected = true;
                    final_pass = false;
                }

                timestep_increment_factor = 1.;
                smooth_timestep_decrease = true;
                counter_smooth_timestep_decrease = 0;

            } // end check source terms too large
        } // end completions
    } // end well loop

    if (final_pass)
    {
        cerr << endl
             << "Well " << well_name << " converged, storing final source terms.";

        double mass_exchange_term(0.), energy_exchange_term_l(0.),
            energy_exchange_term_v(0.), salt_mass_exchange_term(0.), Qr_term(0.),
            lithium_exchange_term(0.), tracer_exchange_term(0.);

        double current_tot_mass_tracer_injected(0.), current_tot_mass_injected(0.),
            old_total_mass_tracer_injected =
            total_tracer_mass_injected; // total mass of tracer injected so far

        for (nit = well_region_ref_.NodesBegin();
             nit != well_region_ref_.NodesEnd(); ++nit)
        {
            (*nit)->Read(completion_key_, completion);
            (*nit)->Read(Qr_key_, Qr);
            (*nit)->Read(pore_volume_key_, pore_volume);
            (*nit)->Read(hCl_key_, enthalpy_content_liquid);
            (*nit)->Read(hCv_key_, enthalpy_content_vapor);
            (*nit)->Read(sh_key_, sh);

            // First we deal with conductive radial heat, which occurs along the whole
            // well
            Qr_term = Qr() * current_dt; // J (radial heat between well and reservoir)
            if (enthalpy_content_liquid() > 0.)
                enthalpy_content_liquid() +=
                    Qr_term / (pore_volume() /* * (1.-sh()) */); // J.m-3;
            else
                enthalpy_content_vapor() +=
                    Qr_term / (pore_volume() /* * (1.-sh()) */); // J.m-3;

            // Now we deal with advective terms (potentially including lithium), which
            // occurs only at completions
            if (completion() > 0)
            { // we limit computations to completions
                (*nit)->Read(mtp_key_, mtp);
                (*nit)->Read(ml_key_, ml);
                (*nit)->Read(mv_key_, mv);
                (*nit)->Read(msp_key_, msp);
                (*nit)->Read(tot_mass_exchanged_key_, tot_mass_exchanged);
                (*nit)->Read(well_mass_transfer_rate_key_, mass_transfer_rate);
                (*nit)->Read(well_energy_transfer_rate_l_key_, energy_transfer_rate_l);
                (*nit)->Read(well_energy_transfer_rate_v_key_, energy_transfer_rate_v);
                (*nit)->Read(well_salt_mass_transfer_rate_key_,
                             salt_mass_transfer_rate);

                mass_exchange_term = -mass_transfer_rate() * current_dt /
                                     (pore_volume() /* * (1.-sh()) */); // kg.m-3
                energy_exchange_term_l = -energy_transfer_rate_l() * current_dt /
                                         (pore_volume() /* * (1.-sh()) */); // J.m-3
                energy_exchange_term_v = -energy_transfer_rate_v() * current_dt /
                                         (pore_volume() /* * (1.-sh()) */); // J.m-3
                salt_mass_exchange_term = -salt_mass_transfer_rate() * current_dt /
                                          (pore_volume() /* * (1.-sh()) */); // kg.m-3

                // Production case for lithium and tracer
                if (with_lithium_)
                {
                    if (mass_transfer_rate() > 0)
                    {
                        (*nit)->Read(lithium_content_liquid_key_, lithium_content_liquid);
                        (*nit)->Read(lithium_content_fluid_key_, lithium_content_fluid);
                        (*nit)->Read(tracer_content_liquid_key_, tracer_content_liquid);
                        (*nit)->Read(mml_key_, mml);
                        (*nit)->Read(mmv_key_, mmv);

                        // LIQUID-BOUND (Benoit 19/08/2026): lithium and tracer are
                        // present in the LIQUID phase only - never in the vapor
                        // (both visitors store their vapor content and vapor
                        // mobility as hard zeros). That is exactly why the total
                        // produced rate must be split by f_l below before either
                        // species is extracted. Lithium additionally has an
                        // immobile solid store ("lithium solid"), which the well
                        // does not touch; the tracer has no solid store at all.

                        // ── PHASE SPLIT (Benoit 19/08/2026) ─────────────────────
                        // mass_transfer_rate() is the TOTAL produced fluid rate
                        // (liquid + vapor) — it is applied to mtp with no phase
                        // split. It was previously paired with a liquid-only
                        // concentration, so a two-phase producer extracted lithium
                        // and tracer using vapor as a carrier, over-extracting by
                        // mt/ml. Exact whenever mv == 0, which is why single-phase
                        // liquid runs never showed it.
                        //
                        // Corrected with the fractional-flow liquid share
                        // f_l = mml/(mml+mmv), the same split
                        // NaClH2OPropertiesVisitorPHX uses to form bfm_l for the
                        // boundary outflow, and consistent with SplitRegionLeakage,
                        // which resolves phase fluxes directly. The well model
                        // resolves phases for ENERGY only (energy_transfer_rate_l/_v),
                        // so no phase-resolved mass rate is available to reuse.
                        //
                        // BASIS: both species now use the LIQUID basis (content/ml)
                        // against a liquid carrier. For lithium the salt cancels
                        // exactly — (mass-salt)*LiCl/(ml-xCl) == mass_l*LiCl/ml — so
                        // salt_mass_transfer_rate is no longer needed here. This
                        // makes the well, the boundary, the FV rebuild and the SB
                        // leakage all structurally identical: carrier * content/ml.
                        //
                        // Guards: ml() > 0. was previously missing entirely, so a
                        // completion that had gone all-vapor divided by zero.
                        const double mm_tot_well = mml() + mmv();
                        const double f_l = (mm_tot_well > 0.) ? (mml() / mm_tot_well) : 0.;
                        // Produced liquid mass rate per unit liquid mass in place.
                        // ml() is folded in HERE so it never appears as a bare
                        // denominator below — guarding only the rate would still
                        // leave 0/0 when a completion has gone all-vapor.
                        const double liquid_rate_over_ml =
                            (ml() > 0.) ? (mass_transfer_rate() * f_l / ml()) : 0.;

                        lithium_exchange_term =
                            -liquid_rate_over_ml * lithium_content_liquid() * current_dt;
                        tracer_exchange_term =
                            -liquid_rate_over_ml * tracer_content_liquid() * current_dt;

                        // ── AVAILABILITY CLAMP (Benoit 19/08/2026) ──────────────
                        // A well must not produce more of a species than the cell
                        // actually holds. Without this the content can be driven
                        // negative and then silently clamped to 0 downstream, which
                        // CREATES mass.
                        //
                        // The failure is reachable because f_l is one step stale:
                        // mml/mmv are only recomputed in FluidRockEquilibration(),
                        // which runs AFTER this (Apply: advection -> well ->
                        // equilibration), so no fresh mobility exists at this point
                        // in the step. ml() and the contents, by contrast, are
                        // post-advection. If a completion dries within a step, the
                        // natural brake (f_l -> 0 as krl -> 0) lags the hazard
                        // (ml -> 0, so content/ml grows), and the extraction can
                        // exceed the holdup. Reordering does not remove this — it
                        // only changes which quantity is stale — so the consequence
                        // is bounded here instead.
                        //
                        // Same guard pattern as the boundary block in LithiumModel /
                        // TracerModel (li_mass_out_l > prev_dissolved) and as
                        // AdjustFluxOut in ExplicitFiniteVolumeTransportPHX.
                        // Exchange terms are in kg; content is kg.m-3 of pore volume.
                        const double lithium_available =
                            lithium_content_liquid() * pore_volume(); // kg
                        const double tracer_available =
                            tracer_content_liquid() * pore_volume();  // kg

                        if (-lithium_exchange_term > lithium_available)
                        {
                            cerr << endl
                                 << "Well " << well_name << ": lithium extraction capped at completion, "
                                 << "requested " << -lithium_exchange_term
                                 << " kg vs available " << lithium_available << " kg";
                            lithium_exchange_term = -lithium_available;
                        }
                        if (-tracer_exchange_term > tracer_available)
                        {
                            cerr << endl
                                 << "Well " << well_name << ": tracer extraction capped at completion, "
                                 << "requested " << -tracer_exchange_term
                                 << " kg vs available " << tracer_available << " kg";
                            tracer_exchange_term = -tracer_available;
                        }

                        // Derived from the CLAMPED term so the reported rate and the
                        // applied source can never disagree. (Benoit 19/08/2026)
                        tracer_mass_transfer_rate =
                            (current_dt > 0.) ? (tracer_exchange_term / current_dt) : 0.;

                        // We update lithium and tracer transport variables
                        lithium_content_liquid() +=
                            lithium_exchange_term /
                            (pore_volume() /* * (1.-sh()) */); // kg.m-3;
                        lithium_content_fluid() +=
                            lithium_exchange_term /
                            (pore_volume() /* * (1.-sh()) )*/); // kg.m-3;
                        tracer_content_liquid() += tracer_exchange_term / (pore_volume());
                        // tracer_content_fluid()      = tracer_content_liquid();

                        // We store the results of advection
                        (*nit)->Store(lithium_content_liquid_key_, lithium_content_liquid);
                        (*nit)->Store(lithium_content_fluid_key_, lithium_content_fluid);
                        (*nit)->Store(tracer_content_liquid_key_, tracer_content_liquid);
                        (*nit)->Store(tracer_content_fluid_key_, tracer_content_liquid);
                        (*nit)->Store(
                            well_tracer_mass_transfer_rate_key_,
                            tracer_mass_transfer_rate); // variable used for tracer
                        // production calculations only
                    }

                    if (mass_transfer_rate() < 0)
                    { // Injection case, here we just compute total injection
                        current_tot_mass_injected -=
                            mass_transfer_rate() * current_dt; // kg

                        // If we did not yet inject "max_mass_tracer" kg of tracer and it is
                        // time to inject, we inject tracer

                        if ((max_mass_tracer < 0. ||
                             old_total_mass_tracer_injected < max_mass_tracer) &&
                            model_time_ >= tracer_start_time)
                        {
                            // Set tracer concentration here

                            // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            //current_tot_mass_tracer_injected -= mass_transfer_rate() * 0.2 *current_dt; // Inject 200g tracer per kg fluid
                            current_tot_mass_tracer_injected -= mass_transfer_rate() * 1. *current_dt;    // Inject 1kg tracer per kg fluid
                            // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        } // end control tracer

                        // Check total tracer mass injected since start and limit to
                        // "max_mass_tracer" kg Calculate total tracer mass injected in this
                        // step
                        total_tracer_mass_injected += current_tot_mass_tracer_injected;

                        if (max_mass_tracer > 0 &&
                            total_tracer_mass_injected > max_mass_tracer)
                        {
                            current_tot_mass_tracer_injected =
                                max_mass_tracer - old_total_mass_tracer_injected;
                            total_tracer_mass_injected = max_mass_tracer;
                        } // end limit tracer
                    } // end if injection
                } // end with lithium aka ~ with_tracer

                if (well_coupling_rejected)
                {
                    cerr << endl << "WE SHOULD NOT HAVE REACHED HERE!!!!!!!!!!!!!!";
                    PAUSE();
                    PAUSE();
                }

                // We update transport variables with source terms
                mtp() += mass_exchange_term;
                enthalpy_content_liquid() += energy_exchange_term_l;
                enthalpy_content_vapor() += energy_exchange_term_v;
                msp() += salt_mass_exchange_term;
                tot_mass_exchanged() -= mass_transfer_rate() * current_dt;

                // We store the results of advection
                (*nit)->Store(mtp_key_, mtp); // Trying to update liq and vap instead
                // What if reservoir reset due to cfl...: not possible, cfl check is earlier
                (*nit)->Store(tot_mass_exchanged_key_, tot_mass_exchanged);
                (*nit)->Store(msp_key_, msp);
            } // end completions

            // We store the results of heat conduction, not just at completions!
            (*nit)->Store(hCl_key_, enthalpy_content_liquid);
            (*nit)->Store(hCv_key_, enthalpy_content_vapor);

        } // end well node loop 1

        cerr << endl
             << "Current mass injected: " << current_tot_mass_injected
             << ", tracer: " << current_tot_mass_tracer_injected;

        // We finish dealing with tracer injection. Also, there is no lithium
        // injection case
        if (with_lithium_)
        {
            // Second loop to inject tracer, distributed based on fractional injection
            // For visualization purposes we output tracer_mass_transfer_rate for
            // production but not injection
            for (nit = well_region_ref_.NodesBegin();
                 nit != well_region_ref_.NodesEnd(); ++nit)
            {
                (*nit)->Read(completion_key_, completion);

                if (completion() > 0)
                { // we limit computations to completions
                    (*nit)->Read(well_mass_transfer_rate_key_, mass_transfer_rate);

                    if (mass_transfer_rate() < 0)
                    { // Injection case
                        (*nit)->Read(pore_volume_key_, pore_volume);
                        (*nit)->Read(tracer_content_liquid_key_, tracer_content_liquid);
                        //(*nit)->Read(tracer_content_fluid_key_, tracer_content_fluid);

                        double fraction_of_total_injection = 0.;
                        if (current_tot_mass_injected > 0.)
                            fraction_of_total_injection = -mass_transfer_rate() * current_dt /
                                                          current_tot_mass_injected;

                        tracer_exchange_term = current_tot_mass_tracer_injected *
                                               fraction_of_total_injection /
                                               pore_volume(); // kg.m-3

                        // cerr<<endl<<"Tracer mass injected =
                        // "<<current_tot_mass_tracer_injected *
                        // fraction_of_total_injection<<
                        //       ", fraction of total fluid mass injection =
                        //       "<<fraction_of_total_injection;

                        tracer_content_liquid() += tracer_exchange_term; // kg.m-3
                        // tracer_content_fluid()  = tracer_content_liquid(); // kg.m-3
                        // We consider that the tracer is only present in the liquid water
                        // phase

                        if(tracer_content_liquid()>mtp()) cerr<<endl<<"More tracer than fluid mass???: tracer_content_liquid(): "<<tracer_content_liquid()<<", mtp(): "<<mtp();

                        // We store results of tracer injection
                        (*nit)->Store(tracer_content_liquid_key_, tracer_content_liquid);
                        (*nit)->Store(tracer_content_fluid_key_, tracer_content_liquid);
                    } // end injection
                } // end completions
            } // end well node loop 2

            // Calculation of total tracer mass in model for info and tracer mobility
            // calculation
            double totalmasstracer(0.); // for monitoring total mass tracer in model
            for (nit = ref.NodesBegin(); nit != ref.NodesEnd(); ++nit)
            {
                (*nit)->Read(pore_volume_key_, pore_volume);
                (*nit)->Read(tracer_content_liquid_key_, tracer_content_liquid);

                totalmasstracer += tracer_content_liquid() * pore_volume();
            } // end model node loop
            cerr << endl << "************Total mass tracer = " << totalmasstracer;
        } // end with_lithium, aka ~with_tracer
        // if(pause_at_end) PAUSE();
    } // end if converged
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::PAUSE() {
    cerr << " press ENTER to continue... ";
    std::cin.get();
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::SetModelTime(double time)
{
    model_time_ = time;
} // end

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ActivateAir(
    double humidity_air_top,    // Constant relative humidity [0–1]
    bool treat_air_and_vapor_as_mixture,
    bool treat_air_and_vapor_as_hydraulic_mixture,
    double rain_mm_per_year)
{
    with_air_ = true;

    // 1. names: add the 6 upwind-control entries (the push_backs gated on with_air)
    names.ActivateAir(thermodynamic_density_in_gravity_term_);

    // 2. upwind_control: allocate air matrices and register keys
    //    names.ActivateAir() must be called FIRST so the name strings exist
    upwind_control.ActivateAir(
        model,
        fv_transport_air,
        names.densities[2],       // "density air" or "density air transport"
        names.relperm_visc[2],    // "relperm viscosity air"
        names.saturations[2],     // "saturation air"
        names.time_control[2],    // "courant air"
        names.velocities[2],      // "velocity air"
        names.pore_velocities[2]  // "pore velocity air"
        );

    // 3. transport: flip the guard bool
    transport.ActivateAir();

    // 4. equilibrator setter
    equilibrator_properties.ActivateAir(humidity_air_top, treat_air_and_vapor_as_mixture, treat_air_and_vapor_as_hydraulic_mixture, rain_mm_per_year);

    // 5. PDE operators
    iface_transferLHS_p.ActivateAir(model); //split boundary

    conductance_p_upwind_air.MultiplyWithTimeIncrement(true);
    grav_air.MultiplyWithTimeIncrement(true);
    grav_air.AddAccumulateLater();
    if (with_gravity_) P_FE.Add(&grav_air);
    P_FE.Add(&conductance_p_upwind_air);

    // cerr << "\n[Scheme] ActivateAirPhase called."
    //      << "\n  treat_as_mixture: " << (treat_air_and_vapor_as_mixture ? "ON" : "OFF") << "\n";
}

template <uint32_t dim>
void CVFEM_PHX_Scheme<dim>::ActivateSplitRegionCoupling(
    bool open_space,
    double fault_perm_anisotropy,
    bool with_air,
    bool with_tracer,
    bool with_lithium,
    bool with_magmatic_fluids,
    bool with_gold)
{
    // NOTE: if air phase is used, ActivateAirPhase() MUST be called before this.
    // with_air_ is read here to configure SplitRegionLeakage with air.

    with_split_region_      = true;
    open_space_             = open_space;
    fault_perm_anisotropy_  = fault_perm_anisotropy;

    // Configure iface_transferLHS_p with the now-known parameters
    iface_transferLHS_p.SetSplitRegionParameters(open_space, fault_perm_anisotropy);

    // PDE operator setup — was in constructor body before
    iface_transferLHS_t.AddAccumulate();
    iface_transferLHS_t.MultiplyWithTimeIncrement(true);
    T_FE.AddSplitBoundaryIntegral(&iface_transferLHS_t);

    iface_transferLHS_p.AddAccumulate();
    iface_transferLHS_p.MultiplyWithTimeIncrement(true);
    P_FE.AddSplitBoundaryIntegral(&iface_transferLHS_p);


    Split_region_leakage = new SplitRegionLeakage<dim>(model, open_space_, fault_perm_anisotropy_);
    Split_region_leakage->WithSpeciesTransport(model, with_air, with_magmatic_fluids, with_tracer, with_gold, with_lithium);

    // Leakage-in-the-loop: hand the leakage to the transport stage. Whether
    // transport actually runs it (end_of_stage / per_substep) or leaves it to
    // this scheme (legacy) is decided by tp_full_subcycling in
    // ThreePhaseTransportPHX.cpp; transport.HandlesLeakage() reports the split.
    transport.SetInterfaceExchange( Split_region_leakage );

    cerr << "\n=== ActivateSplitRegionCoupling ===";
    cerr << "\n  with_lithium:          " << (with_lithium          ? "ON" : "OFF");
    cerr << "\n  with_air:              " << (with_air              ? "ON" : "OFF");
    cerr << "\n  with_tracer:           " << (with_tracer           ? "ON" : "OFF");
    cerr << "\n  with_magmatic_fluids:  " << (with_magmatic_fluids  ? "ON" : "OFF");
    cerr << "\n  with_gold:             " << (with_gold             ? "ON" : "OFF");
    cerr << "\n  open_space:            " << (open_space            ? "ON" : "OFF");
    cerr << "\n  fault_perm_anisotropy: " << fault_perm_anisotropy;
    cerr << "\n===================================\n";
}

template class CVFEM_PHX_Scheme<1U>;
template class CVFEM_PHX_Scheme<2U>;
template class CVFEM_PHX_Scheme<3U>;

} // end namespace csmp
