// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

// =====================================================================
//  Driver for:
//  "Numerical Modeling of Geothermal Heat and Lithium Co-Production
//   in Fault-Hosted Reservoirs"  (Lamy-Chappuis et al., G-Cubed)
//
//  Simulation platform: OpenCSMP++ (Matthai et al., 2026)
//  Method: 3D control-volume finite-element (CVFEM) scheme with the
//          fault represented as a discrete discontinuity (split
//          boundary). The coupled scheme solves fluid mass, salt,
//          energy and lithium transport with a wellbore model, using
//          the H2O-NaCl equation of state for brine properties.
//
//  This minimal driver runs a single calibrated case (a geothermal
//  doublet at Rittershoffen) for a 50-year forecast. Case-defining
//  permeabilities and porosities are set per region in SetPermeability()
//  and SetPorosity(); see Table 2 of the paper for the case values.
//
//  Required input files (same directory):
//    - WellTracesAndFault6C.asc / .dat   (ANSYS-ICEM mesh)
//    - WellTracesAndFault6C-regions.txt (list of regions in mesh)
//    - PhysicalVariablesBenchmarks.txt   (material property table)
//    - WellTracesAndFault6C-wells.txt    (well configuration, one [well] per well)
//    - CVFEM_fault_well_lithium_example.h                        (header)
//
//  Linear solvers: chosen near the top of Run() (search "LINEAR SOLVER
//  SELECTION"), separately for the reservoir and for the well Newton
//  Jacobian. Both are echoed at startup.
// =====================================================================

#include <chrono>
#include <fstream>
#include <map>
#include "CVFEM_fault_well_lithium_example.h"
#include "VTU_Interface.h"
#include "LithiumModel.h"
#include "TracerModel.h"
#include "CVFEM_PHX_Scheme.h"
#include "PoreVolumeVisitor.h"
#include "ModelTime.h"
#include "ModelTimeToInteger.h"
#include "VTK_Interface.h"
#include "PropertyHandle.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "ANSYS_Model3D.h"
#include "time.h"
#include "Boundary.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "CVFEM_SolverChoice.h"
#include "WellConfigurationFile.h"
#include "CSMP_highLevelUtilities.h"   // printRangeOfVariable
#include "SplitBoundaryInterface.h"

#include <iostream>

using namespace std;
namespace csmp {
void CVFEM_fault_well_lithium_example::Specifications() {
    SetTitle("CVFEM Rittershofen-like with fault, wells and lithium");
    SetDifficulty(0);
    SetCategory("CVFEM Examples");
    AddAuthor("BLC");
    AddDescription("source file in: CVFEM_fault_well_lithium_example");
}

// Helper functions are defined below Run().

void CVFEM_fault_well_lithium_example::Run() {
    // Element "thickness" is the out-of-plane extent used to turn an element
    // integral into a volume. Full-dimensional elements get `thickness`; the
    // lower-dimensional fault elements (and, after the split, the middle fault
    // region) get `thickness_LD`, which is the fault aperture — it is what scales
    // their conductance to a real conduit rather than a bare surface.
    double thickness_LD = 1.;   // fault aperture [m], lower-dimensional elements
    double thickness = 1.;      // out-of-plane thickness [m], normal elements

    // ---------------------------------------------------
    // Problem configuration

    bool with_lithium(true);        // solve lithium transport and partitioning
    bool with_well(true);           // couple the INJECTOR/PRODUCER wellbore models
    double heat_flux_bottom(0.11);  // basal heat flux [W/m2], applied at BOTTOM
    double salinity(0.);            // background salt mass fraction [-]; the
    // depth profile is applied later by
    // SpecialSalinity()

    // ---------------------------------------------------
    // ANSYS-ICEM mesh and config files

    // Name of mesh files (without .dat and .asc extensions)
    string icem_mesh_file("CVFEM_fault_well_lithium_example_mesh");

    string output_name = "CVFEM_fault_well_lithium_example_output";

    string region_file("CVFEM_fault_well_lithium_example");

    string phys_var_file("PhysicalVariablesBenchmarks.txt");  // Physical variables file

    // Well configuration, named after the mesh like the regions file:
    //   <mesh>-regions.txt   list of regions
    //   <mesh>-wells.txt     one [well] section per well
    // Generate a commented starting point with
    // WellConfigurationFile::WriteTemplate(well_file).
    string well_file("CVFEM_fault_well_lithium_example-wells.txt");

    bool no_std_output(false);  // suppress the per-step reservoir output
    double &model_time(ModelTime::Instance().modelTime);

    if (with_well) {
        wells_list.push_back("INJECTOR");
        wells_list.push_back("PRODUCER");
    }

    // Every duration below is expressed in `time_unit` and converted to seconds
    // further down, by time_multiplier. Change the unit and they all follow.
    TimeUnit time_unit(TimeUnit::DAYS);
    double max_time(50);              // total simulated time
    long vtk_increment(1);            // reservoir VTU written every ...
    long well_vtk_increment(1);   // well VTU written every ...
    double largest_time_step(0.1);    // upper bound on dt; the scheme cuts below
    // it as CFL and the leakage limiter require
    double initial_time_step(0.01);   // dt for the first step

    // Safety factor on the CFL-limited timestep: 1 takes the full CFL step, less
    // is more conservative.
    double cfl_scaling(0.8);

    bool with_gravity(true);        // gravity in the flow and transport terms

    // TOP boundary conditions, also the starting values for the whole model.
    double temperature_top(15.);    // oC
    double pressure_top(1.01325e5); // Pa (1 atm)

    // Background rock properties. Per-region values are applied afterwards by
    // SetPermeability() and SetPorosity(); these are what any region they do not
    // name will keep.
    double thermal_conductivity(2.);       // W/m/K — the UNIFORM value used by the
    // transient run. The initial geotherm
    // is built on layered values instead;
    // see "INIT TEMPERATURE" below.
    double heat_capacity_rock(830.);       // J/kg/K
    double density_rock(2650.);            // kg/m3
    double compressibility_rock(1.0e-10);  // 1/Pa
    double porosity(0.05);                 // [-], background value

    // Bulk heat capacity of the rock skeleton [J/m3/K].
    double total_heat_capacity(heat_capacity_rock * density_rock * (1. - porosity));

    // ---------------------------------------------------
    // Time variables and objects
    const double year = 31536000.;  //  in sec
    const double day = 86400.;      //  in sec
    const double hour = 3600.;      //  in sec
    const double minute = 60.;      //  in sec

    double time_multiplier = 1.;
    if (time_unit == TimeUnit::YEARS) time_multiplier = year;
    else if (time_unit == TimeUnit::DAYS) time_multiplier = day;
    else if (time_unit == TimeUnit::HOURS) time_multiplier = hour;
    else if (time_unit == TimeUnit::MINUTES) time_multiplier = minute;

    largest_time_step *= time_multiplier;
    initial_time_step *= time_multiplier;
    max_time *= time_multiplier;

    // Everything InitialiseProperties() needs, gathered in one place so the
    // extracted setup does not need a fifteen-argument signature.
    // NOTE: built HERE, after the time-unit scaling above — largest_time_step is
    // in seconds only once it has been multiplied by time_multiplier.
    const RunParameters params{ thickness, thickness_LD, pressure_top, temperature_top,
                               thermal_conductivity, heat_capacity_rock, density_rock,
                               compressibility_rock, porosity, total_heat_capacity,
                               salinity, largest_time_step, with_well, with_lithium };

    double dt(1.);  // Initial time increment
    long model_output(vtk_increment);
    long detailed_output(well_vtk_increment);

    ModelTimeToInteger<double> time_conversion;

    // ---------------------------------------------------
    // Read mesh, verify wells, and set the static initial temperature
    ANSYS_Model3D model(icem_mesh_file.c_str(), region_file.c_str(), phys_var_file.c_str());

    const PropertyDatabase<MODEL_DIM> &p_ref = model.Database();

    for (auto &well_name : wells_list) {
        if (!model.ContainsRegion(well_name)) {
            std::cerr << "Well with the name '" << well_name << "' does not exist." << std::endl;
            throw csmp::Exception(ERROR, "CVFEM_fault_well_lithium_example.", "Well does not exist.");
        }
    }

    model.InstantiateFiniteVolumes();

    // ── Prerequisites of the pre-split conductive geotherm ───────────────────
    // InitialiseProperties() runs much later (it needs the split regions to
    // exist), so the few fields the geotherm reads are set here. All are
    // re-applied there with their final values; only the per-region thermal
    // conductivity is specific to this solve — the uniform value applied later is
    // what the transient run uses.
    model.InputPropertyValue("temperature", ScalarVariable(ANY, temperature_top));
    model.InputPropertyValue("thermal conductivity", ScalarVariable(ANY, thermal_conductivity));
    model.InputPropertyValue("nodal heat flux bottom", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("thickness", ScalarVariable(ANY, thickness));
    if (model.ContainsRegion("FAULT2D_HOST"))
        model.Region("FAULT2D_HOST").InputPropertyValue("thickness", ScalarVariable(ANY, thickness_LD));

    if (model.ContainsRegion("BUNT")) {
        model.Region("SED").InputPropertyValue ("thermal conductivity", ScalarVariable(ANY, 1.25));
        model.Region("KEUP").InputPropertyValue("thermal conductivity", ScalarVariable(ANY, 1.25));
        model.Region("BUNT").InputPropertyValue("thermal conductivity", ScalarVariable(ANY, 9.));
        model.Region("HOST").InputPropertyValue("thermal conductivity", ScalarVariable(ANY, 9.));
        model.Region("BASE").InputPropertyValue("thermal conductivity", ScalarVariable(ANY, thermal_conductivity));
    }

    model.Boundary("TOP").InputPropertyValue("temperature", ScalarVariable(DIRICH, temperature_top));
    HeatFluxBottom(model, heat_flux_bottom);

    // ── LINEAR SOLVER SELECTION ─────────────────────────────────────────────
    // Backends in CVFEM_SolverChoice.h:
    //   SAMG    algebraic multigrid, licensed; tuned in the scheme constructor
    //   PETSc   Krylov + preconditioner, GMRES+ILU by default
    //   Eigen   direct SparseLU; rejected for the reservoir above ~10k DOF
    //
    // RESERVOIR — used by the initialisation solves below and by the transient
    // scheme, so the whole reservoir run shares one backend. SAMG or PETSc only.
    const SolverKind solver_kind = SolverKind::PETSc;

    // WELLS — the per-well Newton Jacobian, independent of the reservoir choice:
    // that system is small and dense-ish, where a direct solve wins and PETSc's
    // per-call setup cost (it rebuilds Mat/Vec/KSP each solve) would dominate.
    const SolverKind well_solver_kind = SolverKind::PETSc;

    // Well DESCRIPTION: geometry, completions, skin, initial state, physics
    // switches — one [well] section per well. Read before the scheme is built so
    // a bad file fails immediately. Keys are documented in the file itself;
    // WellConfigurationFile::WriteTemplate() writes a commented one.
    // (OPERATION — wellhead pressure/temperature, control mode, target rate — is
    // set in code: see StartWell and the time loop.)
    const std::map<std::string, WellConfiguration> well_configs =
        WellConfigurationFile::ReadAll(well_file);

    // ── Split-boundary MIDDLE nodes during initialisation ────────────────────
    // The steady initialisation systems are pure Laplacians: nothing in them
    // couples a MIDDLE node to anything else, because the operator that does
    // that (iface_transferLHS_p) exists only in the transient scheme. Those
    // nodes therefore form a disconnected block with no Dirichlet, and the
    // matrix is singular — a direct factorisation fails outright (PETSc
    // "preonly"+"lu"), while AMG quietly returns one of the valid solutions.
    //
    // Their values are not solved for anyway: CopyInsideSplitBtoMiddle overwrites
    // them with the inside/outside average after every init solve. Flagging them
    // Dirichlet therefore describes the system correctly and makes it solvable by
    // any backend, direct ones included. Released again before the transient loop.
    //
    // Kept on even though AMG copes without it: unpinned, the run depends on
    // solver leniency rather than a well-posed problem. Set false to reproduce
    // the singular behaviour and check a new backend reports it rather than
    // quietly returning garbage.
    const bool pin_split_middle_during_init = true;

    // One bundle for the steady initialisation solves. These integrators use no
    // split-boundary operators, so they stay on the default CompressedRowMatrix
    // (unlike the scheme, which is pinned to SparseMatrix).
    SolverBundle init_bundle(solver_kind);

#if defined(CSMP_WITH_PETSC_SOLVER)
    if (solver_kind == SolverKind::PETSc) {
        // The three initialisation systems are one-off STEADY solves at ~10k DOF:
        // pure Laplacians with permeability spanning ~7 orders of magnitude and no
        // storage term on the diagonal. ILU(0) fails its setup on those
        // (KSP_DIVERGED_PC_FAILED, reported by PETSc_Solver as "error code -11",
        // which is a KSPConvergedReason and not a PetscErrorCode). Verified to be
        // a conditioning problem, not a matrix-format one: the same failure occurs
        // with SparseMatrix and with CompressedRowMatrix.
        //
        // A direct factorisation is affordable at this size and exact. If it turns
        // out too slow, "gmres" + "gamg" is the iterative alternative and a far
        // stronger preconditioner than ILU for this system.
        //
        // The transient systems are unaffected: their capacitance term makes the
        // matrix diagonally dominant, so ILU works there (see the PETSc block in
        // the scheme constructor).

        init_bundle.PETScSettings().SetKSPType("gmres");
        init_bundle.PETScSettings().SetPCType("gamg"); //or hypre possible too
        init_bundle.PETScSettings().SetPrintConvergedReason(true);
    }
#endif
    cerr << "\n\n---- INITIALISATION SOLVES ----";
    cerr << "\n  linear solver: " << init_bundle.Name();
    cerr << "\n  matrix type:   CompressedRowMatrix  (default; these"
            " integrators use no split-boundary operators)";
    cerr << "\n  systems:       initial temperature, lithostatic pressure,"
            " hydrostatic pressure";

    // ════════════════════════════════════════════════════════════════════════
    // INIT TEMPERATURE: conductive geotherm — SOLVED BEFORE THE SPLIT
    //
    // Keep this ahead of the split. Once the fault is split, the two sides share
    // no element and nothing in a steady conduction problem couples INSIDE to
    // OUTSIDE — the operator that does (iface_transferLHS_t) belongs to the
    // transient scheme. Post-split the fault is therefore adiabatic, and with the
    // layered conductivity contrast (SED/KEUP 1.25 vs BUNT/HOST 9) that drives
    // real lateral heat flow, the geotherm comes out wrong. Solving on the
    // unsplit mesh is well posed, and the nodes the split duplicates inherit the
    // field. (The pressure initialisations are unaffected by this and stay after
    // the split: they are gravity-driven and essentially 1D vertical, so an
    // impermeable fault plane does not change the gradient.)
    // ════════════════════════════════════════════════════════════════════════
    cerr << "\n\n---- INIT TEMPERATURE: conductive geotherm (pre-split) ----";

    PDE_Integrator<MODEL_DIM> initial_temperature(init_bundle.Get());
    CVFEM_NumIntegral_dNT_op_dN_dV<MODEL_DIM> conductance(p_ref, "thermal conductivity", "temperature", "temperature", "thickness");
    CVFEM_PointSource_rhsop<MODEL_DIM> heat_bottom(p_ref, "nodal heat flux bottom", "temperature");
    initial_temperature.Add(&conductance);
    initial_temperature.Add(&heat_bottom);
    model.Apply(initial_temperature);
    cerr << endl << "  geotherm solved.";
    printRangeOfVariable(model, "Model", "temperature");
    printRangeOfVariable(model, "Model", "thermal conductivity");

    // Discrete-discontinuity fault. See CreateFaultSplitBoundary() below.
    CreateFaultSplitBoundary(model, params);

    model.CopyRegion("Model", "FLOW REGION");

    // ---------------------------------------------------
    // Create extra properties and input initial values
    PoreVolumeVisitor<MODEL_DIM> Pore_Volume_Visitor(model, "nodal porosity", "bulk volume", "pore volume", "thickness");
    // Properties not declared in PhysicalVariablesBenchmarks.txt, plus the
    // starting value of every field. See InitialiseProperties() below.
    InitialiseProperties(model, params);

    // Per-region values and the geometry-derived fields. See
    // SetUpRegionsAndGeometry() below.
    SetUpRegionsAndGeometry(model, params, heat_flux_bottom, Pore_Volume_Visitor);

    // ---------------------------------------------------
    // Output options

    VTK_Interface<MODEL_DIM> vtk;
    list<string> output_variables;
    list<string> output_variables_well;
    // Reservoir and well VTU field lists. See DefineOutputVariables() below.
    DefineOutputVariables(params, output_variables, output_variables_well);

    // ---------------------------------------------------
    // Visitors instantiation
    LithiumModel<MODEL_DIM> *lithium_visitor = nullptr;
    if (with_lithium) lithium_visitor = new LithiumModel<MODEL_DIM>(model, true);

    TracerModel<MODEL_DIM> *tracer_visitor = nullptr;
    if (with_lithium) tracer_visitor = new TracerModel<MODEL_DIM>(model);

    PropertyHandle<MODEL_DIM> element_id(model, "element ID", SCALAR, ELEMENT);
    PropertyHandle<MODEL_DIM> node_id(model, "node ID", SCALAR, NODE);
    PropertyHandle<MODEL_DIM> dt_handle(model, "dt", SCALAR, NODE);
    csmp::Index n_id_key(p_ref.StorageKey("node ID"));
    csmp::Index e_id_key(p_ref.StorageKey("element ID"));
    ScalarVariable ID;

    for (auto n_it = model.Mesh().NodesBegin(); n_it != model.Mesh().NodesEnd(); n_it++) {
        ID() = n_it->Idx();
        n_it->Store(n_id_key, ID);
    }
    for (auto e_it = model.Mesh().ElementsBegin(); e_it != model.Mesh().ElementsEnd(); e_it++) {
        ID() = e_it->Idx();
        e_it->Store(e_id_key, ID);
    }

    // ---------------------------------------------------
    // CVFEM PHX Scheme
    cerr << endl << "SETTING CVFEM OPTIONS";
    CVFEM_PHX_Scheme<MODEL_DIM> *CVFEM_PHX = nullptr;
    // Last argument selects the linear solver for the transient scheme; the same
    // value is used above for the initialisation solves.
    CVFEM_PHX = new CVFEM_PHX_Scheme<MODEL_DIM>(model, with_gravity, true, wells_list, with_lithium, true,
                                                                       solver_kind, well_solver_kind,
                                                                       well_configs);

    CVFEM_PHX->SetLargestTimeStep(largest_time_step);
    CVFEM_PHX->SetEquilibratorConvergenceSpeedUpTo(true);
    CVFEM_PHX->Adjust_CFL_Criterion(cfl_scaling, true);
    CVFEM_PHX->WithRockLiquidusSolidus(false);
    CVFEM_PHX->SetRockHeatCapacity(heat_capacity_rock);

    CVFEM_PHX->AddFluidContributionToHeatCapacity(true);
    CVFEM_PHX->AttemptToSurviveFluidPropertiesError(true);

    CVFEM_PHX->IncludeWellCalculations(with_well);

    CVFEM_PHX->ActivateSplitRegionCoupling(false, 1.);

    if (with_lithium) {
        CVFEM_PHX->AddAdvectionVariable("lithium content fluid", "lithium content liquid", "liquid lithium mobility", "lithium content vapor", "vapor lithium mobility");
        CVFEM_PHX->AddAdvectionVariable("tracer content fluid", "tracer content liquid", "liquid tracer mobility", "tracer content vapor", "vapor tracer mobility");
    }

    CVFEM_PHX->OpenBoundaries(temperature_top);

    cerr << "\n\n---- CVFEM SCHEME OPTIONS SET ----";

    // ---------------------------------------------------
    // Permeability
    // Background value for any region SetPermeability does not name; without it a
    // missed region keeps k = 0, which makes the steady pressure matrix singular.
    model.InputPropertyValue("permeability", ScalarVariable(ANY, 1.e-17));

    RebuildPermeabilityFields(model);

    model.CopyReplace("permeability", "previous permeability");
    model.CopyReplace("previous permeability", "previous previous permeability");
    CVFEM_PHX->PrepareTransientCalculations();

    // ---------------------------------------------------
    // Calculate initial static pressure and temperature correctly

    // Same bundle as the temperature initialisation above — same backend.
    PDE_Integrator<MODEL_DIM> lithos_pressure(init_bundle.Get());
    PDE_Integrator<MODEL_DIM> initial_pressure_(init_bundle.Get());

    NumIntegral_dNT_op_dN_dV<MODEL_DIM> conductance_p(p_ref, "permeability", "lithostatic pressure", "lithostatic pressure");
    NumIntegral_NT_op_dNi_dV<MODEL_DIM> density_p(p_ref, "nodal density rock", "permeability", "lithostatic pressure");

    NumIntegral_dNT_op_dN_dV<MODEL_DIM> conductance_hp(p_ref, "permeability", "fluid pressure", "fluid pressure");
    NumIntegral_NT_op_dNi_dV<MODEL_DIM> density_hp(p_ref, "fluid density", "permeability", "fluid pressure");

    // ── Pin / release the split-boundary MIDDLE nodes ────────────────────────
    // See the note next to pin_split_middle_during_init. Pinning makes the three
    // steady systems non-singular; releasing restores ordinary unknowns before
    // the transient scheme, which couples the middle nodes through its own
    // split-boundary operator.
    auto set_split_middle_flag = [&](VARIABLE_FLAG flag) {
        if (!pin_split_middle_during_init) return;
        const csmp::Index p_key (model.Database().StorageKey("fluid pressure"));
        const csmp::Index lp_key(model.Database().StorageKey("lithostatic pressure"));
        size_t n_mid(0U);
        for (auto sb = model.SplitBoundariesBegin(); sb != model.SplitBoundariesEnd(); ++sb)
            for (auto sb_c = sb->second.CellsBegin(); sb_c != sb->second.CellsEnd(); ++sb_c)
                for (uint32_t n{0U}; n < (*sb_c)->FE()->Nodes(); ++n) {
                    Node<MODEL_DIM>* in_node  = (*sb_c)->MatchingN(n, INSIDE);
                    Node<MODEL_DIM>* mid_node = (*sb_c)->MatchingN(n, MIDDLE);
                    // Seed from the inside node so the pinned value is sensible;
                    // CopyInsideSplitBtoMiddle replaces it with the in/out average
                    // straight after the solve either way.
                    mid_node->Store(p_key,  makeScalar(flag, in_node->Read(p_key)));
                    mid_node->Store(lp_key, makeScalar(flag, in_node->Read(lp_key)));
                    ++n_mid;
                }
        cerr << endl << "Split-boundary MIDDLE nodes " << (flag == DIRICH ? "pinned (DIRICH)" : "released (ANY)") << ": " << n_mid;
    };

    // Small helper so each initialisation stage announces itself the same way.
    auto stage = [](const char* what) {
        cerr << "\n\n---- " << what << " ----";
    };

    stage("INITIALISATION: pinning split-boundary middle nodes");
    set_split_middle_flag(DIRICH);

    stage("INIT PRESSURE: lithostatic gradient");
    // Calculate the initial lithostatic gradient assuming homogeneous rock density (carefull, it doesn't consider the "bulk" density)
    printRangeOfVariable(model, "Model", "permeability");
    lithos_pressure.Add(&conductance_p);
    lithos_pressure.Add(&density_p);
    model.Apply(lithos_pressure);
    cerr << endl << "  lithostatic pressure solved.";
    printRangeOfVariable(model, "Model", "lithostatic pressure");

    stage("INIT PRESSURE: first hydrostatic solve");
    // Hydrostatic pressure initialization is done further (in equilibration loop)
    // NOTE: "fluid density" is still at its initial value here, so this first
    // solve has no gravity source and mainly establishes the boundary values.
    // The equilibration loop below is what produces the real hydrostatic field.
    initial_pressure_.Add(&conductance_hp);
    initial_pressure_.Add(&density_hp);

    model.Apply(initial_pressure_);
    cerr << endl << "  first hydrostatic solve done.";
    printRangeOfVariable(model, "Model", "fluid pressure");

    model.InterpolateNodeToCellProperty("temperature", "temperature element");

    CopyInsideSplitBtoMiddle(model);

    // Iterate hydrostatic pressure and fluid properties to a consistent
    // initial state (converged when the max nodal pressure change < eps).
    stage("INIT P-T: equilibrating pressure and fluid properties");
    int k = 0;
    bool flag(false);
    // eps is the convergence criterion: the loop stops once the largest nodal
    // pressure change between two passes falls below it. 1000 Pa is ~0.1 m of
    // water column — small against a 23 MPa hydrostatic column, and loose enough
    // that the EOS round-trip noise does not prevent convergence.
    double diff(0.), eps(1000.0);
    while (!flag) {
        cerr << endl << "  iteration " << k << " ...";
        model.CopyReplace("fluid pressure", "test fluid pressure");
        CVFEM_PHX->InitialFluidPropertiesFromPTX();

        model.Apply(initial_pressure_);

        CopyInsideSplitBtoMiddle(model);

        diff = maxDifferenceScalarNodeProperty(model, "fluid pressure", "test fluid pressure");
        cerr << " max pressure change: " << diff << " Pa  (tolerance " << eps << ")";

        if (!(diff == diff)) {   // NaN: the loop would never exit
            cerr << endl << "  NaN in the initial pressure field — aborting equilibration.";
            break;
        }
        if (++k > 100) {
            cerr << endl << "  equilibration did not converge in 100 iterations"
                    " (last change " << diff << " Pa).";
            break;
        }
        if (diff < eps) flag = true;
    }
    cerr << endl << "  equilibration finished after " << k << " iteration(s).";
    printRangeOfVariable(model, "Model", "fluid pressure");
    printRangeOfVariable(model, "Model", "temperature");

    CVFEM_PHX->InitialFluidPropertiesFromPTX();

    cerr << endl << "FINISH initialization P-T";

    // Initialise middle region of splitB

    CopyInsideSplitBtoMiddle(model);

    // Release the MIDDLE nodes: from here on the transient scheme couples them
    // through iface_transferLHS_p, so they must be ordinary unknowns again.
    // CopyInsideSplitBtoMiddle above already stores them as ANY for non-TOP
    // nodes; this makes the release explicit and covers lithostatic pressure,
    // which that function does not touch for the DIRICH/ANY flag.
    stage("INITIALISATION: releasing split-boundary middle nodes");
    set_split_middle_flag(ANY);
    OutputAllVariableRanges(model, output_name, "after_initialisation");

    cerr << "\n\n---- INITIALISATION COMPLETE ----"
            "\n  solver: ";
    cerr << init_bundle.Name()
         << "\n  fields ready: temperature, lithostatic pressure, fluid pressure"
            "\n  handing over to the transient scheme\n";

    model.CopyReplace("fluid pressure", "previous fluid pressure");
    model.CopyReplace("previous fluid pressure", "previous previous fluid pressure");
    model.CopyReplace("previous total enthalpy", "previous previous total enthalpy");

    // ---------------------------------------------------

    // ── WELL OPERATING CONDITIONS ───────────────────────────────────────────
    // STARTING wellhead conditions only. Both wells run on wellhead-pressure
    // control while the fields settle; the time loop then switches them to
    // target-rate control (search "Set_target_rate") and these stop constraining.
    // The injector's wellhead TEMPERATURE is the exception — re-applied every
    // step, since it sets the enthalpy entering the reservoir.

    // Injector: re-injection temperature, and the wellhead pressure it starts at.
    double current_injector_WHT = 80.;    // oC
    double current_injector_WHP = 80.e5;  // Pa (80 bar)

    // Producer: starting wellhead conditions. Both are superseded once the time
    // loop switches this well to target-rate control.
    double current_producer_WHT = 80.;    // oC
    double current_producer_WHP = 80.e5;  // Pa (80 bar)

    // Build each well and put it in its starting operating state. The well
    // DESCRIPTION was read from the configuration file near the top of Run() and
    // is already inside the scheme; what is passed here is only the starting
    // wellhead conditions. See StartWell() below.
    if (with_well) {
        CVFEM_PHX->IncludeWellCalculations(true);   // global: all wells at once

        StartWell(*CVFEM_PHX, "INJECTOR", current_injector_WHP, current_injector_WHT);
        StartWell(*CVFEM_PHX, "PRODUCER", current_producer_WHP, current_producer_WHT);

        //One can set the skin that way, this can be done within the main loop to if skin changes over time
        model.Region("INJECTOR").       InputPropertyValue("well skin",                         ScalarVariable(ANY, 20.));
        model.Region("PRODUCER").       InputPropertyValue("well skin",                         ScalarVariable(ANY, 0.));
    }

    // Uniform thermal conductivity for the transient run (the per-region values
    // applied for the geotherm above are for that solve only).
    model.InputPropertyValue("thermal conductivity", ScalarVariable(ANY, thermal_conductivity));

    if (with_lithium) lithium_visitor->InitializeFromConcentration(&model, 150.e-6);

    if (model.ContainsRegion("BUNT")) {
        model.Region("SED").InputPropertyValue("lithium content liquid", ScalarVariable(ANY, 0.));
        model.Region("SED").InputPropertyValue("lithium concentration liquid", ScalarVariable(ANY, 0.));
        model.Region("SED").InputPropertyValue("lithium content fluid", ScalarVariable(ANY, 0.));
    }
    if (model.ContainsRegion("FAULT_SED")) {
        model.Region("FAULT_SED").InputPropertyValue("lithium content liquid", ScalarVariable(ANY, 0.));
        model.Region("FAULT_SED").InputPropertyValue("lithium concentration liquid", ScalarVariable(ANY, 0.));
        model.Region("FAULT_SED").InputPropertyValue("lithium content fluid", ScalarVariable(ANY, 0.));
    }

    if (with_lithium) {
        lithium_visitor->SetTimeStep(0.);
        model.Accept(*lithium_visitor);

        model.Accept(*tracer_visitor);   // (Benoit 18/08/2026) SetTimeStep removed: TracerModel has no kinetics
    }

    CVFEM_PHX->ChangeTimeStepTo(initial_time_step);

    cerr<<endl<<"OutputToVTU initial";
    OutputToVTU( model, output_name+"_initial", output_variables, 0);
    OutputToVTU_wells(model, output_name+"_initial", output_variables_well, 0);

    if (no_std_output) {
        freopen("/dev/null", "w", stdout);
    }
    // ---------------------------------------------------
    // Time loop: each step applies the coupled CVFEM scheme (returning
    // an adaptive dt), advects lithium, refreshes the permeability field,
    // advances model time, and writes output. Wells are rate-controlled
    // to a 70 kg/s target on the producer (and matching reinjection).
    cerr << "\n\n---- TIME LOOP ----";
    bool stop_sim(false);

    // ── Wall-clock accounting for the time loop ──────────────────────────────
    // Measures the whole iteration (scheme Apply + lithium/tracer + output), not
    // just the linear solve, so it is a like-for-like cost per timestep when
    // comparing solver backends. Steps are counted even when Apply() cuts dt,
    // so "average" is per loop iteration, not per unit of model time.
    using clock_type = std::chrono::steady_clock;
    const auto   run_start = clock_type::now();
    double       accumulated_step_seconds(0.);
    double       slowest_step_seconds(0.);
    size_t       step_count(0U);

    while (model_time <= max_time && !stop_sim) {
        const auto step_start = clock_type::now();
        cerr << defaultfloat;
        cerr.precision(8);

        if (with_well) {
            // Re-applied every step: this sets the enthalpy of the fluid the
            // injector delivers to the reservoir, so it is an operating condition
            // rather than an initial one.
            CVFEM_PHX->Set_well_top_temperature("INJECTOR", current_injector_WHT);

            // Target MASS rate [kg/s], positive produces and negative injects.
            // Setting it puts the well under target-rate control: a controller
            // adjusts the wellhead pressure each step until the rate is matched,
            // so the starting pressures above no longer apply. Change these values
            // — or make them functions of model_time — to run a schedule; passing
            // `false` instead returns a well to wellhead-pressure control.
            CVFEM_PHX->Set_target_rate("PRODUCER", true,  50.);
            CVFEM_PHX->Set_target_rate("INJECTOR", true, -50.);
        }

        // All-in-one T diffusion, P diffusion, advection

        CVFEM_PHX->SetModelTime(model_time);

        dt = CVFEM_PHX->Apply();

        if (with_lithium) {
            lithium_visitor->SetTimeStep(dt);
            model.Accept(*lithium_visitor);
            model.Accept(*tracer_visitor);   // (Benoit 18/08/2026) SetTimeStep removed: TracerModel has no kinetics
        }

        cerr << endl << endl << "////////////////////////";
        cerr << endl << output_name.c_str();
        cerr << endl << "Time in minutes      : " << model_time / minute << ",             in days:  " << model_time / day;
        cerr << endl << "Timestep in minutes  : " << dt / minute << ",             in days:  " << dt / day;
        cerr << endl << "Progress             : " << model_time / max_time * 100 << "%";

        // Close out the timing for this iteration (before the output writing
        // below, which is not part of the physics cost we want to compare).
        {
            const double step_seconds =
                std::chrono::duration<double>(clock_type::now() - step_start).count();
            accumulated_step_seconds += step_seconds;
            if (step_seconds > slowest_step_seconds) slowest_step_seconds = step_seconds;
            ++step_count;

            cerr << endl << "Step wall-clock      : " << step_seconds << " s" << "   (average " << accumulated_step_seconds / step_count << " s over " << step_count << " steps)";
        }

        model.CopyReplace("previous permeability", "previous previous permeability");
        model.CopyReplace("permeability", "previous permeability");

        // additional functionality

        model.InterpolateNodeToCellProperty("temperature", "temperature element");

        RebuildPermeabilityFields(model);

        // ---------------------------------------------------
        // preparing next time step
        model_time += dt;
        model.InputPropertyValue("dt", ScalarVariable(ANY, dt));

        // ---------------------------------------------------
        // output and model saving. See WriteTimestepOutput() below; it advances
        // model_output / detailed_output as it writes.
        WriteTimestepOutput(model, output_name, time_unit, time_conversion, model_time,
                            output_variables, output_variables_well, with_well,
                            model_output, vtk_increment,
                            detailed_output, well_vtk_increment);

        if (dt < 1.) {
            cerr << endl << "Timestep smaller than 1s Aborting ...";
            stop_sim = true;
        }
    }
    // ── Run summary ──────────────────────────────────────────────────────────
    {
        const double total_seconds =
            std::chrono::duration<double>(clock_type::now() - run_start).count();
        const double mean_seconds =
            step_count ? accumulated_step_seconds / step_count : 0.;

        cerr << "\n\n---- RUN SUMMARY ----";
        cerr << "\n  reservoir solver (P and T): " << SolverKindName(solver_kind);
        cerr << "\n  well solver (Newton):       " << SolverKindName(well_solver_kind)
             << (with_well ? "" : "   (wells off)");
        cerr << "\n  timesteps completed:        " << step_count;
        cerr << "\n  model time reached:         " << model_time / day << " days"
             << "   of " << max_time / day << " requested";
        cerr << "\n  mean time per step:         " << mean_seconds << " s";
        cerr << "\n  slowest step:               " << slowest_step_seconds << " s";
        cerr << "\n  total in time loop:         " << accumulated_step_seconds << " s";
        cerr << "\n  total wall clock:           " << total_seconds << " s"
             << "   (includes output writing)";
    }

    cerr << endl << "That's it..." << endl;
    // Returning normally: this run completed. (It used to throw an Exception with
    // "Done", which sent a successful run through the example suite's error path
    // and made every good log end with "ERROR: CVFEM_fault_well_lithium_example".)
}

// Assign permeability per geological region (case-defining values).
// Fault damage zones and reservoir units are set here; the thin fault
// conduit (FRACTURES MID REGIONS) gets the high fault permeability.
void CVFEM_fault_well_lithium_example::SetPermeability(Model<MODEL_DIM> &model) const {
    // Case-defining values; see Table 2 of the paper. Region names come from the
    // ICEM mesh: the geological units (SED, KEUP, BUNT, HOST, BASE) and, for each
    // of the three faults, the damage zone within each unit (FAULT*_SED etc.).
    //
    // The thin fault conduit itself is FRACTURES MID REGIONS — the lower-
    // dimensional region the split boundary creates — and carries fault2d_perm.
    //
    // Any region NOT named here keeps the background value assigned before this
    // is called. That matters: a region left at zero makes the steady pressure
    // matrix singular.
    double damage_perm(2.45e-14);      // m2, fault damage zones
    double fault2d_perm(1.e-12);       // m2, the fault conduit itself
    double overburden_perm(1.e-17);    // m2, low-permeability cover

    if (model.ContainsRegion("FAULT_SED")) model.Region("FAULT_SED").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT_KEUP")) model.Region("FAULT_KEUP").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT_BUNT")) model.Region("FAULT_BUNT").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT_HOST")) model.Region("FAULT_HOST").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT_BASE")) model.Region("FAULT_BASE").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));

    if (model.ContainsRegion("FAULT2_SED")) model.Region("FAULT2_SED").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT2_KEUP")) model.Region("FAULT2_KEUP").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT2_BUNT")) model.Region("FAULT2_BUNT").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT2_HOST")) model.Region("FAULT2_HOST").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT2_BASE")) model.Region("FAULT2_BASE").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));

    if (model.ContainsRegion("FAULT3_SED")) model.Region("FAULT3_SED").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT3_KEUP")) model.Region("FAULT3_KEUP").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("FAULT3_BUNT")) model.Region("FAULT3_BUNT").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT3_HOST")) model.Region("FAULT3_HOST").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));
    if (model.ContainsRegion("FAULT3_BASE")) model.Region("FAULT3_BASE").InputPropertyValue("permeability", ScalarVariable(ANY, damage_perm));

    if (model.ContainsRegion("SED")) model.Region("SED").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("KEUP")) model.Region("KEUP").InputPropertyValue("permeability", ScalarVariable(ANY, overburden_perm));
    if (model.ContainsRegion("BUNT")) model.Region("BUNT").InputPropertyValue("permeability", ScalarVariable(ANY, 1.e-14));
    if (model.ContainsRegion("HOST")) model.Region("HOST").InputPropertyValue("permeability", ScalarVariable(ANY, 1.e-15));
    if (model.ContainsRegion("BASE")) model.Region("BASE").InputPropertyValue("permeability", ScalarVariable(ANY, 1.e-17));
    if (model.ContainsRegion("FRACTURES MID REGIONS")) model.Region("FRACTURES MID REGIONS").InputPropertyValue("permeability", ScalarVariable(ANY, fault2d_perm));
}

// Assign NODAL porosity per geological region. Case-defining values; see Table 2
// of the paper. As with SetPermeability, a region not named here keeps the
// background value assigned before this is called.
//
// Note this sets "nodal porosity"; the element property "porosity" is derived
// from it afterwards by InterpolateNodeToCellProperty.
void CVFEM_fault_well_lithium_example::SetPorosity(Model<MODEL_DIM> &model) const {

    if (model.ContainsRegion("SED")) model.Region("SED").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.01));
    if (model.ContainsRegion("KEUP")) model.Region("KEUP").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.01));
    if (model.ContainsRegion("BUNT")) model.Region("BUNT").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.1));
    if (model.ContainsRegion("HOST")) model.Region("HOST").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("BASE")) model.Region("BASE").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));

    if (model.ContainsRegion("FAULT_SED")) model.Region("FAULT_SED").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT_KEUP")) model.Region("FAULT_KEUP").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT_BUNT")) model.Region("FAULT_BUNT").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT_HOST")) model.Region("FAULT_HOST").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT_BASE")) model.Region("FAULT_BASE").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));

    if (model.ContainsRegion("FAULT2_SED")) model.Region("FAULT2_SED").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT2_KEUP")) model.Region("FAULT2_KEUP").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT2_BUNT")) model.Region("FAULT2_BUNT").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT2_HOST")) model.Region("FAULT2_HOST").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT2_BASE")) model.Region("FAULT2_BASE").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));

    if (model.ContainsRegion("FAULT3_SED")) model.Region("FAULT3_SED").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT3_KEUP")) model.Region("FAULT3_KEUP").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));
    if (model.ContainsRegion("FAULT3_BUNT")) model.Region("FAULT3_BUNT").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT3_HOST")) model.Region("FAULT3_HOST").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.05));
    if (model.ContainsRegion("FAULT3_BASE")) model.Region("FAULT3_BASE").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.02));

    if (model.ContainsRegion("FRACTURES MID REGIONS")) model.Region("FRACTURES MID REGIONS").InputPropertyValue("nodal porosity", ScalarVariable(ANY, 0.04));
}

// Write the requested properties over the whole model to a VTU file.
void CVFEM_fault_well_lithium_example::OutputToVTU(Model<MODEL_DIM> &model, const string &model_name, const list<string> &props, size_t timestep) const {
    static VTU_Interface<MODEL_DIM> vtu(model);
    vtu.OutputDataToVTU((model_name).c_str(), props, model.Region("Model"), timestep);
}

// Write the requested properties over each well region to VTU files.
void CVFEM_fault_well_lithium_example::OutputToVTU_wells(Model<MODEL_DIM> &model, const string &model_name, const list<string> &props, size_t timestep) const {
    static VTU_Interface<MODEL_DIM> vtu(model);
    for (auto &well_name : wells_list) {
        vtu.OutputDataToVTU((model_name).c_str(), props, model.Region(well_name), timestep);
    }
}

// Return the maximum absolute nodal difference between two scalar node
// properties (used as the convergence measure in the P-T init loop).
double CVFEM_fault_well_lithium_example::maxDifferenceScalarNodeProperty(Model<MODEL_DIM> &mdl, const char *snp1, const char *snp2) {
    csmp::Index snp1Key(mdl.Database().StorageKey(snp1));
    csmp::Index snp2Key(mdl.Database().StorageKey(snp2));
    const Region<MODEL_DIM> &mref = mdl.Region("Model");

    double max(0.);

    vector<Node<MODEL_DIM> *>::const_iterator nit;
    for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit) {
        double diff = fabs((*nit)->Read(snp1Key) - (*nit)->Read(snp2Key));

        if (diff > max)
            max = diff;
    }
    return max;
}

// Build the diagonal permeability tensor per element from the separate
// horizontal and vertical permeability scalars.
void CVFEM_fault_well_lithium_example::assignKtensor(Model<MODEL_DIM> &model) {
    csmp::Index vert_perm_key(model.Database().StorageKey("vertical permeability"));
    csmp::Index horiz_perm_key(model.Database().StorageKey("horizontal permeability"));

    csmp::Index perm_tensor_key(model.Database().StorageKey("permeability tensor"));

    ScalarVariable vert_perm(ANY, 0.);
    ScalarVariable horiz_perm(ANY, 0.);
    TensorVariable<MODEL_DIM> perm_tensor(ANY, 0.);

    const Region<MODEL_DIM> &mref = model.Region("Model");

    vector<Element<MODEL_DIM> *>::const_iterator eit;

    for (eit = mref.CellsBegin(); eit != mref.CellsEnd(); ++eit) {
        vert_perm() = (*eit)->Read(vert_perm_key);
        horiz_perm() = (*eit)->Read(horiz_perm_key);

        perm_tensor(0, 0) = horiz_perm();
        perm_tensor(1, 1) = vert_perm();

        if (MODEL_DIM == 3U)
            perm_tensor(2, 2) = horiz_perm();

        (*eit)->Store(perm_tensor_key, perm_tensor);
    }
}

// Assign each node the maximum permeability of its parent (non-line)
// elements, giving a nodal permeability field for the wellbore coupling.

// Create the properties that PhysicalVariablesBenchmarks.txt does not declare (here none),
// and give every field its starting value. Runs once, after the split boundary
// exists, so the split-region and well properties can be created for the regions
// that the split produced.

void CVFEM_fault_well_lithium_example::InitialiseProperties(Model<MODEL_DIM> &model,
                                                      const RunParameters &p) {

    // NOTE (Benoit DD/MM/YYYY): this function creates no properties. They are all
    // declared in PhysicalVariablesBenchmarks.txt, which is the single
    // declaration point — creating them in code lost the unit, the min/max bounds
    // and the short name. Below is only the assignment of starting VALUES; a
    // property missing from the file fails at its InputPropertyValue call.

    //-----
    model.InputPropertyValue("thickness", ScalarVariable(ANY, p.thickness));
    if (model.ContainsRegion("FAULT2D_HOST")) model.Region("FAULT2D_HOST").InputPropertyValue("thickness", ScalarVariable(ANY, p.thickness_LD));
    if (model.ContainsRegion("FRACTURES MID REGIONS")) model.Region("FRACTURES MID REGIONS").InputPropertyValue("thickness", ScalarVariable(ANY, p.thickness_LD));


    model.InputPropertyValue("thermal conductivity", ScalarVariable(ANY, p.thermal_conductivity));
    model.InputPropertyValue("nodal density rock", ScalarVariable(ANY, p.density_rock));
    model.InputPropertyValue("fluid density", ScalarVariable(ANY, 0.));

    model.InputPropertyValue("total heat capacity", ScalarVariable(ANY, p.total_heat_capacity));
    model.InputPropertyValue("nodal heat capacity", ScalarVariable(ANY, p.total_heat_capacity));
    model.InputPropertyValue("nodal heat capacity rock", ScalarVariable(ANY, p.heat_capacity_rock));
    model.InputPropertyValue("nodal heat flux bottom", ScalarVariable(ANY, 0.));

    model.InputPropertyValue("fluid pressure", ScalarVariable(ANY, p.pressure_top));
    model.InputPropertyValue("lithostatic pressure", ScalarVariable(ANY, p.pressure_top));
    model.InputPropertyValue("nodal porosity", ScalarVariable(ANY, p.porosity));
    model.InputPropertyValue("nodal compressibility rock", ScalarVariable(ANY, p.compressibility_rock));
    model.InputPropertyValue("courant liquid", ScalarVariable(ANY, p.largest_time_step));
    model.InputPropertyValue("courant vapor", ScalarVariable(ANY, p.largest_time_step));
    model.InputPropertyValue("salinity", ScalarVariable(ANY, p.salinity));
    model.InputPropertyValue("fluid source rate", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("fluid source h", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("fluid source wt", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("nodal source liquid", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("nodal source vapor", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("upwind control liquid", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("upwind control vapor", ScalarVariable(ANY, 0.));

    // Velocity
    model.InputPropertyValue("velocity", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("pore velocity", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("velocity liquid", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("pore velocity liquid", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("nodal pore velocity", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("velocity vapor", VectorVariable<MODEL_DIM>(ANY, 0.));
    model.InputPropertyValue("pore velocity vapor", VectorVariable<MODEL_DIM>(ANY, 0.));

    // Boundaries
    model.Boundary("TOP").InputPropertyValue("lithostatic pressure", ScalarVariable(DIRICH, p.pressure_top));  // Wrong if TOP is not the surface!!!
    model.Boundary("TOP").InputPropertyValue("temperature", ScalarVariable(DIRICH, p.temperature_top));
    model.Boundary("TOP").InputPropertyValue("fluid pressure", ScalarVariable(DIRICH, p.pressure_top));

    if (p.with_well) {

        //-----
        model.InputPropertyValue("wellhead pressure", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well fluid pressure", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well enthalpy liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well enthalpy vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well enthalpy fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well salt fraction liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well salt fraction vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well salt fraction fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well velocity liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well velocity vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well velocity fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well density liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well density vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well density fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well viscosity fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well Re number", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well saturation liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well saturation vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well saturation halite", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well temperature", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well cross area", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well friction factor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well fluid state", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well completion", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well completion length", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well effective radius", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well injectivity", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well mass transfer rate", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well tracer mass transfer rate", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well energy transfer rate liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well energy transfer rate vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well energy transfer rate", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well salt mass transfer rate", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well index", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well skin", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("well reservoir radial heat", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("total mass exchanged", ScalarVariable(ANY, 0.));
    }

    if (p.with_lithium) {
        model.InputPropertyValue("lithium content fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("lithium content liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("lithium content vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("liquid lithium mobility", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("vapor lithium mobility", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.25));
        model.InputPropertyValue("lithium bulk volumic mass", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("lithium solubility liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("lithium concentration liquid", ScalarVariable(ANY, 0.));

        model.InputPropertyValue("tracer content fluid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("tracer content liquid", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("tracer content vapor", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("liquid tracer mobility", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("vapor tracer mobility", ScalarVariable(ANY, 0.));
    }

    // Split-region coupling fields. Declared in PhysicalVariablesBenchmarks.txt;
    // only their starting values are set here.
    {
        model.InputPropertyValue("split region diffused heat source", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region nodal area", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region mass source", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region energy source", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region salt source", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region hydraulic conductivity", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region gravity mass source", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region potential gradient l", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region potential gradient v", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region relperm l", ScalarVariable(ANY, 0.));
        model.InputPropertyValue("split region relperm v", ScalarVariable(ANY, 0.));
    }

    // Magmatic variables
    model.InputPropertyValue("magmatic fluid mass", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic ratio", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic water ratio", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic fluid mass liquid", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic liquid mass mobility", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic fluid mass vapor", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic vapor mass mobility", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic mass salt", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic salt ratio", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic salt content liquid", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic liquid salt mobility", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic salt content vapor", ScalarVariable(ANY, 0.));
    model.InputPropertyValue("magmatic vapor salt mobility", ScalarVariable(ANY, 0.));

    cerr << "\n\n---- PROPERTIES CREATED AND INITIALISED ----";

}

// Region-by-region setup and the geometry-derived fields: porosity, the lithium
// starting distribution, nodal/barycentre depth and coordinates, the
// depth-dependent salinity profile, the basal heat flux, nodal permeability and
// the pore/bulk volumes.
//
// Runs after InitialiseProperties(), because it overwrites some of the uniform
// starting values with per-region ones. Order matters inside: the coordinate and
// depth assignments must precede SpecialSalinity (which reads nodal depth) and
// the PoreVolumeVisitor (which reads nodal porosity).

void CVFEM_fault_well_lithium_example::SetUpRegionsAndGeometry(Model<MODEL_DIM> &model,
                                                         const RunParameters &p,
                                                         double heat_flux_bottom,
                                                         PoreVolumeVisitor<MODEL_DIM> &Pore_Volume_Visitor) {
    // Special cases init variables

    // Variables setup

    SetPorosity(model);

    if (p.with_lithium) {
        if (model.ContainsRegion("FAULT_SED")) model.Region("FAULT_SED").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("FAULT2_SED")) model.Region("FAULT2_SED").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("FAULT3_SED")) model.Region("FAULT3_SED").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("SED")) model.Region("SED").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("FAULT2D_SED")) model.Region("FAULT2D_SED").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("FAULT2D_SED_REGION0_FAULT_SED_INTERSECTION"))
            model.Region("FAULT2D_SED_REGION0_FAULT_SED_INTERSECTION").InputPropertyValue("lithium solid", ScalarVariable(ANY, 0.));
    }

    cerr << "\n\n---- REGION PROPERTIES AND GEOMETRY SET ----";

    // Final functions

    model.Region("Model").AssignNodeCoordinatesTo("nodal depth", 'y');
    model.Region("Model").AssignNodeCoordinatesTo("coordinate x", 'x');
    model.Region("Model").AssignNodeCoordinatesTo("coordinate y", 'y');
    model.Region("Model").AssignNodeCoordinatesTo("coordinate z", 'z');
    model.InterpolateNodeToCellProperty("nodal depth", "barycenter depth");
    SpecialSalinity(model);
    HeatFluxBottom(model, heat_flux_bottom);
    CalculateNodalPermeability(model);
    model.Accept(Pore_Volume_Visitor);
    model.InterpolateNodeToCellProperty("nodal porosity", "porosity");
}

// Build the two VTU output lists: `output_variables` for the reservoir and
// `output_variables_well` for the well segments. The well list starts as a copy
// of the reservoir list, so anything added before that point appears in both.

void CVFEM_fault_well_lithium_example::DefineOutputVariables(const RunParameters &p,
                                                       std::list<std::string> &output_variables,
                                                       std::list<std::string> &output_variables_well) const {

    if (p.with_lithium) {
        output_variables.push_back("lithium content fluid");
        output_variables.push_back("lithium content liquid");
        output_variables.push_back("lithium content vapor");
        output_variables.push_back("liquid lithium mobility");
        output_variables.push_back("vapor lithium mobility");
        output_variables.push_back("lithium solid");
        output_variables.push_back("lithium bulk volumic mass");
        output_variables.push_back("lithium solubility liquid");
        output_variables.push_back("lithium concentration liquid");

        output_variables.push_back("tracer content fluid");
        output_variables.push_back("tracer content liquid");
        output_variables.push_back("tracer content vapor");
        output_variables.push_back("liquid tracer mobility");
    }

    output_variables.push_back("boundary flow mass");

    output_variables.push_back("nodal depth");
    output_variables.push_back("coordinate x");
    output_variables.push_back("coordinate y");
    if (MODEL_DIM == 3U)
        output_variables.push_back("coordinate z");
    output_variables.push_back("barycenter depth");
    output_variables.push_back("nodal density rock");
    output_variables.push_back("temperature");
    output_variables.push_back("thermal conductivity");
    output_variables.push_back("fluid pressure");
    output_variables.push_back("lithostatic pressure");
    output_variables.push_back("density vapor");
    output_variables.push_back("density liquid");
    output_variables.push_back("volume factor");
    output_variables.push_back("saturation liquid");
    output_variables.push_back("saturation vapor");
    output_variables.push_back("fluid state");
    output_variables.push_back("previous total enthalpy");
    output_variables.push_back("enthalpy content liquid");
    output_variables.push_back("enthalpy content vapor");
    output_variables.push_back("enthalpy liquid");
    output_variables.push_back("enthalpy vapor");
    output_variables.push_back("fluid mass liquid");
    output_variables.push_back("fluid mass vapor");
    output_variables.push_back("fluid density");
    output_variables.push_back("bulk fluid density");
    output_variables.push_back("viscosity liquid");
    output_variables.push_back("viscosity vapor");
    output_variables.push_back("pore velocity liquid");
    output_variables.push_back("pore velocity vapor");
    output_variables.push_back("velocity");
    output_variables.push_back("pore velocity");
    output_variables.push_back("nodal total compressibility");
    output_variables.push_back("compressibility");
    output_variables.push_back("permeability");
    output_variables.push_back("nodal permeability");
    output_variables.push_back("nodal compressibility rock");

    output_variables.push_back("KgradP");
    output_variables.push_back("saturation halite");
    output_variables.push_back("salt content vapor");
    output_variables.push_back("salt content liquid");
    output_variables.push_back("salt fraction vapor");
    output_variables.push_back("salt fraction liquid");
    output_variables.push_back("salt fraction fluid");
    output_variables.push_back("salinity");
    output_variables.push_back("pore volume");
    output_variables.push_back("bulk volume");

    output_variables.push_back("nodal fluid volume source");
    output_variables.push_back("porosity");
    output_variables.push_back("nodal porosity");
    output_variables.push_back("nodal heat capacity rock");
    output_variables.push_back("fluid heat capacity");
    output_variables.push_back("nodal heat capacity");
    output_variables.push_back("nodal heat flux bottom");

    output_variables_well = output_variables;

    if (p.with_well) {
        output_variables_well.push_back("wellhead pressure");
        output_variables_well.push_back("well fluid pressure");
        output_variables_well.push_back("well density liquid");
        output_variables_well.push_back("well density vapor");
        output_variables_well.push_back("well density fluid");
        output_variables_well.push_back("well viscosity fluid");
        output_variables_well.push_back("well Re number");
        output_variables_well.push_back("well saturation liquid");
        output_variables_well.push_back("well saturation vapor");
        output_variables_well.push_back("well saturation halite");
        output_variables_well.push_back("well enthalpy liquid");
        output_variables_well.push_back("well enthalpy vapor");
        output_variables_well.push_back("well enthalpy fluid");
        output_variables_well.push_back("well salt fraction liquid");
        output_variables_well.push_back("well salt fraction vapor");
        output_variables_well.push_back("well salt fraction fluid");
        output_variables_well.push_back("well index");
        output_variables_well.push_back("well reservoir radial heat");
        output_variables_well.push_back("well completion");
        output_variables_well.push_back("well completion length");
        output_variables_well.push_back("well effective radius");
        output_variables_well.push_back("well skin");
        output_variables_well.push_back("well energy transfer rate liquid");
        output_variables_well.push_back("well energy transfer rate vapor");
        output_variables_well.push_back("well injectivity");
        output_variables_well.push_back("well mass transfer rate");
        output_variables_well.push_back("well tracer mass transfer rate");
        output_variables_well.push_back("well salt mass transfer rate");
        output_variables_well.push_back("well velocity liquid");
        output_variables_well.push_back("well velocity vapor");
        output_variables_well.push_back("well velocity fluid");
        output_variables_well.push_back("well cross area");
        output_variables_well.push_back("well friction factor");
        output_variables_well.push_back("well temperature");
        output_variables_well.push_back("total mass exchanged");
    }
}

// Turn the lower-dimensional fault region into a discrete discontinuity: create
// the split boundary from FAULT2D_HOST, insert the middle (fault) region between
// the two sides, and pull them apart by the fault aperture.
//
// The well regions need repairing afterwards: the split duplicates nodes that the
// well line elements referenced, so each well element is re-pointed at the
// surviving node and its subdomain rebuilt. Without that the wells would still be
// attached to nodes that no longer belong to the flow region.
//
// None of the locals here outlive the function — the regions it creates are found
// again later by name ("FRACTURES MID REGIONS").
void CVFEM_fault_well_lithium_example::CreateFaultSplitBoundary(Model<MODEL_DIM> &model,
                                                          const RunParameters &p) {
    // ---------------------------------------------------
    // Discrete-discontinuity fault: create the split boundary from the
    // lower-dimensional fault region, insert the middle (fault) region,
    // and pull the two sides apart by the fault thickness.
    // SplitBoundary<MODEL_DIM> *sb_host = nullptr;
    set<string> new_regions;
    set<string> fractures;
    set<string> fractures_mid_regions;
    vector<string> fractureNames = {"FAULT2D_HOST"};

    // Loop through each fracture name and call CreateSplitBoundaryFrom
    for (const auto &fractureName : fractureNames) {
        //sb_name = *((model.CreateSplitBoundaryFrom(fractureName.c_str(),true)).first.begin());
    }

    pair<set<string>,bool> sb_name  = model.CreateSplitBoundaryFrom( "FAULT2D_HOST", false );
    int32_t material_id{ 1U };
    new_regions = model.InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id );
    SplitBoundary<MODEL_DIM>& s_ref = model.SplitBoundary( (*sb_name.first.begin()) );
    s_ref.PullApartSplitBoundary(p.thickness_LD);

    cerr << endl << "CREATED THE MID REGIONS: ";
    for (auto &item : new_regions) {
        cerr << endl << item;
        fractures_mid_regions.insert(item.c_str());
    }
    model.MergeRegions(fractures_mid_regions, "FRACTURES MID REGIONS");

    if (p.with_well) {
        // Well nodes that intersected the fault were moved to the split-boundary
        // "inside" side; move them back into the middle (fault) region so the
        // wells remain hydraulically connected to the fault conduit.
        Node<MODEL_DIM> *mid_region_node = nullptr;
        Node<MODEL_DIM> *old_well_node = nullptr;

        for (auto sb = model.SplitBoundariesBegin(); sb != model.SplitBoundariesEnd(); ++sb) {
            for (auto sb_c = sb->second.CellsBegin(); sb_c != sb->second.CellsEnd(); ++sb_c) {
                uint32_t sb_nodes = (*sb_c)->FE()->Nodes();

                for (uint32_t n{0U}; n < sb_nodes; ++n) {
                    Node<MODEL_DIM> *in_node = (*sb_c)->MatchingN(n, INSIDE);
                    Node<MODEL_DIM> *mid_node = (*sb_c)->MatchingN(n, MIDDLE);

                    for (auto well_node : model.Region("INJECTOR").NodeVector()) {
                        if (well_node == in_node) {
                            old_well_node = well_node;
                            mid_region_node = mid_node;
                            break;
                        }
                    }
                }
            }
        }

        for (auto &el : model.Region("INJECTOR").CellVector()) {
            for (uint32_t n = 0; n < el->Nodes(); n++) {
                if (el->N(n) == old_well_node)
                    el->Assign(n, mid_region_node);
            }
        }

        model.Region("INJECTOR").RebuildSubDomainAfterChangeOfCellVector();

        mid_region_node = nullptr;
        old_well_node = nullptr;

        for (auto sb = model.SplitBoundariesBegin(); sb != model.SplitBoundariesEnd(); ++sb) {
            for (auto sb_c = sb->second.CellsBegin(); sb_c != sb->second.CellsEnd(); ++sb_c) {
                uint32_t sb_nodes = (*sb_c)->FE()->Nodes();

                for (uint32_t n{0U}; n < sb_nodes; ++n) {
                    Node<MODEL_DIM> *in_node = (*sb_c)->MatchingN(n, INSIDE);
                    Node<MODEL_DIM> *mid_node = (*sb_c)->MatchingN(n, MIDDLE);

                    for (auto well_node : model.Region("PRODUCER").NodeVector()) {
                        if (well_node == in_node) {
                            old_well_node = well_node;
                            mid_region_node = mid_node;
                            break;
                        }
                    }
                }
            }
        }

        for (auto &el : model.Region("PRODUCER").CellVector()) {
            for (uint32_t n = 0; n < el->Nodes(); n++) {
                if (el->N(n) == old_well_node)
                    el->Assign(n, mid_region_node);
            }
        }

        model.Region("PRODUCER").RebuildSubDomainAfterChangeOfCellVector();
    }

    model.InstantiateFiniteVolumes();

    cerr << "\n\n---- FAULT SPLIT BOUNDARY CREATED ----";
}

// Builds ONE well and sets its starting operating conditions. The description
// (geometry, completions, skin, physics switches) came from the configuration
// file via the scheme constructor; this sets only what changes during the run.
//
// Initialize_well() turns that description into a discretised well: resolves the
// radius profile, marks completion nodes, computes completion lengths, seeds the
// skin. It reads reservoir fields, hence called from Run() and not at construction.
//
// Starts under wellhead-pressure control; switch to target-rate or top-injection
// afterwards (see WellControl).

void CVFEM_fault_well_lithium_example::StartWell(CVFEM_PHX_Scheme<MODEL_DIM> &scheme,
                                           const std::string &well_name,
                                           double wellhead_pressure,
                                           double wellhead_temperature) const {

    scheme.Set_well_top_pressure(well_name, wellhead_pressure);
    scheme.Set_well_top_temperature(well_name, wellhead_temperature);
    scheme.Set_well_water_table(well_name, false, 0.);   // displacement off
    scheme.Initialize_well(well_name);

    cerr << endl << "Well " << well_name << " started under wellhead-pressure control: "
         << wellhead_pressure << " Pa, " << wellhead_temperature << " oC"
                                                                    " (geometry and completions from the configuration file)";
}

void CVFEM_fault_well_lithium_example::WriteTimestepOutput(
    Model<MODEL_DIM> &model, const std::string &output_name,
    TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion, double model_time,
    const std::list<std::string> &output_variables,
    const std::list<std::string> &output_variables_well,
    bool with_well,
    long &model_output, long vtk_increment,
    long &detailed_output, long detailed_vtk_increment) const {

    long time_label(0L);
    switch (time_unit) {
    case (TimeUnit::YEARS):
        time_label = time_conversion.ModelTimeInYears(model_time);
        break;
    case (TimeUnit::DAYS):
        time_label = time_conversion.ModelTimeInDays(model_time);
        break;
    case (TimeUnit::HOURS):
        time_label = time_conversion.ModelTimeInHours(model_time);
        break;
    case (TimeUnit::MINUTES):
        time_label = time_conversion.ModelTimeInMinutes(model_time);
        break;
    case (TimeUnit::SECONDS):
        time_label = (long)floor(model_time + 0.5);
        break;
    }

    if (time_label >= model_output) {
        OutputToVTU(model, output_name, output_variables, time_label);
        model_output += vtk_increment;
        cerr << endl << "Writing file, time_label: " << time_label;
    }

    if (time_label >= detailed_output && with_well) {
        OutputToVTU_wells(model, output_name, output_variables_well, time_label);
        detailed_output += detailed_vtk_increment;
        cerr << endl << "Writing detailed output file, time_label: " << time_label;
    }
}

// Rebuild every permeability-derived field from the per-region values.
// Setup and the time loop both call this, so a change to the permeability model
// cannot apply to one and not the other. The wells are overridden AFTER
// SetPermeability because they cut through several regions.
void CVFEM_fault_well_lithium_example::RebuildPermeabilityFields(Model<MODEL_DIM> &model) {
    SetPermeability(model);

    if (!wells_list.empty()) {
        for (const auto &well_name : wells_list)
            model.Region(well_name).InputPropertyValue("permeability", ScalarVariable(ANY, 1.e-17));
    }

    model.CopyReplace("permeability", "vertical permeability");
    model.CopyReplace("permeability", "horizontal permeability");
    assignKtensor(model);
    CalculateNodalPermeability(model);
}

// Dump min/max of every property in the database to a text file. Modelled on
// GeothermalSimulator::OutputAllVariableRangesAfterInitialization: a single
// diff of two of these files tells you exactly which fields two runs disagree on.
void CVFEM_fault_well_lithium_example::OutputAllVariableRanges(Model<MODEL_DIM> &model,
                                                         const std::string &output_name,
                                                         const std::string &tag) const {
    const std::string filename = output_name + "_VariableRanges_" + tag + ".txt";
    std::ofstream out(filename);
    if (!out) {
        cerr << endl << "  could not open " << filename << " for writing";
        return;
    }

    out << "Model:\t" << output_name << "\n";
    out << "Stage:\t" << tag << "\n\n";
    out << "variable\tmin\tmax\n";

    double min_value(0.), max_value(0.);
    size_t n(0U);
    for (auto it = model.Database().Begin(); it != model.Database().End(); ++it) {
        model.MinMaxOf((*it).first.c_str(), min_value, max_value);
        out << (*it).first.c_str() << "\t" << min_value << "\t" << max_value << "\n";
        ++n;
    }
    cerr << endl << "  wrote ranges of " << n << " variables to " << filename;
}

void CVFEM_fault_well_lithium_example::CalculateNodalPermeability(Model<MODEL_DIM> &model) {
    const Region<MODEL_DIM> &mref = model.Region("Model");
    csmp::Index nK_key(model.Database().StorageKey("nodal permeability"));
    csmp::Index K_key(model.Database().StorageKey("permeability"));

    ScalarVariable nK(ANY, 0.);

    int nParents(0);
    double eleK(0.);

    vector<Node<MODEL_DIM> *>::const_iterator nit;
    for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit) {
        double mini_eleK(1000.);
        double max_eleK(0.);
        // number of parent elements of the current well node
        nParents = (*nit)->Parents();

        // for all parent elements
        for (int i = 0; i < nParents; i++) {
            if (!(*nit)->Parent(i)->IsLine()) {
                eleK = (*nit)->Parent(i)->Read(K_key);

                // if( eleK < mini_eleK ) mini_eleK = eleK;
                if (eleK > max_eleK) max_eleK = eleK;
            }
        }

        // nK()=mini_eleK;
        nK() = max_eleK;
        (*nit)->Store(nK_key, nK);
    }
}

// Prescribe the initial salinity profile: linear from 0 at the surface
// to the deep value below 1700 m, then constant with depth.
// Depth-dependent initial salinity: linear from 0 wt% at the surface to 10 wt%
// at 1700 m depth, constant at 10 wt% below that. This stands in for the
// observed brine stratification of the Upper Rhine Graben; it is applied ONCE at
// setup and is then transported like any other solute, so it is an initial
// condition rather than a constraint. The 1700 m break and the 10 wt% plateau are
// case-defining values.
void CVFEM_fault_well_lithium_example::SpecialSalinity(Model<MODEL_DIM> &model) const {
    const Region<MODEL_DIM> &mref = model.Region("Model");

    csmp::Index depth_key(model.Database().StorageKey("nodal depth"));
    csmp::Index S_key(model.Database().StorageKey("salinity"));

    ScalarVariable depth(ANY, 0.);
    ScalarVariable S(ANY, 0.);

    vector<Node<MODEL_DIM> *>::const_iterator nit;
    for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit) {
        depth = (*nit)->Read(depth_key);

        if (depth() > -1700.) {
            S() = -10. / 1700 * depth();
            ;
        } else {
            S() = 10.;
        }
        //cerr << endl << "Depth: " << depth() << ", salinity: " << S();
        (*nit)->Store(S_key, S);
    }
}

// Synchronise pressure and temperature across the split boundary by
// setting the middle (fault) node to the average of the inside/outside
// matrix nodes, preserving the Dirichlet flag on the top boundary.
void CVFEM_fault_well_lithium_example::CopyInsideSplitBtoMiddle(Model<MODEL_DIM> &model) {
    for (auto sb = model.SplitBoundariesBegin(); sb != model.SplitBoundariesEnd(); ++sb) {
        cerr << endl << "Treating SB: " << sb->first;
        for (auto sb_c = sb->second.CellsBegin(); sb_c != sb->second.CellsEnd(); ++sb_c) {
            uint32_t n_nodes = (*sb_c)->FE()->Nodes();
            for (uint32_t n{0U}; n < n_nodes; ++n) {
                Node<MODEL_DIM> *in_node = (*sb_c)->MatchingN(n, INSIDE);
                Node<MODEL_DIM> *out_node = (*sb_c)->MatchingN(n, OUTSIDE);
                Node<MODEL_DIM> *mid_node = (*sb_c)->MatchingN(n, MIDDLE);

                csmp::Index p_key(model.Database().StorageKey("fluid pressure"));
                csmp::Index t_key(model.Database().StorageKey("temperature"));
                csmp::Index lp_key(model.Database().StorageKey("lithostatic pressure"));
                ScalarVariable p, t;

                ScalarVariable p_in, p_mid, p_out;
                p_in = in_node->Read(p_key);
                p_out = out_node->Read(p_key);
                p_mid = mid_node->Read(p_key);

                if (in_node->AtBoundary() == TOP) {
                    cerr << endl << " IN Node at Boundary: " << parseBoundary(in_node->AtBoundary());
                    if (p_in.Flag() == DIRICH) cerr << endl << "Pressure    is DIRICH";
                    else cerr << endl << "Pressure    is not DIRICH!!!!";
                    cerr << endl << " OUT Node at Boundary: " << parseBoundary(out_node->AtBoundary());
                    if (p_out.Flag() == DIRICH) cerr << endl << "Pressure    is DIRICH";
                    else cerr << endl << "Pressure    is not DIRICH!!!!";
                    cerr << endl << " MID Node at Boundary: " << parseBoundary(mid_node->AtBoundary());
                    if (p_mid.Flag() == DIRICH) cerr << endl << "Pressure    is DIRICH";
                    else cerr << endl << "Pressure    is not DIRICH!!!!";
                }

                // new version
                double s_in, s_out, s_mid;
                if (in_node->AtBoundary() == TOP) {
                    s_in = in_node->Read(p_key);
                    s_out = out_node->Read(p_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(p_key, makeScalar(DIRICH, s_in));
                    out_node->Store(p_key, makeScalar(DIRICH, s_out));
                    mid_node->Store(p_key, makeScalar(DIRICH, s_mid));

                    s_in = in_node->Read(t_key);
                    s_out = out_node->Read(t_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(t_key, makeScalar(DIRICH, s_in));
                    out_node->Store(t_key, makeScalar(DIRICH, s_out));
                    mid_node->Store(t_key, makeScalar(DIRICH, s_mid));

                    // Lithostatic pressure: output-only field, but the middle
                    // nodes are otherwise never assigned and show up as spikes
                    // in the VTUs. (Benoit DD/MM/YYYY)
                    s_in = in_node->Read(lp_key);
                    s_out = out_node->Read(lp_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(lp_key, makeScalar(DIRICH, s_in));
                    out_node->Store(lp_key, makeScalar(DIRICH, s_out));
                    mid_node->Store(lp_key, makeScalar(DIRICH, s_mid));

                }

                else {
                    s_in = in_node->Read(p_key);
                    s_out = out_node->Read(p_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(p_key, makeScalar(ANY, s_in));
                    out_node->Store(p_key, makeScalar(ANY, s_out));
                    mid_node->Store(p_key, makeScalar(ANY, s_mid));

                    s_in = in_node->Read(t_key);
                    s_out = out_node->Read(t_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(t_key, makeScalar(ANY, s_in));
                    out_node->Store(t_key, makeScalar(ANY, s_out));
                    mid_node->Store(t_key, makeScalar(ANY, s_mid));

                    s_in = in_node->Read(lp_key);
                    s_out = out_node->Read(lp_key);
                    s_mid = (s_in + s_out) * 0.5;

                    in_node->Store(lp_key, makeScalar(ANY, s_in));
                    out_node->Store(lp_key, makeScalar(ANY, s_out));
                    mid_node->Store(lp_key, makeScalar(ANY, s_mid));
                }
            }
        }
    }
}

// Convert the basal heat-flux boundary condition (W/m^2) into nodal
// point sources on the BOTTOM boundary, weighted by element area.
void CVFEM_fault_well_lithium_example::HeatFluxBottom(Model<MODEL_DIM> &model, double heat_flux_bottom) {
    csmp::Index hfb_key(model.Database().StorageKey("nodal heat flux bottom"));
    csmp::Index thickness_key(model.Database().StorageKey("thickness"));
    ScalarVariable area, thickness, hfb;

    Boundary<MODEL_DIM> &mref = model.Boundary("BOTTOM");
    for (auto eit = mref.CellsBegin(); eit != mref.CellsEnd(); eit++) {
        for (auto i{0U}; i < (*eit)->Nodes(); i++) {
            (*eit)->N(i)->Read(hfb_key, hfb);

            double area = (*eit)->Area();
            area /= (*eit)->Nodes();

            hfb() += heat_flux_bottom * area;

            (*eit)->N(i)->Store(hfb_key, hfb);
        }
    }
}
// namespace csmp

void CVFEM_fault_well_lithium_example::PAUSE() {
    string dummy;
    cerr << " press ENTER to continue... " << endl;
    getline(cin, dummy);
}
}