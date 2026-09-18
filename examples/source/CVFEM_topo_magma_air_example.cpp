// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_topo_magma_air_example.h"

// VTU output
#include "VTU_Interface.h"

// Additional Visitors
#include "MagmaModel.h"
#include "ParmigianiVisitor.h"
#include "RatioVisitor_split2.h"
#include "RatioVisitor_split1.h"
#include "CVFEM_PorosityChangeVisitor.h"
#include "ConservationCheck.h"

// CVFEM transport scheme
#include "CVFEM_PHX_Scheme.h"
#include "PermeabilityVisitor.h"
#include "FailureModeVisitor.h"
#include "PoreVolumeVisitor.h"
#include "ModelTime.h"
#include "ModelTimeToInteger.h"
#include "PropertyHandle.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"

#include "time.h"
#include "Boundary.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "CVFEM_SolverChoice.h"

#include <iostream>
#include <chrono>
#include <fstream>


using namespace std;
namespace csmp {
void CVFEM_topo_magma_air_example::Specifications()
{
    SetTitle("CVFEM 2D with topography, optionally with magma chamber and air phase");
    SetDifficulty(0);
    SetCategory("CVFEM Examples");
    AddAuthor("BLC");
    AddDescription("source file in: CVFEM_topo_magma_air_example");
}

void CVFEM_topo_magma_air_example::Run()
{

    // ════════════════════════════════════════════════════════════════════════
    // 1. WHAT THIS RUN INCLUDES
    //
    // The feature switches. Everything below is parameters for whatever is
    // switched on here.
    // ════════════════════════════════════════════════════════════════════════
    bool        with_air(true);                     // air phase above the water table
    bool        with_gravity(true);                 // gravity in flow and transport
    bool        open_top(true);                     // flow in and out of the TOP boundary
    bool        no_std_output(false);                // suppress the per-step reservoir output

    // Air/vapour treatment. Both feed CVFEM_PHX->ActivateAir() and both enable the
    // gas-phase mixture diagnostics in the VTU output.
    bool        treat_air_and_vapor_as_mixture(false);            // keep false: does not work yet
    bool        treat_air_and_vapor_as_hydraulic_mixture(true);

    // Magma chamber. with_magma_chamber places the intrusion; with_magma_model
    // additionally solves fluid production and crystallisation within it.
    bool        with_magma_chamber(true);

    bool        with_lithostatic_ductile(true);    // fluid pressure = lithostatic in ductile regions

    // Pre-step temperature diffusion, to mimic an intrusion not placed instantaneously.
    bool        with_temperature_halo(true);
    double      initial_T_diffusion_time(3.14e7 * 5000);    // [s]

    // ════════════════════════════════════════════════════════════════════════
    // 2. INPUT AND OUTPUT FILES
    //
    // All expected in the working directory. The mesh name also names the regions
    // file, as <mesh>-regions.txt.
    // ════════════════════════════════════════════════════════════════════════
    string      icem_mesh_file("CVFEM_topo_magma_air_example_mesh");    // ANSYS-ICEM mesh, without .asc/.dat
    string      region_file("CVFEM_topo_magma_air_example");       // region list

    string      phys_var_file("PhysicalVariablesBenchmarks.txt");

    string      output_name("CVFEM_topo_magma_air_example_output");       // prefix for every output file
    string      restart_file("saved-CVFEM_topo_magma_air_example");     // read when restarting
    string      save_file("saved-CVFEM_topo_magma_air_example");        // written for a later restart

    // ════════════════════════════════════════════════════════════════════════
    // 3. LINEAR SOLVER
    //
    // One choice for every steady solve here — the conductive geotherm, the
    // transient temperature diffusion and the two initial pressure systems.
    // Backends are in CVFEM_SolverChoice.h; all derive from csmp::Solver.
    //
    //   SolverKind::SAMG   algebraic multigrid. Needs a licence and the build
    //                      flag CSMP_WITH_SAMG_SOLVER.
    //   SolverKind::PETSc  Krylov + preconditioner; needs CSMP_WITH_PETSC_SOLVER.
    //                      Its GMRES+ILU default is weak for the steady Laplacians
    //                      here, so ApplySolverSettings() switches it to AMG.
    // ════════════════════════════════════════════════════════════════════════
    const SolverKind solver_kind = SolverKind::PETSc;

    // Applied to every bundle, so the solves cannot drift apart.
    auto ApplySolverSettings = [&](SolverBundle &bundle)
    {
#ifdef CSMP_WITH_SAMG_SOLVER
        if (bundle.Kind() == SolverKind::SAMG)
        {
            SAMG_Settings &s = bundle.SAMGSettings();
            s.Set_nred(1);
            s.Set_a_cmplx(3.);
            s.Set_g_cmplx(2.);
            s.Set_w_avrge(2.);
            s.Set_idmp(-1);
            s.Set_iout1(-1);
            s.Set_iout2(-1);
        }
#endif
#if defined(CSMP_WITH_PETSC_SOLVER)
        if (bundle.Kind() == SolverKind::PETSc)
        {
            // Steady Laplacians with a wide permeability contrast and no storage
            // term on the diagonal; ILU fails its setup on those, so use AMG.
            // ("gamg" is always available; "hypre" is usually faster where built.)
            bundle.PETScSettings().SetKSPType("gmres");
            bundle.PETScSettings().SetPCType("gamg");
            bundle.PETScSettings().SetRelativeTolerance(1.e-12);
            bundle.PETScSettings().SetPrintConvergedReason(true);
        }
#endif
    };

    // ════════════════════════════════════════════════════════════════════════
    // 4. TIME STEPPING
    //
    // Every duration here is in `time_unit` and is converted to seconds further
    // down by time_multiplier. Change the unit and they all follow.
    // ════════════════════════════════════════════════════════════════════════
    TimeUnit    time_unit(TimeUnit::YEARS);
    double      max_time(5000);              // total simulated time
    long        vtk_increment(100);           // VTU written every ...
    double      largest_time_step(0.1);     // upper bound on dt
    double      initial_time_step(0.001);   // dt for the first step
    double      cfl_scaling(0.5);           // safety factor on the CFL-limited dt

    // ════════════════════════════════════════════════════════════════════════
    // 5. BOUNDARY CONDITIONS
    // ════════════════════════════════════════════════════════════════════════
    double      temperature_top(15.);       // oC, TOP boundary
    double      pressure_top(1.01325e5);    // Pa, TOP boundary (1 atm)

    int         heat_flux_bottom_(100);                             // mW/m2, as quoted
    double      heat_flux_bottom(heat_flux_bottom_ / 1000.);        // W/m2, as used

    double      humidity_air = 0.;          // [-], constant relative humidity at TOP
    double      rain_mm_per_year = 100.;    // rain infiltration at TOP

    // Salinity top and model initial
    double      salinity_top(0.);
    double      salinity(0.);


    // ════════════════════════════════════════════════════════════════════════
    // 6. ROCK PROPERTIES
    // ════════════════════════════════════════════════════════════════════════
    double      perm(pow(10., -15));        // m2, baseline permeability
    double      max_perm(1.e-14);           // m2, upper bound enforced by the visitor
    double      min_perm(1.e-22);           // m2, lower bound
    double      log_perm(log10(perm));
    double      log_max_perm(log10(max_perm));
    double      log_min_perm(log10(min_perm));

    double      porosity(0.05);             // [-], element property
    double      thermal_conductivity(2.);   // W/m/K
    double      heat_capacity_rock(830.);   // J/kg/K
    double      density_rock(2650.);        // kg/m3 — CHECK CONSISTENCY IN THE MAGMA REGION
    double      compressibility_rock(1.0e-10);   // 1/Pa

    // Bulk heat capacity of the rock skeleton [J/m3/K].
    double      total_heat_capacity(heat_capacity_rock * density_rock * (1. - porosity));

    // ════════════════════════════════════════════════════════════════════════
    // 7. MAGMA PROPERTIES
    //
    // Used when with_magma_chamber is on; the crystallisation parameters matter
    // only when with_magma_model is also on.
    // ════════════════════════════════════════════════════════════════════════
    double      starting_temperature(950.);      // oC — must be >= the liquidus
    double      solidus_temperature(700.);       // oC
    double      liquidus_temperature(900.);      // oC
    double      crystal_density(2700.);          // kg/m3
    double      compressibility_magma(1.0e-8);   // 1/Pa

    double      initial_wtpercent_volatile_in_magma(0.);          // wt%, capped at saturation
    double      initial_volatile_vol_frac_chamber(porosity);      // [-], exsolved volatiles

    // Crystallisation curve. One id drives both the rock/melt latent-heat curve
    // used by the CVFEM scheme and the curve used inside MagmaModel, so the two
    // cannot drift apart.
    //   0 = power law         (uses crystallization_curve_exponent below)
    //   1 = error function    (uses sigma1_ below)
    //   2 = marxer ulmer
    uint32_t    crystallization_curve_id       = 0;
    double      crystallization_curve_exponent = 1.;    // b, power-law exponent
    double      sigma1_                        = 1.;    // erf law; also 2/3, 4/3

    std::string crystallization_curve;
    if      (crystallization_curve_id == 0) crystallization_curve = "power law";
    else if (crystallization_curve_id == 1) crystallization_curve = "error function";
    else if (crystallization_curve_id == 2) crystallization_curve = "marxer ulmer";
    else throw Exception(FATAL_ERROR, "CVFEM_topo_magma_air_example",
                        "Invalid crystallization curve: " + std::to_string(crystallization_curve_id) +
                            "\nAvailable options: 0 = power law; 1 = error function; 2 = marxer ulmer");

    // Magmatic fluid salinity [wt%] and Cl-partitioning switch — new MagmaModel ctor args.
    double      salinity_magmatic_fluids       = 0.;    // wt% NaCl in exsolved fluids
    bool        with_Cl_partitioning           = false;

    // ════════════════════════════════════════════════════════════════════════
    // 8. PERMEABILITY AND FAILURE VISITORS
    // ════════════════════════════════════════════════════════════════════════
    bool        depth_dependent(true);          // depth-dependent permeability
    double      depth_dependence_factor(1.);
    bool        temperature_dependent(true);     // temperature-dependent permeability
    bool        immediate_closure(true);         // close cracks as soon as fluid pressure drops
    bool        hydrofracturing(false);          // SET TO TRUE FOR FRACTURING
    bool        change_brittle_ductile_transtion_T(true);

    int         bdt_s(400);                      // oC, brittle-ductile transition start
    int         bdt_e(540);                      // oC, brittle-ductile transition end
    double      brittle_ductile_transition_start(1. * bdt_s);
    double      brittle_ductile_transition_end(1. * bdt_e);
    double      brittle_ductile_transition_ductile(0.5 * (bdt_s + bdt_e)); // midpoint
    double      log_k_start(-14.);      // background log10(k)
    double      log_k_duct(-17.);       // ductile-regime log10(k)
    double      log_k_end(-22.);        // minimum log10(k)

    double      C(10.e6);                        // Pa, rock cohesive strength
    double      diff_stress(0.);                 // Pa, set by the failure visitor
    bool        near_critically_pressured(true);
    bool        temperature_dependent_differential_stress(true);

    // ════════════════════════════════════════════════════════════════════════
    // 9. DERIVED TIME VARIABLES
    //
    // From here on every duration is in SECONDS.
    // ════════════════════════════════════════════════════════════════════════
    const double  year   = 31536000.;     //  in sec
    const double  day    = 86400.;        //  in sec
    const double  hour   = 3600.;         //  in sec
    const double  minute = 60.;           //  in sec

    double                          time_multiplier = 1.;
    if (time_unit == TimeUnit::YEARS)    time_multiplier = year;
    else if (time_unit == TimeUnit::DAYS)     time_multiplier = day;
    else if (time_unit == TimeUnit::HOURS)    time_multiplier = hour;
    else if (time_unit == TimeUnit::MINUTES)  time_multiplier = minute;

    largest_time_step   *= time_multiplier;
    initial_time_step   *= time_multiplier;
    max_time            *= time_multiplier;

    double     &model_time(ModelTime::Instance().modelTime);
    double      dt(1.);                          // initial time increment
    long        model_output(vtk_increment);

    ModelTimeToInteger<double> time_conversion;
    // ---------------------------------------------------
    // Further initialization of problem

    cerr<<endl<<"Read ansys mesh";

    if (MODEL_DIM == 2U)
    {
        ANSYS_Model2D model(icem_mesh_file.c_str(), region_file.c_str(), phys_var_file.c_str());
        model.OutputToBinaryFile(restart_file.c_str());
    }

    if (MODEL_DIM == 3U)
    {
        ANSYS_Model3D model(icem_mesh_file.c_str(), region_file.c_str(), phys_var_file.c_str());
        model.OutputToBinaryFile(restart_file.c_str());
    }

    cerr<<endl<<"Read ansys mesh done";
    Model<MODEL_DIM> model(restart_file);
    cerr<<endl<<"Created model";
    const PropertyDatabase<MODEL_DIM>    &p_ref = model.Database();

    cerr<<endl<<"InstantiateFiniteVolumes()";
    model.InstantiateFiniteVolumes();
    cerr<<endl<<"InstantiateFiniteVolumes() done";

    //---------------------------------------------------------


    model.                      InputPropertyValue("temperature",                              ScalarVariable(ANY, temperature_top));
    model.                      InputPropertyValue("thermal conductivity",                     ScalarVariable(ANY, thermal_conductivity));
    model.                      InputPropertyValue("nodal heat flux bottom",                   ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("thickness",                                ScalarVariable(ANY, 1.));

    //model.                      InputPropertyValue("total heat capacity",                      ScalarVariable(ANY, total_heat_capacity));
    model.                      InputPropertyValue("nodal heat capacity",                      ScalarVariable(ANY, total_heat_capacity));
    model.                      InputPropertyValue("nodal heat capacity rock",                 ScalarVariable(ANY, heat_capacity_rock));


    // Boundaries
    model.                      Boundary("TOP").InputPropertyValue("temperature",           ScalarVariable(DIRICH, temperature_top));
    HeatFluxBottom(model, heat_flux_bottom);
    cerr<<endl<<"HeatFluxBottom done";

    // Linear solver for the temperature systems; backend chosen by
    // solver_kind above. See CVFEM_SolverChoice.h.
    SolverBundle T_bundle(solver_kind);
    ApplySolverSettings(T_bundle);

    cerr << endl << "INIT TEMPERATURE";
    PDE_Integrator<MODEL_DIM>                    initial_temperature(T_bundle.Get());
    CVFEM_NumIntegral_dNT_op_dN_dV<MODEL_DIM>    conductance(p_ref, "thermal conductivity", "temperature", "temperature","thickness");//Benoit add thickness
    CVFEM_PointSource_rhsop<MODEL_DIM>           heat_bottom(p_ref, "nodal heat flux bottom", "temperature");

    // Calculate the initial geothermal gradient
    initial_temperature.Add(&conductance);
    initial_temperature.Add(&heat_bottom);
    model.Apply(initial_temperature);

    // setup to compute a temperature halo arounf magma chamber if with_temperature_halo=true
    PDE_Integrator<MODEL_DIM>              initial_temperature_halo(T_bundle.Get());
    NumIntegral_dNT_op_dN_dV<MODEL_DIM>    t_conductance( p_ref,"thermal conductivity","temperature","temperature" );
    NumIntegral_NT_lhsop_N_dV<MODEL_DIM>   t_capacitance_lhs( p_ref,"nodal heat capacity","temperature","temperature" );
    NumIntegral_NT_rhsop_N_dV<MODEL_DIM>   t_capacitance_rhs( p_ref,"nodal heat capacity","temperature" );
    CVFEM_PointSource_rhsop<MODEL_DIM>     t_heat_bottom(p_ref,"nodal heat flux bottom","temperature")   ;

    t_capacitance_lhs.LumpedFormulation(true);
    t_capacitance_rhs.LumpedFormulation(true);
    t_conductance.MultiplyWithTimeIncrement(true);
    t_heat_bottom.AddAccumulateLater();

    initial_temperature_halo.Add(&t_conductance);
    initial_temperature_halo.Add(&t_capacitance_lhs);
    initial_temperature_halo.Add(&t_capacitance_rhs);
    initial_temperature_halo.Add(&t_heat_bottom);
    initial_temperature_halo.TimeIncrement(initial_T_diffusion_time);
    initial_temperature_halo.Verbose(false);



    cerr<<endl<<"Initial temperature range: ";
    printRangeOfVariable(model, "temperature");


    if (with_magma_chamber)
    {
        // Magma chamber based on coordinates, quick test
        Point<2U> p1( 3000., 1000. );
        Point<2U> p2( 5000., 2000. );


        //model.FormRectangularRegion( "CHAMBER", p1, p2 );

        model.                      CreateProperty("porous flag element",                      "PFlagE",       "none", SCALAR, ELEMENT);
        model.                      CreateProperty("porous flag node",                         "PFlagN",       "none", SCALAR, NODE);
        model.                      CreateProperty("volatile dissolved fraction",            "volatiledissFrac",   "none", SCALAR, NODE);

        if (model.ContainsRegion("CHAMBER"))  model.Region("CHAMBER"). InputPropertyValue("porous flag element", ScalarVariable(ANY, 0.));
        if (model.ContainsRegion("CHAMBER2")) model.Region("CHAMBER2").InputPropertyValue("porous flag element", ScalarVariable(ANY, 0.));

        // Form regions based on the "porous flag element" variable
        model.FormRegionFrom("pluton",     "porous flag element", 0, 0.5);
        model.FormRegionFrom("cst_pluton", "porous flag element", 0, 0.5);
        model.FormRegionFrom("host rock",  "porous flag element", 0.5, 1.0);

        model.ExtrapolateCellToNodeProperty("porous flag element", "porous flag node");
    }

    if (model.ContainsRegion("pluton"))
    {
        model.Region("pluton").     InputPropertyValue("temperature", ScalarVariable(ANY, starting_temperature));

        if (with_temperature_halo)
        {
            model.Region("pluton").     ChangePropertyStatus("temperature", DIRICH);

            model.Apply( initial_temperature_halo );

            model.Region("pluton").     ChangePropertyStatus("temperature", ANY);
            model.Region("pluton").     InputPropertyValue("temperature", ScalarVariable(ANY, starting_temperature));

            cerr<<endl<<"Initial temperature range after halo: ";
            printRangeOfVariable(model, "temperature");
        }
    }

    model.                      InterpolateNodeToCellProperty("temperature", "temperature element");
    model.                      FormRegionFrom("high temperature region", "temperature element", brittle_ductile_transition_end, 1100.);
    model.                      FormRegionFrom("low temperature region", "temperature element", 0., brittle_ductile_transition_end);

    //---------------------------------------------------------


    model.InstantiateFiniteVolumes(); // New in open-csmp!
    model.CopyRegion("Model", "FLOW REGION");//the whole model is a flow region


    long time_label = TimeLabel(time_unit, time_conversion, model_time);

    // ---------------------------------------------------
    // Create extra properties and input initial values
    PoreVolumeVisitor<MODEL_DIM> Pore_Volume_Visitor(model, "nodal porosity", "bulk volume", "pore volume", "thickness");



    model.                      CreateProperty("nfvs scaling factor",                   "nfvs scaling factor",  "none", SCALAR, NODE);
    model.                      InputPropertyValue("nfvs scaling factor",               ScalarVariable(ANY, 1.));

    //-----
    model.                      InputPropertyValue("thickness",                                ScalarVariable(ANY, 1.));



    // fracturing reference: 0 or 1 per element, gating whether PermeabilityVisitor
    // applies the normal laws (>= 1) or forces log_k_min (< 1). 1 globally here,
    // overridden to 0 in the pluton below. The old ComputeSimpleFracturableFlag
    // wrote 0/1 per element based on solidus temperature, but it hard-codes
    // model.Region("CHAMBER"), which does not exist here (this example uses
    // pluton / cst_pluton), so it would throw. Region-based masking is what
    // GeothermalSimulator does.
    model.                      InputPropertyValue("fracturing reference",                     ScalarVariable(ANY, 1.));
    model.                      InputPropertyValue("failure pressure",                         ScalarVariable(ANY, pressure_top));
    model.                      InputPropertyValue("permeability",                             ScalarVariable(ANY, perm));

    // gradP scaling: read in PermeabilityVisitor::CalculateGradPScaling and used
    // as a multiplier on k() when 0 < gradP_sc < 1. The property default in
    // PhysicalVariablesBenchmarks.txt is 1e-10, which would silently drop
    // permeability by ten orders of magnitude; initialise it to 1 (no correction).
    model.                      InputPropertyValue("gradP scaling",                             ScalarVariable(ANY, 1.));

    // Keys the new PermeabilityVisitor reads. Without "BDT start permeability
    // element" the visitor reads 0 and log10(0) = -inf, silently killing the
    // depth-dependent branch. "permeability anisotropy factor" is read only on
    // the anisotropic path, but initialise it so it's never NaN.
    model.                      InputPropertyValue("BDT start permeability element",            ScalarVariable(ANY, perm));
    model.                      InputPropertyValue("permeability anisotropy factor",            ScalarVariable(ANY, 1.));

    model.CopyReplace("permeability","vertical permeability");
    model.CopyReplace("permeability","horizontal permeability");

    PermeabilityVisitor<MODEL_DIM> permeability_visitor(model);
    if (change_brittle_ductile_transtion_T)
        permeability_visitor.ChangeBrittleDuctileTransitionTemperature(
            true,
            brittle_ductile_transition_start,
            brittle_ductile_transition_ductile,
            brittle_ductile_transition_end,
            log_k_start,
            log_k_duct,
            log_k_end);
    if (depth_dependent)
        permeability_visitor.DepthDependent(true);
    if (temperature_dependent)
        permeability_visitor.TemperatureDependent(true);

    permeability_visitor.AssignPermeabilityTensor(model);

    model.                      InputPropertyValue("nodal density rock",                       ScalarVariable(ANY, density_rock));
    model.                      InputPropertyValue("fluid density",                            ScalarVariable(ANY, 0.));

    model.                      InputPropertyValue("fluid pressure",                           ScalarVariable(ANY, pressure_top));
    model.                      InputPropertyValue("lithostatic pressure",                     ScalarVariable(ANY, pressure_top));
    model.                      InputPropertyValue("nodal porosity",                           ScalarVariable(ANY, porosity));
    model.                      InputPropertyValue("nodal compressibility rock",               ScalarVariable(ANY, compressibility_rock));
    model.                      InputPropertyValue("courant liquid",                           ScalarVariable(ANY, largest_time_step));
    model.                      InputPropertyValue("courant vapor",                            ScalarVariable(ANY, largest_time_step));
    model.                      InputPropertyValue("courant air",                              ScalarVariable(ANY, largest_time_step));
    model.                      InputPropertyValue("salinity",                                 ScalarVariable(ANY, salinity));
    model.                      InputPropertyValue("fluid source rate",                        ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fluid source h",                           ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fluid source wt",                          ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("nodal source liquid",                      ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("nodal source vapor",                       ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("nodal source air",                         ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("upwind control liquid",                    ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("upwind control vapor",                     ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("upwind control air",                       ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fracturing events",                        ScalarVariable(ANY, 0.));

    // Magmatic variables
    model.                      InputPropertyValue("magmatic fluid mass",                  ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic ratio",                       ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic water ratio",                 ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic fluid mass liquid",           ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic liquid mass mobility",        ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic fluid mass vapor",            ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic vapor mass mobility",         ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic mass salt",                   ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic salt ratio",                  ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic salt content liquid",         ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic liquid salt mobility",        ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic salt content vapor",          ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic vapor salt mobility",         ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("NaCl KCl ratio",                       ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("mass KCl",                             ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("salt content liquid KCl",              ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("liquid salt mobility KCl",             ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("salt content vapor KCl",               ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("vapor salt mobility KCl",              ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("copper content fluid",                 ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("copper fraction fluid",                ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("copper content liquid",                ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("copper content vapor",                 ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("liquid copper mobility",               ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("vapor copper mobility",                ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("copper precipitation",                 ScalarVariable(ANY, 0.));

    // Fluid trackers
    model.                      InputPropertyValue("salt flux integral",                   ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("energy flux integral",                 ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fluid flux integral",                  ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic fluid flux integral",         ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("shell magmatic fluid flux integral",   ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("magmatic production rate",             ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fluid flux per second",                ScalarVariable(ANY, 0.));
    model.                      InputPropertyValue("fluid flux with direction integral",   VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("fluid flux in",                ScalarVariable(ANY, 0.));


    // Velocity
    model.                      InputPropertyValue("velocity",                             VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("pore velocity",                        VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("velocity liquid",                      VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("pore velocity liquid",                 VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("nodal pore velocity",                  VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("velocity vapor",                       VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("pore velocity vapor",                  VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("velocity air",                       VectorVariable<MODEL_DIM>(ANY, 0.));
    model.                      InputPropertyValue("pore velocity air",                  VectorVariable<MODEL_DIM>(ANY, 0.));

    // Boundaries
    model.                      Boundary("TOP").InputPropertyValue("lithostatic pressure",  ScalarVariable(DIRICH, pressure_top));//Wrong if TOP is not the surface!!!
    model.                      Boundary("TOP").InputPropertyValue("fluid pressure",        ScalarVariable(DIRICH, pressure_top));
    //model.                      Boundary("BOTTOM").InputPropertyValue("fluid pressure",     ScalarVariable(DIRICH, pressure_top*100));

    model.                       Boundary("TOP").InputPropertyValue("salinity",              ScalarVariable(ANY, salinity_top));



    // Pluton overrides. NOTE: this block must come AFTER the global
    // InputPropertyValue("fracturing reference", 1.) above — otherwise the
    // global 1 would clobber the pluton's 0 and the intrusion would be treated
    // as fracturable, defeating the whole point.
    if (model.ContainsRegion("pluton"))
    {
        model.Region("pluton").     InputPropertyValue("fracturing reference",        ScalarVariable(ANY, 0.));
        model.Region("pluton").     InputPropertyValue("volatile dissolved fraction", ScalarVariable(ANY, initial_wtpercent_volatile_in_magma));
        model.Region("pluton").     InputPropertyValue("nodal porosity",              ScalarVariable(ANY, initial_volatile_vol_frac_chamber));
        model.Region("pluton").     InputPropertyValue("porosity",                    ScalarVariable(ANY, initial_volatile_vol_frac_chamber));
        model.Region("pluton").     InputPropertyValue("salinity",                    ScalarVariable(ANY, salinity));
    }

    model.Region("Model").      AssignNodeCoordinatesTo("nodal depth",  'y');
    model.Region("Model").      AssignNodeCoordinatesTo("coordinate x", 'x');
    model.Region("Model").      AssignNodeCoordinatesTo("coordinate y", 'y');
    if (MODEL_DIM == 3U)
        model.Region("Model").  AssignNodeCoordinatesTo("coordinate z", 'z');
    model.InterpolateNodeToCellProperty("nodal depth", "barycenter depth");


    // Calculate liquidus and solidus temperatures (depth-dependent)
    cerr << endl << "ComputeLiquidusSolidus";
    ComputeLiquidusSolidus(model, "liquidus temperature", "solidus temperature", "nodal depth", liquidus_temperature, solidus_temperature);
    cerr << endl << "HeatFluxBottom";
    HeatFluxBottom(model, heat_flux_bottom);

    cerr << endl << "PORE VOLUME VISITOR";
    model.                          Accept(Pore_Volume_Visitor);
    model.                          InterpolateNodeToCellProperty("nodal porosity", "porosity");



    // ---------------------------------------------------
    // Output options
    list<string> output_variables;

    // VTU field list. See DefineOutputVariables() below.
    DefineOutputVariables({ with_air, with_magma_chamber,
                           treat_air_and_vapor_as_mixture,
                           treat_air_and_vapor_as_hydraulic_mixture },
                          output_variables);
    // ---------------------------------------------------



    // ---------------------------------------------------
    // Visitors instantiation

    CVFEM_PorosityChangeVisitor<MODEL_DIM>        *d_phi_visitor = nullptr;
    RatioVisitor_split1<MODEL_DIM>                *magmatic_ratio = nullptr;
    RatioVisitor_split2<MODEL_DIM>                *magmatic_ratio2 = nullptr;

    PropertyHandle<MODEL_DIM> element_id(model, "element ID", SCALAR, ELEMENT);
    PropertyHandle<MODEL_DIM> node_id(model, "node ID", SCALAR, NODE);
    PropertyHandle<MODEL_DIM> dt_handle(model, "dt", SCALAR, NODE);
    csmp::Index   n_id_key(p_ref.StorageKey("node ID"));
    csmp::Index   e_id_key(p_ref.StorageKey("element ID"));
    ScalarVariable ID;

    for (auto n_it = model.Mesh().NodesBegin(); n_it != model.Mesh().NodesEnd(); n_it++)
    {
        ID() = n_it->Idx();
        n_it->Store(n_id_key, ID);
    }
    for (auto e_it = model.Mesh().ElementsBegin(); e_it != model.Mesh().ElementsEnd(); e_it++)
    {
        ID() = e_it->Idx();
        e_it->Store(e_id_key, ID);
    }



    // ---------------------------------------------------
    // CVFEM PHX Scheme
    cerr << endl << "SETTING CVFEM OPTIONS";

    // No wells in this example: the empty well list leaves the scheme's wellbore
    // machinery dormant. Lithium and zinc transport are likewise off.
    CVFEM_PHX_Scheme<MODEL_DIM> *CVFEM_PHX = nullptr;
    CVFEM_PHX = new CVFEM_PHX_Scheme<MODEL_DIM>(
        model,
        with_gravity,
        true,        // thermodynamic_density_in_gravity_term
        {},          // wells: none
        false,       // with_zinc
        false,       // with_lithium
        solver_kind);


    CVFEM_PHX->  SetLargestTimeStep(largest_time_step);
    CVFEM_PHX->  SetCutLogLabel(output_name); // names the dt-cut CSV: dt_cuts_<output_name>_<timestamp>.csv
    CVFEM_PHX->  SetEquilibratorConvergenceSpeedUpTo(true);
    CVFEM_PHX->  Adjust_CFL_Criterion(cfl_scaling, true); //true= pore velocity based//////////////////////////

    if(with_air) CVFEM_PHX->ActivateAir(humidity_air, treat_air_and_vapor_as_mixture, treat_air_and_vapor_as_hydraulic_mixture, rain_mm_per_year);//Needs to come before openboundaries!

    CVFEM_PHX->  WithRockLiquidusSolidus();

    CVFEM_PHX->  SetRockHeatCapacity(heat_capacity_rock);
    CVFEM_PHX->  SetRockCrystallizationCurve(1.78, sigma1_, 300000., crystallization_curve_exponent, crystallization_curve);
    CVFEM_PHX->  AddFluidContributionToHeatCapacity(true);
    CVFEM_PHX->  AttemptToSurviveFluidPropertiesError(true);

    if (open_top)
    {
        CVFEM_PHX->OpenBoundaries(salinity_top);
    }

    cerr << endl << "END SETTING CVFEM OPTIONS";

    // Linear solver for the steady pressure systems.
    SolverBundle P_bundle(solver_kind);
    ApplySolverSettings(P_bundle);

    PDE_Integrator<MODEL_DIM>  lithos_pressure(P_bundle.Get());
    PDE_Integrator<MODEL_DIM>  initial_pressure_(P_bundle.Get());
    NumIntegral_dNT_op_dN_dV<MODEL_DIM>    conductance_p(p_ref,"permeability","lithostatic pressure","lithostatic pressure");
    NumIntegral_NT_op_dNi_dV<MODEL_DIM>    density_p(p_ref,"nodal density rock","permeability","lithostatic pressure");

    NumIntegral_dNT_op_dN_dV<MODEL_DIM>    conductance_hp(p_ref,"permeability","fluid pressure","fluid pressure");
    NumIntegral_NT_op_dNi_dV<MODEL_DIM>    density_hp(p_ref,"fluid density","permeability","fluid pressure");


    cerr << endl << "INIT PRESSURE";
    // Calculate the initial lithostatic gradient assuming homogeneous rock density (carefull, it doesn't consider the "bulk" density)
    lithos_pressure.Add(&conductance_p);
    lithos_pressure.Add(&density_p);
    model.Apply(lithos_pressure);

    // Hydrostatic pressure initialization is done further (in equilibration loop)
    initial_pressure_.Add(&conductance_hp);
    initial_pressure_.Add(&density_hp);
    model.Apply(initial_pressure_);

    model.Accept(permeability_visitor);
    model.CopyReplace("permeability","vertical permeability");
    model.CopyReplace("permeability","horizontal permeability");
    permeability_visitor.AssignPermeabilityTensor(model);

    int k = 0;
    bool flag(false);
    double diff(0.), eps(1000.0);
    while (!flag)
    {
        cerr << endl << "Initialization P-T, iteration: " << k;
        printRangeOfVariable(model, "temperature");
        model.CopyReplace("fluid pressure", "test fluid pressure");
        CVFEM_PHX->InitialFluidPropertiesFromPTX();

        model.Apply(initial_pressure_);

        diff = maxDifferenceScalarNodeProperty(model, "fluid pressure", "test fluid pressure");
        cerr << endl << "*********************************Max diff pressure: " << diff;
        if (diff < eps) flag = true;
        k++;
    }

    model.CopyReplace("fluid pressure", "reference pressure");

    k = 0;
    if (with_lithostatic_ductile)
    {
        flag = false;
        diff = 0.;
        eps = 1000.0;
        while (!flag)
        {
            cerr << endl << "Initialization P-T with lithostatic ductile, iteration: " << k;
            printRangeOfVariable(model, "temperature");
            model.CopyReplace("fluid pressure", "test fluid pressure");
            CVFEM_PHX->InitialFluidPropertiesFromPTX();

            model.Apply(initial_pressure_);

            if (model.ContainsRegion("high temperature region"))
            {
                model.Region("high temperature region").CopyReplace("lithostatic pressure", "fluid pressure");
            }


            diff = maxDifferenceScalarNodeProperty(model, "fluid pressure", "test fluid pressure");
            cerr << endl << "*********************************Max diff pressure: " << diff;
            if (diff < eps) flag = true;
            k++;
        }
    }

    CVFEM_PHX->InitialFluidPropertiesFromPTX();
    cerr << endl << "FINISH initialization P-T";

    model.CopyReplace("fluid pressure", "previous fluid pressure");
    model.CopyReplace("previous fluid pressure", "previous previous fluid pressure");
    model.CopyReplace("previous total enthalpy", "previous previous total enthalpy");

    if (model.ContainsRegion("high temperature region"))
    {
        model.RemoveRegion("high temperature region", false);
    }

    if (model.ContainsRegion("low temperature region"))
        model.RemoveRegion("low temperature region", false);


    if (depth_dependent||temperature_dependent)
    {
        model.Accept( permeability_visitor );
        model.CopyReplace("permeability","vertical permeability");
        model.CopyReplace("permeability","horizontal permeability");
        permeability_visitor.AssignPermeabilityTensor(model);
    }

    if (hydrofracturing)    permeability_visitor.Hydrofracturing(true, log_max_perm, log_min_perm);
    if (immediate_closure)  permeability_visitor.SetImmediateClosureTo(immediate_closure);

    model.                      CopyReplace( "permeability"             ,   "previous permeability" );
    model.                      CopyReplace( "previous permeability"    ,   "previous previous permeability" );
    cerr<<endl<<"CVFEM_PHX->PrepareTransientCalculations();";
    CVFEM_PHX->PrepareTransientCalculations();

    // ---------------------------------------------------
    // calculate failure pressure
    cerr<<endl<<"FailureModeVisitor";
    FailureModeVisitor<MODEL_DIM> failure_mode(model, C, diff_stress, near_critically_pressured);
    if (temperature_dependent_differential_stress)
        failure_mode.SetTemperatureRelaxationTo(temperature_dependent_differential_stress);
    if (change_brittle_ductile_transtion_T)
        failure_mode.ChangeBrittleDuctileTransitionTemperature(brittle_ductile_transition_start, brittle_ductile_transition_end);
    model.                      Accept(failure_mode);

    cerr<<endl<<"OutputToVTU initial";
    OutputToVTU( model, output_name+"_initial", output_variables, 0);

    CVFEM_PHX->ChangeTimeStepTo(initial_time_step);

    double start_time_minute = model_time / minute;
    if (no_std_output)
    {
        freopen("/dev/null", "w", stdout);
    }

    ConservationCheck<MODEL_DIM> cons(model, "TestConsVisitNEW", with_air);
    cons.BeginPass(); model.Accept(cons); cons.EndPass(0., 0.);   // baseline, once

    OutputAllVariableRanges(model, output_name, "after_initialisation");

    cerr << "\n\n---- TIME LOOP ----";
    bool stop_sim(false);

    // ── Wall-clock accounting ───────────────────────────────────────────────
    // Covers the whole iteration, so it is a like-for-like cost per timestep when
    // comparing solver backends. Steps are counted per loop iteration, so a run
    // that cuts dt repeatedly does many cheap steps and gets a flattering mean —
    // read it together with "model time reached" below.
    using clock_type = std::chrono::steady_clock;
    const auto run_start = clock_type::now();
    double accumulated_step_seconds(0.), slowest_step_seconds(0.);
    size_t step_count(0U);

    while (model_time <= max_time && !stop_sim)
    {
        const auto step_start = clock_type::now();


        cerr << defaultfloat;
        cerr.precision(8);

        // All-in-one T diffusion, P diffusion, advection

        CVFEM_PHX->SetModelTime(model_time);
        dt = CVFEM_PHX->Apply();

        cons.BeginPass(); model.Accept(cons); cons.EndPass(model_time, dt);


        cerr << endl << endl << "////////////////////////";
        cerr << endl << output_name.c_str() ;
        cerr << endl << "Start time           : " << start_time_minute;
        cerr << endl << "Time in days      : " << model_time / day << ",             in years:  " << model_time / year;
        cerr << endl << "Timestep in minutes  : " << dt / minute << ",             in days:  " << dt / day;
        cerr << endl << "Progress             : " << model_time / max_time * 100 << "%";


        model.                      CopyReplace("previous permeability", "previous previous permeability");
        model.                      CopyReplace("permeability", "previous permeability");


        // additional functionality
        // permeability stuff
        model.                      Accept(failure_mode);
        model.InterpolateNodeToCellProperty("temperature", "temperature element");
        // (SetTimeIncrement(dt) was called here. SimplePermeabilityVisitor stored the
        //  value and never read it, so it did nothing; PermeabilityVisitor has no
        //  such method.)
        model.                      Accept(permeability_visitor);


        model.CopyReplace("permeability", "vertical permeability");
        model.CopyReplace("permeability", "horizontal permeability");
        permeability_visitor.AssignPermeabilityTensor(model);

        // ---------------------------------------------------
        // preparing next time step
        model_time += dt;
        model.InputPropertyValue("dt", ScalarVariable(ANY, dt));

        // ---------------------------------------------------
        // output and model saving. WriteTimestepOutput advances model_output.
        WriteTimestepOutput(model, output_name, time_unit, time_conversion, model_time,
                            output_variables, model_output, vtk_increment);

        if (dt < 1.)
        {
            cerr << endl << "Timestep smaller than 1 s Aborting ...";
            stop_sim = true;
        }

        {
            const double step_seconds =
                std::chrono::duration<double>(clock_type::now() - step_start).count();
            accumulated_step_seconds += step_seconds;
            if (step_seconds > slowest_step_seconds) slowest_step_seconds = step_seconds;
            ++step_count;
            cerr << endl << "Step wall-clock      : " << step_seconds << " s"
                 << "   (average " << accumulated_step_seconds / step_count
                 << " s over " << step_count << " steps)";
        }
    }

    // ── Run summary ─────────────────────────────────────────────────────────
    {
        const double total_seconds =
            std::chrono::duration<double>(clock_type::now() - run_start).count();
        const double mean_seconds =
            step_count ? accumulated_step_seconds / step_count : 0.;

        cerr << "\n\n---- RUN SUMMARY ----";
        cerr << "\n  linear solver:        " << SolverKindName(solver_kind);
        cerr << "\n  magma chamber:        " << (with_magma_chamber ? "ON" : "off");
        cerr << "\n  air phase:            " << (with_air ? "ON" : "off");
        cerr << "\n  timesteps completed:  " << step_count;
        cerr << "\n  model time reached:   " << model_time / year << " years"
             << "   of " << max_time / year << " requested";
        cerr << "\n  mean time per step:   " << mean_seconds << " s";
        cerr << "\n  slowest step:         " << slowest_step_seconds << " s";
        cerr << "\n  total in time loop:   " << accumulated_step_seconds << " s";
        cerr << "\n  total wall clock:     " << total_seconds << " s";
    }

    OutputAllVariableRanges(model, output_name, "end_of_run");

    cerr << endl << "That's it..." << endl;
    throw Exception(ERROR, "CVFEM_topo_magma_air_example", "Done");
}


// Build the VTU field list. Which fields appear depends on the features in play,
// so the flags come in as a small struct rather than three separate arguments.

void CVFEM_topo_magma_air_example::DefineOutputVariables(const OutputFlags &p,
                                                         std::list<std::string> &output_variables) const {

    if (p.with_air)
    {
        output_variables.push_back("density air");
        output_variables.push_back("enthalpy air");
        output_variables.push_back("boundary flow air");
        output_variables.push_back("fluid mass air");
        output_variables.push_back("enthalpy content air");
        output_variables.push_back("saturation air");
        output_variables.push_back("viscosity air");
        output_variables.push_back("volumetric enthalpy air");
        output_variables.push_back("density air transport");
        output_variables.push_back("nodal compressibility air");
        output_variables.push_back("air mass mobility");
        output_variables.push_back("air mass mobility density");
        output_variables.push_back("velocity air");
        output_variables.push_back("pore velocity air");
        output_variables.push_back("courant air");
        output_variables.push_back("relperm viscosity air");

        // ── Boundary / rain variables ───────────────────────────────────────

        output_variables.push_back("specific enthalpy vapor top");
        output_variables.push_back("specific enthalpy air top");
        output_variables.push_back("mass fraction inflow vapor in air top");
        output_variables.push_back("mass fraction inflow air top");
        output_variables.push_back("rain recharge top");
        //output_variables.push_back("reference specific enthalpy top");

        if (p.treat_air_and_vapor_as_mixture or p.treat_air_and_vapor_as_hydraulic_mixture)
        {
            // ── Gas-phase mixture diagnostics (Dalton's law mode) ──────────
            output_variables.push_back("gas mix molar fraction air");
            output_variables.push_back("gas mix molar fraction vapor");
            output_variables.push_back("partial pressure air");
            output_variables.push_back("partial pressure vapor");
            output_variables.push_back("gas mix density");
            output_variables.push_back("gas mix viscosity");
            output_variables.push_back("gas mix enthalpy");
            output_variables.push_back("gas mix compressibility");
        }
    }

    output_variables.push_back("specific enthalpy liquid top");
    output_variables.push_back("courant liquid");
    output_variables.push_back("courant vapor");
    output_variables.push_back("relperm viscosity liquid");
    output_variables.push_back("relperm viscosity vapor");
    output_variables.push_back("node ID");
    output_variables.push_back("element ID");
    output_variables.push_back("thickness");
    output_variables.push_back("permeability tensor");
    output_variables.push_back("nodal depth");
    output_variables.push_back("barycenter depth");
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
    output_variables.push_back("velocity");
    output_variables.push_back("pore velocity");
    output_variables.push_back("pore velocity liquid");
    output_variables.push_back("pore velocity vapor");
    output_variables.push_back("permeability");
    output_variables.push_back("nodal permeability");
    output_variables.push_back("nodal total compressibility");
    output_variables.push_back("compressibility");
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
    output_variables.push_back("nodal fluid volume source");
    output_variables.push_back("nfvs scaling factor");
    output_variables.push_back("nodal porosity");
    output_variables.push_back("nodal heat capacity rock");
    output_variables.push_back("fluid heat capacity");
    output_variables.push_back("nodal heat capacity");
    output_variables.push_back("coordinate x");
    output_variables.push_back("coordinate y");

}

long CVFEM_topo_magma_air_example::TimeLabel(TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion,
                                             double model_time) const
{
    switch (time_unit)
    {
    case TimeUnit::YEARS:   return time_conversion.ModelTimeInYears(model_time);
    case TimeUnit::DAYS:    return time_conversion.ModelTimeInDays(model_time);
    case TimeUnit::HOURS:   return time_conversion.ModelTimeInHours(model_time);
    case TimeUnit::MINUTES: return time_conversion.ModelTimeInMinutes(model_time);
    case TimeUnit::SECONDS: return (long)floor(model_time + 0.5);
    }
    return 0L;
}

void CVFEM_topo_magma_air_example::WriteTimestepOutput(Model<MODEL_DIM> &model, const string &output_name,
                                                       TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion,
                                                       double model_time,
                                                       const list<string> &output_variables,
                                                       long &model_output, long vtk_increment) const
{
    const long time_label = TimeLabel(time_unit, time_conversion, model_time);

    if (time_label >= model_output)
    {
        OutputToVTU(model, output_name, output_variables, time_label);
        model_output += vtk_increment;
        cerr << endl << "Writing file, time_label: " << time_label;
    }
}

// Dump min/max of every property in the database. Two of these files diffed
// against each other show exactly which fields two runs disagree on — far faster
// than opening VTUs.
void CVFEM_topo_magma_air_example::OutputAllVariableRanges(Model<MODEL_DIM> &model,
                                                           const string &output_name,
                                                           const string &tag) const
{
    const string filename = output_name + "_VariableRanges_" + tag + ".txt";
    std::ofstream out(filename);
    if (!out)
    {
        cerr << endl << "  could not open " << filename << " for writing";
        return;
    }

    out << "Model:\t" << output_name << "\n";
    out << "Stage:\t" << tag << "\n\n";
    out << "variable\tmin\tmax\n";

    double min_value(0.), max_value(0.);
    size_t n(0U);
    for (auto it = model.Database().Begin(); it != model.Database().End(); ++it)
    {
        model.MinMaxOf((*it).first.c_str(), min_value, max_value);
        out << (*it).first.c_str() << "\t" << min_value << "\t" << max_value << "\n";
        ++n;
    }
    cerr << endl << "  wrote ranges of " << n << " variables to " << filename;
}

void CVFEM_topo_magma_air_example::OutputToVTU(Model<MODEL_DIM> &model, const string &model_name, const list<string> &props, size_t timestep) const
{
    static VTU_Interface<MODEL_DIM> vtu(model);
    vtu.OutputDataToVTU((model_name).c_str(), props, model.Region("Model"), timestep);
}



double CVFEM_topo_magma_air_example::maxDifferenceScalarNodeProperty(Model<MODEL_DIM> &mdl, const char *snp1, const char *snp2)
{
    csmp::Index snp1Key(mdl.Database().StorageKey(snp1));
    csmp::Index snp2Key(mdl.Database().StorageKey(snp2));
    const Region<MODEL_DIM>   &mref = mdl.Region("Model");

    double max(0.);

    vector<Node<MODEL_DIM>*>::const_iterator nit;
    for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit)
    {
        double diff = fabs((*nit)->Read(snp1Key) - (*nit)->Read(snp2Key));

        if (diff > max)
            max = diff;
    }
    return max;
}

void CVFEM_topo_magma_air_example::ComputeLiquidusSolidus(Model<MODEL_DIM> &model, const char *liquidus, const char *solidus, const char *depth, double TL, double TS)
{
    const Region<MODEL_DIM>   &mref = model.Region("Model");
    csmp::Index liquidusKey(model.Database().StorageKey(liquidus));
    csmp::Index solidusKey(model.Database().StorageKey(solidus));
    csmp::Index depthKey(model.Database().StorageKey(depth));

    ScalarVariable tl(ANY, 0.);
    ScalarVariable ts(ANY, 0.);
    ScalarVariable ndepth(ANY, 0.);

    vector<Node<MODEL_DIM>*>::const_iterator nit;
    for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit)
    {
        ndepth = (*nit)->Read(depthKey);
        tl() = TL;
        ts() = TS;
        (*nit)->Store(liquidusKey, tl);
        (*nit)->Store(solidusKey, ts);
    }
}

// NOTE: the local assignKtensor member function has been removed. All call sites
// now use permeability_visitor.AssignPermeabilityTensor(model), which is the canonical
// implementation (identical logic) and avoids duplicating it here.

void CVFEM_topo_magma_air_example::HeatFluxBottom(Model<MODEL_DIM> &model, double heat_flux_bottom)
{

    csmp::Index         hfb_key(model.Database().StorageKey("nodal heat flux bottom"));
    csmp::Index         thickness_key(model.Database().StorageKey("thickness"));
    ScalarVariable      area, thickness, hfb;

    Boundary<MODEL_DIM>&   mref = model.Boundary("BOTTOM");
    for ( auto eit = mref.CellsBegin(); eit != mref.CellsEnd(); eit++ )
    {

        for ( auto i{0U}; i<(*eit)->Nodes(); i++ )
        {
            (*eit)->N( i )->Read( hfb_key, hfb);

            double area = (*eit)->Area();
            area /= (*eit)->Nodes();

            hfb() += heat_flux_bottom * area;

            (*eit)->N( i )->Store( hfb_key, hfb );
        }
    }
}
} //csmp