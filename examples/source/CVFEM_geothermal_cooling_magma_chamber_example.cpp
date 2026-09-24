// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

// =====================================================================

// =====================================================================

#include "CVFEM_geothermal_cooling_magma_chamber_example.h"


#ifdef CSMP_WITH_PETSC_SOLVER
#include <petscsys.h>
#endif

using namespace std;
namespace csmp {

// ── Dimension-specific mesh construction ────────────────────────────────────
// The only thing that actually differs between 2D and 3D. Specialized once
// per dim; the constructor itself (below) is written generically and calls
// this, rather than existing as two independently-maintained copies.
template<>
Model<2U>* CVFEM_geothermal_cooling_magma_chamber_example<2U>::CreateANSYSModel(
    const std::string& mesh, const std::string& region, const std::string& vars)
{
    return new ANSYS_Model2D(mesh.c_str(), region.c_str(), vars.c_str());
}

template<>
Model<3U>* CVFEM_geothermal_cooling_magma_chamber_example<3U>::CreateANSYSModel(
    const std::string& mesh, const std::string& region, const std::string& vars)
{
    return new ANSYS_Model3D(mesh.c_str(), region.c_str(), vars.c_str());
}

// ── Constructor ──────────────────────────────────────────────────────────────
// argv[1] = mesh file basename, argv[2] = region file basename,
// argv[3] = output path, argv[4] = output name.
template<uint32_t dim>
CVFEM_geothermal_cooling_magma_chamber_example<dim>::CVFEM_geothermal_cooling_magma_chamber_example(int argc, char** argv)
    : geometry_name_(argc > 1 ? argv[1] : "2D-intrusion"),
    config_file_name_  (argc > 2 ? argv[2] : "CVFEM_geothermal"),
    output_path_  (argc > 3 ? argv[3] : "./"),
    output_name_  (argc > 4 ? argv[4] : "CVFEM_geothermal_example"),
    vars_name_("PhysicalVariables.txt"), // created by the ModelBuilder

    model(nullptr),
    pd_ref(nullptr),
    CVFEM_PHX(nullptr),
    permeability_visitor(nullptr),
    Pore_Volume_Visitor(nullptr),
    failure_mode(nullptr),
    vtu(nullptr)
{
#ifdef CSMP_WITH_PETSC_SOLVER
    PetscInitializeNoArguments();
#endif

    solver_kind      = SolverKind::PETSc; // change solver here
    well_solver_kind = SolverKind::PETSc; // change solver here

} // end constructor

template<uint32_t dim>
CVFEM_geothermal_cooling_magma_chamber_example<dim>::~CVFEM_geothermal_cooling_magma_chamber_example()
{
    delete permeability_visitor;
    delete Pore_Volume_Visitor;
    delete failure_mode;
    delete CVFEM_PHX;
    delete model;

#ifdef CSMP_WITH_PETSC_SOLVER
    PetscFinalize();
#endif
}



template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::Specifications() {
    SetTitle("Magmatic heat-driven hydrothermal convection");
    SetDifficulty(0);
    SetCategory("CVFEM Examples");
    AddAuthor("JK");
    AddDescription("source file in: /examples/example_inputs/CVFEM_examples");
}




template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::Run() {

    //! Create Model and read *-configuration.txt file
    // CreateModel();

    //! Time variables and objects
    double& model_time  = ModelTime::Instance().modelTime;

    SetTimeVariables();

    //! Problem configuration
    bool no_std_output = model->Read( pd_ref->StorageKey( "no standard output" ) );

    //! Output variables
    DefineOutputVariables();
    output_regions_list = {}; // list of specific regions to output additional .vtu-files for

    ///////////////////////////////////////////////////////
    //! Define extra properties and input initial values
    InitializeParameters();

    //! Calculate pore-volume
    Pore_Volume_Visitor = new PoreVolumeVisitor<dim> (*model, "nodal porosity", "bulk volume", "pore volume", "thickness");
    model->Accept(*Pore_Volume_Visitor);
    model->CopyReplace("pore volume", "previous pore volume");
    model->InterpolateNodeToCellProperty("nodal porosity", "porosity"); // Note that this interpolation creats a minor transition between two regions and no sharp boundaries

    //! CVFEM PHX Scheme instantiation
    InstatiateCVFEM();

    ///////////////////////////////////////////////////////
    //! Define initial permeability
    // Permeability calculation based on defined criterions
    // e.g., depth-dependent, temperature-dependent, etc.
    InstantiatePermeabilityVisitor();
    if (hydrofracturing)
        InstantiateFailureModeVisitor();

    //! Initialize background temperature and hydrostatic pressure
    // Temperature defined by assigned bottom heat flux
    ApplyBottomHeatFlux();
    InitializeTemperature();

    // Calculate permeability
    model->Region("Model").Accept( *permeability_visitor );
    permeability_visitor->AssignPermeabilityTensor(*model);

    // Hydrostatic pressure
    InitializeHydrostaticPressure();

    //! Emplace magmatic intrusion
    // Set intrusion temperature;
    // optional: calculate initial temperature diffusion to mimic heat diffusion during magma emplacement
    AddMagmaticIntrusion();

    // Calculate permeability
    model->Region("Model").Accept( *permeability_visitor );
    permeability_visitor->AssignPermeabilityTensor(*model);

    // Recalculate pressure and apply lithostatic pressure in impermeable regions
    InitializePressure();

    //! Activate hydrofracturing and calculate failure pressure
    ActivateFracturing(hydrofracturing);

    //! Prepare Transient Calculations
    CVFEM_PHX->PrepareTransientCalculations();

    //! Create initial output file
    long time_label;
    time_label = model_time / time_multiplier;

    OutputToVTU(output_path_ + output_name_, output_variables, 0, output_regions_list);


    ///////////////////////////////////////////////////////
    ///
    //! Main loop
    ///
    // Set initial timestep
    CVFEM_PHX->ChangeTimeStepTo(initial_time_step);

    std::cerr << "\n\n*** Starting time-stepping ***" << std::endl << std::endl;
    if ( no_std_output )
    {
        std::cerr << "-------------------------" << std::endl;
        std::cerr << " Stopping all std output " << std::endl;
        std::cerr << "-------------------------" << std::endl;

        // Reopen standard stream, and associate /dev/null
        freopen("/dev/null", "w", stdout);
    }


    // track duration of iteration
    typedef std::chrono::high_resolution_clock Clock;

    bool stop = false;
    while ( model_time <= max_time && !stop )
    {
        auto start_loop = std::chrono::high_resolution_clock::now();

        PrintProgressToScreen();

        //! Apply CVFEM scheme
        CVFEM_PHX->SetModelTime(model_time);
        dt = CVFEM_PHX->Apply();
        FluxPerSecond(dt);

        //! Apply Failure Mode Visitor
        if (hydrofracturing)
            model->Accept( *failure_mode );

        //! Permeability visitors and hydrofracturing
        if (depth_dependent || temperature_dependent || pore_fluid_factor_dependent || hydrofracturing)
        {
            model->CopyReplace("previous permeability", "previous previous permeability");
            model->CopyReplace("permeability", "previous permeability");
            model->InterpolateNodeToCellProperty("temperature", "temperature element"); // needed to calculate temperature-dependent pemermeability

            permeability_visitor->ComputeSimpleFracturableFlag(); // set fracturing reference to 1 in case T < solidus temperature

            model->Accept( *permeability_visitor );
            permeability_visitor->AssignPermeabilityTensor(*model);
        }


        //! Preparing next time step and output + save
        model_time += dt;
        model->InputPropertyValue("dt", ScalarVariable(ANY, dt));
        model->InputPropertyValue("model time", ScalarVariable(ANY, model_time));

        //! output vtu file
        CreateOutputFile(output_regions_list);

        // Abort simulation if timestep is to small
        if (dt < 1.e-4)
        {
            cerr << "Timestep smaller than 0.1 ms Aborting ..." << std::endl;
            stop = true;
        }

        //! Duration of iteration
        auto stop_loop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = (stop_loop - start_loop);
        std::cerr << fixed << setprecision(2) << "\nIteration took:\t" << time.count() << " s " << std::endl;
        cerr << "////////////////////////" << endl << endl;
    } // end main while loop


    std::cerr << "\nmain: That's it. Time for a coffee... \n\n" << std::endl;

} // end CVFEM_geothermal_cooling_magma_chamber_example::Run()



// Helper functions are defined below Run().

/// Create model

template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::CreateModel()
{
    //! Get setter boolean variables to define key model options
    open_top                        = GetBoolFromConfigFile("open top");
    with_magma_chamber              = GetBoolFromConfigFile("with magma chamber");

    // permeability options
    depth_dependent                 = GetBoolFromConfigFile("depth dependent permeability");
    temperature_dependent           = GetBoolFromConfigFile("temperature dependent permeability");
    pore_fluid_factor_dependent     = GetBoolFromConfigFile("pore fluid factor dependent permeability");
    hydrofracturing                 = GetBoolFromConfigFile("hydrofracturing");
    anisotropic_k                   = GetBoolFromConfigFile("with anisotropic permeability");
    T_dependent_differential_stress = GetBoolFromConfigFile("temperature dependent differential stress flag");

    //! Create model and reference to property database
    // Create PhysicalVariables file
    ModelBuilder<dim> builder;
    builder.CreateVariablesFile(false, false, false, false, false, false);
    builder.CreateLimiterLogFile();

    // Create Model
    model  = CreateANSYSModel(geometry_name_, geometry_name_, vars_name_);
    pd_ref = &(model->Database());
    model->InstantiateFiniteVolumes();

    //! Read config file
    InputDataManager<dim> model_configuration;
    model_configuration.ConfigureFromFile(*model, config_file_name_.c_str(),
                                          false,    ///< regionname from parameter range
                                          true,     ///< default property values
                                          true,     ///< regional property values
                                          false,    ///< boundary conditions for box-shaped model
                                          false,    ///< regional property conditions
                                          true,     ///< boundary conditions for arbitrary-shaped model
                                          run_settings );
} // end CreateModel()


/// Used to read 'with'-booleans (e.g. with magma model, with gold, ...) before initialising the model based on the config-file.
/// These booleans are used in the ModelBuilder (creates PhysicalVariables.txt) to only create variables actually used.
template<uint32_t dim>
bool CVFEM_geothermal_cooling_magma_chamber_example<dim>::GetBoolFromConfigFile(const std::string& var)
{
    std::ifstream file(std::string(config_file_name_) + "-configuration.txt");
    if (!file.is_open())
        throw std::runtime_error("Cannot open config file: " + std::string(config_file_name_));

    std::string line;
    while (std::getline(file, line))
    {
        // Skip comment and blank lines
        if (line.empty() || line[0] == '#')
            continue;

        // Strip inline comments
        std::size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);

        // Check the line starts with the variable name
        if (line.rfind(var, 0) != 0)  // must match at position 0
            continue;

        // Ensure the character after the name is whitespace (avoid partial matches)
        if (line.size() <= var.size() || !std::isspace(line[var.size()]))
            continue;

        // Tokenize the remainder by tabs/spaces to find the value
        std::istringstream iss(line.substr(var.size()));
        std::string token;
        while (iss >> token)  // >> skips all whitespace including tabs
        {
            try
            {
                int value = std::stoi(token);
                if (value == 0) return false;
                if (value == 1) return true;
                throw std::runtime_error("Expected 0 or 1 for bool variable: " + var);
            }
            catch (const std::invalid_argument&)
            {
                continue; // token wasn't a number, try the next
            }
        }
    }

    throw std::runtime_error("Bool variable not found in config: " + var);
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::SetTimeVariables()
{
    //! Time variables and objects
    //! Time variables defined in the config file are in unit time_unit defined below
    double& model_time  = ModelTime::Instance().modelTime;
    max_time            = run_settings.Duration();      // Max simulation time
    largest_time_step   = run_settings.TimeIncrement(); // Maximum timestep
    initial_time_step   = largest_time_step / 10.;      // Initial timestep
    output_increment    = run_settings.NearestOutputTime(0.0); // Output interval
    model_output        = output_increment;

    // Set time unites
    TimeUnit time_unit  = YEARS; // could be made a variable in config file
    if      (time_unit == YEARS)    time_multiplier = year;
    else if (time_unit == DAYS)     time_multiplier = day;
    else if (time_unit == HOURS)    time_multiplier = hour;
    else if (time_unit == MINUTES)  time_multiplier = minute;

    largest_time_step           *= time_multiplier;
    initial_time_step           *= time_multiplier;
    max_time                    *= time_multiplier;
    dt                           = initial_time_step;
    model_time                   = 0.;
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InstatiateCVFEM()
{
    //! CVFEM
    bool with_gravity            = model->Read(pd_ref->StorageKey("gravity flag"));
    bool with_custom_latent_heat = model->Read(pd_ref->StorageKey("with custom latent heat"));
    double cfl_scaling           = model->Read(pd_ref->StorageKey("cfl scaling"));
    double heat_capacity_rock    = model->Region("Model").NodeVector().front()->Read(pd_ref->StorageKey( "nodal heat capacity rock" ) );

    // Latent heat of fusion
    uint32_t crystallization_curve_         = model->Read(pd_ref->StorageKey("crystallization curve"));           // 0 = power law, 1 = error function, 2 = marxer ulmer
    double   crystallization_curve_exponent = model->Read(pd_ref->StorageKey("crystallization curve exponent"));  // If "power law" crystallization is used, b is the "power law" exponent
    double   sigma1                         = model->Read(pd_ref->StorageKey("sigma1 MagmaModel"));               // If "error function" crystallization is used, other possible values: 2/3, 1, 4/3
    double   heat_capacity_melt             = model->Read(pd_ref->StorageKey("heat capacity melt"));
    double   latent_heat_fusion             = model->Read(pd_ref->StorageKey("latent heat of fusion"));

    // Define crystallization curve
    std::string crystallization_curve;
    if      (crystallization_curve_ == 0)   crystallization_curve = "power law";
    else if (crystallization_curve_ == 1)   crystallization_curve = "error function";
    else if (crystallization_curve_ == 2)   crystallization_curve = "marxer ulmer";
    else
        throw csmp::Exception( FATAL_ERROR, "GeothermalSimulator<dim>::InstatiateCVFEM",
                              "Invalid crystallization curve: " + std::to_string(crystallization_curve_) +
                                  "\nAvailable options: 0 = power law; 1 = error function; 2 = marxer ulmer" );

    // Top boundary
    bool open_top       = model->Read( pd_ref->StorageKey( "open top" ) );
    double salinity_top = model->Read( pd_ref->StorageKey( "salinity top" ) );


    // Well geometry, completions, initial state and physics switches: one
    // [well] section per well, read from <mesh>-wells.txt before the scheme is
    // built so a bad file fails immediately. The values that used to be passed
    // to Initialize_well() (feedzone depths, initial temperature mode, initial
    // bulk steam mass fraction) live in this file now.
    // WellConfigurationFile::WriteTemplate(well_file) writes a commented
    // template with the defaults.
    // No wells used in this example.
    const string well_file("");
    const std::map<std::string, WellConfiguration> well_configs =
        wells_list.empty() ? std::map<std::string, WellConfiguration>()
                           : WellConfigurationFile::ReadAll(well_file);

    //! Instantiate CVFEM Scheme
    CVFEM_PHX = new CVFEM_PHX_Scheme<dim>(*model, with_gravity, true,
                                          wells_list,
                                          false, false,      // with_zinc, with_lithium
                                          solver_kind,       // reservoir pressure/temperature systems
                                          well_solver_kind,  // well Newton Jacobian
                                          well_configs);

    CVFEM_PHX->  SetEquilibratorConvergenceSpeedUpTo(true);
    CVFEM_PHX->  Adjust_CFL_Criterion(cfl_scaling, true); //true = pore velocity based; false = velocity based
    CVFEM_PHX->  WithRockLiquidusSolidus(true);
    CVFEM_PHX->  SetRockHeatCapacity(heat_capacity_rock);
    CVFEM_PHX->  AddFluidContributionToHeatCapacity(true);
    CVFEM_PHX->  AttemptToSurviveFluidPropertiesError(true);
    CVFEM_PHX->  SetLargestTimeStep(largest_time_step);

    // latent heat
    if (!with_custom_latent_heat)
        CVFEM_PHX->SetRockCrystallizationCurve(1.78, sigma1, 300000., crystallization_curve_exponent, crystallization_curve);
    else
        CVFEM_PHX->SetRockCrystallizationCurve(heat_capacity_melt / heat_capacity_rock, sigma1, latent_heat_fusion, crystallization_curve_exponent, crystallization_curve);

    // well model
    CVFEM_PHX->  IncludeWellCalculations(false);

    // top boundary
    if (open_top)
        CVFEM_PHX->OpenBoundaries(salinity_top);
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InstantiatePermeabilityVisitor()
{
    ///< Permeability visitor
    permeability_visitor = new PermeabilityVisitor<dim>(*model);

    //! Background permeability calculation based on defined criterions
    //! e.g., depth-dependent, temperature-dependent, etc.
    PropertyHandle<dim> prev_k(*model, "previous permeability", SCALAR, ELEMENT );
    PropertyHandle<dim> prev_prev_k(*model, "previous previous permeability", SCALAR, ELEMENT );

    // read variables defined in config file
    bool    with_reference_depth = model->Read( pd_ref->StorageKey( "with reference depth" ) );                  // Option to use a mesh where the surface is not at y=0
    double  reference_depth      = model->Read( pd_ref->StorageKey( "reference depth" ) );                       // y-coordinate of the surface

    bool    change_brittle_ductile_transtion_T = model->Read( pd_ref->StorageKey( "change BDT temperature" ) ); // update BDT temperature

    double  brittle_ductile_transition_start   = model->Read( pd_ref->StorageKey( "BDT start temperature" ) );     // Start to become ductile at this temperature
    double  brittle_ductile_transition_ductile = model->Read( pd_ref->StorageKey( "BDT ductile temperature" ) );   // Ductile at this temperature
    double  brittle_ductile_transition_end     = model->Read( pd_ref->StorageKey( "BDT end temperature" ) );       // T at which rock is impermeable

    double  log_k_start     = model->Read( pd_ref->StorageKey( "BDT start log permeability" ) );       // permeability at T_bdt_s
    double  log_k_duct      = model->Read( pd_ref->StorageKey( "BDT ductile log permeability" ) );     // permeability at T_bdt_ductile
    double  log_k_end       = model->Read( pd_ref->StorageKey( "BDT end log permeability" ) );         // permeability at T_bdt_e

    // region ID of intrusion.
    // Define here to avoid applying CalculateGradPScaling() to the cooled intrusion.
    // Makes permeability within intrusion static after T < brittle_ductile_transition_start.
    uint32_t intrusion_ID = 1;

    // Set permeability options
    permeability_visitor->DepthDependent( depth_dependent );
    permeability_visitor->SetReferenceDepthTo( with_reference_depth, reference_depth );

    permeability_visitor->AnisotropicPermeabilityTensor( anisotropic_k );
    permeability_visitor->TemperatureDependent( temperature_dependent );
    permeability_visitor->PoreFluidFactorDependent(pore_fluid_factor_dependent, false); // true = average nodal Pff; false = maximum nodal Pff
    permeability_visitor->ComputeSimpleFracturableFlag(); // set fracturing reference to 1 in case T < solidus temperature

    permeability_visitor->ChangeBrittleDuctileTransitionTemperature(change_brittle_ductile_transtion_T, // bool
                                                                    brittle_ductile_transition_start,   // default: 360 degC
                                                                    brittle_ductile_transition_ductile, // default: 400 degC
                                                                    brittle_ductile_transition_end,     // default: 500 degC
                                                                    log_k_start,                        // default: -15 m2
                                                                    log_k_duct,                         // default: -17 m2
                                                                    log_k_end);                         // default: -22 m2
    permeability_visitor->SetIntrusionRegionID( intrusion_ID );

    model->InterpolateNodeToCellProperty("temperature", "temperature element"); // used for temperature-dependent permeability

} // end InstantiatePermeabilityPorosityVisitors


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InstantiateFailureModeVisitor()
{
    // Read data from config file
    bool    near_critically_pressured           = model->Read( pd_ref->StorageKey( "near critically pressured" ) );
    bool    change_brittle_ductile_transtion_T  = model->Read( pd_ref->StorageKey( "change BDT temperature" ) );
    double  C                                   = model->Read( pd_ref->StorageKey( "cohesion" ) );              // Rock cohesive strength
    double  diff_stress                         = model->Read( pd_ref->StorageKey( "differential stress" ) );   // used in case NOT "near critically pressured"!!! otherwise calculated in FailureModeVisitor
    double  brittle_ductile_transition_start    = model->Read( pd_ref->StorageKey( "BDT start temperature" ) ); // Start to become ductile at this temperature
    double  brittle_ductile_transition_end      = model->Read( pd_ref->StorageKey( "BDT end temperature" ) );   // T at which rock is impermeable

    // Instantiate Visitor and define settings
    failure_mode = new FailureModeVisitor<dim>(*model, C, diff_stress, near_critically_pressured);

    failure_mode->SetTemperatureRelaxationTo(T_dependent_differential_stress);

    if (change_brittle_ductile_transtion_T)
        failure_mode->ChangeBrittleDuctileTransitionTemperature( brittle_ductile_transition_start, brittle_ductile_transition_end );
}

template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::ActivateFracturing(bool hydrofracturing)
{
    bool    immediate_closure   = model->Read( pd_ref->StorageKey( "immediate fracture closure" ) );        // Immediate fracture closure if fluid pressure decreases
    double  log_max_perm        = model->Read( pd_ref->StorageKey( "maximum log permeability model" ) );    // Maximum allowed permeability
    double  log_min_perm        = model->Read( pd_ref->StorageKey( "BDT end log permeability" ) );
    permeability_visitor->Hydrofracturing(hydrofracturing, log_max_perm, log_min_perm);
    permeability_visitor->SetImmediateClosureTo(immediate_closure);
    if (hydrofracturing)
        model->Accept(*failure_mode);
}

template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializeTemperature()
{
    SolverBundle init_bundle(solver_kind);
#if defined(CSMP_WITH_PETSC_SOLVER)
    if (solver_kind == SolverKind::PETSc) {
        init_bundle.PETScSettings().SetKSPType("gmres");
        init_bundle.PETScSettings().SetPCType("gamg"); //or hypre possible too
        init_bundle.PETScSettings().SetPrintConvergedReason(true);
    }
#endif

    //! Read variables
    PDE_Integrator<dim>                 initial_temperature(init_bundle.Get());
    CVFEM_NumIntegral_dNT_op_dN_dV<dim> conductance(*pd_ref, "thermal conductivity", "temperature", "temperature", "thickness");
    CVFEM_PointSource_rhsop<dim>        heat_bottom(*pd_ref, "nodal heat flux bottom", "temperature");

    initial_temperature.Add(&conductance);
    initial_temperature.Add(&heat_bottom);

    //! Calculate initial temperature
    model->Apply(initial_temperature);
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::AddMagmaticIntrusion()
{
    SolverBundle init_bundle(solver_kind);
#if defined(CSMP_WITH_PETSC_SOLVER)
    if (solver_kind == SolverKind::PETSc) {
        init_bundle.PETScSettings().SetKSPType("gmres");
        init_bundle.PETScSettings().SetPCType("gamg"); //or hypre possible too
        init_bundle.PETScSettings().SetPrintConvergedReason(true);
    }
#endif


    //! Read variables
    double  starting_temperature        = model->Read( model->Database().StorageKey( "initial intrusion temperature" ) );
    bool    with_temperature_halo       = model->Read( pd_ref->StorageKey( "with temperature halo" ) );
    double  initial_T_diffusion_time    = model->Read( pd_ref->StorageKey( "initial T diffusion time" ) );
    initial_T_diffusion_time            *= time_multiplier;

    PDE_Integrator<dim>                 temperature_diffusion(init_bundle.Get());

    NumIntegral_dNT_op_dN_dV<dim>       t_conductance(*pd_ref,"thermal conductivity","temperature","temperature" );
    NumIntegral_NT_lhsop_N_dV<dim>      t_capacitance_lhs(*pd_ref,"nodal heat capacity","temperature","temperature" );
    NumIntegral_NT_rhsop_N_dV<dim>      t_capacitance_rhs(*pd_ref,"nodal heat capacity","temperature" );
    CVFEM_PointSource_rhsop<dim>        t_heat_bottom(*pd_ref,"nodal heat flux bottom","temperature");

    t_capacitance_lhs.LumpedFormulation(true);
    t_capacitance_rhs.LumpedFormulation(true);
    t_conductance.MultiplyWithTimeIncrement(true);
    t_heat_bottom.AddAccumulateLater();

    temperature_diffusion.Add(&t_conductance);
    temperature_diffusion.Add(&t_capacitance_lhs);
    temperature_diffusion.Add(&t_capacitance_rhs);
    temperature_diffusion.Add(&t_heat_bottom);
    temperature_diffusion.TimeIncrement(initial_T_diffusion_time);
    temperature_diffusion.Verbose(false);


    //! Add magmatic intrusion
    if (with_magma_chamber && model->ContainsRegion("CHAMBER"))
    {
        model->Region("CHAMBER").InputPropertyValue("temperature", ScalarVariable(ANY, starting_temperature));

        //! Calculate initial temperature diffusion to mimic heat diffusion during magma emplacement
        if (with_temperature_halo)
        {
            model->Region("CHAMBER").ChangePropertyStatus("temperature", DIRICH);

            temperature_diffusion.IntegrateOver(model->Region("Model"));

            model->Region("CHAMBER").ChangePropertyStatus("temperature", ANY);
            model->Region("CHAMBER").InputPropertyValue("temperature", ScalarVariable(ANY, starting_temperature));

            std::cerr<<endl<<"Initial temperature range after halo: " << printRangeOfVariable(*model, "temperature") << std::endl;
        }
    }

    model->InterpolateNodeToCellProperty("temperature", "temperature element");

    double  brittle_ductile_transition_end = model->Read( pd_ref->StorageKey( "BDT end temperature" ) );            // T at which rock is impermeable
    model->FormRegionFrom("high temperature region", "temperature element", brittle_ductile_transition_end, 1100.); // used to assign lithostatic pressure
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializeHydrostaticPressure()
{
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

    PDE_Integrator<dim>             initial_pressure_(init_bundle.Get());
    NumIntegral_dNT_op_dN_dV<dim>   conductance_hp(*pd_ref,"permeability","fluid pressure","fluid pressure");
    NumIntegral_NT_op_dNi_dV<dim>   density_hp(*pd_ref,"fluid density","permeability","fluid pressure");

    // Hydrostatic pressure initialization is done further (in equilibration loop)
    initial_pressure_.Add(&conductance_hp);
    initial_pressure_.Add(&density_hp);
    model->Apply(initial_pressure_);

    //! Calculate reference pressure
    //! Temperature is only defined by bottom heat flow (without intrusion)
    int k = 0;
    bool flag(false);
    double diff(0.), eps(1000.0);

    while (!flag) {
        cerr << endl << "Initialization P-T, iteration: " << k;
        model->CopyReplace("fluid pressure", "test fluid pressure");
        CVFEM_PHX->InitialFluidPropertiesFromPTX();

        model->Apply(initial_pressure_);

        diff = maxDifferenceScalarNodeProperty("fluid pressure", "test fluid pressure");
        cerr << endl << "*********************************Max diff pressure: " << diff;
        if (diff < eps) flag = true;
        k++;
    }

    // Store reference pressure and reset temperature to T with intrusion before calculating hydrostatic pressure
    model->CopyReplace("fluid pressure", "reference pressure");

} // end InitializeHydrostaticPressure



template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializePressure()
{
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



    //! Calculate the lithostatic gradient assuming homogeneous rock density (carefull, it doesn't consider the "bulk" density)
    PDE_Integrator<dim>             lithos_pressure(init_bundle.Get());
    NumIntegral_dNT_op_dN_dV<dim>   conductance_p(*pd_ref,"permeability","lithostatic pressure","lithostatic pressure");
    NumIntegral_NT_op_dNi_dV<dim>   density_p(*pd_ref,"nodal density rock","permeability","lithostatic pressure");

    lithos_pressure.Add(&conductance_p);
    lithos_pressure.Add(&density_p);
    model->Apply(lithos_pressure);

    //! Hydrostatic pressure
    PDE_Integrator<dim>             initial_pressure_(init_bundle.Get());
    NumIntegral_dNT_op_dN_dV<dim>   conductance_hp(*pd_ref,"permeability","fluid pressure","fluid pressure");
    NumIntegral_NT_op_dNi_dV<dim>   density_hp(*pd_ref,"fluid density","permeability","fluid pressure");
    initial_pressure_.Add(&conductance_hp);
    initial_pressure_.Add(&density_hp);
    model->Apply(initial_pressure_);

    //! Calculate hydrostatic pressure
    int k = 0;
    bool flag(false);
    double diff(0.), eps(1000.0);

    while (!flag) {
        cerr << endl << "Initialization P-T, iteration: " << k;
        model->CopyReplace("fluid pressure", "test fluid pressure");
        CVFEM_PHX->InitialFluidPropertiesFromPTX();

        model->Apply(initial_pressure_);

        diff = maxDifferenceScalarNodeProperty("fluid pressure", "test fluid pressure");
        cerr << endl << "*********************************Max diff pressure: " << diff;
        if (diff < eps) flag = true;
        k++;
    }


    //! Set lithostatic pressure at T>= solidus temperature
    bool with_lithostatic_pressure = model->Read( pd_ref->StorageKey( "with lithostatic pressure" ) );
    if (with_lithostatic_pressure) {
        flag = false; k = 0; diff = 0.; eps = 1000.0;

        while (!flag) {
            cerr << endl << "Initialization P-T with lithostatic ductile, iteration: " << k;
            model->CopyReplace("fluid pressure", "test fluid pressure");
            CVFEM_PHX->InitialFluidPropertiesFromPTX();

            model->Apply(initial_pressure_);

            if (model->ContainsRegion("high temperature region"))
                model->Region("high temperature region").CopyReplace("lithostatic pressure", "fluid pressure");

            diff = maxDifferenceScalarNodeProperty("fluid pressure", "test fluid pressure");
            cerr << endl << "*********************************Max diff pressure: " << diff;
            if (diff < eps) flag = true;
            k++;
        }
    }

    // remove temporary regions
    if (model->ContainsRegion("high temperature region"))
        model->RemoveRegion("high temperature region", false);


    //! Calculate failure pressure
    if (hydrofracturing)
        model->Accept(*failure_mode);

} // end InitializePressure



template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializeTotalHeatCapacity()
{
    csmp::Index phi_key( model->Database().StorageKey("nodal porosity") );
    csmp::Index rho_key( model->Database().StorageKey("nodal density rock") );
    csmp::Index Cpr_key( model->Database().StorageKey("nodal heat capacity rock") );
    csmp::Index Cp_key ( model->Database().StorageKey("nodal heat capacity") );

    ScalarVariable phi, rho, Cpr, Cp;

    typename vector<Node<dim>*>::const_iterator nit;
    for ( nit = model->Region("Model").NodesBegin(); nit != model->Region("Model").NodesEnd(); ++nit )
    {
        // read vars
        (*nit)->Read(phi_key, phi);
        (*nit)->Read(rho_key, rho);
        (*nit)->Read(Cpr_key, Cpr);
        Cp() = 0;

        // calculate total heat capacity
        Cp() = Cpr() * rho() * (1-phi());

        // store nodal value
        (*nit)->Store(Cp_key, Cp);
    }
    model->InterpolateNodeToCellProperty("nodal heat capacity", "total heat capacity");
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializeAllNaNsToZero()
{

    double min(0.), max(0.);

    std::map<std::string, csmp::Parameter>::const_iterator it;

    std::cerr<<"Initializing NaNs to Zero!" << std::endl << std::endl;
    for (it = model->Database().Begin(); it != model->Database().End(); it++)
    {
        model->MinMaxOf( (*it).first.c_str(), min, max );

        if (std::isnan(min) || std::isnan(max))
        {
            std::cerr<<"Setting to Zero: " << (*it).first.c_str() <<  std::endl;

            if ((*it).second.key.type == SCALAR)
                model->InputPropertyValue((*it).first.c_str(), ScalarVariable(ANY, 0.0));
            else if ((*it).second.key.type == VECTOR)
                model->InputPropertyValue((*it).first.c_str(), VectorVariable<dim>(ANY, 0.0));
            else if ((*it).second.key.type == TENSOR)
                model->InputPropertyValue((*it).first.c_str(), TensorVariable<dim>(ANY, 0.0));
        }
    }
    std::cerr<<"Initializing NaNs to Zero: Done!" << std::endl << std::endl;
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::InitializeParameters ()
{
    double    largest_time_step    = run_settings.TimeIncrement() * YEARS; // largest timestep in seconds

    // Input
    // coordinates
    model->Region("Model"). AssignNodeCoordinatesTo( "nodal depth", 'y' );
    model->Region("Model"). AssignNodeCoordinatesTo( "coordinate x", 'x' );
    model->Region("Model"). AssignNodeCoordinatesTo( "coordinate y", 'y' );

    if (dim == 3U)
        model->Region("Model"). AssignNodeCoordinatesTo( "coordinate z", 'z' );

    model->InterpolateNodeToCellProperty("nodal depth", "barycenter depth");

    InitializeTotalHeatCapacity(); // needs correct nodal porosity for initial temperature diffusion
    model-> InputPropertyValue( "thickness", ScalarVariable(ANY, 1.0));

    model-> InputPropertyValue( "courant liquid", ScalarVariable(ANY, largest_time_step));
    model-> InputPropertyValue( "courant vapor", ScalarVariable(ANY, largest_time_step));

    model-> InputPropertyValue( "liquid mass mobility density", ScalarVariable(ANY, 1.e-10));
    model-> InputPropertyValue( "vapor mass mobility density", ScalarVariable(ANY, 1.e-10));
    model-> InputPropertyValue( "density liquid transport", ScalarVariable(ANY, 1.0));
    model-> InputPropertyValue( "density vapor transport", ScalarVariable(ANY, 1.0));

    model-> InputPropertyValue( "liquid enthalpy mobility", ScalarVariable(ANY, 1.e-10));
    model-> InputPropertyValue( "vapor enthalpy mobility", ScalarVariable(ANY, 1.e-10));
    model-> InputPropertyValue( "liquid enthalpy mobility density", ScalarVariable(ANY, 1.e-10));
    model-> InputPropertyValue( "vapor enthalpy mobility density", ScalarVariable(ANY, 1.e-10));

    model-> InputPropertyValue( "density liquid", ScalarVariable(ANY, 1.));
    model-> InputPropertyValue( "density vapor", ScalarVariable(ANY, 1.));
    model-> InputPropertyValue( "viscosity liquid", ScalarVariable(ANY, 1.e-4));
    model-> InputPropertyValue( "viscosity vapor", ScalarVariable(ANY, 1.e-5));
    model-> InputPropertyValue( "enthalpy liquid", ScalarVariable(ANY, 100.));
    model-> InputPropertyValue( "enthalpy vapor", ScalarVariable(ANY, 100.));
    model-> InputPropertyValue( "fluid heat capacity", ScalarVariable(ANY, 1000.));
    model-> InputPropertyValue( "bulk fluid density", ScalarVariable(ANY, 1000.));

    model-> InputPropertyValue( "relperm viscosity liquid", ScalarVariable(ANY, 1.e-5));
    model-> InputPropertyValue( "relperm viscosity vapor", ScalarVariable(ANY, 1.e-5));

    model-> InputPropertyValue( "density halite", ScalarVariable(ANY, 1700.));
    model-> InputPropertyValue( "enthalpy halite", ScalarVariable(ANY, 200.));

    model-> InputPropertyValue( "previous compressibility", ScalarVariable(ANY, 0.0001));
    model-> InputPropertyValue( "reference compressibility", ScalarVariable(ANY, 0.0001));

    model-> InputPropertyValue( "conductivity", ScalarVariable(ANY, 1.));
    model-> InputPropertyValue( "mass conductivity", ScalarVariable(ANY, 1.));
    model-> InputPropertyValue( "fluid viscosity", ScalarVariable(ANY, 1.e-4));

    model-> InputPropertyValue( "fracturing reference", ScalarVariable(ANY, 30.0));

    // Permeability
    model-> CopyReplace("permeability","BDT start permeability element");        // Define BDT start peremability as initial rock permeability
    model-> ExtrapolateCellToNodeProperty("permeability","nodal permeability");  // Could cause conflicts in case of anisotropic permeability ?

    //! Set all variables that are NaNs to zero
    InitializeAllNaNsToZero();
}


template<uint32_t dim>
double CVFEM_geothermal_cooling_magma_chamber_example<dim>::maxDifferenceScalarNodeProperty(const char *snp1, const char *snp2)
{
    csmp::Index snp1Key( model->Database().StorageKey(snp1) );
    csmp::Index snp2Key( model->Database().StorageKey(snp2) );

    double diff(0.), max(0.);

    typename vector<Node<dim>*>::const_iterator nit;

    for ( nit = model->Region("Model").NodesBegin(); nit != model->Region("Model").NodesEnd(); ++nit )
    {
        diff = fabs( (*nit)->Read(snp1Key) - (*nit)->Read(snp2Key) );

        if (diff > max)
            max = diff;
    }

    return max;
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::ApplyBottomHeatFlux ()
{
    csmp::Index       hfb_key (model->Database().StorageKey("nodal heat flux bottom"));
    ScalarVariable    heat_flux, hfb;
    double area;
    uint32_t nnodes;

    const Boundary<dim>   &bref(model->Boundary("BOTTOM"));
    typename vector<Node< dim>*>::const_iterator nit( bref.NodesBegin() );
    heat_flux = (*nit)->Read(hfb_key);

    model->InputPropertyValue( "nodal heat flux bottom", ScalarVariable(ANY, 0.) );

    if (heat_flux() != 0.)
    {
        const typename vector<Face< dim>*>::const_iterator boundaryCellsEnd(
            bref.CellsEnd() );

        for ( typename vector<Face<dim>*>::const_iterator it( bref.CellsBegin() );
             it != boundaryCellsEnd; ++it )
        {
            nnodes = (*it)->Nodes();
            area = (*it)->Area();

            for (uint32_t i = 0U; i < nnodes; i++)
            {
                hfb()  = (*it)->N(i)->Read(hfb_key);
                hfb() += heat_flux() * area / nnodes;
                (*it)->N(i)->Store(hfb_key, hfb);
            }
        }
    }
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>:: FluxPerSecond (double dt)
{
    const Region<dim> &mref = model->Region("Model");
    csmp::Index ff_key( model->Database().StorageKey("fluid flux") );
    csmp::Index ffps_key( model->Database().StorageKey("fluid flux per second") );
    csmp::Index ef_key( model->Database().StorageKey("energy flux") );
    csmp::Index efps_key( model->Database().StorageKey("energy flux per second") );

    ScalarVariable ff (ANY, 0.0);
    ScalarVariable ffps (ANY, 0.0);
    ScalarVariable ef (ANY, 0.0);
    ScalarVariable efps (ANY, 0.0);

    typename vector<Node<dim>*>::const_iterator nit;

    for ( nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit )
    {
        ffps = (*nit)->Read(ff_key) / dt;
        (*nit)->Store(ffps_key, ffps);

        efps = (*nit)->Read(ef_key) / dt;
        (*nit)->Store(efps_key, efps);
    }
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::PrintProgressToScreen()
{
    double& model_time  = ModelTime::Instance().modelTime;

    std::cerr << std::endl << std::endl << "////////////////////////";
    std::cerr << std::endl << "////////////////////////";
    std::cerr << std::endl << output_name_.c_str() << std::endl;
    std::cerr << "\nModel time:\t" << model_time / year << "  years" << std::endl;

    if (dt/day >= 1.)  std::cerr << fixed << setprecision(2) << "Timestep:\t" << dt / day << "  days" << std::endl;
    else               std::cerr << fixed << setprecision(2) << "Timestep:\t" << dt / hour << "  hours" << std::endl;

    std::cerr << fixed << setprecision(2) << "Progress:\t" << model_time / max_time * 100 << "  %" << std::endl << std::endl;

    std::cerr << defaultfloat;
    std::cerr.precision(4);
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::DefineOutputVariables ()
{
    // 1. THERMODYNAMIC / FLUID STATE
    // -------------------------------
    output_variables.push_back("fluid state");
    output_variables.push_back("temperature");
    output_variables.push_back("fluid pressure");
    output_variables.push_back("reference pressure");
    output_variables.push_back("lithostatic pressure");
    output_variables.push_back("density liquid");
    output_variables.push_back("density vapor");
    output_variables.push_back("fluid density");
    output_variables.push_back("bulk fluid density");
    output_variables.push_back("saturation liquid");
    output_variables.push_back("saturation vapor");
    output_variables.push_back("viscosity liquid");
    output_variables.push_back("viscosity vapor");

    // 2. ENTHALPY / HEAT CONTENT
    // ---------------------------
    output_variables.push_back("fluid enthalpy");
    output_variables.push_back("enthalpy liquid");
    output_variables.push_back("enthalpy vapor");
    output_variables.push_back("enthalpy content liquid");
    output_variables.push_back("enthalpy content vapor");
    output_variables.push_back("previous total enthalpy");
    output_variables.push_back("specific enthalpy liquid top");

    // 3. THERMAL PROPERTIES
    // ----------------------
    output_variables.push_back("thermal conductivity");
    output_variables.push_back("nodal heat capacity");
    output_variables.push_back("nodal heat capacity rock");
    output_variables.push_back("fluid heat capacity");
    output_variables.push_back("nodal heat flux bottom");

    // 4. ROCK / MEDIUM PROPERTIES
    // -----------------------------
    output_variables.push_back("porosity");
    output_variables.push_back("nodal porosity");
    output_variables.push_back("permeability tensor");
    output_variables.push_back("permeability");
    output_variables.push_back("permeability ID");
    output_variables.push_back("BDT start permeability element");
    output_variables.push_back("region ID");

    // 5. FLOW & TRANSPORT
    // ---------------------
    output_variables.push_back("velocity");
    output_variables.push_back("velocity liquid");
    output_variables.push_back("velocity vapor");
    output_variables.push_back("boundary flow mass");
    output_variables.push_back("energy flux");
    output_variables.push_back("fluid flux");
    output_variables.push_back("KgradP");

    // 6. MASS & VOLUME
    // ------------------
    output_variables.push_back("fluid mass liquid");
    output_variables.push_back("fluid mass vapor");
    output_variables.push_back("pore volume");
    output_variables.push_back("bulk volume");

    // 7. CHEMISTRY / COMPOSITION
    // -----------------------------
    output_variables.push_back("salinity");

    // 8. SIMULATION CONTROL
    // ------------------------
    output_variables.push_back("dt");


    // optional output variables depending on model setup
    if (anisotropic_k)
        output_variables.push_back("permeability anisotropy factor");
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::CreateOutputFile(
    std::vector<std::string> output_regions_list)
{
    //! Create vtu output file
    double& model_time  = ModelTime::Instance().modelTime;
    long time_label = model_time / time_multiplier;
    if ( time_label >= model_output )
    {
        OutputToVTU(output_path_ + output_name_, output_variables, time_label, output_regions_list);
        model_output += output_increment;
        std::cerr << "vtu file created for year: " << time_label << std::endl;
    }
}


template<uint32_t dim>
void CVFEM_geothermal_cooling_magma_chamber_example<dim>::OutputToVTU(string model_name,
                                                                      const list<string> &props,
                                                                      size_t timestep,
                                                                      std::vector<std::string> additional_region) const
{
    static VTU_Interface<dim> vtu(*model);
    vtu.OmitZeroInFileName(false);
    vtu.OutputDataToVTU( (model_name).c_str(), props, model->Region("Model"), timestep);

    // output additional vtu files for individual regions
    if (!additional_region.empty())
    {
        for (const auto &region : additional_region)
            if (model->ContainsRegion(region))
                vtu.OutputDataToVTU( (model_name).c_str(), props, model->Region(region), timestep);
    }
}


template class CVFEM_geothermal_cooling_magma_chamber_example<2U>;
template class CVFEM_geothermal_cooling_magma_chamber_example<3U>;
} // end namespace csmp