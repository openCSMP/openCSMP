#include "CVFEM_2D_VVCase.h"
#include <sys/stat.h>


using namespace std;



/** 2D Benchmarks by Jonas Köpping (based on BLC's BasicTask6)
   This code is experimental...

   Config File:  test-201a, test-201b, test-202a, test-202b, test-203a, test-203b
   Regions File: rectangle_coarse
   Var File:     Created by ModelBuilder

**/
namespace csmp
{

// config_file_name_ (passed in as configName below) is already the full path to this
// test case's own benchmark subfolder — see getTestDataPath_2D() in CVFEM_tests_main.cpp,
// which returns e.g. ".../2D_Benchmarks/test-202a/". No separate root/directory lookup
// is needed here; this just appends the file suffix onto that already-complete path.
// e.g. benchmarkFilePath(".../2D_Benchmarks/test-202a/", "_Initial_Properties")
//      -> ".../2D_Benchmarks/test-202a/_Initial_Properties"
namespace {
std::string benchmarkFilePath(const std::string& configName, const std::string& suffix)
{
    return configName + suffix;
}
} // anonymous namespace


// Extracts the short test name from a full benchmark path, e.g.
//   ".../2D_Benchmarks/test-202a/" -> "test-202a"
namespace {
std::string extractTestName(const std::string& path)
{
    std::string trimmed = path;
    while (!trimmed.empty() && trimmed.back() == '/')
        trimmed.pop_back();

    const auto slash = trimmed.find_last_of('/');
    if (slash != std::string::npos)
        return trimmed.substr(slash + 1);
    return trimmed;
}
}// anonymous namespace


// input arguments: mesh, PhysicalVariables, config-file, output-path, output-name
CVFEM_2D_VVCase::CVFEM_2D_VVCase(int argc, char **argv, const uint32_t test_type)
    : geometry_name_(argv[1]),
    config_file_name_(argv[2]),
    output_name_(argv[2]),
    test_type(test_type),
    vars_name_("PhysicalVariables.txt")
{
    //! Create PhysicalVariables file
    builder.CreateVariablesFile(false, false, false, false, false, false);
    builder.CreateLimiterLogFile();

    //! Create model and reference to property database
    model = new ANSYS_Model2D(geometry_name_, geometry_name_, vars_name_);
    model->InstantiateFiniteVolumes();
    pd_ref = &(model->Database());

    //! Read config file
    InputDataManager<dim>           model_configuration;
    model_configuration.ConfigureFromFile(*model, config_file_name_,
                                          false,   // groupname from parameter range
                                          true,    // default property values
                                          true,    // regional property values
                                          true,    // boundary conditions for box-shaped model //JK: Set to true?!?
                                          true,    // essential conditions for groups
                                          true,    // csmp::Boundary properties
                                          run_settings );

    //! assign the gravity flag from the config file
    csmp::Index   g_key(model->Database().StorageKey("gravity flag"));
    with_gravity = bool(model->Read(g_key));

    //! create new CVFEM PHX Scheme object
    // ── LINEAR SOLVER SELECTION ─────────────────────────────────────────────
    // Backends in CVFEM_SolverChoice.h:
    //   SAMG    algebraic multigrid, licensed; tuned in the scheme constructor
    //   PETSc   Krylov + preconditioner, GMRES+ILU by default
    //   Eigen   direct SparseLU; rejected for the reservoir above ~10k DOF
    //
    // RESERVOIR — used by the initialisation solves below and by the transient
    // scheme, so the whole reservoir run shares one backend. SAMG or PETSc only.
    const SolverKind solver_kind = SolverKind::PETSc;

    // One bundle for the steady initialisation solves. These integrators use no
    // split-boundary operators, so they stay on the default CompressedRowMatrix
    // (unlike the scheme, which is pinned to SparseMatrix).
    SolverBundle init_bundle (solver_kind);

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
    // Well geometry, completions, initial state and physics switches: one
    // [well] section per well, read from <mesh>-wells.txt before the scheme is
    // built so a bad file fails immediately. The values that used to be passed
    // to Initialize_well() (feedzone depths, initial temperature mode, initial
    // bulk steam mass fraction) live in this file now.
    // WellConfigurationFile::WriteTemplate(well_file) writes a commented
    // template with the defaults.
    std::vector<std::string> wells_list = {};
    const string well_file("");
    const std::map<std::string, WellConfiguration> well_configs =
        wells_list.empty() ? std::map<std::string, WellConfiguration>()
                           : WellConfigurationFile::ReadAll(well_file);

    //! Instantiate CVFEM Scheme
    CVFEM_PHX = new CVFEM_PHX_Scheme<dim>(*model, with_gravity, true,
                                                                 wells_list,
                                                                 false, false,      // with_zinc, with_lithium
                                                                 solver_kind,       // reservoir pressure/temperature systems
                                                                 solver_kind,       // well Newton Jacobian
                                                                 well_configs);

    CVFEM_PHX->ActivatePointInjectionBenchmarks(extractTestName(config_file_name_).c_str()); // Activate point source for specific 2D Benchmarks

    //! create new permeability visitor and assign permeability tensor
    permeability_visitor = new PermeabilityVisitor<dim>(*model);
    permeability_visitor->AssignPermeabilityTensor(*model);

    //! calculate pore volume
    PoreVolumeVisitor<dim> Pore_Volume_Visitor(*model, "nodal porosity", "bulk volume", "pore volume", "thickness");
    model->Accept(Pore_Volume_Visitor);

    //! vtu interface fro output
    vtu = new VTU_Interface<dim> (*model);
    vtu->OmitZeroInFileName( false );

    //! output main variables
    output_props_.push_back("boundary flow mass");
    output_props_.push_back("nodal porosity");
    output_props_.push_back("pore volume");
    output_props_.push_back("porosity");
    output_props_.push_back("fluid pressure");
    output_props_.push_back("velocity");
    output_props_.push_back("velocity liquid");
    output_props_.push_back("velocity vapor");
    output_props_.push_back("temperature");
    output_props_.push_back("salinity");
    output_props_.push_back("fluid state");
    output_props_.push_back("fluid density");
    output_props_.push_back("fluid enthalpy");
    output_props_.push_back("saturation liquid");
    output_props_.push_back("saturation vapor");
    output_props_.push_back("saturation halite");
    output_props_.push_back("fluid mass liquid");
    output_props_.push_back("fluid mass vapor");
    output_props_.push_back("permeability tensor");
    output_props_.push_back("nodal heat flux bottom");
}

CVFEM_2D_VVCase::~CVFEM_2D_VVCase()
{
    //! clean memory
    delete CVFEM_PHX;

    if (permeability_visitor != NULL)
        delete permeability_visitor;

    delete model;
    delete vtu;
}


double CVFEM_2D_VVCase::maxDifferenceScalarNodeProperty( Model<dim>& mdl, const char* snp1, const char* snp2 )
{
    csmp::Index snp1Key( mdl.Database().StorageKey(snp1) );
    csmp::Index snp2Key( mdl.Database().StorageKey(snp2) );
    const Region<dim>&   mref = mdl.Region("Model");

    double diff(0.), max(0.);

    typename vector<Node<dim>*>::const_iterator nit;
    for( nit = mref.NodesBegin(); nit != mref.NodesEnd(); ++nit )
    {
        diff = fabs( (*nit)->Read(snp1Key) - (*nit)->Read(snp2Key) );

        if (diff > max)
            max = diff;
    }
    return max;
}

void CVFEM_2D_VVCase::InitializeHydrostaticPressure()
{
    const SolverKind solver_kind = SolverKind::PETSc;
    const SolverKind well_solver_kind = SolverKind::PETSc;
    const bool pin_split_middle_during_init = true;
    SolverBundle init_bundle (solver_kind);

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


    PDE_Integrator<dim>  initial_pressure_(init_bundle.Get());

    NumIntegral_dNT_op_dN_dV<dim> conductance_hp(*pd_ref, "permeability", "fluid pressure", "fluid pressure");
    NumIntegral_NT_op_dNi_dV<dim> density_hp(*pd_ref, "fluid density", "permeability", "fluid pressure");

    initial_pressure_.Add(&conductance_hp);
    initial_pressure_.Add(&density_hp);
    model->Apply(initial_pressure_);

    int k = 0; bool flag(false); double diff(0.), eps(1000.0);
    while (!flag)
    {
        cerr << endl << "Initialization P-T, iteration: " << k;
        model->CopyReplace("fluid pressure", "test fluid pressure");
        CVFEM_PHX->InitialFluidPropertiesFromPTX();

        model->Apply(initial_pressure_);

        diff = maxDifferenceScalarNodeProperty(*model, "fluid pressure", "test fluid pressure");
        cerr << endl << "*********************************Max diff pressure: " << diff;
        if (diff < eps) flag = true;
        k++;
    }
}

void CVFEM_2D_VVCase::InitializeAllNaNsToZero()
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

void CVFEM_2D_VVCase::ApplyBottomHeatFlux()
{
    csmp::Index hfb_key (pd_ref->StorageKey("nodal heat flux bottom"));
    ScalarVariable    heat_flux, heat_flux1, hfb;
    double area;
    uint32_t nnodes;


    const Boundary<dim>   &bref(model->Boundary("BOTTOM"));
    typename vector<Node< dim>*>::const_iterator nit( bref.NodesBegin() );
    heat_flux  = (*nit)->Read(hfb_key);
    std::cerr<<"heat flux: " << heat_flux() << std::endl;

    if (model->ContainsRegion("HOT_BOTTOM"))
    {
        const Region<dim> &rref(model->Region("HOT_BOTTOM"));
        typename vector<Node< dim>*>::const_iterator nit( rref.NodesBegin() );
        heat_flux1 = (*nit)->Read(hfb_key);
        std::cerr<<"heat flux HOT_BOTTOM: " << heat_flux1() << std::endl;
    }

    model->InputPropertyValue( "nodal heat flux bottom", ScalarVariable(ANY, 0.) );

    if (heat_flux() != 0.)
    {
        const typename vector<Face< dim>*>::const_iterator boundaryCellsEnd( bref.CellsEnd() );
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


    if (model->ContainsRegion("HOT_BOTTOM") && heat_flux1() != 0.)
    {
        const Region<dim> &rref(model->Region("HOT_BOTTOM"));
        vector<Element<dim>*>::const_iterator eit;
        for ( eit = rref.CellsBegin(); eit != rref.CellsEnd(); ++eit )
        {
            nnodes = (*eit)->Nodes();
            area = (*eit)->Volume(); // okay to use in case of a 1D bottom boundary ( volume = length = area in case thickness = 1 m)

            for (uint32_t i = 0U; i < nnodes; i++)
            {
                hfb()  = (*eit)->N(i)->Read(hfb_key);
                hfb() += heat_flux1() * area / nnodes;
                (*eit)->N(i)->Store(hfb_key, hfb);
            }
        }
    }


}

void CVFEM_2D_VVCase::InitializeModel()
{
    //! Apply bottom heat flux
    ApplyBottomHeatFlux();

    //! read time increment from config file
    // Time variables
    if      (time_unit == YEARS)
        time_multiplier = year;
    else if (time_unit == DAYS)
        time_multiplier = day;
    else if (time_unit == HOURS)
        time_multiplier = hour;
    else if (time_unit == MINUTES)
        time_multiplier = minute;

    max_time                    = run_settings.Duration();
    double& model_time          = ModelTime::Instance().modelTime;
    model_time                  = 0.0; // reset to zero needed here in Benchmarks only to allow for running multiple tests in a row
    output_increment            = run_settings.NearestOutputTime(0.0);
    model_output                = output_increment;
    double largest_time_step    = run_settings.TimeIncrement();

    largest_time_step           *= time_multiplier;
    max_time                    *= time_multiplier;

    //! Set all unassigned values to 0
    InitializeAllNaNsToZero();

    //! Set up CVFEM
    double cfl_scaling = model->Read(pd_ref->StorageKey("cfl scaling")); // cfl = 0.1 used in Weis et al., 2014

    CVFEM_PHX-> Adjust_CFL_Criterion(cfl_scaling, true); //true = pore velocity based; false = velocity based
    CVFEM_PHX-> SetLargestTimeStep( largest_time_step );


    InitializeHydrostaticPressure();

    //! Create initial output file
    vtu->OutputDataToVTU( string(output_name_).c_str(), output_props_, model->Region("Model"), 0);
    // model->OutputToBinaryFile(( string(config_file_name_) + "_Initial_Properties" ).c_str());

    //! Compare initial output with benchmark simulations
    /*
    set<string> all_variables; // empty set prompts model to read all the variables contained in the binary
    string initial_props = benchmarkFilePath( config_file_name_, "_Initial_Properties" );
    Model<dim>* model_initial = new Model<dim>(initial_props, all_variables);

    // Create property maps
    std::cerr << "\nCOMPARING INITIAL MODEL STATE:"<<std::endl;
    for (const std::string& prop : output_props_ )
    {
        std::map<Point<dim>, std::vector<double>> property_map_bench, property_map_current;

        // TO DO: distinguish between point and cell data!
        property_map_bench   = BuildPropertyMap( *model_initial, prop );
        property_map_current = BuildPropertyMap( *model, prop );

        // Compare property maps
        bool matching_maps(false);
        std::string error_message;
        error_message ="";

        matching_maps = ComparePropertyMaps(property_map_bench, property_map_current, prop, error_message);
        (matching_maps) ? std::cerr << "Matching property maps: " << prop << std::endl
                        : std::cerr << "Property maps not matching: " << prop << std::endl << error_message << std::endl;

        do_test( matching_maps, "Initial state: property \"" + prop + "\" matches benchmark", __FILE__, __LINE__ );
    }

    delete model_initial;
    // */
    //////////////////////////////// end Compare initial output with benschmark simulations

    CVFEM_PHX->PrepareTransientCalculations();
    CVFEM_PHX->ChangeTimeStepTo(largest_time_step/10.);

    //! Avoid std output
    bool no_std_output = model->Read( pd_ref->StorageKey( "no standard output" ) );
    if ( no_std_output )
    {
        std::cerr << "-------------------------" << std::endl;
        std::cerr << " Stopping all std output " << std::endl;
        std::cerr << "-------------------------" << std::endl;

        // Reopen standard stream, and associate /dev/null
        freopen("/dev/null", "w", stdout);
    }
}

void CVFEM_2D_VVCase::OutputAndBookKeeping()
{
    //! Preparing next time step
    double& model_time  = ModelTime::Instance().modelTime;
    model_time += dt;
    model->InputPropertyValue("dt",ScalarVariable(PLAIN,dt));
    timestep++;

    //! Output *.vtu file
    time_label = model_time / time_multiplier;
    if ( time_label >= model_output )
    {
        vtu->OutputDataToVTU( string(output_name_).c_str(), output_props_, model->Region("Model"), time_label);
        model_output += output_increment;
    }

    //! Cancel in case dt is too small
    if (dt < 1.e-4)
    {
        cerr << "Timestep smaller than 0.1 ms Aborting ..." << endl;
        stop = true;
    }
}

void CVFEM_2D_VVCase::PrintProgressToScreen()
{
    double& model_time  = ModelTime::Instance().modelTime;

    std::cerr << std::endl << std::endl << "////////////////////////";
    std::cerr << std::endl << "////////////////////////";
    std::cerr << std::endl << output_name_.c_str() << std::endl;
    std::cerr << fixed << setprecision(0) << "Model time \t" << "in years: " << model_time / year << std::endl;
    if (dt/day >= 1.)  std::cerr << fixed << setprecision(2) << "Timestep in days:  " << dt / day << std::endl;
    else               std::cerr << fixed << setprecision(2) << "Timestep in hours:  " << dt / hour << std::endl;
    std::cerr << fixed << setprecision(2) << "Progress:\t" << model_time / max_time * 100 << "%" << std::endl << std::endl;

    std::cerr << defaultfloat;
    std::cerr.precision(4);
}


//! Main loop to run CVFEM Scheme
void CVFEM_2D_VVCase::run()
{
    double& model_time  = ModelTime::Instance().modelTime;

    InitializeModel();

    time_label = 0;
    timestep = 0;
    stop = false;

    //! Main loop
    while ( model_time <= max_time && !stop )
    {
        PrintProgressToScreen();

        //! Apply CVFEM scheme
        dt = CVFEM_PHX->Apply();

        //! Update dt and output vtu file
        OutputAndBookKeeping();


        //! Compare first iteration with benchmark simulations
        if (timestep==1)
        {
            //! save model after first iteration for vset comparison
            // model->OutputToBinaryFile((string(output_name_) + "_First-iteration_Properties").c_str());
            // model->Region("Model").UpdateMemberIndexes(); // needed because node and element IDs are somehow modified when writing the BinaryFile

            //! Compare model to benchmark vset file
            /*
            const double saved_model_time = model_time;   // protect the live simulation's time
            set<string> all_variables; // empty set prompts model to read all the variables contained in the binary
            string first_it_props  = benchmarkFilePath( config_file_name_, "_First-iteration_Properties" );
            Model<dim>* model_first_it = new Model<dim>(first_it_props, all_variables);

            // Create property maps
            std::cerr << "\nCOMPARING FIRST ITERATION:"<<std::endl;
            for (const std::string& prop : output_props_ )
            {
                std::map<Point<dim>, std::vector<double>> property_map_bench, property_map_current;
                property_map_bench   = BuildPropertyMap( *model_first_it, prop );
                property_map_current = BuildPropertyMap( *model, prop );

                // Compare property maps
                bool matching_maps(false);
                std::string error_message;
                error_message ="";

                matching_maps = ComparePropertyMaps(property_map_bench, property_map_current, prop, error_message);
                (matching_maps) ? std::cerr << "Matching property maps: " << prop << std::endl
                                : std::cerr << "Property maps not matching: " << prop << std::endl << error_message << std::endl;

                do_test( matching_maps, "First iteration: property \"" + prop + "\" matches benchmark", __FILE__, __LINE__ );
            }

            delete model_first_it;

            // restore before continuing the simulation
            model_time = saved_model_time;

            // */

            // Exit while loop in case only initialisation and first iteration are tested
            if (test_type==0)
                break;
        }

    } // end main while loop

    //! save final model state
    // model->OutputToBinaryFile(( string(config_file_name_) + "_Final_Properties" ).c_str());
    std::cout << "Finished Running CVFEM_2D_VVCase " << std::endl;

    //! Compare final state with benchmark simulations
    /*
    if (timestep>1)
    {
        set<string> all_variables; // empty set prompts model to read all the variables contained in the binary
        string final_props  = benchmarkFilePath( config_file_name_, "_Final_Properties" );
        Model<dim>* model_final = new Model<dim>(final_props, all_variables);

        // Create property maps
        std::cerr << "\nCOMPARING FINAL MODEL STATE:"<<std::endl;
        for (const std::string& prop : output_props_ )
        {
            std::map<Point<dim>, std::vector<double>> property_map_bench, property_map_current; property_map_current;
            property_map_bench   = BuildPropertyMap( *model_final, prop );
            property_map_current = BuildPropertyMap( *model, prop );

            // Compare property maps
            bool matching_maps(false);
            std::string error_message;
            error_message ="";

            matching_maps = ComparePropertyMaps(property_map_bench, property_map_current, prop, error_message);
            (matching_maps) ? std::cerr << "Matching property maps: " << prop << std::endl
                            : std::cerr << "Property maps not matching: " << prop << std::endl << error_message << std::endl;

            do_test( matching_maps, "Final state: property \"" + prop + "\" matches benchmark", __FILE__, __LINE__ );
        }

        delete model_final;
    }
    //*/

} // end run CVFEM




//! Create property maps
std::map<Point<2U>, std::vector<double>> CVFEM_2D_VVCase::BuildPropertyMap( const Model<dim>& sg, const std::string& property_name )
{
    std::map<Point<dim>, std::vector<double>> property_map;
    csmp::Index prop_key = sg.Database().StorageKey(property_name.c_str());

    // Reads whatever is stored at prop_key into a flat list of component values —
    // 1 entry for a scalar, N entries for a vector/tensor's N components.
    auto readComponents = [&prop_key]( auto& entity ) -> std::vector<double>
    {
        std::vector<double> values;
        if ( prop_key.type == SCALAR )
        {
            values.push_back( entity->Read(prop_key) );
        }
        else if ( prop_key.type == VECTOR )
        {
            VectorVariable<dim> var;
            entity->Read(prop_key, var);
            for ( size_t i = 0; i < var.Size(); ++i )
                values.push_back( var[i] );
        }
        else if ( prop_key.type == TENSOR )
        {
            TensorVariable<dim> var;
            entity->Read(prop_key, var);
            for ( uint32_t i = 0; i < var.Size(); ++i )
                values.push_back( var.Component(i) );
        }
        return values;
    };

    if ( prop_key.place == NODE)
    {
        // Loop over all nodes and create property map
        for ( typename std::vector<Node<dim>*>::const_iterator
                 nit = sg.Region("Model").NodesBegin(); nit != sg.Region("Model").NodesEnd(); ++nit )
        {
            Point<dim> coord = (*nit)->Coordinate();
            property_map[coord] = readComponents(*nit);
        }
    }

    if ( prop_key.place == ELEMENT)
    {
        // Loop over all cells and create property map
        const typename std::vector<Element<dim>*>::const_iterator CellsEnd( model->Region("Model").CellsEnd() );
        for ( typename std::vector<Element<dim>*>::const_iterator eit( model->Region("Model").CellsBegin() );
             eit != CellsEnd; ++eit )
        {
            Point<dim> coord = (*eit)->BaryCenter();
            property_map[coord] = readComponents(*eit);
        }
    }

    return property_map;
}


//! Compare property maps
bool CVFEM_2D_VVCase::ComparePropertyMaps(
    const std::map<Point<dim>, std::vector<double>>& map_test,
    const std::map<Point<dim>, std::vector<double>>& map_benchmark,
    const std::string& property_name,
    std::string& error_message )
{
    std::ostringstream oss;
    bool all_match = true;

    if ( map_test.size() != map_benchmark.size() )
    {
        oss << "Map size mismatch: test has " << map_test.size()
        << " nodes, benchmark has " << map_benchmark.size() << " nodes.\n";
        all_match = false;
    }

    // catch points that exist only in benchmark
    for ( const auto& [coord, values] : map_benchmark )
    {
        if ( map_test.find(coord) == map_test.end() )
        {
            oss << "Point " << coord << " present in benchmark map but not in test map.\n";
            all_match = false;
        }
    }

    for ( const auto& [coord, values_test] : map_test )
    {
        // catch points that exist only in current simulation
        auto it = map_benchmark.find(coord);
        if ( it == map_benchmark.end() )
        {
            oss << "Point " << coord << " present in test map but not in benchmark map.\n";
            all_match = false;
            continue;
        }


        const std::vector<double>& values_bench = it->second;
        if ( values_test.size() != values_bench.size() )
        {
            oss << "Component count mismatch at " << coord
                << ": test has " << values_test.size()
                << " components, benchmark has " << values_bench.size() << ".\n";
            all_match = false;
            continue;
        }

        for ( size_t c = 0; c < values_test.size(); ++c )
        {
            if(!essentiallyEqual(values_test[c], values_bench[c], numeric_limits<double>::epsilon()) )
            {
                double diff= std::abs(values_test[c] - values_bench[c]);
                oss <<"Mismatch at " << coord <<", property " << property_name
                    <<": test = " << values_test[c]
                    <<", benchmark = " << values_bench[c]
                    <<", diff = " << diff;

                if(definitelyLessThan(diff, values_bench[c]*1.e-6, numeric_limits<double>::epsilon()))
                    oss <<", diff < 0.00001 %\n";
                else
                    oss  <<"\n";

                all_match=false;
            }
        }
    }

    error_message = oss.str();
    return all_match;
}

}//end csmp
