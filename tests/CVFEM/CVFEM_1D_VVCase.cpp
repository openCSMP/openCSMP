#include "CVFEM_1D_VVCase.h"
#include <sys/stat.h>

using namespace std;



/** Philipp Weis tests ported from Alina Yapparova's Geyser Simulator by Edoardo Pezzulli fixed by Jonas Köpping (just to add more names to the list)
   This code is experimental...

   Config File:  Test101,
                 Test101-wg
   Regions File: Test101,
                 Test101-wg
   Var File:     Created by ModelBuilder



**/
namespace csmp
{

namespace {
std::string benchmarkFilePath(const std::string& configName, const std::string& suffix)
{
    return configName + suffix;
}
} // anonymous namespace



CVFEM_1D_VVCase::CVFEM_1D_VVCase(const char* config_file, const uint32_t test_type )
    : vars_name_("PhysicalVariables.txt"),
    config_file_name_(config_file),
    test_type(test_type),
    builder()
{
    std::cout << "configuration file: " << config_file_name_ << std::endl;
    // cin.get();

    //Dimensions and resolution of 1D mesh
    double len = 2000; // [m]
    int      n = 200; // number of elements

    double length(len);
    size_t n_elements (n);

    if (length == 0.0 || n_elements == 0.0)
        throw csmp::Exception(FATAL_ERROR,
                              "CVFEM_1D_VVCase<dim>::CVFEM_1D_VVCase()",
                              "invalid input parameters: model length or number of elements");

    //! Create PhysicalVariables file
    builder.CreateVariablesFile(false, false, false, false, false, false);
    builder.CreateLimiterLogFile();

    //! Create a 1D model
    model = new Model1D<dim>("Line", vars_name_, length, n_elements);
    model->InstantiateFiniteVolumes();

    pd_ref = &(model->Database()); //reference to the models property database.

    InputDataManager<dim>  model_configuration;
    model_configuration.ConfigureFromFile(*model, config_file_name_,
                                          false,    // groupname from parameter range
                                          true,    // default property values
                                          true,    // regional property values
                                          true,    // boundary conditions for box-shaped model
                                          true,    // essential conditions for groups
                                          true,    // csmp::Boundary properties
                                          run_settings );


    //! assign the gravity flag from the config file
    bool with_gravity = false;
    csmp::Index   g_key(model->Database().StorageKey("gravity flag"));
    with_gravity = bool(model->Read(g_key));

    //! read time increment from config file
    double largest_time_step = run_settings.TimeIncrement();

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
    CVFEM_PHX->SetLargestTimeStep(largest_time_step);

    //! create new permeability visitor and assign permeability tensor
    permeability_visitor = new PermeabilityVisitor<dim>(*model);
    permeability_visitor->AssignPermeabilityTensor(*model);

    //! Set halite saturation (Fig. 6c-d; Weis et al., 2014)
    model->InputPropertyValue( "saturation halite", ScalarVariable(ANY, 0.1));

    //! calculate pore volume
    PoreVolumeVisitor<dim> Pore_Volume_Visitor(*model, "nodal porosity", "bulk volume", "pore volume", "thickness");
    model->Accept(Pore_Volume_Visitor);

    //! vtu interface fro output
    vtu = new VTU_Interface<dim> (*model);
    vtu->OmitZeroInFileName( false );

    //! output main variables
    output_props_.push_back("nodal porosity");
    output_props_.push_back("pore volume");
    output_props_.push_back("porosity");
    output_props_.push_back("fluid pressure");
    output_props_.push_back("velocity");
    output_props_.push_back("velocity liquid");
    output_props_.push_back("velocity vapor");
    output_props_.push_back("temperature");
    output_props_.push_back("salinity");
    output_props_.push_back("fluid density");
    output_props_.push_back("fluid enthalpy");
    output_props_.push_back("saturation liquid");
    output_props_.push_back("saturation vapor");
    output_props_.push_back("saturation halite");
    output_props_.push_back("fluid mass liquid");
    output_props_.push_back("fluid mass vapor");
    output_props_.push_back("permeability tensor");

    //! Standard output in terminal: yay or nay
    bool no_std_output(true);
    if ( no_std_output )
    {
        std::cerr << "-------------------------" << std::endl;
        std::cerr << " Stopping all std output " << std::endl;
        std::cerr << "-------------------------" << std::endl;

        // Reopen standard stream, and associate /dev/null
        freopen("/dev/null", "w", stdout);
    }
}

CVFEM_1D_VVCase::~CVFEM_1D_VVCase()
{
    delete CVFEM_PHX;

    delete model;
    delete vtu;
}


int CVFEM_1D_VVCase::InitializeFluidPropertiesLinearPressure()
{
    //! calculate and assign initial pressure gradient
    double p_value;
    ScalarVariable P_value;
    double P_top, P_bottom, T_top, T_bottom;

    csmp::Index   p_key(model->Database().StorageKey("fluid pressure"));
    csmp::Index   t_key(model->Database().StorageKey("temperature"));
    csmp::Index   wt_key(model->Database().StorageKey("salinity"));

    double p_bottom, p_top, p_diff, ymax;
    p_bottom = p_top = p_diff = ymax = std::numeric_limits<double>::quiet_NaN();

    const typename vector<Node<dim>*>::const_iterator modelNodesEnd( model->Region("Model").NodesEnd() );
    for( typename vector<Node<dim>*>::const_iterator it( model->Region("Model").NodesBegin() );
         it != modelNodesEnd; ++it )
    {
        if ((*it)->AtBoundary() == CNR1) // boundary LEFT
        {
            (*it)->Read(p_key, P_value);
            p_top = P_value();
        }
        if ((*it)->AtBoundary() == CNR2) // boundary RIGHT
        {
            (*it)->Read(p_key, P_value);
            p_bottom = P_value();
            ymax = (*it)->y();
        }
    }

    p_diff = p_bottom - p_top;

    for( typename vector<Node<dim>*>::const_iterator it( model->Region("Model").NodesBegin() );
         it != modelNodesEnd; ++it )
    {
        if ((*it)->Status(p_key) != DIRICH)
        {
            p_value = p_top + (*it)->y() / ymax * p_diff;
            (*it)->Store(p_key, makeScalar(PLAIN, p_value));
        }
    }


    //! Initialize fluid properties from PTX
    CVFEM_PHX->InitialFluidPropertiesFromPTX();

    return 0;
}

void CVFEM_1D_VVCase::InitializeAllNaNsToZero()
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


void CVFEM_1D_VVCase::PrintProgressToScreen()
{
    double& model_time  = ModelTime::Instance().modelTime;

    std::cerr << std::endl << std::endl << "////////////////////////";
    std::cerr << std::endl << "////////////////////////";
    std::cerr << std::endl << string(config_file_name_).c_str() << std::endl;
    std::cerr << fixed << setprecision(0) << "Model time \t" << "in years: " << model_time / year << std::endl;
    if (dt/day >= 1.)  std::cerr << fixed << setprecision(2) << "Timestep in days:  " << dt / day << std::endl;
    else               std::cerr << fixed << setprecision(2) << "Timestep in hours:  " << dt / hour << std::endl;
    std::cerr << fixed << setprecision(2) << "Progress:\t" << model_time / max_time * 100 << "%" << std::endl << std::endl;

    std::cerr << defaultfloat;
    std::cerr.precision(4);
}

//This runs the CVFEM scheme
void CVFEM_1D_VVCase::run()
{
    double& model_time  = ModelTime::Instance().modelTime;
    model_time          = 0.; // reset to zero needed here in Benchmarks only to allow for running multiple tests in a row

    //! Calculate linear pressure gradient, initialize fluid properties from PTX, initial output to vtu
    InitializeAllNaNsToZero();
    InitializeFluidPropertiesLinearPressure();

    permeability_visitor->AssignPermeabilityTensor(*model);
    // */

    //! initial output
    cout<<"\n Output Initial Properties";
    vtu->OutputDataToVTU( ( string(config_file_name_) + "_Initial_Properties" ).c_str(), output_props_, model->Region("Model"), 0);
    // model->OutputToBinaryFile(( string(config_file_name_) + "_Initial_Properties" ).c_str());

    //! Compare initial output with benchmark simulations
    /*
    set<string> all_variables; // empty set prompts model to read all the variables contained in the binary
    string initial_props = benchmarkFilePath( config_file_name_, "_Initial_Properties" );
    Model<1U>* model_initial = new Model<1U>(initial_props, all_variables);

    // Create property maps
    std::cerr << "\nCOMPARING INITIAL MODEL STATE:"<<std::endl;
    for (const std::string& prop : output_props_ )
    {
        std::map<Point<1U>, std::vector<double>> property_map_bench, property_map_current;

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
    //*/



    ///////////////////////////////////////////////////////
    //! setup time variables
    vtu_dt = run_settings.NearestOutputTime(0.0);
    max_time = run_settings.Duration();
    timestep = 0;
    output_counter = 0;

    typedef std::chrono::high_resolution_clock Clock; // track duration of iteration
    //! main time loop
    while (model_time <= max_time)
    {
        //! Start timer
        auto start_loop = std::chrono::high_resolution_clock::now();

        PrintProgressToScreen();

        dt = CVFEM_PHX->Apply();

        permeability_visitor->AssignPermeabilityTensor(*model);

        model_time += dt;
        cout<<"\n model time = "<<model_time<<" seconds"<<endl;

        //! output to vtk
        if (model_time >= (output_counter + 1)*vtu_dt)
        {
            vtu->OutputDataToVTU( ( string(config_file_name_) + "_Properties" ).c_str(), output_props_, model->Region("Model"), output_counter);
            output_counter++;
        }

        timestep++;

        //! Stop timer
        auto stop_loop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = (stop_loop - start_loop);
        std::cerr << fixed << setprecision(4) << "Iteration took:\t" << time.count() << " s " << std::endl;

        //! save model after first iteration for vset comparison
        if (timestep==1)
        {
            // model->OutputToBinaryFile(( string(config_file_name_) + "_First-iteration_Properties" ).c_str());
            // model->Region("Model").UpdateMemberIndexes(); // needed because node and element IDs are somehow modified when writing the BinaryFile

            //! Compare first iteration with benchmark simulations
            /*
            set<string> all_variables; // empty set prompts model to read all the variables contained in the binary
            string first_it_props  = benchmarkFilePath( config_file_name_, "_First-iteration_Properties" );
            Model<1U>* model_first_it = new Model<1U>(first_it_props, all_variables);

            // Create property maps
            std::cerr << "\nCOMPARING FIRST ITERATION:"<<std::endl;
            for (const std::string& prop : output_props_ )
            {
                std::map<Point<1U>, std::vector<double>> property_map_bench, property_map_current;
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
            // */


            // Exit while loop in case only initialisation and first iteration are tested
            if (test_type==0)
                break;
        }

    }

    //! save final model state
    // model->OutputToBinaryFile(( string(config_file_name_) + "_Final_Properties" ).c_str());
    std::cout << "Finished Running CVFEM_1D_VVCase " << std::endl;

    //! Compare final state with benchmark simulations
    /*
    if (timestep>1)
    {
        string final_props  = benchmarkFilePath( config_file_name_, "_Final_Properties" );
        Model<1U>* model_final = new Model<1U>(final_props, all_variables);

        // Create property maps
        std::cerr << "\nCOMPARING FINAL MODEL STATE:"<<std::endl;
        for (const std::string& prop : output_props_ )
        {
            std::map<Point<1U>, std::vector<double>> property_map_bench, property_map_current; property_map_current;
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

}


//! Create property maps
std::map<Point<1U>, std::vector<double>> CVFEM_1D_VVCase::BuildPropertyMap( const Model<1U>& sg, const std::string& property_name )
{
    std::map<Point<1U>, std::vector<double>> property_map;
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
            VectorVariable<1U> var;
            entity->Read(prop_key, var);
            for ( size_t i = 0; i < var.Size(); ++i )
                values.push_back( var[i] );
        }
        else if ( prop_key.type == TENSOR )
        {
            TensorVariable<1U> var;
            entity->Read(prop_key, var);
            for ( uint32_t i = 0; i < var.Size(); ++i )
                values.push_back( var.Component(i) );
        }
        return values;
    };

    if ( prop_key.place == NODE)
    {
        // Loop over all nodes and create property map
        for ( typename std::vector<Node<1U>*>::const_iterator
                 nit = sg.Region("Model").NodesBegin(); nit != sg.Region("Model").NodesEnd(); ++nit )
        {
            Point<1U> coord = (*nit)->Coordinate();
            property_map[coord] = readComponents(*nit);
        }
    }

    if ( prop_key.place == ELEMENT)
    {
        // Loop over all cells and create property map
        const typename std::vector<Element<1U>*>::const_iterator CellsEnd( model->Region("Model").CellsEnd() );
        for ( typename std::vector<Element<1U>*>::const_iterator eit( model->Region("Model").CellsBegin() );
             eit != CellsEnd; ++eit )
        {
            Point<1U> coord = (*eit)->BaryCenter();
            property_map[coord] = readComponents(*eit);
        }
    }

    return property_map;
}


//! Compare property maps
bool CVFEM_1D_VVCase::ComparePropertyMaps(
    const std::map<Point<1U>, std::vector<double>>& map_test,
    const std::map<Point<1U>, std::vector<double>>& map_benchmark,
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
