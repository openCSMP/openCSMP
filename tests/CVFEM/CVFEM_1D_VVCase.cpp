#include "CVFEM_1D_VVCase.h"
#include <sys/stat.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>

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
    builder(),
    conduction_erfc_test_(std::string(config_file).find("test-105") != std::string::npos),
    fluid_heat_capacity_(std::string(config_file).find("test-105b") != std::string::npos
                         || std::string(config_file).find("test-105c") != std::string::npos
                         || std::string(config_file).find("test-106") != std::string::npos
                         || std::string(config_file).find("test-107") != std::string::npos
                         || std::string(config_file).find("test-108") != std::string::npos),
    stefan_test_(std::string(config_file).find("test-106") != std::string::npos),
    dike_test_(std::string(config_file).find("test-107") != std::string::npos),
    wide_interval_test_(std::string(config_file).find("test-108") != std::string::npos),
    solidus_(0.),
    liquidus_(0.),
    latent_heat_(0.),
    exponent_(1.),
    cold_temperature_(0.),
    diffusivity_(0.),
    front_parameter_(0.),
    boundary_heat_(0.),
    initial_energy_(0.),
    next_record_time_(record_interval_),
    next_profile_time_(0.),
    front_rows_(),
    profile_rows_(),
    profile_times_(),
    profile_temperatures_()
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

    //! S. Scott, 2026-11-18: rock heat capacity from the configuration file (the scheme's Rock object otherwise uses its
    //! default, 880 J kg-1 K-1, which is the value of tests 101-104); tests 105b/c also put the fluid in the conduction
    //! heat capacity
    CVFEM_PHX->SetRockHeatCapacity( model->Region("Model").NodeVector().front()->Read( pd_ref->StorageKey("nodal heat capacity rock") ) );
    CVFEM_PHX->AddFluidContributionToHeatCapacity( fluid_heat_capacity_ );
    if (stefan_test_ || dike_test_ || wide_interval_test_)
        ConfigureLatentHeat();

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


/** S. Scott, 2026-11-18: largest deviation of the nodal temperature from the error-function solution for a suddenly
    heated end, T(x, t) = T_right + (T_left - T_right) erfc( x / (2 sqrt(kappa t)) ) with x the distance from the LEFT end
    and kappa = lambda / ((1 - phi) rho_r c_r + phi rho_w c_w) from the model's properties (Carslaw & Jaeger 1959). Valid
    while the warming has not reached the RIGHT end (tests 105a-c stop when 2 sqrt(kappa t) is a fifth of the length). */
double CVFEM_1D_VVCase::ConductionErfcDeviation() const
{
    const double model_time( ModelTime::Instance().modelTime );
    const Index t_key( pd_ref->StorageKey("temperature") );
    const Index phi_key( pd_ref->StorageKey("nodal porosity") );
    const Index rho_r_key( pd_ref->StorageKey("nodal density rock") );
    const Index cp_r_key( pd_ref->StorageKey("nodal heat capacity rock") );
    const Index rho_w_key( pd_ref->StorageKey("density liquid") );
    const Index cp_w_key( pd_ref->StorageKey("fluid heat capacity") );
    const Index lambda_key( pd_ref->StorageKey("thermal conductivity") );
    double t_left(0.), t_right(0.), x_left(0.), heat_capacity(0.), nodes(0.), conductivity(0.), cells(0.);
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        if ((*it)->AtBoundary() == CNR1) { t_left = (*it)->Read(t_key); x_left = (*it)->y(); }
        if ((*it)->AtBoundary() == CNR2) t_right = (*it)->Read(t_key);
        const double phi( (*it)->Read(phi_key) );
        heat_capacity += (1. - phi) * (*it)->Read(rho_r_key) * (*it)->Read(cp_r_key)
                         + phi * (*it)->Read(rho_w_key) * (*it)->Read(cp_w_key);
        nodes += 1.;
    }
    for ( auto eit = model->Region("Model").CellsBegin(); eit != model->Region("Model").CellsEnd(); ++eit )
    {
        conductivity += (*eit)->Read(lambda_key);
        cells += 1.;
    }
    const double kappa( (conductivity / cells) / (heat_capacity / nodes) );
    double deviation(0.);
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        const double x( std::fabs((*it)->y() - x_left) );
        const double exact( t_right + (t_left - t_right) * std::erfc( x / (2. * std::sqrt(kappa * model_time)) ) );
        deviation = std::max( deviation, std::fabs((*it)->Read(t_key) - exact) );
    }
    return deviation;
}

/** S. Scott, 2026-11-20: tests 106, 107 and 108. Switches the rock liquidus/solidus on (nodal "liquidus temperature"
    and "solidus temperature" from the configuration file) with the power-law crystallisation curve, the
    configuration's "latent heat of fusion", "crystallization curve exponent" and melt-to-rock heat capacity ratio
    ("heat capacity melt" / "nodal heat capacity rock", which must be 1 as in the analytical solutions). */
void CVFEM_1D_VVCase::ConfigureLatentHeat()
{
    const Node<dim>* node( model->Region("Model").NodeVector().front() );
    liquidus_    = node->Read( pd_ref->StorageKey("liquidus temperature") );
    solidus_     = node->Read( pd_ref->StorageKey("solidus temperature") );
    latent_heat_ = model->Read( pd_ref->StorageKey("latent heat of fusion") );
    exponent_    = model->Read( pd_ref->StorageKey("crystallization curve exponent") );
    const double heat_capacity_ratio( model->Read( pd_ref->StorageKey("heat capacity melt") )
                                      / node->Read( pd_ref->StorageKey("nodal heat capacity rock") ) );
    if (!(liquidus_ > solidus_) || std::fabs(heat_capacity_ratio - 1.) > 1.e-12)
        throw csmp::Exception( FATAL_ERROR, "CVFEM_1D_VVCase::ConfigureLatentHeat",
                               "tests 106, 107 and 108 need a liquidus above the solidus and a melt heat capacity equal to the rock heat capacity" );
    CVFEM_PHX->WithRockLiquidusSolidus( true );
    CVFEM_PHX->SetRockCrystallizationCurve( heat_capacity_ratio, 1., latent_heat_, exponent_, "power law" );
    std::cerr << "\nCVFEM_1D_VVCase::ConfigureLatentHeat: solidus " << solidus_ << " C, liquidus " << liquidus_
              << " C, latent heat " << latent_heat_ << " J kg-1, exponent " << exponent_ << std::endl;
}

void CVFEM_1D_VVCase::InitializeDikeTemperature()
{
    const Index t_key( pd_ref->StorageKey("temperature") );
    double y_left(0.);
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
        if ((*it)->AtBoundary() == CNR1) y_left = (*it)->y();
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        ScalarVariable temperature;
        (*it)->Read(t_key, temperature);
        if (std::fabs((*it)->y() - y_left) < dike_half_width_ && temperature.Flag() != DIRICH)
            (*it)->Store(t_key, ScalarVariable(temperature.Flag(), liquidus_));
    }
}

void CVFEM_1D_VVCase::ComputeAnalyticalParameters()
{
    const Profile profile( SortedProfile() );
    cold_temperature_ = dike_test_ ? profile.temperature.back() : profile.temperature.front();
    const double heat_capacity( SensibleHeatCapacity() );
    diffusivity_ = MeanConductivity() / heat_capacity;
    // L sqrt(pi) / (c (T_m - T_0)) of T&S eqs. 4-141 and 4-149, with heat capacity and latent heat per bulk volume
    const double target( RockMassPerVolume() * latent_heat_ * std::sqrt(std::numbers::pi)
                         / (heat_capacity * (liquidus_ - cold_temperature_)) );
    if (stefan_test_)
        front_parameter_ = Root( [](double x) { return std::exp(-x * x) / (x * std::erf(x)); }, target );
    if (dike_test_)
        front_parameter_ = Root( [](double x) { return std::exp(-x * x) / (x * (1. + std::erf(x))); }, target );
}

template <typename Function>
double CVFEM_1D_VVCase::Root(Function f, double target)
{
    double low{1.e-9};
    double high{5.};
    for (unsigned iteration = 0U; iteration < 200U; ++iteration)
    {
        const double middle( 0.5 * (low + high) );
        if (f(middle) > target) low = middle;
        else                    high = middle;
    }
    return 0.5 * (low + high);
}

double CVFEM_1D_VVCase::MeltFraction(double temperature) const
{
    if (temperature <= solidus_) return 0.;
    if (temperature >= liquidus_) return 1.;
    return std::pow( (temperature - solidus_) / (liquidus_ - solidus_), exponent_ );
}

double CVFEM_1D_VVCase::RockEnthalpy(double temperature) const
{
    const double heat_capacity_rock( model->Region("Model").NodeVector().front()->Read( pd_ref->StorageKey("nodal heat capacity rock") ) );
    return heat_capacity_rock * temperature + latent_heat_ * MeltFraction(temperature);
}

double CVFEM_1D_VVCase::SensibleHeatCapacity() const
{
    const Index phi_key( pd_ref->StorageKey("nodal porosity") );
    const Index rho_r_key( pd_ref->StorageKey("nodal density rock") );
    const Index cp_r_key( pd_ref->StorageKey("nodal heat capacity rock") );
    const Index rho_f_key( pd_ref->StorageKey("fluid density") );
    const Index cp_f_key( pd_ref->StorageKey("fluid heat capacity") );
    double heat_capacity{0.};
    double nodes{0.};
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        const double phi( (*it)->Read(phi_key) );
        heat_capacity += (1. - phi) * (*it)->Read(rho_r_key) * (*it)->Read(cp_r_key)
                         + phi * (*it)->Read(rho_f_key) * (*it)->Read(cp_f_key);
        nodes += 1.;
    }
    return heat_capacity / nodes;
}

double CVFEM_1D_VVCase::MeanConductivity() const
{
    const Index lambda_key( pd_ref->StorageKey("thermal conductivity") );
    double conductivity{0.};
    double cells{0.};
    for ( auto eit = model->Region("Model").CellsBegin(); eit != model->Region("Model").CellsEnd(); ++eit )
    {
        conductivity += (*eit)->Read(lambda_key);
        cells += 1.;
    }
    return conductivity / cells;
}

double CVFEM_1D_VVCase::RockMassPerVolume() const
{
    const Index phi_key( pd_ref->StorageKey("nodal porosity") );
    const Index rho_r_key( pd_ref->StorageKey("nodal density rock") );
    double mass{0.};
    double nodes{0.};
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        mass += (1. - (*it)->Read(phi_key)) * (*it)->Read(rho_r_key);
        nodes += 1.;
    }
    return mass / nodes;
}

double CVFEM_1D_VVCase::StoredEnergy() const
{
    const Index t_key( pd_ref->StorageKey("temperature") );
    const Index phi_key( pd_ref->StorageKey("nodal porosity") );
    const Index rho_r_key( pd_ref->StorageKey("nodal density rock") );
    const Index volume_key( pd_ref->StorageKey("bulk volume") );
    const Index rho_f_key( pd_ref->StorageKey("fluid density") );
    const Index h_f_key( pd_ref->StorageKey("fluid enthalpy") );
    double energy{0.};
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
    {
        const double phi( (*it)->Read(phi_key) );
        energy += (*it)->Read(volume_key) * ( (1. - phi) * (*it)->Read(rho_r_key) * RockEnthalpy( (*it)->Read(t_key) )
                                              + phi * (*it)->Read(rho_f_key) * (*it)->Read(h_f_key) );
    }
    return energy;
}

CVFEM_1D_VVCase::Profile CVFEM_1D_VVCase::SortedProfile() const
{
    const Index t_key( pd_ref->StorageKey("temperature") );
    double y_left(0.);
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
        if ((*it)->AtBoundary() == CNR1) y_left = (*it)->y();
    std::vector<std::pair<double, double>> nodes;
    for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
        nodes.emplace_back( std::fabs((*it)->y() - y_left), (*it)->Read(t_key) );
    std::sort( nodes.begin(), nodes.end() );
    Profile profile;
    for (const auto& [distance, temperature] : nodes)
    {
        profile.distance.push_back(distance);
        profile.temperature.push_back(temperature);
    }
    return profile;
}

double CVFEM_1D_VVCase::Interpolate(const Profile& profile, double distance)
{
    const auto upper( std::upper_bound(profile.distance.begin(), profile.distance.end(), distance) );
    const std::size_t i( std::clamp<std::size_t>( static_cast<std::size_t>(upper - profile.distance.begin()), 1U,
                                                  profile.distance.size() - 1U ) );
    const double weight( (distance - profile.distance[i - 1U]) / (profile.distance[i] - profile.distance[i - 1U]) );
    return (1. - weight) * profile.temperature[i - 1U] + weight * profile.temperature[i];
}

double CVFEM_1D_VVCase::BoundaryHeatFlow() const
{
    const Profile profile( SortedProfile() );
    const std::size_t last( profile.distance.size() - 1U );
    const double conductivity( MeanConductivity() );
    double flow( conductivity * (profile.temperature[last - 1U] - profile.temperature[last])
                 / (profile.distance[last] - profile.distance[last - 1U]) );
    if (!dike_test_)   // the LEFT end is insulated in test-107
        flow += conductivity * (profile.temperature[1] - profile.temperature[0]) / (profile.distance[1] - profile.distance[0]);
    return flow;
}

double CVFEM_1D_VVCase::FrontPosition(const Profile& profile) const
{
    // first crossing of the solidus, walking from the cold end (106 and 108: LEFT, 107: RIGHT) towards the melt
    const std::size_t nodes( profile.distance.size() );
    for (std::size_t k = 1U; k < nodes; ++k)
    {
        const std::size_t solid( dike_test_ ? nodes - k : k - 1U );
        const std::size_t melt( dike_test_ ? nodes - k - 1U : k );
        if (profile.temperature[melt] >= solidus_)
        {
            const double weight( (solidus_ - profile.temperature[solid])
                                 / (profile.temperature[melt] - profile.temperature[solid]) );
            return profile.distance[solid] + weight * (profile.distance[melt] - profile.distance[solid]);
        }
    }
    return dike_test_ ? 0. : profile.distance.back();
}

double CVFEM_1D_VVCase::AnalyticalTemperature(double distance, double time, bool with_latent_heat) const
{
    const double length( 2. * std::sqrt(diffusivity_ * time) );
    const double step( liquidus_ - cold_temperature_ );
    if (dike_test_)
    {
        // eta from the contact, positive into the country rock
        const double eta( (distance - dike_half_width_) / length );
        if (!with_latent_heat)   // two half-spaces at T_m and T_0 in contact
            return cold_temperature_ + step * 0.5 * std::erfc(eta);
        if (eta <= -front_parameter_) return liquidus_;
        return cold_temperature_ + step * std::erfc(eta) / (1. + std::erf(front_parameter_));
    }
    const double eta( distance / length );
    if (!with_latent_heat || !stefan_test_)   // half-space at T_m, end held at T_0
        return cold_temperature_ + step * std::erf(eta);
    if (eta >= front_parameter_) return liquidus_;
    return cold_temperature_ + step * std::erf(eta) / std::erf(front_parameter_);
}

double CVFEM_1D_VVCase::AnalyticalFront(double time, bool with_latent_heat) const
{
    const double length( 2. * std::sqrt(diffusivity_ * time) );
    const double fraction( (solidus_ - cold_temperature_) / (liquidus_ - cold_temperature_) );
    if (dike_test_)
    {
        if (with_latent_heat) return dike_half_width_ - front_parameter_ * length;
        // 0.5 erfc(eta) = fraction, eta < 0
        return dike_half_width_ - length * Root( [](double x) { return -0.5 * std::erfc(-x); }, -fraction );
    }
    if (with_latent_heat) return front_parameter_ * length;
    return length * Root( [](double x) { return -std::erf(x); }, -fraction );
}

void CVFEM_1D_VVCase::RecordLatentHeatState()
{
    const double model_time( ModelTime::Instance().modelTime );
    const Profile profile( SortedProfile() );
    if (model_time >= next_record_time_)
    {
        front_rows_ << model_time / year;
        if (wide_interval_test_)
            front_rows_ << "," << boundary_heat_ << "," << StoredEnergy() - initial_energy_;
        else
        {
            front_rows_ << "," << FrontPosition(profile) << "," << AnalyticalFront(model_time, true) << ","
                        << AnalyticalFront(model_time, false);
            if (dike_test_)
                front_rows_ << "," << Interpolate(profile, dike_half_width_) << ","
                            << AnalyticalTemperature(dike_half_width_, model_time, true) << ","
                            << AnalyticalTemperature(dike_half_width_, model_time, false);
        }
        front_rows_ << "\n";
        next_record_time_ += record_interval_;
    }
    if (model_time >= next_profile_time_)
    {
        if (wide_interval_test_)
        {
            profile_times_.push_back(model_time);
            profile_temperatures_.push_back(profile.temperature);
        }
        else
            for (std::size_t i = 0U; i < profile.distance.size(); ++i)
                profile_rows_ << model_time / year << "," << profile.distance[i] << "," << profile.temperature[i] << ","
                              << AnalyticalTemperature(profile.distance[i], model_time, true) << ","
                              << AnalyticalTemperature(profile.distance[i], model_time, false) << "\n";
        next_profile_time_ += 0.25 * max_time;
    }
}

std::vector<std::vector<double>> CVFEM_1D_VVCase::EnthalpyMethodReference(const std::vector<double>& distances,
                                                                          const std::vector<double>& times) const
{
    // fixed end temperatures (the initial state of test-108: melt at the RIGHT end temperature, LEFT end cold)
    const std::size_t intervals( (distances.size() - 1U) * refinement_ );
    const double dx( distances.back() / static_cast<double>(intervals) );
    const double heat_capacity( SensibleHeatCapacity() );
    const double latent_heat_per_volume( RockMassPerVolume() * latent_heat_ );
    const double conductivity( MeanConductivity() );
    const double stable_time_step( 0.25 * dx * dx * heat_capacity / conductivity );   // explicit Euler: < 0.5 dx^2 / kappa

    // volumetric enthalpy E(T) = C T + (1 - phi) rho_r L X(T) and its inverse
    auto enthalpy = [&](double temperature) { return heat_capacity * temperature + latent_heat_per_volume * MeltFraction(temperature); };
    auto temperature_of = [&](double energy)
    {
        if (energy <= enthalpy(solidus_)) return energy / heat_capacity;
        if (energy >= enthalpy(liquidus_)) return (energy - latent_heat_per_volume) / heat_capacity;
        double low(solidus_), high(liquidus_);
        for (unsigned iteration = 0U; iteration < 60U; ++iteration)
        {
            const double middle( 0.5 * (low + high) );
            if (enthalpy(middle) < energy) low = middle;
            else                           high = middle;
        }
        return 0.5 * (low + high);
    };

    std::vector<double> temperature( intervals + 1U, liquidus_ );
    temperature.front() = cold_temperature_;
    std::vector<double> energy( intervals + 1U );
    for (std::size_t i = 0U; i <= intervals; ++i) energy[i] = enthalpy(temperature[i]);
    std::vector<std::vector<double>> profiles;
    double time{0.};
    for (const double target : times)
    {
        const std::size_t steps( static_cast<std::size_t>(std::ceil((target - time) / stable_time_step)) );
        const double factor( (target - time) / static_cast<double>(steps) * conductivity / (dx * dx) );
        for (std::size_t step = 0U; step < steps; ++step)
        {
            for (std::size_t i = 1U; i < intervals; ++i)
                energy[i] += factor * (temperature[i - 1U] - 2. * temperature[i] + temperature[i + 1U]);
            for (std::size_t i = 1U; i < intervals; ++i)
                temperature[i] = temperature_of(energy[i]);
        }
        time = target;
        std::vector<double> interpolated;
        for (const double distance : distances)
        {
            const std::size_t i( std::min(static_cast<std::size_t>(distance / dx), intervals - 1U) );
            const double weight( (distance - static_cast<double>(i) * dx) / dx );
            interpolated.push_back( (1. - weight) * temperature[i] + weight * temperature[i + 1U] );
        }
        profiles.push_back(interpolated);
    }
    return profiles;
}

void CVFEM_1D_VVCase::WriteLatentHeatOutput()
{
    const std::string base( config_file_name_ );
    const std::string front_name( base + (wide_interval_test_ ? "_energy.csv" : "_front.csv") );
    std::ofstream front( front_name );
    if (wide_interval_test_)
        front << "time_yr,heat_conducted_out_J_m2,stored_energy_change_J_m2\n";
    else if (dike_test_)
        front << "time_yr,y_front_model_m,y_front_latent_m,y_front_no_latent_m,"
                 "T_contact_model_C,T_contact_latent_C,T_contact_no_latent_C\n";
    else
        front << "time_yr,y_front_model_m,y_front_latent_m,y_front_no_latent_m\n";
    front << front_rows_.str();

    std::ofstream profiles( base + "_profiles.csv" );
    profiles << (wide_interval_test_ ? "time_yr,y_m,T_model_C,T_reference_latent_C,T_no_latent_C\n"
                                     : "time_yr,y_m,T_model_C,T_analytical_latent_C,T_analytical_no_latent_C\n");
    if (wide_interval_test_)
    {
        const Profile profile( SortedProfile() );
        const std::vector<std::vector<double>> reference( EnthalpyMethodReference(profile.distance, profile_times_) );
        for (std::size_t k = 0U; k < profile_times_.size(); ++k)
            for (std::size_t i = 0U; i < profile.distance.size(); ++i)
                profiles << profile_times_[k] / year << "," << profile.distance[i] << "," << profile_temperatures_[k][i]
                         << "," << reference[k][i] << ","
                         << AnalyticalTemperature(profile.distance[i], profile_times_[k], false) << "\n";
    }
    else
        profiles << profile_rows_.str();
    std::cerr << "CVFEM_1D_VVCase::WriteLatentHeatOutput: " << front_name << ", " << base << "_profiles.csv" << std::endl;
}

void CVFEM_1D_VVCase::CheckStefanSolidification()
{
    const double model_time( ModelTime::Instance().modelTime );
    const Profile profile( SortedProfile() );
    const double front_latent( AnalyticalFront(model_time, true) );
    const double front_model( FrontPosition(profile) );
    const double deviation( front_model / front_latent - 1. );
    // solid layer away from the front (within 80 % of y_m): the front cell holds the narrow interval
    double temperature_deviation{0.};
    for (std::size_t i = 0U; i < profile.distance.size(); ++i)
        if (profile.distance[i] <= 0.8 * front_latent)
            temperature_deviation = std::max( temperature_deviation, std::fabs(profile.temperature[i]
                                              - AnalyticalTemperature(profile.distance[i], model_time, true)) );
    const double temperature_limit( temperature_tolerance_ * (liquidus_ - cold_temperature_) );
    std::cerr << "\nCVFEM_1D_VVCase::run: " << config_file_name_ << ": after " << model_time / year << " years, front (solidus) at "
              << front_model << " m; T&S eq. 4-136 (latent heat, lambda1 " << front_parameter_ << ") " << front_latent
              << " m, without latent heat " << AnalyticalFront(model_time, false) << " m; relative deviation " << deviation
              << " (limit " << front_tolerance_ << "); solid-layer temperature deviation from eq. 4-137 " << temperature_deviation
              << " K (limit " << temperature_limit << " K)" << std::endl;
    do_test( std::fabs(deviation) < front_tolerance_, "Stefan problem: solidification front follows T&S eq. 4-136", __FILE__, __LINE__ );
    do_test( temperature_deviation < temperature_limit, "Stefan problem: solid-layer temperature follows T&S eq. 4-137", __FILE__, __LINE__ );
}

void CVFEM_1D_VVCase::CheckDikeSolidification()
{
    const double model_time( ModelTime::Instance().modelTime );
    const Profile profile( SortedProfile() );
    const double travel_latent( dike_half_width_ - AnalyticalFront(model_time, true) );
    const double travel_model( dike_half_width_ - FrontPosition(profile) );
    const double deviation( travel_model / travel_latent - 1. );
    const double contact_model( Interpolate(profile, dike_half_width_) );
    const double contact_latent( AnalyticalTemperature(dike_half_width_, model_time, true) );
    const double temperature_limit( temperature_tolerance_ * (liquidus_ - cold_temperature_) );
    std::cerr << "\nCVFEM_1D_VVCase::run: " << config_file_name_ << ": after " << model_time / year << " years, front moved "
              << travel_model << " m into the dike; T&S eq. 4-142 (latent heat, lambda2 " << front_parameter_ << ") "
              << travel_latent << " m, without latent heat " << dike_half_width_ - AnalyticalFront(model_time, false)
              << " m; relative deviation " << deviation << " (limit " << front_tolerance_ << "); contact temperature "
              << contact_model << " C, eq. 4-146 " << contact_latent << " C, without latent heat "
              << AnalyticalTemperature(dike_half_width_, model_time, false) << " C (limit " << temperature_limit << " K)" << std::endl;
    do_test( std::fabs(deviation) < front_tolerance_, "Dike solidification: front follows T&S eq. 4-142", __FILE__, __LINE__ );
    do_test( std::fabs(contact_model - contact_latent) < temperature_limit,
             "Dike solidification: contact temperature follows T&S eq. 4-146", __FILE__, __LINE__ );
}

void CVFEM_1D_VVCase::CheckWideIntervalCooling()
{
    const double model_time( ModelTime::Instance().modelTime );
    const double energy_change( StoredEnergy() - initial_energy_ );
    const double residual( (energy_change + boundary_heat_) / boundary_heat_ );
    const Profile profile( SortedProfile() );
    const std::vector<double> reference( EnthalpyMethodReference(profile.distance, {model_time}).front() );
    double deviation{0.};
    for (std::size_t i = 0U; i < profile.distance.size(); ++i)
        deviation = std::max( deviation, std::fabs(profile.temperature[i] - reference[i]) );
    const double temperature_limit( temperature_tolerance_ * (liquidus_ - cold_temperature_) );
    std::cerr << "\nCVFEM_1D_VVCase::run: " << config_file_name_ << ": after " << model_time / year << " years, heat conducted out "
              << boundary_heat_ << " J m-2, change of stored energy (rock incl. latent heat + fluid) " << energy_change
              << " J m-2, residual " << residual << " of the heat conducted out (limit " << energy_tolerance_
              << "); largest deviation from the enthalpy-method reference " << deviation << " K (limit "
              << temperature_limit << " K)" << std::endl;
    do_test( std::fabs(residual) < energy_tolerance_, "Wide melting interval: energy balance with latent heat", __FILE__, __LINE__ );
    do_test( deviation < temperature_limit, "Wide melting interval: temperatures match the enthalpy-method reference", __FILE__, __LINE__ );
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
    const bool latent_heat_test( stefan_test_ || dike_test_ || wide_interval_test_ );   // S. Scott, 2026-11-20
    if (dike_test_)
        InitializeDikeTemperature();
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
    //! S. Scott, 2026-11-18: prepare the transient calculation as CVFEM_2D_VVCase does (previous-step quantities), first
    //! time step a tenth of the largest; only for tests 105a-c, so that tests 101-104 are unchanged (S. Scott,
    //! 2026-11-20: also for tests 106, 107 and 108, which in addition record the initial state)
    const bool analytical_test( conduction_erfc_test_ || latent_heat_test );
    if (analytical_test)
    {
        CVFEM_PHX->PrepareTransientCalculations();
        CVFEM_PHX->ChangeTimeStepTo(run_settings.TimeIncrement() / 10.);
    }
    vtu_dt = run_settings.NearestOutputTime(0.0);
    max_time = run_settings.Duration();
    timestep = 0;
    output_counter = 0;
    if (latent_heat_test)
    {
        ComputeAnalyticalParameters();
        initial_energy_ = StoredEnergy();
        boundary_heat_ = 0.;
        next_profile_time_ = 0.25 * max_time;
    }

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
        if (latent_heat_test)   // S. Scott, 2026-11-20: heat conducted out through the fixed-temperature ends, CSV rows
        {
            boundary_heat_ += dt * BoundaryHeatFlow();
            RecordLatentHeatState();
        }

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
            if (test_type==0 && !analytical_test)   // S. Scott, 2026-11-18: tests 105a-c (and 106-108) always run to the end
                break;
        }

    }

    //! save final model state
    // model->OutputToBinaryFile(( string(config_file_name_) + "_Final_Properties" ).c_str());
    std::cout << "Finished Running CVFEM_1D_VVCase " << std::endl;

    //! S. Scott, 2026-11-18: tests 105a-c, heat conduction from the heated LEFT end against the error-function solution
    if (conduction_erfc_test_)
    {
        double t_left(0.), t_right(0.);
        const Index t_key( pd_ref->StorageKey("temperature") );
        for ( auto it = model->Region("Model").NodesBegin(); it != model->Region("Model").NodesEnd(); ++it )
        {
            if ((*it)->AtBoundary() == CNR1) t_left  = (*it)->Read(t_key);
            if ((*it)->AtBoundary() == CNR2) t_right = (*it)->Read(t_key);
        }
        const double deviation( ConductionErfcDeviation() );
        const double limit( 0.01 * std::fabs(t_left - t_right) );   // 1 % of the temperature step
        std::cerr << "\nCVFEM_1D_VVCase::run: " << config_file_name_ << ": largest deviation from the error-function solution "
                  << deviation << " K (limit " << limit << " K) after " << model_time / year << " years" << std::endl;
        do_test( deviation < limit, "Heat conduction from the heated end matches the error-function solution", __FILE__, __LINE__ );
    }
    if (latent_heat_test)
        WriteLatentHeatOutput();
    if (stefan_test_)
        CheckStefanSolidification();
    if (dike_test_)
        CheckDikeSolidification();
    if (wide_interval_test_)
        CheckWideIntervalCooling();

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
