#include "CVFEM_1D_VVCase.h"

using namespace std;



/** Phillip Weiss tests ported from Alina Yapparova's Geyser Simulator by Edoardo Pezzulli
   This code is experimental...

   Config File:  Test101,
                 Test101-wg
   Regions File: Test101,
                 Test101-wg
   Var File:     PhysicalVariablesBenchmark.txt



**/
namespace csmp
{

CVFEM_1D_VVCase::CVFEM_1D_VVCase(const char* config_file )
  : vars_name_("PhysicalVariablesBenchmarks.txt"),
    config_file_name_(config_file)
  {

    enum{dim=1}; //if you change here you must change everywhere

    //The output directory needs to exist for the files to be saved!
    std::string output_Directory = "../Output/Workshop/";


    //Old parameters passed in Geyser command line
    double arg3 = 2000;
    int arg4 = 200;

    double length(arg3),  dx;
    size_t n_elements (arg4);

    if (length == 0.0 || n_elements == 0.0)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_1D_VVCase<dim>::CVFEM_1D_VVCase()",
          "invalid input parameters: model length or number of elements");

    //! Create a 1D model
    dx = length/(double)n_elements;

    model = new Model1D<dim>("Line", vars_name_, length, n_elements);
    model->InstantiateFiniteVolumes();

    pd_ref = &(model->Database()); //reference to the models property database.

    InputDataManager<dim>  model_configuration;
    model_configuration.ConfigureFromFile(*model, config_file_name_, false,    // groupname from parameter range
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
    CVFEM_PHX = new CVFEM_PHX_Scheme<dim>(*model, with_gravity);
    CVFEM_PHX->SetLargestTimeStep(largest_time_step);

    //! vtu interface fro output
    vtu = new VTU_Interface<dim> (*model);

    //! output main variables
    output_props_.push_back("nodal porosity");
    output_props_.push_back("porosity");
    output_props_.push_back("fluid pressure");
    output_props_.push_back("velocity");
    output_props_.push_back("velocity liquid");
    output_props_.push_back("velocity vapor");
    output_props_.push_back("temperature");
    output_props_.push_back("salinity");
    output_props_.push_back("effective diffusivity");
    output_props_.push_back("fluid density");
    //output_props_.push_back("fluid viscosity");
    output_props_.push_back("fluid enthalpy");
    output_props_.push_back("saturation liquid");
    output_props_.push_back("fluid mass liquid");
    output_props_.push_back("fluid mass vapor");

  }

CVFEM_1D_VVCase::~CVFEM_1D_VVCase()
  {
      delete CVFEM_PHX;

      delete model;
      delete vtu;
  }

int CVFEM_1D_VVCase::CalculateSteadyStatePressure ()
{
  enum{dim=1}; //if you change here you must change everywhere


  Index lambda_key          = model->Database ().StorageKey ("conductivity");
  Index lambda_tot_key      = model->Database ().StorageKey ("mass conductivity");

  Index gravity_key         = model->Database().StorageKey ("gravity term");
  Index mass_gravity_key    = model->Database().StorageKey ("mass gravity term");

  Index k_key               = model->Database ().StorageKey ("permeability");
  Index rho_key             = model->Database ().StorageKey ("fluid density");
  Index mu_key              = model->Database ().StorageKey ("fluid viscosity");

  ScalarVariable k, rho, mu, lambda, lambda_tot;
  VectorVariable<dim> gravityVector;

  //CalculateMassConductivity();
  const typename vector<Element< dim>*>::const_iterator modelElementsEnd( model->Region("Model").CellsEnd() );
  for( typename vector<Element<dim>*>::const_iterator it( model->Region("Model").CellsBegin() );
       it != modelElementsEnd; ++it )
  {
    k() = (*it)->Read(k_key);
    for ( uint32_t ip{0u}; ip<(*it)->IntegrationPoints(); ++ip )
    {
      (*it)->PropertyValueAtIntegrationPoint( rho_key, ip, rho );
      (*it)->PropertyValueAtIntegrationPoint( mu_key, ip, mu );
      lambda() = k()/mu();
      lambda_tot() = rho()*k()/mu();
      (*it)->Store( ip, lambda_key, lambda );
      (*it)->Store( ip, lambda_tot_key, lambda_tot );
    }
  }
  //CalculateMassGravityTerm();
  ComputeGravityTermVisitor<dim> gravity_visitor( *model, "gravity term", "permeability",
                                                          "fluid viscosity", "fluid density" );
  (*model).Accept (gravity_visitor);

  for( typename vector<Element<dim>*>::const_iterator it( model->Region("Model").CellsBegin() );
       it != modelElementsEnd; ++it )
  {
    for ( uint32_t ip{0u}; ip<(*it)->IntegrationPoints(); ++ip )
    {
      (*it)->Read(ip, gravity_key, gravityVector);
      (*it)->PropertyValueAtIntegrationPoint( rho_key, ip, rho );
      (*it)->Store( ip, mass_gravity_key, gravityVector*rho() );
    }
  }

  //! steady state pressure
  //SAMG_Solver *p_solver;
  //p_solver = new SAMG_Solver();
  LUdcmp_Solver p_solver;

  PDE_Integrator<dim, Element>  steady_state_pressure( p_solver );

  NumIntegral_dNT_op_dN_dV <dim>  p_conductance( *pd_ref, "mass conductivity",
                                                                        "fluid pressure", "fluid pressure" );
  NumIntegral_dNT_op_dV <dim> p_gravity( *pd_ref, "mass gravity term", "fluid pressure" );
  VelocityAndVolumeFlux<dim>  velocity( *model, "conductivity", "nodal porosity",
                                                              "fluid pressure", "fluid density", false );

  steady_state_pressure.Add( &p_conductance );
  steady_state_pressure.Add( &p_gravity );
  steady_state_pressure.AddPostProcess( &velocity );
  (*model).Apply (steady_state_pressure);

  return 0;
}

int CVFEM_1D_VVCase::CalculateSteadyStateTemperature()
{
  enum{dim=1}; //if you change here you must change everywhere

  LUdcmp_Solver t_solver;

  PDE_Integrator<dim, Element>  steady_state_temperature( t_solver );
  NumIntegral_dNT_op_dN_dV<dim>   t_conductance( *pd_ref,
                  "thermal conductivity",
                  "temperature",
                  "temperature" );
  NumIntegral_SetRHS_to_Zero<dim> t_rhs(*pd_ref, "temperature");

  steady_state_temperature.Add(&t_conductance);
  steady_state_temperature.Add(&t_rhs);

  (*model).Apply (steady_state_temperature);

  return 0;

}

int CVFEM_1D_VVCase::InitializeFluidPropertiesLinearPressure()
{
  enum{dim=1}; //if you change here you must change everywhere

  //! calculate and assign initial pressure gradient
  double p_value;
  ScalarVariable P_value;

  csmp::Index   p_key(model->Database().StorageKey("fluid pressure"));
  csmp::Index   wt_key(model->Database().StorageKey("salinity"));

  double p_bottom, p_top, p_diff, ymax;
  p_bottom = p_top = p_diff = ymax = std::numeric_limits<double>::quiet_NaN();

  const typename vector<Node<dim>*>::const_iterator modelNodesEnd( model->Region("Model").NodesEnd() );
  for( typename vector<Node<dim>*>::const_iterator it( model->Region("Model").NodesBegin() );
       it != modelNodesEnd; ++it )
  {
      if ((*it)->AtBoundary() == CNR1)
        {
          (*it)->Read(p_key, P_value);
          p_top = P_value();
        }
      if ((*it)->AtBoundary() == CNR2)
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
          p_value = p_top + (*it)->y() / ymax * p_diff; // FIXME: different from philipp's code
          (*it)->Store(p_key, makeScalar(PLAIN, p_value));
      }

  }

  //CalculateSteadyStateTemperature();

  //! Initialize fluid properties from PTX
CVFEM_PHX->InitialFluidPropertiesFromPTX();

  //! initial output
  //cout<<"\n Output Initial Properties";
  //vtu->OutputDataToVTU( ( string(config_file_name_) + "_Initial_Properties" ).c_str(),
  //                          output_props_, model->Region("Model"), 0);

  return 0;
}

int CVFEM_1D_VVCase::InitializeFluidProperties()
{
  enum{dim=1}; //if you change here you must change everywhere

  for (int i=0; i<5; i++)
  {
    //! Calculate steady state pressure with gravity
    CalculateSteadyStatePressure();
    //CalculateSteadyStateTemperature();

    //! Initialize fluid properties from PTX
    CVFEM_PHX->InitialFluidPropertiesFromPTX();

    vtu->OutputDataToVTU( ( string(config_file_name_) + "_Initial_Properties_loop" ).c_str(),
                            output_props_, model->Region("Model"), i);
  }

  //! initial output
  cout<<"\n Output Initial Properties";
  vtu->OutputDataToVTU( ( string(config_file_name_) + "_Initial_Properties" ).c_str(),
                            output_props_, model->Region("Model"), 0);

  return 0;
}

//This runs the CVFEM scheme
void CVFEM_1D_VVCase::run()
{
  enum{dim=1}; //if you change here you must change everywhere

  //! Calculate linear pressure gradient, initialize fluid properties from PTX, initial output to vtu
  InitializeFluidPropertiesLinearPressure();
  //cin.get();

  //! setup time variables
  model_time = 0.;
  vtu_dt = run_settings.NearestOutputTime(0.0);
  max_time = run_settings.Duration();
  timestep = 0;
  output_counter = 0;

  //! main time loop
  while (model_time <= max_time)
  {
    dt = CVFEM_PHX->Apply();
    cout<<"\ndt = "<<dt;

    model_time += dt;
    cout<<"\n model time = "<<model_time<<" seconds"<<endl;

    //! output to vtk
    if (model_time >= (output_counter + 1)*vtu_dt)
    {
      vtu->OutputDataToVTU( ( string(config_file_name_) + "_Properties" ).c_str(),
                                   output_props_, model->Region("Model"), output_counter);
      output_counter++;
      //cin.get();
    }

  timestep++;
  }

  std::cout << "Finished Running CVFEM_1D_VVCase " << std::endl;

}


}//end csmp
