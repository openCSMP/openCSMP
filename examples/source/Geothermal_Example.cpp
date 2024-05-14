#include "Geothermal_Example.h"

#include "ModelTime.h"
#include "GlobalVerbose.h"
#include "Model.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "ANSYS_Model2D.h"
#include "SAMG_Solver.h"
#include "CSMP_highLevelUtilities.h"
#include "LUdcmp_Solver.h"
#include "ComputationalSettings.h"

// finite volumes 
#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"

// finite elements
#include "PDE_Integrator.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "PointSource_rhsop.h"
#include "VelocityAndVolumeFlux.h"

// visitors
//#include "ComputeGravityTermVisitor.h"
#include "ThermalVisitor.h"
#include "SourceVisitor.h"
#include "ConductivityVisitor.h"  /// @todo change name to something more specific
#include "PropertyAtPointVisitor.h"

// utilities and monitoring
#include "InputDataManager.h"


using namespace std;

namespace csmp {

void Geothermal_Example::Specifications()
{
  SetTitle( "Geothermal circulation in box-shaped model" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "Alina Yapparova" );
  AddDescription( "Heat transport with pressure-enthalpy-temperature scheme, see Yapparova et al. (2014). Realistic simulation of an aquifer thermal \
energy storage: Effects of injection temperature, well placement and groundwater flow. \
Energy, 76, 1011-1018." );
  AddDescription( "Comparison with TOUGH" );
  AddDescription( "source in: Geothermal_Example.cpp" );
  AddRequirement( "input file set: 2000x1000_mesh(csmp binary files)" );
  AddRequirement( "variable file (Geothermal_pseudo1D_VVCase.txt), Geothermal_Example-configuration.txt" );
  AddRequirement( "H2OPropertiesLookupTableSinglePhase.bin");
}


/**
@todo: SKM: fix configuration file, it requests face properties that are not present in variables file
@todo: SKM: add gravity to the calculations

*/
void Geothermal_Example::Run()
{
  bool& globalVerbose( GlobalVerbose::Instance().globalVerbose );

  //! set as desired
  globalVerbose = true;

  /*
  //! -------------------------------------------------------
  //! 1. Model construction - modify to use alternative model
  //! -------------------------------------------------------
  const std::string geometry_name( "2000x1000_mesh" );
  const std::string regions_name( "2000x1000_mesh" );
  const std::string config_name( "Geothermal_Example" );
  const std::string vars_name( "Geothermal_Example-variables.txt" );

  ANSYS_Model2D model( geometry_name.c_str(), regions_name.c_str(), vars_name.c_str(), false, true, true );
   */


  // ------------------------------------------------------------
  // 1. Load CSMP native format model
  // ------------------------------------------------------------
  string model_name;
  cout<< "\nPlease enter the name of input model, or press ENTER to use the default model '2000x1000_mesh':"<<endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "2000x1000_mesh";

  //find the name of current example source file
  string file_name = GetExampleFileName(__FILE__);
  string variable_file = "Geothermal_pseudo1D_VVCase.txt";
  string config_name ="Geothermal_pseudo1D_VVCase";
  //create a working directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_name);
  //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model<2U>  model(model_name, variable_file);
  model.InstantiateFiniteVolumes();

  const PropertyDatabase<DIM>&  pd_ref( model.Database() ); //reference to the models property database.
  printModelDimensions<DIM>( model, true );

  //copy in H2O Properties Lookup Table (H2OPropertiesLookupTableSinglePhase.bin) as it takes too long to compute
  string path = "../../example_inputs/variables_and_configuration_files/";
  string name = "H2OPropertiesLookupTableSinglePhase.bin";
  file_name = path + name;
  if (filesystem::exists(file_name)) filesystem::copy(file_name, "./");
  else {
    string error_message = "\n\nError: file '";
    string input_directory = (filesystem::current_path().parent_path().parent_path()).string();
    input_directory += "/example_inputs/variables_and_configuration_files/";
    error_message += (name + "' does not exist in directory "  + input_directory);
    error_message += (", example cannot run, please copy this file into this directory\n");
    throw std::runtime_error(error_message);
  }


  /*
  // Finite Element Mesh Construction
  string name = "Geothermal_pseudo1D_VVCase";
  std::string geometry_name ("2000x1000_mesh");
  std::string regions_name (name);
  std::string config_name (name);
  std::string vars_name (name +".txt");

  ANSYS_Model2D model(geometry_name.c_str(), regions_name.c_str(), vars_name.c_str() , false, true );
  model.InstantiateFiniteVolumes();
  const PropertyDatabase<DIM>&  pd_ref(model.Database()); //reference to the models property database.

  printModelDimensions<DIM>(model, true );
  */

  //! Assignment of material properties, initial conditions, and boundary
  InputDataManager<DIM>  model_configuration;
  ComputationalSettings  run_settings;

  // boolean flags set reading to: 2) default prop.values, 3) group prop.values, 4) essential conditions for box-shaped model
  model_configuration.ConfigureFromFile( model, config_name.c_str(), false,    // groupname from parameter range
                                         true,    // default property values
                                         true,    // regional property values
                                         true,    // boundary conditions for box-shaped model
                                         true,    // essential conditions for groups
                                         false,    // csmp::Boundary properties
                                         run_settings );

  //! node visitor for computation of equation of state and transport properties
  ThermalVisitor<DIM>   thermal_equilibrator( model );
  //! element visitor to calculate thermal contraction or expansion related source/sink terms
  SourceVisitor<DIM>    source_calculator( model );

  //! transform rock properties from the nodes to the elements (why is this needed?)
  model.ExtrapolateCellToNodeProperty( "porosity", "nodal porosity" );
  model.ExtrapolateCellToNodeProperty( "density rock", "nodal density rock" );
  model.ExtrapolateCellToNodeProperty( "heat capacity rock", "nodal heat capacity rock" );
  model.ExtrapolateCellToNodeProperty( "compressibility rock", "nodal compressibility rock" );

  //! to quickly replace SAMG with LU solver when necessary, uncomment
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver solver;
#else
  EigenSolver solver;
#endif

  //! --------------------------------------------------------
  //! 2. Setting up finite element algorithms for PT diffusion
  //! --------------------------------------------------------
  //! computation of initial steady-state fluid pressure
  PDE_Integrator<DIM,Element>  steady_state_pressure( solver );

  NumIntegral_dNT_op_dN_dV<DIM>    p_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
  // replace if you want to work with fluid sources and sinks
  NumIntegral_SetRHS_to_Zero<DIM>  zero_fluid_src( pd_ref, "fluid pressure" );

  VelocityAndVolumeFlux<DIM>  velocity( model, "conductivity", "porosity", "fluid pressure", false );

  steady_state_pressure.Add( &p_conductance );
  steady_state_pressure.Add( &zero_fluid_src );
  steady_state_pressure.AddPostProcess( &velocity );


  //! finite element computation of transient pressure
  PDE_Integrator<DIM,Element>  transient_pressure( solver );

  NumIntegral_dNT_op_dN_dV<DIM>  pt_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
  pt_conductance.MultiplyWithTimeIncrement( true );

  NumIntegral_NT_lhsop_N_dV<DIM> pt_capacitance_lhs( pd_ref, "total compressibility", "fluid pressure", "fluid pressure" );
  pt_capacitance_lhs.LumpedFormulation( true );

  NumIntegral_NT_op_N_dV<DIM>    pt_capacitance_rhs( pd_ref, "total compressibility", "fluid pressure" );
  pt_capacitance_rhs.LumpedFormulation( true );

  //! nodal fluid volume source for mass correction term
  PointSource_rhsop <DIM> fluid_src( pd_ref, "nodal fluid volume source", "fluid pressure" );
  fluid_src.MultiplyWithTimeIncrement( false );
  fluid_src.LumpedFormulation( true );
  fluid_src.AddAccumulateLater();

  VelocityAndVolumeFlux<DIM>  t_velocity( model, "conductivity", "porosity", "fluid pressure", false );

  transient_pressure.Add( &pt_conductance );
  transient_pressure.Add( &pt_capacitance_lhs );
  transient_pressure.Add( &pt_capacitance_rhs );
  transient_pressure.Add( &fluid_src );
  transient_pressure.AddPostProcess( &t_velocity );

  //! transient thermal diffusion
  PDE_Integrator<DIM,Element>  temperature_diffusion( solver );

  NumIntegral_dNT_op_dN_dV<DIM>   t_conductance( pd_ref, "thermal conductivity", "temperature", "temperature" );

  NumIntegral_NT_lhsop_N_dV<DIM>  t_capacitance_lhs( pd_ref, "total heat capacity", "temperature", "temperature" );
  t_capacitance_lhs.LumpedFormulation( true );
  t_capacitance_lhs.MultiplyWithTimeIncrement( true );

  NumIntegral_NT_op_N_dV<DIM>  t_capacitance_rhs( pd_ref, "total heat capacity", "temperature" );
  t_capacitance_rhs.LumpedFormulation(true); //added
  t_capacitance_rhs.MultiplyWithTimeIncrement( true );

  temperature_diffusion.Add( &t_conductance );
  temperature_diffusion.Add( &t_capacitance_lhs );
  temperature_diffusion.Add( &t_capacitance_rhs );


  //! ----------------------------------------------------------------------
  //! 3. Setting up finite volume algorithms for mass and enthalpy advection
  //! ----------------------------------------------------------------------
  ExplicitMassBasedTransport<DIM, MassBasedStencilProcessor> mass_advection( "Model", model,
                                                                             "porosity",          //porosity    
                                                                             "fluid density",     //advected property lhs 
                                                                             "density liquid",    //advected property rhs 
                                                                             "velocity",          //transport velocity
                                                                             "nodal heat source", //nodal source
                                                                             false );              //second order in space

  ExplicitMassBasedTransport<DIM, MassBasedStencilProcessor> enthalpy_advection( "Model", model,
                                                                                 "porosity",                    //porosity    
                                                                                 "enthalpy content liquid",     //advected property lhs  
                                                                                 "volumetric enthalpy liquid",  //advected property  rhs
                                                                                 "velocity",                    //transport velocity
                                                                                 "nodal heat source",           //nodal source
                                                                                 false );                        //second order in space


  //! ------------------------------------------------------------------
  //! 4. Pressure, temperature, mass and enthalpy transfer computation
  //! ------------------------------------------------------------------
  //! time
  double max_time = 3.15e9 * 30; //3000 years
  double time_increment = 3.15e7; //1 year
  size_t time_step = 0;

  //! definition of output properties - comment out if not wanted
  list<string> output_props;
  list<string> region_names;

  output_props.push_back( "fluid pressure" );
  output_props.push_back( "velocity" );
  output_props.push_back( "temperature" );
  output_props.push_back( "fluid density" );
  output_props.push_back( "density liquid" );
  output_props.push_back( "volumetric enthalpy liquid" );
  output_props.push_back( "enthalpy content liquid" );
  output_props.push_back( "nodal fluid volume source" );
  output_props.push_back( "total compressibility" );
  output_props.push_back( "total heat capacity" );

  region_names.push_back( "Model" );


  //! Model-wide node-property initialisation including H2O-EOS calculation of fluid density etc.
  thermal_equilibrator.SetInitialProperties( &model );
  //! extrapolation of nodal heat capacity, compressibility and density to the element
  model.InterpolateNodeToCellProperty( "nodal total heat capacity", "total heat capacity" );
  model.InterpolateNodeToCellProperty( "nodal total compressibility", "total compressibility" );
  model.InterpolateNodeToCellProperty( "density liquid", "density liquid element" );

  //! computation of the hydraulic conductivity K=rho k / mu
  ComputeMassConductivity( model );

  printRangeOfVariable( model, "conductivity" );
  printRangeOfVariable( model, "mass conductivity" );
  printRangeOfVariable( model, "porosity" );
  printRangeOfVariable( model, "fluid pressure" );

  //! initializing pressure in the model followed by a repeated fluid and rock property calculation
  model.Apply( steady_state_pressure );

  thermal_equilibrator.SetInitialProperties( &model );
  model.InterpolateNodeToCellProperty( "nodal total heat capacity", "total heat capacity" );
  model.InterpolateNodeToCellProperty( "nodal total compressibility", "total compressibility" );
  model.InterpolateNodeToCellProperty( "density liquid", "density liquid element" );

  //! setting mass transfer and enthalphy (advected properties) to Dirich at those boundaries where p, t are Dirichlet
  model.Boundary( "LEFT" ).ChangePropertyStatus( "fluid density", DIRICH );
  model.Boundary( "RIGHT" ).ChangePropertyStatus( "fluid density", DIRICH );
  model.Boundary( "LEFT" ).ChangePropertyStatus( "enthalpy content liquid", DIRICH );
  model.Boundary( "RIGHT" ).ChangePropertyStatus( "enthalpy content liquid", DIRICH );

  //! initial output of model state into timestep 0
  OutputToVTU( model, config_name, output_props, time_step );
  time_step++;

  //! ------------------------------------------------------------------
  //! 5. PT-equlibration loop
  //! ------------------------------------------------------------------
  double& global_time( ModelTime::Instance().modelTime );
  global_time = 0.;

  while ( global_time <= max_time )
  {
    //! thermal diffusion
    temperature_diffusion.TimeIncrement( 1. / time_increment );
    model.Apply( temperature_diffusion );

    //! enthalphy transport
    double time_advection( 0. );
    while ( time_advection < time_increment )
    {
      double time_increment_advection( mass_advection.AnisotropicCourantIncrement() );
      if ( time_advection + time_increment_advection > time_increment )
        time_increment_advection = time_increment - time_advection;

      mass_advection.AdvectVariable( time_increment_advection );
      enthalpy_advection.AdvectVariable( time_increment_advection );

      time_advection += time_increment_advection;
    }

    //! update of PT dependent properties
    model.Accept( thermal_equilibrator );
    model.Accept( source_calculator );
    model.InterpolateNodeToCellProperty( "nodal total heat capacity", "total heat capacity" );
    model.InterpolateNodeToCellProperty( "nodal total compressibility", "total compressibility" );
    model.InterpolateNodeToCellProperty( "density liquid", "density liquid element" );
    ComputeMassConductivity( model );

    //! pressure diffusion
    transient_pressure.TimeIncrement( time_increment );
    model.Apply( transient_pressure );

    //! output of results
    if ( time_step % 100 == 0 )
      OutputToVTU( model, config_name, output_props, time_step );

    if ( time_step == 1600 ) {
      if ( Compare( model, "Geothermal_Example_Comparison_TOUGH_data_1600yr.txt" ) )
        cout << "\n Comparison OK, timestep =  " << time_step;
    }
    else if ( time_step == 1900 ) {
      if ( Compare( model, "Geothermal_Example_Comparison_TOUGH_data_1900yr.txt" ) )
        cout << "\n Comparison OK, timestep =  " << time_step;
    }
    else if ( time_step == 3000 ) {
      if ( Compare( model, "Geothermal_Example_Comparison_TOUGH_data_3000yr.txt" ) )
        cout << "\n Comparison OK, timestep =  " << time_step;
    }

    global_time += time_increment;
    time_step++;
  }

  filesystem::current_path("../../example_inputs/");

} // end run



/**
computes hydraulic conductivity - density product, K = rho A k / mu  (kg s-1),
where mu is the dynamic viscosity of water and, A and k and the flow-cross-sectional area and the
intrinsic permeability, respectively.
*/
void Geothermal_Example::ComputeMassConductivity( Model<DIM>& model )
{
  PropertyHandle<DIM> rho_ph( model, "density liquid", SCALAR, NODE );
  PropertyHandle<DIM> kappa_ph( model, "conductivity", SCALAR, ELEMENT );
  PropertyHandle<DIM> lambda_ph( model, "mass conductivity", SCALAR, ELEMENT );

  ConductivityVisitor<DIM> conductivity_visitor( model, "conductivity", "permeability", "fluid viscosity" );
  model.Accept( conductivity_visitor );

  lambda_ph = 0.;
  lambda_ph += rho_ph;
  lambda_ph *= kappa_ph;
}



/**
shorthand for the output of multiple variables to VTU
*/
void Geothermal_Example::OutputToVTU( Model<DIM>& model, string model_name, const list<string>& props, size_t timestep ) const
{
  static VTU_Interface<DIM> vtu( model );
  vtu.OutputDataToVTU( (model_name + "_Properties").c_str(), props, model.Region( "Model" ), timestep );
}



/**
Comparison between CSMP and TOUGH results.

@todo replace dark ages C code and add a bit more flexibility
*/
bool Geothermal_Example::Compare( Model<DIM>& model, const string& file ) const
{
  FILE* fin;
  double x, y, val, res, tol;
  size_t ind( 0 );
  size_t flag( 0 );
  size_t nPoints;
  ScalarVariable result( PLAIN, 0.0 );
  map<size_t, std::vector<double> > points;
  vector<double>  values;

  fin = fopen( file.c_str(), "r" );
  if ( fin == NULL ) {
    printf( "\n File could not be opened" );
    return false;
  }

  // TOUGH simulation results
  // reads (x_coord, y_coord, temperature) data from a text file
  while ( flag != EOF ) {
    flag = fscanf( fin, "%lf\t%lf\t%lf\n", &x, &y, &val );
    points[ind].push_back( x );
    points[ind].push_back( y );
    values.push_back( val );
    ind++;
  }

  fclose( fin );

  nPoints = points.size();

  // finds CSMP computed temperature values at given points
  PropertyAtPointVisitor<DIM> pAt( model, points, "temperature" );
  model.Accept( pAt );

  res = 0.;
  tol = 1.;

  // computes C-norm of the difference (sum of the absolute values)
  for ( size_t ind2 = 0; ind2<nPoints; ind2++ ) {
    pAt.PropertyValueAt( ind2, result );
    res += fabs( result() - values[ind2] );
  }

  // divided by the number of points
  res /= nPoints;

  //printf("\n res = %f", res);

  // if residual is smaller than tolerance, test passed
  if ( res < tol )
    return true;
  else
    return false;
}


} //end csmp

