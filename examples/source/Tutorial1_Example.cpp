#include "Tutorial1_Example.h"

// the CSMP model
#include "Model.h"
#include "ModelTime.h"

// a simple FE mesh generator
#include "Triangulator.h"

// the FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
#include "Integral_NT_op_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

// output interfaces
#include "VTK_Interface.h"
#include "MatlabInterface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "ConstantFactor.h"

using namespace std;

namespace csmp{

void Tutorial1_Example::Specifications()
{
  SetTitle( "Tutorial 1: Transient fluid-pressure diffusion" );
  SetDifficulty( 3 );
  SetCategory( "Tutorials (composite functionality)" );
  AddAuthor( "Sebastian Geiger" );
  AddDescription( "A transient pressure diffusion equation solved using a fully" );
  AddDescription( "implicit FE discretisation. A simple triangular FE mesh is generated automatically in CSMP" );
  AddDescription( "from which the Model is built. Boundary and initial conditions are applied and the" );
  AddDescription( "hydraulic conductivity is computed as a function of the permeability and viscosity. Output");
  AddDescription( "is written to VTK and Matlab files.");
  AddRequirement( "tutorial1_input, tutorial1_variables.txt");
} // Initialize()

// **********************************************************************************************
//
// A simple CSMP file that solves the transient pressure diffusion equation using a fully
// implicit FE discretisation. A simple triangular FE mesh is generated automatically in CSMP
// from which the Model is built. Boundary and initial conditions are applied and the
// hydraulic conductivity is computed as a function of the permeability and viscosity. Output
// is written to VTK and Matlab files.
//
// Tasks and exercises:
//
// 1. Generate different text-input files with varying permeabilities
//    (see User's Guide Page 45 "Regular Meshes from Pixel Data")
//    NB: Uncomment the line where a fixed, uniform permeability is set!
// 2. Visualise results and observe how pressure diffusion varies for models with
//    and without permeability contrast
// 3. Vary boundary and initial conditions (e.g., fixed pressure on two sides vs. one side,
//    change compressibility, etc.)
// 4. Visualise results and observe what happens
// 5. Generate a simple mesh that is only 1 FE high, use a uniform permeability and apply
//    boundary conditions such that you can compare numerical results with your analytical
//    solution for the diffusion equation
// 6. Write a new function that reads in porosity and computes the compressibility as
//    a function of the user-specified fluid- and rock-compressibilities, i.e. solves
//    c = phi * c_f + (1-phi) * c_r for each FE with phi = porosity, c_f = fluid compressibility
//    and c_r = rock compressibility. Look at other classes in /source_code/interrelations for
//    inspiration.
//
// **********************************************************************************************

void Tutorial1_Example::Run()
{
    double& model_time( ModelTime::Instance().modelTime );
    model_time =  0.; // time

    /*
    // -----------------------------------------------------------------------
    // 1.0 Build the CSMP Model from a simple color-coded text input file
    //     that generates a uniform triangular FE mesh. The function
    //     readTextPixelData() reads in the file, generates the FE mesh
    //     (of uniform triangles) and returns it as a VSet
    // -----------------------------------------------------------------------
  
    // txt file defines the physical variables to be used in the simulation
    VSet<2U>   vset=readTextPixelData();
    Model<2U>  model( vset, "tutorial1_variables.txt" );
    */

    // ------------------------------------------------------------
    // 1.0 Load CSMP native format model
    // ------------------------------------------------------------
    string model_name;
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'tutorial1_input':"<<endl;
    cin.ignore();
    getline(cin, model_name);
    if (model_name.length() == 0) model_name = "tutorial1_input";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "tutorial1_variables.txt";
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    Model<2U>  model(model_name, variable_file);

    const PropertyDatabase<2>& p_ref(model.Database());  // constant reference to the property database

    // give the model dimensions
    printModelDimensions( model, true );

    // -----------------------------------------------------------------------
    // 2.0 Now we apply boundary and initial conditions (this can also be done,
    //     more conveniently, in a configuration file for more realistic runs)
    // -----------------------------------------------------------------------

    // assigning material properties
    model.InputPropertyValue( "porosity",         makeScalar(PLAIN,0.1) );     // always as a fraction
    model.InputPropertyValue( "permeability",     makeScalar(PLAIN,1.0e-15) ); // always in m2 (comment out if heterogeneous k-field is used in input file)
    model.InputPropertyValue( "compressibility",  makeScalar(PLAIN,5.0e-10) ); // for fluid and rock, in Pa-1

    // assigning initial conditions
    model.InputPropertyValue( "fluid pressure",      makeScalar(PLAIN,1.0e+07) );  // always in Pascal
    model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.0) );      // no sources/sinks (units m3 m-2 s-1)

    // assigning boundary conditions for fluid pressure at the LEFT and RIGHT model boundaries
    // such that a pressure wave travels from left to right through the model
    model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,3.0e+07) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.0e+07) );


    // ----------------------------------------------------------------------------------------------
    // 3.0 Now we use an Interrelation (ConstantFactor, inherited from base class Interrelation)
    //     to compute the hydraulic conductivity K = k/mu (k = permeability, mu = viscosity) at each
    //     finite element
    // ----------------------------------------------------------------------------------------------
    ConstantFactor<2U,divides>  conductivity( p_ref, "conductivity", "permeability", 0.001 ); // viscosity 1 cp = 0.001 Pa s

    // the Model applies the object "conductivity", which is instantiated from class ConstantFactor
    // the result variable "conductivity" is computed automatically and its range is checked
    model.Apply( conductivity );

    // output the range of the result variable
    printRangeOfVariable( model, "conductivity" );

    // ------------------------------------------------------------------------------------------
    // 4.0 Setting up an FE algorithm to solve the diffusion equation c dp/dt = div(K grad p) + S
    //     p = fluid pressure
    //     c = compressibility (fluid and rock)
    //     K = k/mu = hydraulic conductivity (from above)
    //     S = volumetric source term
    //
    //     We solve the discretised equation full implict as
    //
    //     ([c]/dt + [K]){p}t+dt = {c}/dt{p}t + {S}t+dt
    //
    //     Note: [] denotes a matrix, {} a vector
    //
    //     This results in the linear system [A] * {x} = {b}
    //     where [A] is the discretisation of div(K grad p) and c dp/dt
    //     {b} contains the known pressure at time t and the unknown source at
    //     time t+dt; {x} is the unknown pressure at time t+dt that we are solving for
    //
    // ------------------------------------------------------------------------------------------

    // create the CSMP FE Algorithm with SAMG solver to invert linear system
#ifdef USE_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    PDE_Integrator<2U,Region>  fluid_pressure(samg_solver);
#else
    EigenSolver  linear_solver;
    PDE_Integrator<2U,Region>  fluid_pressure(linear_solver);
#endif

    // LHS stiffness matrix                              operand         basis function    test function
    Integral_dNT_op_dN_dV<2U,Element<2U> >  stiffness_matrix( p_ref, "conductivity", "fluid pressure", "fluid pressure" );

    // LHS mass matrix
    Integral_NT_lhsop_N_dV<2U,Element<2U> > mass_matrix_lhs( p_ref, "compressibility", "fluid pressure", "fluid pressure" );

    // RHS mass vector
    Integral_NT_op_N_dV<2U,Element<2U> >    mass_matrix_rhs( p_ref, "compressibility", "fluid pressure" );

    // RHS mass vector for integrating source term
    Integral_NT_op_N_dV<2U,Element<2U> >    source_term( p_ref, "fluid volume source", "fluid pressure" );

    // mass matrices for dp/dt term must be divided by time increment
    mass_matrix_lhs.MultiplyWithTimeIncrement(true);
    mass_matrix_rhs.MultiplyWithTimeIncrement(true);

    // use lumped formulation for all mass matrices (i.e., diagonalise matrices)
    mass_matrix_lhs.LumpedFormulation(true);
    mass_matrix_rhs.LumpedFormulation(true);
    source_term.LumpedFormulation(true);

    // evalute source term last
    source_term.AddAccumulateLater();

    // define a post-processing step that computes the velocity in each finite element by solving Darcy's law
    VelocityAndVolumeFlux<2U,Element<2U> >     velo( model, "conductivity", "porosity", "fluid pressure", true ); // true = extrapolate element velocities to nodes

    // now add each FE operation (i.e., PDE Operator) to the FE algorithm
    fluid_pressure.Add( &stiffness_matrix );
    fluid_pressure.Add( &source_term );
    fluid_pressure.Add( &mass_matrix_lhs );
    fluid_pressure.Add( &mass_matrix_rhs );
    fluid_pressure.AddPostProcess( &velo );


    // -----------------------
    // 5.0 Time Loop Variables
    // -----------------------
    // write output in VTK format and for Matlab
    VTK_Interface<2U>  vtk_output;
    MatlabInterface    matlab;

    vtk_output.OutputDataToVTK( model, "fluid_pressure", "fluid pressure", 0 );
    matlab.Write2DMatlabFile(   model, "fluid_pressure", "fluid pressure", 0 );

    // define some constant variables
    const double     hour(3600.0);
    const double     max_time(240.0*hour); // run for 10 days
    double           time_increment(2.0*hour); // timestep 2 hours
    const long         save_frequency(6); // write results to file every 12 hours
    size_t	           save_counter(1), time;

    // set the time increment for the FE algorithm
    fluid_pressure.TimeIncrement( 1.0/time_increment );


    // -----------------------
    // 6.0 Transient loop
    // -----------------------
    while ( model_time < max_time )
      {
         // pressure diffusion is computed as the Model applies the FE algorithm
         model.Apply( fluid_pressure );

         // show results
         printRangeOfVariable( model, "fluid pressure" );
         printRangeOfVariable( model, "velocity" );

         // increment time
         model_time += time_increment;

         // output variables
         if ( save_counter == save_frequency ) {
              time = static_cast<long>(model_time/hour);
              // to VTK files
              vtk_output.OutputDataToVTK( model, "fluid_pressure", "fluid pressure", time );
              vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       time );
              vtk_output.OutputDataToVTK( model, "volume_flux",    "volume flux",    time );
              // to Matlab files
              matlab.Write2DMatlabFile(  model, "fluid_pressure", "fluid pressure",    time );
              matlab.Write2DMatlabFile(  model, "volume_flux",    "nodal volume flux", time );
              save_counter = 0;
          }
         save_counter++;

         // runtime info
         cout <<"\n\nmain: RUNTIME (HRS): "<< model_time/hour << endl << endl;
      }

    // final output
    // VTK
    time = static_cast<long>(model_time/hour);
    vtk_output.OutputDataToVTK( model, "fluid_pressure", "fluid pressure", time );
    vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       time );
    vtk_output.OutputDataToVTK( model, "flux",           "volume flux",    time );
    // Matlab
    matlab.Write2DMatlabFile(  model, "fluid_pressure", "fluid pressure",    time );
    matlab.Write2DMatlabFile(  model, "volume_flux",    "nodal volume flux", time );

    // terminate
    cout <<"\nmain: That's it..."<< endl;

    fs::current_path("../../example_inputs/");

} // Run

} // csmp
