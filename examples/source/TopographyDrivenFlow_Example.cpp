#include "TopographyDrivenFlow_Example.h"

#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"

#include "TRIANGLE_Interface.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "LinearSolver.h"

#include "PropertyHandle.h"

// Interrelations
#include "Concatenate.h"
#include "ConstantFactor.h"
#include "HydraulicHead.h"
#include "GroundwaterDarcyVelocity.h"

// Data output to Visualization Toolkit
#include "TextInterface.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp{

void TopographyDrivenFlow_Example::Specifications()
{
  SetTitle( "Topography-driven flow" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "Topography-driven flow in a cross-sectional model for a steady-state fluid pressure distribution" );
  AddDescription( "source in: TopographyDrivenFlow_Example.cpp, example16.txt(variable file)" );
  AddRequirement( "file set: 'topo.1'");
}

// TODO: use correct terminology (hydraulic head etc.) instead of fluid pressure
void TopographyDrivenFlow_Example::Run()
{
  /* **********************************************************************

     Use 'topo.1' file set as input geometry.

     Topography driven flow steady-state pressure distribution.

     - use 'topo.1' file set as input model. It constitutes a cross-section
       through a valley with a fault zone outcropping at the valley floor
       which is offsetting an aquifer horizon.

  *********************************************************************** */
    /*
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;

  // ------------------------------------------------------------
  // 0. building Region object
  // ------------------------------------------------------------
    char  file_name[200];
    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
    cin >> file_name;
    mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );
    Model<2U>  model( mesh_container, "example16.txt" );
    */

    // ------------------------------------------------------------
    // 0. Load CSMP native format model
    // ------------------------------------------------------------
    string model_name;
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'topo.1':"<<endl;
    cin.ignore();
    getline(cin, model_name);
    if (model_name.length() == 0) model_name = "topo.1";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "example16.txt";
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    Model<2U>  model(model_name, variable_file);

    printModelDimensions( model, true );

  // ---------------------------------------------------------------
  // 1. assign initial values & Dirichlet boundary conditions
  // ---------------------------------------------------------------
    model.InputPropertyValue ( "fluid density",       makeScalar(PLAIN,1000.) );
    model.InputPropertyValue ( "fluid volume source", makeScalar(PLAIN,0.) );

    // atmospheric fluid pressure is assigned to earth surface
    const double  p_atm(101325.);
    model.InputBoundaryValue( TOP, "absolute fluid pressure", makeScalar(DIRICH,p_atm) );


  // --------------------------------------------------------------
  // 2. calculate hydraulic conductivity from permeability
  // --------------------------------------------------------------
    const double  fluid_viscosity(1.6e-3); // Pa s-1
    ConstantFactor<2U,divides>  conductivity( model.Database(),
                                             "conductivity", "permeability",
                                              fluid_viscosity );
    model.Apply( conductivity );
    printRangeOfVariable( model, "conductivity" );


  // -----------------------------------------------------------------------------
  // 3. compute absolute fluid pressure taking into account topography
  // -----------------------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    PDE_Integrator<2U,Element>  steady_state_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<2U,Element>  steady_state_pressure(linear_solver);
#endif

    // conductance matrix [K] on the lefthand side
    NumIntegral_dNT_lhsop_dN_dV<2U>    conductance( model.Database(), "conductivity", "absolute fluid pressure",  "absolute fluid pressure" );
    // source vector {Q} on the righthand side
    NumIntegral_NT_rhsop_N_dV<2U>   source( model.Database(),  "fluid volume source", "absolute fluid pressure" );
    // gravity term on righthand side
    NumIntegral_NT_op_dNi_dV<2U>    gravity( model.Database(),  "fluid density",  "conductivity", "absolute fluid pressure" );

    // add PDE_Operators to the FE Algorithm
    steady_state_pressure.Add( &conductance );
    steady_state_pressure.Add( &source );
    steady_state_pressure.Add( &gravity );

    model.Apply( steady_state_pressure );


  // -----------------------------------------------------------------------------
  // 4. assigning the Y-coordinate value to the variable 'elevation'
  // -----------------------------------------------------------------------------
    model.AssignNodeCoordinatesTo( "elevation", 'y' );


  // ------------------------------------------------------------
  // 5. calculate hydraulic head, Darcy velocity, and scalar flux
  // ------------------------------------------------------------
    HydraulicHead<2U>  head( model.Database() );
    model.Apply( head );
    model.CopyGradientOfProperty_A_To_B( "hydraulic head", "hydraulic head gradient" );

    GroundwaterDarcyVelocity<2U>  velo( model.Database() );
    model.Apply( velo );


  // ------------------------------------------------------------
  // 6. Output of results
  // ------------------------------------------------------------
    TextInterface  text_output;

    // outputting the logarithm of conductivity (this modifies the original variable)
    PropertyHandle<2U>  K( model, "conductivity" );
    K.Log10();

    text_output.OutputDataAsTextColumns( model, "velocity",       "velocity",       1 );
    text_output.OutputDataAsTextColumns( model, "conductivity",   "conductivity",   1 );
    text_output.OutputDataAsTextColumns( model, "hydraulic-head", "hydraulic head",   1 );
    text_output.OutputDataAsTextColumns( model, "absolute-fluid-pressure", "absolute fluid pressure",   1 );


    VTK_Interface<2U>  vtk_output;

    vtk_output.OutputDataToVTK( model, "conductivity",            "conductivity", 1 );
    vtk_output.OutputDataToVTK( model, "velocity",                "velocity", 1 );
    vtk_output.OutputDataToVTK( model, "hydraulic-head",          "hydraulic head", 1 );
    vtk_output.OutputDataToVTK( model, "absolute-fluid-pressure", "absolute fluid pressure", 1 );
  
    // nodal velocities for streamlines (however inaccurate at material interfaces)
    model.ExtrapolateCellToNodeProperty( "velocity", "nodal velocity" );
    vtk_output.OutputDataToVTK( model, "nodal-velocity", "nodal velocity", 1 );

    cout <<"\nmain: That's it..."<< endl;

    filesystem::current_path("../../example_inputs/");

} // Run()

} // csmp
