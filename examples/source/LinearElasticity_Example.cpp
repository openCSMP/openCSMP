#include "LinearElasticity_Example.h"

#include "Region.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "meshManagementUtilities.h"
#include "CSMP_definitions.h"

// File I/O and Initialization
#include "TRIANGLE_Interface.h"
#include "VTK_Interface.h"
#include "FiniteElementManager.h"
#include "InputDataManager.h"
#include "Standard_IO_Handler.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "PropertyHandle.h"

// Interrelations
#include "ExtractVectorVariableComponent.h"
#include "ExtractTensorVariableComponent.h"
#include "ConstantFactor.h"

#include "PDE_Integrator.h"
#include "SteadyStateDiffusor.h"

// PDE Operators
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "PT_op.h"
#include "StressesAndStrains2.h"

#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "VelocityAndVolumeFlux.h"

// comment-in this preprocessor directive if you want to run directly with the Triangle mesher input file deck
// #define RUN_DIRECTLY_WITH_TRIANGLE_INPUT_FILE

using namespace std;

namespace csmp {

void LinearElasticity_Example::Specifications()
{
  SetTitle( "Linear Elasticity" );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "Linear elasticity with fractures and regional computations" );
  AddDescription( "source in: LinearElasticity_Example.cpp" );
  AddRequirement( "input file set: 'blunt30deg.1', blunt30deg.1-configuration.txt, LinearElasticity_Example-variables.txt");
}



/** *****************************************************************************************

   Linear elasticity with option of restriction to region 'rock' 
   and (optionally) taking into account the effect of pore fluid pressure assuming a
   Biot coefficient of 1.
 
   use 'blunt30deg.1'  an inclined blunt-tip fracture for example calculation

   for this model you have the option to include or exclude the fracture interior from the calculation.
   If you would like to specify fluid pressure in the fracture (as computed from the
   boundary conditions) you need to include it in the computations.

   @attention 2-dimensional example but code with DIM=3 should work also for a 3-dimensional model

   Use Tcl/Tk ELASTICITY_GROUP_viewer.tcl as VTK visualisation script. Open it and turn entries
   on and off to visualise different properties.

   *****************************************************************************************

*/
void LinearElasticity_Example::Run()
{
#ifdef RUN_DIRECTLY_WITH_TRIANGLE_INPUT_FILE
    // ---------------------------------------------------------------------------------------
    // 0. Import model from Shewchuk's Triangle mesher
    // ---------------------------------------------------------------------------------------
    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;
    string              file_name;
    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
    cin >> file_name;
    const bool isoparametric{true};
    mesh_interface.ReadTriangle2DMesh( file_name.c_str(), mesh_container, isoparametric );
    VSetConverter<2U>  converter;
    converter.ConvertLinearToQuadraticTriangles( mesh_container );

    // 'LinearElasticity_Example-variables.txt' is the text file that defines the variables used in this example
    Model<2U>  model( mesh_container, "LinearElasticity_Example-variables.txt" );
    string     config_file{ file_name };
    // give this model a name as none is supplied to the constructor
    model.Name( file_name.c_str() );
    mesh_container.Erase();
    PrintModelProperties( model );
#else
    // ---------------------------------------------------------------------------------------
    // 1. Load CSMP native format model
    // ---------------------------------------------------------------------------------------
    string model_name;
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'blunt30deg.1':"<<endl;
    cin.ignore();
    getline(cin, model_name);
    if (model_name.length() == 0) model_name = "blunt30deg.1";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "LinearElasticity_Example-variables.txt";
    string config_file = model_name;
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    Model<2U>  model(model_name, variable_file);
#endif
    printModelDimensions( model, true );

    // assigns the element area to a distributed variable called 'area'
    model.AssignCellCharacteristicsTo( "area", "area" );


  // ---------------------------------------------------------------------------------------
  // 2. Making group of entire model to monitor volume change
  // ---------------------------------------------------------------------------------------
    Region<2U>& model_domain(model.Region("Model"));
    double  volume = model_domain.Volume();
    cout <<"\nmain: The model has a volume of: "<< volume <<" m^3."<< endl;

  // ---------------------------------------------------------------------------------------
  // 3. Material properties, groups etc. & initial & essential conditions
  // ---------------------------------------------------------------------------------------
    Standard_IO_Handler  stdio;
    bool  restricted_to_rock      =  stdio.RecordLogicalChoice("Would you like to restrict computation to group 'rock'");
    bool  with_plane_stress       =  stdio.RecordLogicalChoice("Would you like to model 'plane stress'");
    bool  with_boundary_stresses  =  stdio.RecordLogicalChoice("Would you like to model 'boundary stresses'");
    bool  with_body_forces        =  stdio.RecordLogicalChoice("Would you like to model 'body forces'");
    bool  with_pore_pressure(false);
    if ( !restricted_to_rock ) with_pore_pressure = stdio.RecordLogicalChoice("Would you like to model 'pore pressure'");
    bool  with_volume_strains(false);
    if ( with_pore_pressure ) with_volume_strains =  stdio.RecordLogicalChoice("Would you like to model 'volume strains'");

    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile( model, config_file.c_str() );

    // writing user-choices to file
    stdio.Out();
    VTK_Interface<2U>  vtk_output;

    vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );



  // ---------------------------------------------------------------------------------------
  // 5. Computing fluid pressure
  // ---------------------------------------------------------------------------------------
  // interrelation demontrating the use of STL unary functions
    if ( with_pore_pressure ) {
         const double fluid_viscosity(1.0e-03);
         ConstantFactor<2U,divides>  conductivity( model.Database(),
                                                  "conductivity", "permeability",
                                                   fluid_viscosity );
         model.Apply( conductivity );
         // input variable values
         printRangeOfVariable( model, "conductivity" );
         vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
         // fluid pressure computation
         SteadyStatePressure( model );
         vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
      }


  // ---------------------------------------------------------------------------------------
  // 4. Defining and computing linear elastic response
  // ---------------------------------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach(2); // this is important because it sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
                               // which is needed for deformation simulations
    SAMG_Solver  solver(&settings);
    PDE_Integrator<2U,Element>  deformation( solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Element>  deformation( solver );
#endif

    PT_op<2U>                  bforces( model.Database(), "force", "displacement" );
    NumIntegral_BT_D_B_dV<2U>  stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement");
    if ( with_plane_stress )   stiffness.PlaneStress();
    NumIntegral_PT_op_dS<2U>   bstresses( model.Database(), "Neumann stress", "displacement" );
    NumIntegral_PT_op_dV<2U>   bodyforce( model.Database(), "gravity force", "displacement");
    NumIntegral_BT_D_op_dV<2U> volstrain( model.Database(), "dilatation", "Young's modulus", "Poisson's ratio", "displacement");
    NumIntegral_BT_op_dV<2U>   porepressure( model.Database(), "fluid pressure", "displacement");

    deformation.Add( &stiffness );
    deformation.Add( &bforces );  // force vector must always be there so that Dirichlet conditions are accumulated
    if ( with_body_forces )    deformation.Add( &bodyforce );
    if ( with_volume_strains ) deformation.Add( &volstrain );
    if ( with_pore_pressure )  deformation.Add( &porepressure );
    if ( with_boundary_stresses ) deformation.AddBoundaryIntegral( &bstresses );

    const bool plane_strain(true);
    const bool Eigen_vectors(true);
    StressesAndStrains<2U>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement",
                                     plane_strain, Eigen_vectors );

    if ( with_plane_stress ) postpro.PlaneStress();
     deformation.AddPostProcess( &postpro );

    if ( !restricted_to_rock ) {
        // printing the variable values that enter the equations for the whole model domain
        printRangeOfVariable( model, "Young's modulus" );
        printRangeOfVariable( model, "Poisson's ratio" );
        if ( with_body_forces ) printRangeOfVariable( model, "gravity force" );
        if ( with_volume_strains ) printRangeOfVariable( model, "dilatation" );
        printRangeOfVariable( model, "displacement" );
        vtk_output.OutputDataToVTK( model, "displacement", "displacement", 0, true );
        if ( with_pore_pressure ) printRangeOfVariable( model, "fluid pressure" );
        // solving the mechanics problem on whole model
        deformation.IntegrateOver( model_domain );
      }
    // note! - the model must be supplied here so that the algorithm can search for boundaries that touch the computational domain
    else {
        // printing the variable values that enter the equations for the whole model domain
        printRangeOfVariable( model, "rock", "Young's modulus" );
        printRangeOfVariable( model, "rock", "Poisson's ratio" );
        if ( with_body_forces ) printRangeOfVariable( model, "rock", "gravity force" );
        if ( with_volume_strains ) printRangeOfVariable( model, "rock", "dilatation" );
        printRangeOfVariable( model, "rock", "displacement" );
        vtk_output.OutputDataToVTK( model, "rock", "displacement", "displacement", 0, true );
        if ( with_pore_pressure ) printRangeOfVariable( model, "rock", "fluid pressure" );
        // solving the mechanics problem only for the region 'rock'
        deformation.IntegrateOver( model, model.Region("rock") );
      }


  // ---------------------------------------------------------------------------------------
  // 6. Extracting the vertical stress component
  // ---------------------------------------------------------------------------------------
    ExtractTensorVariableComponent<2U>  ystress(   model.Database(), "stress", "stress-y",  1,1 );
    ExtractTensorVariableComponent<2U>  stress_xy( model.Database(), "stress", "stress-xy", 0,1 );
    if ( restricted_to_rock ) {
         model.Region("rock").Apply( ystress );
         model.Region("rock").Apply( stress_xy );
      }
    else {
         model.Apply( ystress );
         model.Apply( stress_xy );
      }


  // ---------------------------------------------------------------------------------------
  // 7. applying DISPLACEMENT TO MESH
  // ---------------------------------------------------------------------------------------
    model.MoveNodeCoordinatesBy("displacement");

  // ---------------------------------------------------------------------------------------
  // 8. output of deformed mesh etc.
  // ---------------------------------------------------------------------------------------
    // output of the region 'rock' only
    if ( restricted_to_rock ) {
      vtk_output.OutputDataToVTK( model, "rock", "displacement", "displacement", 1, true );
      vtk_output.OutputDataToVTK( model, "rock", "strain",       "strain",       1, true );
      vtk_output.OutputDataToVTK( model, "rock", "stress",       "stress",       1, true );
      vtk_output.OutputDataToVTK( model, "rock", "mean-stress",  "mean stress",  1, true );
      vtk_output.OutputDataToVTK( model, "rock", "dilatation",   "dilatation",   1, true );
      vtk_output.OutputDataToVTK( model, "rock", "sigma1_",      "sigma1",       1, true );
      vtk_output.OutputDataToVTK( model, "rock", "sigma2_",      "sigma2",       1, true );
      vtk_output.OutputDataToVTK( model, "rock", "strain1_",     "strain1",      1, true );
      vtk_output.OutputDataToVTK( model, "rock", "strain2_",     "strain2",      1, true );
      vtk_output.OutputDataToVTK( model, "rock", "stress-y",     "stress-y",     1, true );
      vtk_output.OutputDataToVTK( model, "rock", "stress-xy",    "stress-xy",    1, true );
      }
    // output of the entire model
    else
      {
      vtk_output.OutputDataToVTK( model, "displacement", "displacement", 1, true );
      vtk_output.OutputDataToVTK( model, "strain",       "strain",       1, true );
      vtk_output.OutputDataToVTK( model, "stress",       "stress",       1, true );
      vtk_output.OutputDataToVTK( model, "mean-stress",  "mean stress",  1, true );
      vtk_output.OutputDataToVTK( model, "dilatation",   "dilatation",   1, true );
      vtk_output.OutputDataToVTK( model, "sigma1_",      "sigma1",       1, true );
      vtk_output.OutputDataToVTK( model, "sigma2_",      "sigma2",       1, true );
      vtk_output.OutputDataToVTK( model, "strain1_",     "strain1",      1, true );
      vtk_output.OutputDataToVTK( model, "strain2_",     "strain2",      1, true );
      vtk_output.OutputDataToVTK( model, "stress-y",     "stress-y",     1, true );
      vtk_output.OutputDataToVTK( model, "stress-xy",    "stress-xy",    1, true );

      // fluid pressure and Darcy velocity from post-processing
        if ( with_pore_pressure ) {
           vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
           vtk_output.OutputDataToVTK( model, "velocity",       "velocity", 1 );
        }
      }

    cout <<"\nmain: That's it..."<< endl;

#ifndef RUN_DIRECTLY_WITH_TRIANGLE_INPUT_FILE
    filesystem::current_path("../../example_inputs/");
#endif
} // end Run





/**
   For the pressure applied at the model TOP and BOTTOM boundary, this
   method computes the pore pressure in the region of interest.
   Subsequently, the Darcy velocity is computed in a post-processing step.
*/
void LinearElasticity_Example::SteadyStatePressure( Model<2U>& model )
 {
   // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
   //ostream &cout = *GetStream();
   
   SteadyStateDiffusor<2U,Element>  pressure( model, "conductivity", "fluid pressure", "fluid volume source" );

   VelocityAndVolumeFlux<2U>  postpro( model, "conductivity", "porosity", "fluid pressure" );

   pressure.AddPostProcess( &postpro );
   
#if defined CSMP_WITH_SAMG_SOLVER && defined SAMG_MULTIPLE_INSTANCES
   pressure.GetSolverSettings().SetSolverInstance(2);
#endif
   pressure.ComputeSteadyState( model.Region("Model") );

   double  fmin, fmax;
   model.MinMaxOf( "velocity", fmin, fmax );
   cout <<"\nsteadyStatePressure: computed 'Darcy velocity' range (min/max, m s-1): "<< fmin <<", "<< fmax << endl;
   model.MinMaxOf( "fluid pressure", fmin, fmax );
   cout <<"\nsteadyStatePressure: computed 'fluid pressure' range (min/max, Pa): "<< fmin <<", "<< fmax << endl;

} // end SteadyStatePressure



void LinearElasticity_Example::PrintModelProperties( const Model<2U>& model )
 {
    cout << "\n==================================================================================================";
    cout << "\nModel '"<< model.Name() <<"' has been established successfully ";
    if ( model.Mesh().Elements() > 0 ) {
         size_t volume_elmts{0U}, surface_elmts{0U}, line_elmts{0U};
         cout <<"(total cells "<< currentCellTypes( model.Mesh(), ELEMENT, volume_elmts, surface_elmts, line_elmts );
         cout <<", nodes "<< model.Mesh().Nodes() <<")";
         cout <<"\n\t\t\t("<< model.Mesh().Elements() <<" elements: volumes "<< volume_elmts <<", surfaces "<< surface_elmts <<", lines "<< line_elmts <<")";
      }
    if ( model.Mesh().Faces() > 0 ) {
         size_t volume_faces{0U}, surface_faces{0U}, line_faces{0U};
         currentCellTypes( model.Mesh(), FACE, volume_faces, surface_faces, line_faces );
         assert( volume_faces == 0U );
         cout <<"\n\t\t\t("<< model.Mesh().Faces() <<" faces: surfaces "<< surface_faces <<", lines "<< line_faces <<")";
      }
    if ( model.Mesh().Interfaces() > 0 ) {
         size_t volume_ifaces{0U}, surface_ifaces{0U}, line_ifaces{0U};
         currentCellTypes( model.Mesh(), INTER_FACE, volume_ifaces, surface_ifaces, line_ifaces );
         assert( volume_ifaces == 0U );
         cout <<"\n\t\t\t("<< model.Mesh().Interfaces() <<" interfaces: surfaces "<< surface_ifaces <<", lines "<< line_ifaces <<")";
      }
    cout << "\n==================================================================================================";
    cout << endl;
    
} // end PrintModelProperties

} // csmp
