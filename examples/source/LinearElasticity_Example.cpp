#include "LinearElasticity_Example.h"

#include "Region.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
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
    /*
    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;
    char                file_name[200];
    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
    cin >> file_name;
    mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );
    VSetConverter<2U>  converter;
    converter.ConvertLinearToQuadraticTriangles( mesh_container );

    // 'example12.txt' is the text file that defines the variables used in this example
    Model<2U>  model( mesh_container, "LinearElasticity_Example-variables.txt" );
    mesh_container.Erase();
    */

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
  // NB: interrelations are deprecated now; nonetheless this demontrates use of STL unary functions
    if ( with_pore_pressure ) {
         const double fluid_viscosity(1.0e-03);
         ConstantFactor<2U,divides>  conductivity( model.Database(),
                                                  "conductivity", "permeability",
                                                   fluid_viscosity );
         model.Apply( conductivity );
         SteadyStatePressure( model );
         vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
      }


  // ---------------------------------------------------------------------------------------
  // 4. Defining and computing linear elastic response
  // ---------------------------------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach(2); // this is important because it sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
                               // which is needed for deformation simulations
    SAMG_Solver  solver(&settings);
    PDE_Integrator<2U,Region>  deformation( solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Region>  deformation( solver );
#endif

    PT_op<2U,Element<2U> >     bforces( model.Database(), "force", "displacement" );
    NumIntegral_BT_D_B_dV<2U>  stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement");
    if ( with_plane_stress ) stiffness.PlaneStress();
    NumIntegral_PT_op_dS<2U>   bstresses( model.Database(), "Neumann stress", "displacement" );
    NumIntegral_PT_op_dV<2U>   bodyforce( model.Database(), "gravity force", "displacement");
    NumIntegral_BT_D_op_dV<2U> volstrain( model.Database(), "dilatation", "Young's modulus", "Poisson's ratio", "displacement");
    NumIntegral_BT_op_dV<2U>   porepressure( model.Database(), "fluid pressure", "displacement");

    deformation.Add( &stiffness );
    deformation.Add( &bforces );  // force vector must always be there so that Dirichlet conditions are accumulated
    if ( with_body_forces )       deformation.Add( &bodyforce );
    if ( with_volume_strains )    deformation.Add( &volstrain );
    if ( with_pore_pressure )     deformation.Add( &porepressure );

#ifdef CSMP_WITH_SAMG_SOLVER
    if ( with_boundary_stresses ) deformation.AddBoundaryIntegral( &bstresses );
#endif

    const bool plane_strain(true);
    const bool principal_vectors(true);
    StressesAndStrains<2U>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", plane_strain, principal_vectors );

    if ( with_plane_stress ) postpro.PlaneStress();
     deformation.AddPostProcess( &postpro );

    if ( !restricted_to_rock ) deformation.IntegrateOver( model_domain );
    // note! - the model must be supplied here so that the algorithm can search for boundaries that touch the computational domain
    else deformation.IntegrateOver( model, model.Region("rock") );


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

    fs::current_path("../../example_inputs/");
  
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
   
   SteadyStateDiffusor<2U,Region>  pressure( model, "conductivity", "fluid pressure", "fluid volume source" );

   VelocityAndVolumeFlux<2U,Element<2U> >  postpro( model, "conductivity", "porosity", "fluid pressure" );

   pressure.AddPostProcess( &postpro );
   
#ifdef SAMG_MULTIPLE_INSTANCES
   pressure.GetSolverSettings().SetSolverInstance(2);
#endif
   pressure.ComputeSteadyState( model.Region("Model") );

   double  fmin, fmax;
   model.MinMaxOf( "velocity", fmin, fmax );
   cout <<"\nsteadyStatePressure: computed 'Darcy velocity' range (min/max, m s-1): "<< fmin <<", "<< fmax << endl;
   model.MinMaxOf( "fluid pressure", fmin, fmax );
   cout <<"\nsteadyStatePressure: computed 'fluid pressure' range (min/max, Pa): "<< fmin <<", "<< fmax << endl;

} // end SteadyStatePressure



} // csmp
