// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PETSc_Example.h"

#include "Region.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "meshManagementUtilities.h"
#include "CSMP_definitions.h"

// File I/O and Initialization
#include "TRIANGLE_Interface.h"
#include "VTU_Interface.h"
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

// uncomment to compare with EigenSolver
/*
#if defined(CSMP_WITH_PETSC_SOLVER)
#undef CSMP_WITH_PETSC_SOLVER
#endif
*/

#if defined(CSMP_WITH_PETSC_SOLVER)
#include "PETSc_Solver.h"
#include "PETSc_Settings.h"
#else
#include "LinearSolver.h"
#endif

// PDE Operators
// mechanics
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "PT_op.h"
//#include "StressesAndStrains2.h"
#include "StressesAndStrains.h"
// fluid flow
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "VelocityAndVolumeFlux.h"

// comment-in this preprocessor directive if you want to run directly with the Triangle mesher input file deck
// #define RUN_DIRECTLY_WITH_TRIANGLE_INPUT_FILE

using namespace std;

namespace csmp {

void PETSc_Example::Specifications()
{
  SetTitle( "PETSc_Example" );
  SetDifficulty( 3 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "SKM" );
  AddDescription( "Linear elasticity in model with fracture void; demonstration of PETSc_Solver application." );
  AddDescription( "Example shows both, the solution of a vector (displacement) and a scalar (fluid pressure) variable problem." );
  AddDescription( "source in: PETSc_Example.cpp" );
  AddRequirement( "default input file set (Triangle mesher): 'blunt30deg.1', blunt30deg.1-configuration.txt, PETSc_Example-variables.txt");
}



/** *****************************************************************************************

   Solves 'fluid pressure' everywhere and linear elasticity problem within the model subregion region 'rock'
   Pore pressure is taken into account assuming a  Biot coupling coefficient of 1.
   
   The PETSc_Solver  is demonstrated here for separate scalar and vector primary (solution) variables and combined boundary conditions
 
   Example uses model 'blunt30deg.1'  an inclined blunt-tip fracture in the calculation

   *****************************************************************************************

 @author S.K. Matthaei
 
*/
void PETSc_Example::Run()
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

    // 'PETSc_Example-variables.txt' is the text file that defines the variables used in this example
    Model<2U>  model( mesh_container, "PETSc_Example-variables.txt" );
    string     config_file{ file_name };
    // give this model a name as none is supplied to the constructor
    model.Name( file_name.c_str() );
    mesh_container.Erase();
    PrintModelProperties( model );
#else
    // ---------------------------------------------------------------------------------------
    // 1. Load CSMP native format model
    // ---------------------------------------------------------------------------------------
    const string model_name("blunt30deg.1");
    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "PETSc_Example-variables.txt";
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
    double  initial_area = model_domain.Volume();
    cout <<"\nmain: The model has an initial area of: "<< initial_area<<" m^2."<< endl;

  // ---------------------------------------------------------------------------------------
  // 3. Material properties, groups etc. & initial & essential conditions
  // ---------------------------------------------------------------------------------------
    Standard_IO_Handler  stdio;
    bool  restricted_to_rock      =  false; // stdio.RecordLogicalChoice("Would you like to restrict computation to group 'rock'");
    bool  with_plane_stress       =  true; // stdio.RecordLogicalChoice("Would you like to model 'plane stress'");
    bool  with_boundary_stresses  =  true; // (load on top boundary) stdio.RecordLogicalChoice("Would you like to model 'boundary stresses'");
    bool  with_body_forces        =  true; // (action of gravity on rock mass) stdio.RecordLogicalChoice("Would you like to model 'body forces'");
    bool  with_pore_pressure(true);
    bool  with_volume_strains(true);

    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile( model, config_file.c_str() );

    // writing user-choices to file
    stdio.Out();


  // ---------------------------------------------------------------------------------------
  // 5. Computing fluid pressure
  // ---------------------------------------------------------------------------------------
  // interrelation demontrating the use of STL unary functions
    if ( with_pore_pressure ) {
         constexpr double fluid_viscosity(1.0e-03);
         ConstantFactor<2U,divides>  conductivity( model.Database(),
                                                  "conductivity", "permeability",
                                                   fluid_viscosity );
         model.Apply( conductivity );
         // input variable values
         printRangeOfVariable( model, "conductivity" );
         // fluid pressure computation
         SteadyStatePressure( model );
      }


  // ---------------------------------------------------------------------------------------
  // 4. Defining and computing linear elastic response
  //    (result is the vector variable "displacement" which requires specific solution approach)
  // ---------------------------------------------------------------------------------------
#if defined(CSMP_WITH_PETSC_SOLVER)
    PETSc_Settings settings;
    
    // linear elasticity solver profile
    settings.SetKSPType("cg");
    settings.SetPCType("gamg");
    settings.SetRelativeTolerance(1.0e-8);
    settings.SetMaximumIterations(10000);
    
    PETSc_Solver   solver(&settings);
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
        if ( with_pore_pressure ) printRangeOfVariable( model, "fluid pressure" );
        
        // solving the mechanics problem on whole model
        // ============================================
        deformation.IntegrateOver( model_domain );
      }
    else {
        // printing the variable values that enter the equations for the whole model domain
        printRangeOfVariable( model, "rock", "Young's modulus" );
        printRangeOfVariable( model, "rock", "Poisson's ratio" );
        if ( with_body_forces ) printRangeOfVariable( model, "rock", "gravity force" );
        if ( with_volume_strains ) printRangeOfVariable( model, "rock", "dilatation" );
        printRangeOfVariable( model, "rock", "displacement" );
        if ( with_pore_pressure ) printRangeOfVariable( model, "rock", "fluid pressure" );
        
        // solving the mechanics problem only for the region 'rock'
        // =======================================================
        // note: Model must be supplied to the PDE_Integrator here so that
        // it can search for potential Neumann boundaries that touch the computational domain
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
    double  final_area = model_domain.Volume();
    cout <<"\nmain: The model now has an area of: "<< final_area <<" m^2."<< endl;
    cout <<" Its area changed by "<< final_area - initial_area <<" m2"<< endl;


  // ---------------------------------------------------------------------------------------
  // 8. output of deformed mesh etc.
  // ---------------------------------------------------------------------------------------
  list<string> output_props{ "displacement", "strain", "stress", "mean stress", "stress-xy",
                             "dilatation", "fluid pressure", "area", "velocity" };
                             
  VTU_Interface<2> vtu_quadratic( model );
  vtu_quadratic.OutputDataToVTU( model_name, output_props, string("rock"), static_cast<long>(1) );
  //vtu_quadratic.OutputDataToVTU( model_name, output_props, string("Model"), static_cast<long>(1) );

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
void PETSc_Example::SteadyStatePressure( Model<2U>& model )
 {
   // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
   //ostream &cout = *GetStream();

#if defined(CSMP_WITH_PETSC_SOLVER)
    PETSc_Settings settings;
    
    // steady state pressure solver profile
    settings.SetKSPType("cg");
    settings.SetPCType("gamg");
    settings.SetRelativeTolerance(1.0e-8);
    settings.SetMaximumIterations(10000);
    
    PETSc_Solver   solver(&settings);
    
    // using PETSc's defaults to solve this scalar equation
    PDE_Integrator<2U,Element>  pressure( solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Element>  pressure( solver );
#endif
    NumIntegral_dNT_lhsop_dN_dV<2U> conductance( model.Database(), "conductivity", "fluid pressure",  "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<2U>   source( model.Database(),  "fluid volume source", "fluid pressure" );
    VelocityAndVolumeFlux<2U>       postpro( model, "conductivity", "porosity", "fluid pressure" );

    pressure.Add( &conductance );        // conductance matrix
    pressure.Add( &source );             // mass matrix
    pressure.AddPostProcess( &postpro ); // velocity and volume flux post-processing

    model.Apply( pressure );

    double  fmin, fmax;
    model.MinMaxOf( "fluid pressure", fmin, fmax );
    cout <<"\nsteadyStatePressure: computed 'fluid pressure' range (min/max, Pa): "<< fmin <<", "<< fmax << endl;
    model.MinMaxOf( "velocity", fmin, fmax );
    cout <<"\nsteadyStatePressure: computed 'Darcy velocity' range (min/max, m s-1): "<< fmin <<", "<< fmax << endl;

} // end SteadyStatePressure




void PETSc_Example::PrintModelProperties( const Model<2U>& model )
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
