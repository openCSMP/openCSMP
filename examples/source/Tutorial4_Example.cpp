// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Tutorial4_Example.h"

// CSMP model
#include "ANSYS_Model2D.h" // TODO: replace and use CSMP native model instead
#include "Boundary.h"

// FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_dNi_dV_sc.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dS.h"

#include "LinearSolver.h"

// interfaces
#include "VTU_Interface.h"
#include "VTK_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"

// visitors
#include "StabilizationParameterVisitor.h"

using namespace std;

namespace csmp {

void Tutorial4_Example::Specifications()
{
  SetTitle( "Tutorial 4: Stokes lubrication equation" );
  SetDifficulty( 3 );
  SetCategory( "Tutorials (composite functionality)" );
  AddAuthor( "Sebastian Geiger. Edited by Luat and Hani" );
  AddDescription( "A CSMP main file that uses an solves the coupled Stokes equation in a 2D geometry representing" );
  AddDescription( "pores and grains in a carbonate rock. The solution of the Stokes equation provides the pressure" );
  AddDescription( "and velocity fields inside the connected pore scale. Note that special SAMG settings are needed" );
  AddDescription( "for solving the Stokes equation using a standard FE method. Furthermore, a stabilisation parameter" );
  AddDescription( "is introduced to allow the discretisation of velocity and pressure in the Stokes equation using" );
  AddDescription( "the same FE basis functions." );
  AddRequirement( "source code in: 'Tutorial4_Example.cpp'" );
  AddRequirement( "model 'pores' binary files, stokes_variables.txt, no need of configuration" );
}


/** 

Stokes (lubrication) equation (no inertia, incompressible, viscous flow)

A CSMP main file that uses an solves the coupled Stokes equation in a 2D geometry representing
pores and grains in a carbonate rock. The solution of the Stokes equation provides the pressure
and velocity fields inside the connected pore scale. Note that special SAMG settings are needed
for solving the Stokes equation using a standard FE method. Furthermore, a stabilisation parameter
is introduced to allow the discretisation of velocity and pressure in the Stokes equation using
the same FE basis functions.

For details on the numerical scheme see Zaretskiy, Geiger, Sorbie and Foerster
Advances in Water Resources (2010), doi:10.1016/j.advwatres.2010.08.008

*/
void Tutorial4_Example::Run()
{
    // -----------------------------------------
    // 1.0 Build the CSMP Model from ANSYS files
    // -----------------------------------------
    string  input_file;
#ifdef BUILD_MODEL_FROM_ANSYS_FILE_SET
    cout << "\nTutorial4_Example: Please enter the name of input mesh ( default: pores ): ";
    cout.flush();
    cin >> input_file;
    ANSYS_Model2D model( input_file.c_str(), "stokes_variables.txt" );
#else
    // ------------------------------------------------------------
    // 1.1 Load CSMP native format model
    // ------------------------------------------------------------
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'pores':"<<endl;
    cin.ignore();
    getline(cin, input_file);
    if (input_file.length() == 0) input_file = "pores";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "stokes_variables.txt";
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, input_file, variable_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    Model<2U>  model(input_file, variable_file);
#endif
    // scaling model to sub-mm scale
    scaleRegion( model, 20000.0 ); // scale the model size by 1/value
    printModelDimensions( model, true );


  // ------------------------------------
  // 2.0 Boundary and initial conditions
  // ------------------------------------
  // input material and fluid properties
  const double viscosity{ 0.001 }; // in Pa sec
  model.InputPropertyValue( "viscosity", makeScalar( PLAIN, viscosity ) );
  model.InputPropertyValue( "porosity", makeScalar( PLAIN, 1.0 ) );
  model.InputPropertyValue( "zero", makeScalar( PLAIN, 0.0 ) ); // dummy variable needed to close coupled FE algorithm

  // assigning initial conditions
  model.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( "nodal velocity x", makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( "nodal velocity y", makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( "Neumann traction x", makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( "Neumann traction y", makeScalar( PLAIN, 0.0 ) );

  // Dirichlet (no-slip) conditions for velocity at horizontal model boundaries
  model.InputBoundaryValue( TOP, "nodal velocity x", makeScalar( DIRICH, 0.0 ) );
  model.InputBoundaryValue( TOP, "nodal velocity y", makeScalar( DIRICH, 0.0 ) );
  model.InputBoundaryValue( BOTTOM, "nodal velocity x", makeScalar( DIRICH, 0.0 ) );
  model.InputBoundaryValue( BOTTOM, "nodal velocity y", makeScalar( DIRICH, 0.0 ) );

  // Dirichlet (no-slip) conditions for velocity at boundaries of mineral grains
  model.Region( "GRAINS" ).ChangePropertyStatus( "nodal velocity x", DIRICH );
  model.Region( "GRAINS" ).ChangePropertyStatus( "nodal velocity y", DIRICH );

  // Neumann condition for traction term
  model.InputBoundaryValue( LEFT, "Neumann traction x", makeScalar( NEUMANN, 200.0 ) );
  model.InputBoundaryValue( LEFT, "Neumann traction y", makeScalar( NEUMANN, 0.0 ) );
//  model.InputBoundaryValue( RIGHT, "Neumann traction x", makeScalar( NEUMANN, -200.0 ) ); FAIL
//  model.InputBoundaryValue( RIGHT, "Neumann traction y", makeScalar( NEUMANN, 0.0 ) ); FAIL
//  model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar( DIRICH, 1.0325e+5 ) );
//  model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar( DIRICH, 1.0325e+5 ) );

  // visitor for calculation of the stabilization parameter, which is needed to compute
  // Stokes flow using FEs which have the same basis functions for velocity and pressure
  StabilizationParameterVisitor<2U> stabpar_visitor( model, viscosity, 2.0, "stabilization parameter" );
  model.Accept( stabpar_visitor );
  printRangeOfVariable( model, "stabilization parameter" );

  // ------------------------------------------------------------
  // 3.0 Coupled FE algorithm for solving Stokes equation
  // ------------------------------------------------------------

#ifdef CSMP_WITH_SAMG_SOLVER
  // custom SAMG settings for multi-variable solution as identified by Malte Foerster
  SAMG_Settings settings;

  // nsolve
  settings.Set_napproach( 2 ); // coupled system approach, set nsys>1: interpolation separate for each unknown
  settings.Set_nxtyp( 9 );     // Uzawa smoothing is best option for saddle-point problems like this Stokes p,u one
  settings.Set_igam(1);     // 1 = standard V-cycle approach
  // ncycle
  settings.Set_ncgrad( 3 ); // for GMRES
  settings.Set_nkdim( 9 );
  settings.Set_ncycle( 50 );
 
  // FE algorithm with specialised SAMG settings
  SAMG_Solver  solver( &settings );
  PDE_Integrator<2U,Element>  stokes_flow( solver );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  solver;
  PDE_Integrator<2U,Element>  stokes_flow( solver );
#endif

  // Stokes lubrication equation
  const PropertyDatabase<2U>& p_ref(model.Database());

  // LHS operators                                        operand       basic function      test function
  NumIntegral_dNT_lhsop_dN_dV<2U>  viscosity_matr_x( p_ref, "viscosity", "nodal velocity x", "nodal velocity x" );
  NumIntegral_NT_dNi_dV_sc<2U>  gradient_x( p_ref, "zero", "fluid pressure", "nodal velocity x" );
  gradient_x.SpatialDerivative( X_DIRECTION );
  gradient_x.MultiplyBy( -1.0 );

  NumIntegral_dNT_lhsop_dN_dV<2U>  viscosity_matr_y( p_ref, "viscosity", "nodal velocity y", "nodal velocity y" );
  NumIntegral_NT_dNi_dV_sc<2U>  gradient_y( p_ref, "zero", "fluid pressure", "nodal velocity y" );
  gradient_y.SpatialDerivative( Y_DIRECTION );
  gradient_y.MultiplyBy( -1.0 );

  // LHS operators                                    operand  basic function      test function
  NumIntegral_NT_dNi_dV_sc<2U>  divergence_x( p_ref, "zero", "nodal velocity x", "fluid pressure" );
  divergence_x.SpatialDerivative( X_DIRECTION );
  divergence_x.Transposed();

  NumIntegral_NT_dNi_dV_sc<2U>  divergence_y( p_ref, "zero", "nodal velocity y", "fluid pressure" );
  divergence_y.SpatialDerivative( Y_DIRECTION );
  divergence_y.Transposed();

  NumIntegral_dNT_lhsop_dN_dV<2U>  stab_matrix( p_ref, "stabilization parameter", "fluid pressure", "fluid pressure" );

  // RHS operators                                  operand               test function
  NumIntegral_NT_op_N_dS<2U> extpressure_x( p_ref, "Neumann traction x", "nodal velocity x" );
  NumIntegral_NT_op_N_dS<2U> extpressure_y( p_ref, "Neumann traction y", "nodal velocity y" );

  // RHS operators                            operand   test function
  NumIntegral_NT_rhsop_N_dV<2U>  dummy( p_ref, "zero", "fluid pressure" );

  // add each PDE Operator to the FE algorithm

  // Stokes lubrication equation
  stokes_flow.Add( &viscosity_matr_x );
  stokes_flow.Add( &gradient_x );
  stokes_flow.Add( &viscosity_matr_y );
  stokes_flow.Add( &gradient_y );
  stokes_flow.AddBoundaryIntegral( &extpressure_x );
  stokes_flow.AddBoundaryIntegral( &extpressure_y );

  // Continuity equation
  stokes_flow.Add( &divergence_x );
  stokes_flow.Add( &divergence_y );
  stokes_flow.Add( &stab_matrix );
  stokes_flow.Add( &dummy );

  // ----------------------------------------------------------------------------------
  // 3.1 Record CPU time required for simulation
  // ----------------------------------------------------------------------------------
  clock_t start( clock() );

  // solve the Stokes equation and release memory (model is needed tgo get boundary integral right)
  stokes_flow.IntegrateOver( model, model.Region("PORES") );
  stokes_flow.Reset();

  // record timing of the simulation
  clock_t end( clock() );
  cerr << "\n\nmain: CPU time was " << static_cast<double>((end - start) / CLOCKS_PER_SEC) << " seconds " << endl << endl;

  // output resulting variable ranges
  printRangeOfVariable( model, "nodal velocity x" );
  printRangeOfVariable( model, "nodal velocity y" );
  printRangeOfVariable( model, "fluid pressure" );

  // construct nodal velocity vector
  constructVelocityVector( model );
  printRangeOfVariable( model, "nodal velocity" );

  // output results to VTU
  VTU_Interface<2U>   vtu( model );
  std::list<std::string> outputProps;
  outputProps.push_back( "fluid pressure" );
  outputProps.push_back( "nodal velocity" );

  // output pressure and (nodal) velocity to the same VTU file
  vtu.OutputDataToVTU( "PoreScaleVariables", outputProps, "PORES", static_cast<int>(0) );
  
  // terminate
  cerr << "\n\nmain: That's it, run completed successfully..." << endl;

} // Run()




  /**
  Quick hack that loops over a 2D boundary, reads in a nodal source term and scales it by 0.5 of the length
  of the edges of the connected FEs wich lie at the model boundary (here it is the LEFT boundary)
  */
void Tutorial4_Example::assignFluxToPointSource( Model<2U>& mdl, const char* flux )
{
  // backup the fluxes that are orginally assigned in a temporary variable
  const string temp_flux( "flux" );

  if ( !mdl.Database().IsDefined( temp_flux.c_str() ) )
    mdl.CreateProperty( temp_flux.c_str(), "Q", "m3/s", SCALAR, NODE );
  mdl.CopyReplace( flux, temp_flux.c_str() );

  // define variables needed to compute the length of the FE edges that lie at the boundary
  // and read the temporary flux and store the final flux
  const csmp::Index tf_key( mdl.Database().StorageKey( temp_flux.c_str() ) ), // nodal source
                    f_key( mdl.Database().StorageKey( flux ) );

  cout << "\nassignFluxToPointSource: Translating '" << flux << "' into a nodal point source" << endl;

  // loop over all finite elements and identify elements that lie at the model boundary of interest (here LEFT)
  Boundary<2U>&   left = mdl.Boundary( "LEFT" );
  for ( auto eit = left.CellsBegin(); eit != left.CellsEnd(); eit++ )
    {
       // getting the area of the face
        double area = (*eit)->Area();

        // second loop to calculate and scale nodal flux
        for ( uint32_t i{0u}; i<(*eit)->Nodes(); i++ ) {
            // read existing flux at node i
            double tf = (*eit)->N(i)->Read( tf_key );
            // read existing source at node i
            double f = (*eit)->N(i)->Read( f_key );
            // scale the influx with area / nodes
            f += tf * (area / static_cast<double>((*eit)->Nodes()));
            // store it back to node
// SKM FIX            (*eit)->N(i)->Store( f_key, makeScalar(NEUMANN,f) );
            (*eit)->N(i)->Store( f_key, makeScalar(DIRICH,f) );
          }
    }

  printRangeOfVariable( mdl, flux );
  cout << "\nassignFluxToPointSource: Successfully assigned '" << flux << "' to the nodes... " << endl;
}



/// simple combine x and y component of computed velocities into a vector and saves it to the grod
void Tutorial4_Example::constructVelocityVector( Model<2U>& mdl )
{
  const csmp::Index vx_key( mdl.Database().StorageKey( "nodal velocity x" ) ),
                    vy_key( mdl.Database().StorageKey( "nodal velocity y" ) ),
                     v_key( mdl.Database().StorageKey( "nodal velocity" ) );

  Region<2U>& mref = mdl.Region( "Model" );
  for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
    (*nit)->Store( v_key, makeVector(PLAIN, PLAIN, (*nit)->Read( vx_key ), (*nit)->Read( vy_key ) ) );
}



/// scale the size of the CSMP model (NB divides by the provided factor!)
void Tutorial4_Example::scaleRegion( Model<2U>& mdl, double scale_factor )
{
  Region<2U>& mref = mdl.Region( "Model" );
  for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
    {
      double x = (*nit)->x();
      (*nit)->x( x / scale_factor );
      double y = (*nit)->y();
      (*nit)->y( y / scale_factor );
    }
}

} // csmp
