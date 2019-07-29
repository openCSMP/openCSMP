#include "DESAdvectionDiffusion2D_Example.h"

// the CSMP model
#include "Model.h"

// the FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

// DES algorithms
#include "DESAdvectionDiffusion.h"

// output interfaces
#include "VTK_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "ConstantFactor.h"
#include "Standard_IO_Handler.h"

// FE grid generation
#include "Quadrilaterator.h"

using namespace std;

namespace csmp {

void DESAdvectionDiffusion2D_Example::Specifications()
{
  SetTitle( "2D Passive advection of concentration using DES (discrete event simulation)" );
  SetDifficulty( 3 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "Qi Shao" );
  AddDescription( "A steady state pressure diffusion equation solved using a fully" );
  AddDescription( "implicit FE discretisation and computes the velocity field averwards. This velocity field is" );
  AddDescription( "then used to compute the advection of a non-reacting chemical species with both the " );
  AddDescription( "discrete event simulation (DES) and the time-driven simulation (TDS)." );
  AddDescription( "Their outputs are written to VTK files and and their efficiency are compared." );
  AddDescription( "A simple quadrilateral FE mesh is generated automatically in CSMP from which the Model is built." );
  AddRequirement( " tutorial2_input_20x20 (20 x 20m), DES_variables.txt");
} // Initialize()

// **********************************************************************************************
//
// A CSMP main file that first solves the steady state pressure diffusion equation using a fully
// implicit FE discretisation and computes the velocity field averwards. This velocity field is
// then used to compute the advection of a non-reacting chemical species with both the
// discrete event simulation (DES) and the time-driven simulation (TDS).
// Their outputs are written to VTK files and and their efficiency are compared.
// A simple quadrilateral FE mesh is generated automatically in CSMP from
// which the Model is built.
//
// **********************************************************************************************

void DESAdvectionDiffusion2D_Example::Run()
{
    // -----------------------------------------------------------------------
    // 1.0 Generate a simple quadrilateral FE mesh from color-coded input file
    // -----------------------------------------------------------------------
    Quadrilaterator    quadrilaterator; // simple FE mesher
    VSet<2U>           mesh_container;  // container to store the input mesh
    string             file_name;
    double64           x, y;
    cerr << "\nmain: Enter the pixel-based input geometry for the quadrilaterator: " << endl;
    cin >> file_name;
    cerr << "\nmain: The x- and y-dimensions of your model (in m): " << endl;
    cin >> x;
    cin >> y;

    // check that model dimensions are ok
    if ( x <= 0.0 or y <= 0.0 ) {
        cerr << "\nmain: Non-physical model dimensions in x- and/or y-direction, terminating... " << endl;
        return;
      }

    // read in file and generate mesh
    quadrilaterator.QuadrilateralsFromRegularGrid( mesh_container, file_name.c_str(), x, y );

    // --------------------------------------------
    // 2.0 Create Model with isoparametric FEs
    // --------------------------------------------
    Model<2U>                  model( mesh_container, "DES_variables.txt", true ); // true = isoparametric FEs
    const PropertyDatabase<2>& p_ref = model.Database();

    // give the model dimensions
    printModelDimensions( model, true );

    // --------------------------------------------
    // 3.0 Applying boundary and initial conditions
    // --------------------------------------------
    // assigning material properties
    model.InputPropertyValue( "porosity",         makeScalar(PLAIN,0.1) );     // always as a fraction
    model.InputPropertyValue( "diffusivity",      makeScalar(PLAIN,1.0e-06) ); // solute diffusivity (units m2 s-1)


    // assigning initial conditions
    model.InputPropertyValue( "fluid pressure",        makeScalar(PLAIN,1.0e+07) );  // always in Pascal
    model.InputPropertyValue( "fluid volume source",   makeScalar(PLAIN,0.0) );      // no sources/sinks (units m3 m-2 s-1)
    model.InputPropertyValue( "concentration",         makeScalar(PLAIN,1.0) );      // initially one (units kg m-3)
    model.InputPropertyValue( "new concentration",     makeScalar(PLAIN,0.0) );
    model.InputPropertyValue( "nodal concentration source",  makeScalar(PLAIN,0.0) );      // no sources/sinks for solute (units kg m-3 s-1)


    // assigning boundary conditions such that flow is from LEFT to RIGHT and solute
    // transport occurs in the same direction
    model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,3.0e+07) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.0e+07) );
    model.InputBoundaryValue( LEFT,  "concentration",  makeScalar(DIRICH,10.0) );
    model.InputBoundaryValue( RIGHT, "concentration",  makeScalar(DIRICH,1.0) );

    // ------------------------------------------------------------------------------------
    // 4.0 Form a new model sub-region based on the permeability of the input mesh, change
    //     some values within this region
    // ------------------------------------------------------------------------------------
    string region_name("CENTRAL_REGION");
    cout << "\nmain: Forming the Region '" << region_name << "': " << endl;

    // all elements in permeability range 1.0e-21 to 1.0e-19 (Color value 0 and 236 to 255) are collected in the Region CENTRAL REGION
    // NB: Other color values of 128 corresponds to permeability of 1.0e-15
    model.FormRegionFrom( region_name.c_str(), "permeability", 1.0e-21, 1.0e-19 );

    // Change the porosity to 0.3 and permeability 1.0e-13 to everywhere in the Region
    model.Region(region_name.c_str()).InputPropertyValue( "porosity", makeScalar(PLAIN,0.3) );
    model.Region(region_name.c_str()).InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-13) );

    // -----------------------------------------------------------------------------------------
    // 5.0 Interrelation to compute hydraulic conductivity from permeability and fixed viscosity
    // -----------------------------------------------------------------------------------------
    ConstantFactor<2U,divides>  conductivity( p_ref, "conductivity", "permeability", 0.001 ); // viscosity 1 cp = 0.001 Pa s

    // calculate values
    model.Apply( conductivity );

    // output the range of the result variable
    printRangeOfVariable( model, "conductivity" );

    // ------------------------------------------------------------------------------------------
    // 6.0 Setting up an FE algorithm to solve the diffusion equation 0 = div(K grad p) + S
    //     p = fluid pressure
    //     K = k/mu = hydraulic conductivity (from above)
    //     S = volumetric source term
    //
    //     We solve the discretised equation full implict as
    //
    //     [K]{p} = {S}
    //
    // ------------------------------------------------------------------------------------------

    // create the CSMP FE Algorithm with SAMG solver
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    PDE_Integrator<2U,Region>  fluid_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<2U,Region>  fluid_pressure(linear_solver);
#endif

    // LHS stiffness matrix                                 operand         basis function    test function
    NumIntegral_dNT_op_dN_dV<2U,Element<2U> >  stiffness_matrix( p_ref, "conductivity", "fluid pressure", "fluid pressure" );

    // RHS mass matrix for integrating source term
    NumIntegral_NT_op_N_dV<2U,Element<2U> >    source_term( p_ref, "fluid volume source", "fluid pressure" );

    // use lumped formulation for all mass matrices (i.e., diagonalise matrices)
    source_term.LumpedFormulation(true);

    // define a post-processing step that computes the velocity in each finite element by solving Darcy's law
    VelocityAndVolumeFlux<2U,Element<2U> >     velo( model, "conductivity", "porosity", "fluid pressure", true );

    // now add each FE operation (i.e., PDE Operator) to the FE algorithm
    fluid_pressure.Add( &stiffness_matrix );
    fluid_pressure.Add( &source_term );
    fluid_pressure.AddPostProcess( &velo );

    // compute fluid pressure field
    model.Apply( fluid_pressure );

    // reset algorithm and return memory because it is no longer needed (important for large models)
    fluid_pressure.Reset();

    // show results
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );

    // -----------------------------
    // 7.0 Output initial conditions
    // -----------------------------
    // output to vtk and matlab
    VTK_Interface<2U>     vtk_output;


    // to VTK files
    vtk_output.OutputDataToVTK( model, "concentration",  "concentration",  0 );
    vtk_output.OutputDataToVTK( model, "fluid_pressure", "fluid pressure", 0 );
    vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       0 );
    vtk_output.OutputDataToVTK( model, "volume_flux",    "volume flux",    0 );

    // -------------------------------------------------------------
    // 8.0 Construct the finite volume grid and transport algorithms
    // -------------------------------------------------------------
    Standard_IO_Handler  stdio;
    double64 cfl_multiplier; // for explicit transport, CFL should be mulitplied by 0.5
    cerr <<"\nEnter CFL multiplier: (suggested value 0.3 ~ 0.7)";
    cin  >> cfl_multiplier;    

    // -----------------------
    // 9.0 Time Loop Variables
    // -----------------------
    // define some constant variables
    const double64    hour(3600.0);
    double64	hours = 720.0;
    double64	max_time=hours*hour;    
    double64	time_increment(2.0*hour); // timestep 2 hours
    const long	save_frequency(hours/40);  // write results to file
    
    double64	model_time(0.);	// global time for simulated runtime
    size_t	save_counter(1), time;     
    double64	PEP_factor=1.0;
    size_t	n_threads=1;    

    bool  DES = stdio.YesNo("Do you want to solve the advection equation with DES (y=DES, n=TDS)");

    if (DES)
    {
    // ----------------------------------------
    // 10.0 Transient loop using DES simulation
    // ----------------------------------------
    cerr <<"\n\nmain: Starting DES simulation: "<<endl;
    clock_t T_begin= clock();    
    DESAdvectionDiffusion<2U>* DES_transport = new DESAdvectionDiffusion<2U>(model, "Model", cfl_multiplier, PEP_factor, false);   
   
    while ( model_time < max_time )
    {
         // compute advection of solute with DES
         DES_transport->AdvectVariable_DES( model_time+time_increment, n_threads);   
         //DES_transport->AdvectVariable_TDS( time_increment, model_time+time_increment, cfl_multiplier, PEP_factor);
                 
         // increment time
         model_time += time_increment;

         // output variables
         if ( save_counter == save_frequency ) {
              time = static_cast<long>(model_time/hour);
              // to VTK files
              vtk_output.OutputDataToVTK( model, "DES_concentration", "concentration", time );
              vtk_output.OutputDataToVTK( model, "DES_update_count", "update count", time );
              vtk_output.OutputDataToVTK( model, "DES_variation_rate_count", "rate count", time );
              vtk_output.OutputDataToVTK( model, "DES_schedule_count", "schedule count", time );
              vtk_output.OutputDataToVTK( model, "DES_synchronization_count", "synchronize count", time );
                          
              save_counter = 0;
          }
         save_counter++;
         
         // runtime info
         cout <<"\n\nmain: RUNTIME (HRS): "<< model_time/hour << endl << endl;
    }
    cerr <<"\nmain: Finished DES simulation, using time (sec): "<< (clock() - T_begin)/double64(CLOCKS_PER_SEC) << endl;
    }

    else
    {
    // ----------------------------------------
    // 11.0 Transient loop using TDS simulation
    // ----------------------------------------
    cerr <<"\n\nmain: Starting TDS simulation: "<<endl;
    clock_t	T_begin= clock();   
    DESAdvectionDiffusion<2U>* TDS_transport = new DESAdvectionDiffusion<2U>(model, "Model", cfl_multiplier, PEP_factor, false);     
    
    while ( model_time < max_time )
    {
         // compute advection of solute with DES
         TDS_transport->AdvectVariable_TDS( time_increment, n_threads);
                 
         // increment time
         model_time += time_increment;

         // output variables
         if ( save_counter == save_frequency ) {
              time = static_cast<long>(model_time/hour);
              // to VTK files
              vtk_output.OutputDataToVTK( model, "TDS_concentration", "concentration", time );
                          
              save_counter = 0;
          }
         save_counter++;
         
         // runtime info
         cout <<"\n\nmain: RUNTIME (HRS): "<< model_time/hour << endl << endl;
    }
    cerr <<"\nmain: Finished TDS simulation, using time (sec): "<< (clock() - T_begin)/double64(CLOCKS_PER_SEC) << endl;
    }
    
    // terminate
    cerr <<"\nmain: That's it..."<< endl;

} // Run()

} // csmp
