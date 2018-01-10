#include "TopographyDrivenFlow_Example.h"

#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"

#include "TRIANGLE_Interface.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"

#include "Integral_NT_op_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "Integral_NT_op_dNi_dV.h"

#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"

#include "PropertyHandle.h"

// transport calculation
#include "Standard_IO_Handler.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"

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
  AddDescription( "source in: TopographyDrivenFlow_Example.cpp" );
  AddRequirement( "file set: 'topo.1'");
}

void computeNodalVelocityAndVolumeFlux( Model<2U>& );
void simulateTracerTransport( Model<2U>&  model, VTK_Interface<2U>& vtk_output, ComputationalSettings& settings );
 

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
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  // ------------------------------------------------------------
  // 0. building Region object
  // ------------------------------------------------------------
  /*
    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;
    char  file_name[200];
    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
    cin >> file_name;
    mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );
    Model<2U>  model( mesh_container, "example16.txt" );
    printModelDimensions( model, true );
  */

    char  file_name[200];
    cout <<"\nmain: Enter name of 'ANSYS' binary input file set: ";
    cin >> file_name;
    ANSYS_Model2D  model( file_name, "example16.txt" );

    // give the model dimensions
    printModelDimensions( model, true );



  // ---------------------------------------------------------------
  // 1. assign initial values & Dirichlet boundary conditions
  // ---------------------------------------------------------------
   InputDataManager<2U>  model_configuration;
   ComputationalSettings settings;

   model_configuration.ConfigureFromFile( model, file_name,
                                          false,           ///< region name from parameter range
                                          true,            ///< default property values
                                          true,            ///< regional property values
                                          true,            ///< boundary conditions for box-shaped model
                                          true,           ///< essential conditions for regions
                                          true,           ///< boundary conditions for arbitrary-shaped model
                                          settings );

  // --------------------------------------------------------------
  // 2. calculate hydraulic conductivity from permeability
  // --------------------------------------------------------------
    printRangeOfVariable( model, "permeability" );
    const double64  fluid_viscosity(1.6e-3); // Pa s-1
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
    PDE_Integrator<2U,Region>  steady_state_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<2U,Region>  steady_state_pressure(linear_solver);
#endif

    // conductance matrix [K] on the lefthand side
    NumIntegral_dNT_op_dN_dV<2U,Element<2U> >    conductance( model.Database(), "conductivity", "absolute fluid pressure",  "absolute fluid pressure" );
    // source vector {Q} on the righthand side
    NumIntegral_NT_op_N_dV<2U,Element<2U> >      source( model.Database(),  "fluid volume source", "absolute fluid pressure" );
    // gravity term on righthand side
    NumIntegral_NT_op_dNi_dV<2U,Element<2U> >    gravity( model.Database(),  "fluid density",  "conductivity", "absolute fluid pressure" );

    // add PDE_Operators to the FE Algorithm
    steady_state_pressure.Add( &conductance );
    steady_state_pressure.Add( &source );
    steady_state_pressure.Add( &gravity );

    model.Apply( steady_state_pressure );
    printRangeOfVariable( model, "absolute fluid pressure" );


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
    computeNodalVelocityAndVolumeFlux( model );


  // ------------------------------------------------------------
  // 6. Output of results
  // ------------------------------------------------------------
    TextInterface  text_output;

    // outputting the logarithm of conductivity (this modifies the original variable)
    PropertyHandle<2U>  K( model, "conductivity" );
    K.Log10();

    text_output.OutputDataAsTextColumns( model, "velocity",       "velocity",       1 );
    text_output.OutputDataAsTextColumns( model, "log10-conductivity",   "conductivity",   1 );
    text_output.OutputDataAsTextColumns( model, "hydraulic-head", "hydraulic head",   1 );
    text_output.OutputDataAsTextColumns( model, "absolute-fluid-pressure", "absolute fluid pressure", 1 );


    VTK_Interface<2U>  vtk_output;

    vtk_output.OutputDataToVTK( model, "log10-hydraulic-conductivity",  "conductivity", 1 );
    vtk_output.OutputDataToVTK( model, "porosity", "porosity", 1 );
    vtk_output.OutputDataToVTK( model, "velocity",                "velocity", 1 );
    vtk_output.OutputDataToVTK( model, "hydraulic-head",          "hydraulic head", 1 );
    vtk_output.OutputDataToVTK( model, "absolute-fluid-pressure", "absolute fluid pressure", 1 );

    vtk_output.OutputDataToVTK( model, "volume-flux", "volume flux", 1 );
    vtk_output.OutputDataToVTK( model, "transport-velocity", "transport velocity", 1 );
    vtk_output.OutputDataToVTK( model, "nodal-transport-velocity", "nodal transport velocity", 1 );
  
  
  // ------------------------------------------------------------
  // 7. Tracer transport simulation
  // ------------------------------------------------------------
    simulateTracerTransport( model, vtk_output, settings );
    printRangeOfVariable( model, "concentration" );

    cout <<"\nmain: That's it..."<< endl;

} // Run()





/**
    Post-processing for streamline routing and visualisation
*/
void computeNodalVelocityAndVolumeFlux( Model<2U>& model )
 {
    const csmp::Index  vD_key = model.Database().StorageKey("velocity");
    const csmp::Index  vf_key = model.Database().StorageKey("volume flux");
    const csmp::Index  po_key = model.Database().StorageKey("porosity");
    const csmp::Index  vn_key = model.Database().StorageKey("nodal transport velocity");
    const csmp::Index  vt_key = model.Database().StorageKey("transport velocity");
    VectorVariable<2U> velo;
   
    Region<2U>& model_domain(model.Region("Model"));
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         (*it)->Read( vD_key, velo );
         // volume flux
         (*it)->Store( vf_key, makeScalar(PLAIN,velo.Length()) );
         // transport velocity
         double64 porosity = (*it)->Read( po_key );
         velo /= porosity;
         (*it)->Store( vt_key, velo );
      
      }
   
    model.ExtrapolateElementToNodeProperty( "transport velocity", "nodal transport velocity" );
   
 } // end computeNodalVelocityAndVolumeFlux






/**
    Carry out a tracer transport calculation.
*/
void simulateTracerTransport( Model<2U>&  model, VTK_Interface<2U>& vtk_output, ComputationalSettings& settings )
 {
    double64 pf_multiplier, adv_multiplier;
    settings.EstablishMultipliers( pf_multiplier, adv_multiplier );
 
    // -------------------------------------------------------------
    // 8.0 Construct the finite volume grid and transport algorithms
    // -------------------------------------------------------------
    Standard_IO_Handler  stdio;
    double64 cfl_multiplier(0.5); // for explicit transport, CFL should be mulitplied by 0.5

    // query user if implicit or explicit FV scheme should be used
    cerr << "\nsimulateTracerTransport: Hit enter to continue..." << endl;
    bool  implicit = stdio.YesNo("simulateTracerTransport: Do you want to solve the advection equation implicitly (n=explicitly)");
   
    model.Region("Model").UpdateMemberIndexes();

    // NULL pointer to FV transport algorithm
    NodeCenteredFiniteVolumeTransport<2U>*  transport(NULL);

    if ( implicit ) {
        // implicit finite volume scheme
        transport = new NodeCenteredFiniteVolumeTransport<2U>( "Model",                         // default region is named 'Model'
                                                               model,
                                                               "porosity",                      // porosity multiplier
                                                               "concentration",                 // advected variable
                                                               "velocity",                      // advecting variable
                                                               "concentration source",          // source terms
                                                               false,                           // second-order accuracy
                                                               false );                         // higher-order temporal approximation should not be used

        cfl_multiplier *= adv_multiplier; // overstep the CFL critertion by a factor 10000 in implict method
      }
    else {
        // query user if higher order FV scheme with slope reconstruction should be used
        bool second_order = stdio.YesNo("simulateTracerTransport: Do you want to use a second order accurate FV advection algorithm in space (n=first order accuracy)");

        // explicit finite volume scheme
        transport = new ExplicitNodeCenteredFiniteVolumeTransport<2U,ExplicitStencilProcessor>
                                                                  ( "Model",                          // default region is named 'Model'
                                                                     model,
                                                                     "porosity",                      // porosity multiplier
                                                                     "concentration",                 // advected variable
                                                                     "transport velocity",            // advecting variable
                                                                     "concentration source",          // source terms
                                                                     second_order );                  // second-order accuracy
      }

    // -----------------------
    // 9.0 Time Loop Variables
    // -----------------------
    // define some constant variables
    const double64    hour(3600.0), time_tolerance(100.);
    double64          time_increment(8640000.); // 100 days
    settings.TimeIncrement(time_increment);
    size_t	          time;

    // -----------------------
    // 10.0 Transient loop
    // -----------------------
    // define the variables used throughout the simulation
    double64  model_time(0.);                                 // global time for simulated runtime

    while ( model_time <= settings.Duration() )
      {
         // getting the timing right so that output time are honoured
         time_increment = settings.TimeIncrement();
         if ( (model_time + time_increment) > settings.Duration() ) time_increment = settings.Duration() - model_time;

         time_increment = (settings.TimeToNearestOutputTime(model_time) < time_increment) ?
                           settings.TimeToNearestOutputTime(model_time) : time_increment;

         // compute advection of solute
         // first boolean: check and correct for divergence of flow field
         // second boolean: update velocity field (set to true if it changes in time)
         transport->AdvectVariable( time_increment, cfl_multiplier, true, false );

         // increment time
         model_time += time_increment;

         // output variables
         if ( settings.IsOutputTime( model_time, time_tolerance ) ) {
              time = static_cast<long>(model_time/hour);
              // to VTK files
              printRangeOfVariable( model, "concentration" );
              vtk_output.OutputDataToVTK( model, "concentration", "concentration", time );
          }

         // runtime info
         cout <<"\n\nsimulateTracerTransport: RUNTIME (HRS): "<< model_time/hour << endl << endl;
      }

    // final output
    // to VTK files
    time = static_cast<long>(model_time/hour);
    vtk_output.OutputDataToVTK( model, "concentration", "concentration", time );

} // end transport simulation



} // csmp















