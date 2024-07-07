#include "FiniteVolumeStencilSpeed_Test.h"

#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"

#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "VelocityAndVolumeFlux.h"

#include "finiteVolumeAuxiliaryFunctions.h"
#include "ExplicitTransport.h"

#include "VTK_Interface.h"

#include "Timer.hpp"

using namespace std;

namespace csmp {

// helper
static bool isThereFileCalled(const string& name ) {
    ifstream f(name.c_str());
    return f.good();
}

void FiniteVolumeStencilSpeed_Test::run()
  {
      Timer timer;
      double ansys_build_time, native_build_time;

//      const string model_name{"Clair"};    // 799206 elmts, 145178 nodes
      const string model_name{"prism_test"}; //  53821 elements, 10831 nodes, 3040 faces
      
      // checking whether there already is a native csmp file set in place
      if ( !isThereFileCalled( model_name + "_variables.dat") )
        {
          // creating the model from ANSYS .asc and .dat files and converting it to CSMP native binary fileset
          ANSYS_Model3D ansys_model( model_name.c_str(), model_name.c_str(), "VariableSet_TracerTransfer-variables.txt", true );
          cout << "\n"<<"FiniteVolumeStencilSpeed_Test::run: time taking to build model from ANSYS: " << (ansys_build_time=timer.StopClock()) << "\n\n\n";
          printModelDimensions( ansys_model );
          // saving model to disk
          ansys_model.OutputToBinaryFile(model_name.c_str());
        }
      // reading from binary file
      timer.Start();
      Model<3U> model( model_name );
      cout << "\nTime taking to build model '"<< model_name <<"' from CSMP binary file set: " << (native_build_time=timer.StopClock()) << "\n";
      //model.Database().IndexTrackerOut();
      cout <<"\n"<<"datadepth of current indices in model:"<< endl;
      const IndexTracker& indices = model.Database().VariableIndexes();
      LocalVariables::int_type max_index{0ul};
      for ( auto it=indices.IndicesBegin(); it!=indices.IndicesEnd(); ++it )
        max_index = std::max( max_index, static_cast<LocalVariables::int_type>((*it).first->dataDepth + (*it).first->dataOffset) );
        cout <<"\t\t"<< max_index << endl;

      // compute volume, pore volume and prescribed 'total velocity'
      // -----------------------------------------------------------
      Region<3U>&  model_domain = model.Region("Model");
      const bool initialize_flux( true ); // prescibed 'total velocity'
      initializeFiniteVolumeProperties( model, model_domain, initialize_flux );
      AssignFlowProperties( model );

      // computing flux balance for computed divergence free 'total velocity' field
      // --------------------------------------------------------------------------
      // computing divergence free 'total velocity' field and 'facet flux'
      DivergenceFreeTotalVelocityField( model, 2.0e5 ); // delta p = 2 bars over 10 meters
      {
        ExplicitTransport<3U>  transport( model, "Model" );
        _equal( transport.IncomingVolumetricFlow(),
                transport.OutgoingVolumetricFlow(), numeric_limits<double>::epsilon() * transport.IncomingVolumetricFlow() * 1e5 );
      }
      // speed test 1: flow through model with TVD concentration
      // -------------------------------------------------------
      // prescribed would set velocity vector to 1,0,0
      const bool prescribed_velocity(false);
      timer.Start();
      TestFlowThroughModel( model, prescribed_velocity );
      cout << "\n"<<"FiniteVolumeStencilSpeed_Test::run: "<< timer.StopClock() << endl;
      
    } // end run






    /// thickness, permeability, porosity, total velocity
void  FiniteVolumeStencilSpeed_Test::AssignFlowProperties( Model<3U>& model )
 {
     // property assigment
    model.InputPropertyValue( "thickness",     makeScalar(ANY,1.0) );
    model.InputPropertyValue( "permeability",  makeScalar(ANY,1.0e-13) );
    model.InputPropertyValue( "conductivity",  makeScalar(ANY,1.0e-10) ); // dyn visc = 1.0e-3
    model.InputPropertyValue( "porosity",      makeScalar(ANY,0.25) );
    model.InputPropertyValue( "concentration", makeScalar(ANY,0.) );
    model.InputPropertyValue( "fluid volume source", makeScalar(ANY,0.) );
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(ANY,0.) );
    model.InputPropertyValue( "fluid viscosity", makeScalar(ANY,1.0e-3) );
    
    if ( string(model.Name()) == "prism_test" ) {
         Region<3U>& fracdomain = model.Region("FRAC_VOLUMES");
         fracdomain.InputPropertyValue( "permeability",  makeScalar(ANY,1.0e-10) );
         fracdomain.InputPropertyValue( "conductivity",  makeScalar(ANY,1.0e-7) ); // dyn visc = 1.0e-3
         fracdomain.InputPropertyValue( "porosity",      makeScalar(ANY,1.) );
     }
    
    if ( verbose_ ) {
         VTK_Interface<3U> vtk_output;
         vtk_output.OutputDataToVTK( model, "test_output", "porosity", 0 );
      }
    // call of testee
    // here the facet fluxes as precomputed are used
    VectorVariable<3U> velo(ANY,ANY,ANY, 1., 1., 1. );
    velo.EuclideanNormalize(); // to 1.
    model.InputPropertyValue( "total velocity", velo );
    if ( verbose_ ) {
         cout <<"\nExplicitTransport_Test::AssignFlowProperties: 'total velocity' magnitude: "<< velo.Length() <<"\n";
         printRangeOfVariable( model, "permeability" );
         printRangeOfVariable( model, "conductivity" );
         printRangeOfVariable( model, "porosity" );
      }
      
} // end AssignFlowProperties




/**
    Applies uniform pressures on left and right side of the model,
    and solves for the steady-state fluid pressure distribution in the
    absence of fluid sources and sinks.
*/
void  FiniteVolumeStencilSpeed_Test::DivergenceFreeTotalVelocityField( Model<3U>& model, double delta_pf )
 {
    // 1.  Building the steady-state FE Algorithm "fluid_pressure"
    // -----------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps(0.);
    SAMG_Solver    solver(&settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
    PDE_Integrator<3U,Element>  fluid_pressure( solver );

    NumIntegral_dNT_op_dN_dV<3U> conductance( model.Database(), "conductivity", "fluid pressure",  "fluid pressure" );
    NumIntegral_NT_op_N_dV<3U>   source( model.Database(),  "fluid volume source", "fluid pressure" );
  /// @todo influx surface integral: NumIntegral_NT_op_N_dS<2U>   influx( model.Database(), "influx", "fluid pressure" );
    VelocityAndVolumeFlux<3U>    velocity( model,  "conductivity", "porosity", "fluid pressure", false );

    // add PDE_Operators and post-processor to the FE Algorithm
    fluid_pressure.Add( &conductance );
    fluid_pressure.Add( &source );
    //  fluid_pressure.Add( &influx );
    fluid_pressure.AddPostProcess( &velocity );


    // 2. Assign boundary conditions
    // -----------------------------
    // pressure range between 1 bar and (1 bar + delta_pf)
    const double bar(100325.);
    model.InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH,bar+delta_pf) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,bar) );


    // 3.  Pass the FE algorithm to the Model and reset it after the solution
    // ----------------------------------------------------------------------
    fluid_pressure.IntegrateOver( model.Region("Model") );
    model.CopyReplace( "velocity", "total velocity" );


    // 4.  Output the initial range of the variables
    // ---------------------------------------------
    if ( verbose_ ) {
          printRangeOfVariable( model, "fluid pressure" );
          printRangeOfVariable( model, "velocity" );
       }

    // 5.  Output the initial conditions to VTK
    // -----------------------------------------
#ifdef DEBUG
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0, true );
    vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       0, true );
#endif

 } // end






/**
    A single pressure solve followed by TDS  passive tracer advection simulation; no special flow laws.
*/
void  FiniteVolumeStencilSpeed_Test::TestFlowThroughModel( Model<3U>& model, bool prescribed_velocity )
 {
    AssignFlowProperties( model );
 
    // constant velocity field, left-to-right, velocity = 1m/s
    double  velo_magnitude(1.);
    if ( prescribed_velocity ) {
          VectorVariable<3U>  vc1(ANY,ANY,ANY,1.,0.,0.);
          model.InputPropertyValue( "velocity", vc1 ); // NB: transport scheme uses 'velocity'
      }
    const bool print_vmaximum{false};
    velo_magnitude = printRangeOfVariable( model, "velocity", print_vmaximum );
    
    // initial and boundary conditions for tracer transport
    // - concentration
    // - constraints at boundary
    // - initial amount of tracer in the system
    model.InputPropertyValue( "concentration", makeScalar(ANY,0.) );
    const double inlet_concentration(3.);
    model.InputBoundaryValue( LEFT, "concentration", makeScalar(DIRICH,inlet_concentration) );
    
#ifdef DEBUG
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "concentration", "concentration", 0, true );
#endif

    // uses the current velocity field, to initialise facet fluxes in construction
    ExplicitTransport<3U>  transport( model, "Model" );
 
     if ( verbose_ ) {
         cout <<"\nExplicitTransport_Test::TestFlowThroughModel: key variable ranges in current model '"<< model.Name() <<"'\n";
         printRangeOfVariable( model, "FV pore volume" );
         printRangeOfVariable( model, "concentration" );
      }
 
    // assumes flow and model long axis are aligned with the X-axis
    Point<3U> xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max ); 
    const double model_length(xyz_max[0]-xyz_min[0]), xsect_area((xyz_max[1]-xyz_min[1]) * (xyz_max[2]-xyz_min[2]));
    const double time_interval( (model_length/velo_magnitude) * 1e-4 ); // ~10-m travel distance (4 O(M) is difference between fastest and slowest flow speed
    double       duration(0.); // calculated from velocity and model length

    // 0. testing whether inflow and outflow from the model have the expected values
    // -----------------------------------------------------------------------------
    if ( prescribed_velocity ) {
          cout <<"\nrun: prescribed_velocity 'velocity' magnitude: "<< velo_magnitude << endl;
          const double expected_volume_flux(xsect_area * velo_magnitude);
          _equal( transport.IncomingVolumetricFlow(), expected_volume_flux, numeric_limits<double>::epsilon() * expected_volume_flux );
          _equal( transport.OutgoingVolumetricFlow(), expected_volume_flux, numeric_limits<double>::epsilon() * expected_volume_flux );
      }

    // 1. tracer tranport and conservation tests
    // -----------------------------------------
    // 1.1 getting some tracer into model
    transport.AdvectVariable( time_interval );
    duration += time_interval;
#ifdef DEBUG
    vtk_output.OutputDataToVTK( model, "concentration", "concentration", 1, true );
#endif

    // 1.2 switching supply off and transporting more
    model.InputBoundaryValue( LEFT, "concentration", makeScalar(DIRICH,0.) );
    transport.AdvectVariable( time_interval );
    duration += time_interval;
#ifdef DEBUG
    vtk_output.OutputDataToVTK( model, "concentration", "concentration", 2, true );
#endif

    // 1.3 transporting for trice the time
    transport.AdvectVariable( time_interval );
    duration += time_interval;
 #ifdef DEBUG
   vtk_output.OutputDataToVTK( model, "concentration", "concentration", 3, true );
 #endif
    const bool print_maximum{true};
    const double max_concentration = printRangeOfVariable( model, "concentration", print_maximum );
    _test( max_concentration <= inlet_concentration );

    // 1.4 transporting tracer across outflow boundary, verifying that there is no build up
    transport.AdvectVariable( time_interval * 2. );
    duration += time_interval * 2.;
    _test( printRangeOfVariable( model, "concentration", print_maximum ) <= inlet_concentration );
#ifdef DEBUG
    vtk_output.OutputDataToVTK( model, "concentration", "concentration", 4, true );
#endif

 } // end TestFlowThroughModel




  } // end csmp
