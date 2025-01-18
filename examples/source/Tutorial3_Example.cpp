#include "Tutorial3_Example.h"

// the CSMP model
#include "ANSYS_Model2D.h"

// FE algorithm
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// FV algorithms
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "TwoPhaseExplicitNodeCenteredFVTransport.h"
#include "ExplicitStencilProcessor.h"

// relative permeability calculations
#include "BrooksCorey.h"

// monitoring individual regions
#include "RegionMonitor.h"

// interfaces
#include "InputDataManager.h"
#include "MatlabInterface.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"

using namespace std;

namespace csmp{


void Tutorial3_Example::Specifications()
{
  SetTitle( "Tutorial 3: Incompressible 2-phase flow" );
  SetDifficulty( 3 );
  SetCategory( "Tutorials (composite functionality)" );
  AddAuthor( "Sebastian Geiger" );
  AddDescription( "A CSMP main file that uses an ANSYS-ICEM FE mesh to build the Model and simulate" );
  AddDescription( "incompressible two-phase flow using an IMPES (Implicit Pressure Explicit Saturation) or" );
  AddDescription( "IMPIS (Implicit Pressure Implicit Saturation) approach based on the FE-FV scheme of CSMP." );
  AddDescription( "Output is written to Matlab files and individual regions are monitored");
  AddRequirement( "box2d_fault, tutorial3_variables.txt");
}

/** **********************************************************************************************

Example program that uses a mesh created with ANSYS-ICEM FE  to build the Model and simulate
incompressible two-phase flow using an IMPES (Implicit Pressure Explicit Saturation) or
IMPIS (Implicit Pressure Implicit Saturation) approach based on the FE-FV scheme of CSMP.
Output is written to Matlab files and individual groups are monitored

Tasks and exercises:

1. Generate a 2D ANSYS mesh with a new geometry, check comphg.wikidot.com for a small how-to

2. Compare results for implicit and explicit FV schemes and different timestepping.

3. Compare results for homogeneous porous media to the analytical solution for the
   Buckley-Leverett problem
   
4. Allow the user to define, at runtime (!), to choose between a different relperm models
   available in the two_phase directory
   
5. Design and program a Visitor that computes the total mobility and phase mobilities
   using the kr-model that the user defines. Get together as a team to compare your visitors
   and generate a final version for the SVN (NB: The visitor should check that the results
   are within range and give error statements if the wrong variable placements are used or
   variables are not defined)
   
6. Program a function that computes the capillary pressure at the element center as well as
   the derivative dpc/dS
   
7. Program a diffusion algorithm, using the FE method, to solve the non-linear diffusion
   equation for the capillary pressure, i.e. phi dS/dt = div * ( lambda_star dpc/dS grad S )
   
8. Compare simulations with and without capillary effects for explicit/implicit FV schemes
   and different time-steps. Analyse how the FV simulation of capillary spreading compares
   to your FE algorithm

**********************************************************************************************/
void Tutorial3_Example::Run()
{
    // -------------------------------------
    // 0.0 Set variables used throughout the simulation
    // -------------------------------------
    // clocking the runtime
    clock_t start(clock());

    double model_time(0.); // time

    /*
    // ---------------------------------------------------
    // 1.0 Create Model directly from ANSYS-ICEM mesh
    // ---------------------------------------------------
    string input_file, yesno;
    bool   implicit;

    cout<< "\nPlease enter the name of input mesh ( default: box2d_fault ):"<< endl;
    cin >> input_file;

    cout<< "\nPlease choose the type of FV numerical scheme ( Y/y - Implicit, N/n - Explicit ):"<< endl;
    cin >> yesno;

    if ( yesno == "Y" or yesno == "y" ) implicit = true;
    else                                implicit = false;

    ANSYS_Model2D                model( input_file.c_str(), "tutorial3_variables.txt" );
    */


    // ------------------------------------------------------------
    // 1.0 Load CSMP native format model
    // ------------------------------------------------------------
    string input_file;
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'box2d_fault':"<<endl;
    cin.ignore();
    getline(cin, input_file);
    if (input_file.length() == 0) input_file = "box2d_fault";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "tutorial3_variables.txt";
    string config_file = input_file;
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, input_file, variable_file, config_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    Model<2U>  model(input_file, variable_file);

    const PropertyDatabase<2>&   p_ref = model.Database();

    // give the model dimensions
    printModelDimensions( model, true );

    Standard_IO_Handler  stdio;
    bool  implicit = stdio.YesNo("\nPlease choose the type of FV numerical scheme ( Y/y - Implicit, N/n - Explicit ):");

    // --------------------------------------------
    // 3.0 Configure the simulation from a file
    // --------------------------------------------
    InputDataManager<2U>  model_configuration;
    model_configuration.Configure_ANSYS_ModelFromFile( model, input_file.c_str() );

    // ---------------------------------------------------------------------
    // 4.0 Use the Brooks-Corey model to compute the relative permeabilities
    // ---------------------------------------------------------------------
    BrooksCorey<2U> relperm_model( p_ref, "brooks corey parameter", "entry pressure" );

    computeTotalMobility( model, relperm_model );

    // output the range of the result variable
    printRangeOfVariable( model, "total mobility" );

    // ------------------------------------------------------------------------------------------
    // 5.0 Setting up an FE algorithm to solve the diffusion equation 0 = div(lambda_t grad p)
    //
    //     lambda_t = total mobility
    //     p = fluid pressure
    // ------------------------------------------------------------------------------------------

    // create a steady-state CSMP FE Algorithm using a high-level class
    SteadyStateDiffusor<2U,Element> fluid_pressure( model, "total mobility",
                                                          "fluid pressure",
                                                          "fluid volume source" );

    // operation to compute velocity
    VelocityAndVolumeFlux<2U>  velo( model,
                                                 "total mobility",
                                                 "porosity",
                                                 "fluid pressure", true );

    // add velocity calculation as post process
    fluid_pressure.AddPostProcess( &velo );

    // calculate initial pressure distribution
    fluid_pressure.ComputeSteadyState( model.Region("Model") );

    // show results
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );

    // -----------------------------
    // 6.0 Output conditions
    // -----------------------------
    MatlabInterface    matlab;
    VTU_Interface<2U>  vtu(model);

    // -------------------------------------------------------------
    // 7.0 Construct the finite volume grid and transport algorithms
    // -------------------------------------------------------------
    // NULL pointer to FV transport algorithm
    model.InstantiateFiniteVolumes();
    NodeCenteredFiniteVolumeTransport<2U>*  transport(NULL);

    if ( implicit ) {
        // implicit finite volume scheme (NB function arguments spelled out for completeness
        // constructor also has default arguments!)
        transport = new TwoPhaseImplicitNodeCenteredFVTransport<2U,StencilProcessor>( "Model", model,
                                                                                      "porosity",
                                                                                      "viscosity oil",
                                                                                      "viscosity water",
                                                                                      "density oil",
                                                                                      "density water",
                                                                                      "saturation water",
                                                                                      "saturation oil",
                                                                                      "velocity",
                                                                                      "nodal fluid volume source" );

        // overstep CFL: as a consequence, the BL front will be retarded because the transport is not linearized to the shock speed
        transport->CFL_Multiplier(10.);
     }
    else {
        // explicit finite volume scheme with 1st order accuracy in space (NB: constructor has default arguments as well!)
        // note the difference between the implicit case where the first constructor argument must be provided
        // while here the 1st argument needs only to be provided if a model sub-region is specified
        transport = new TwoPhaseExplicitNodeCenteredFVTransport<2U,ExplicitStencilProcessor>( "Model", model,
                                                                                              "porosity",
                                                                                              "diffusivity",
                                                                                              "saturation water",
                                                                                              "saturation oil",
                                                                                              "velocity",
                                                                                              "nodal fluid volume source",
                                                                                              false );

        dynamic_cast<TwoPhaseExplicitNodeCenteredFVTransport<2U,ExplicitStencilProcessor>*>(transport)->DisableCapillarySpreading(); // simulate via FE method
      }

    model.InputPropertyValue( "nodal fluid volume source", makeScalar( PLAIN,0.0) ); // no FV sources/sinks
    model.InputPropertyValue( "diffusivity", makeScalar(PLAIN,1.0e-25) ); // very small value


    // ----------------------------------------------------
    // 12. Monitor the properties of the different groups
    // ----------------------------------------------------
    // output names
    string ranges(input_file);
    ranges += "_range_properties";

    // statistics for integrated properties and min-max values
    list<std::string> integral_properties;  integral_properties.push_back("saturation oil");
    list<std::string> range_properties;     range_properties.push_back("total mobility");
                                            range_properties.push_back("velocity");

    // initalise monitoring class and record values at t = 0
    RegionMonitor<2U>  monitor( model, integral_properties, range_properties );
    monitor.ScalarPropertyRanges( model, model_time );
    monitor.ScalarPropertyIntegrals( model, model_time );

    // -----------------------
    // 8.0 Time Loop Variables
    // -----------------------

    // define some constant variables
    const double    day(86400.0);
    const double    max_time (100. * day);   // run for 100 days
    double          time_increment(day);     // timestep 1 day
    const long        save_frequency(5);       // write results to file every 10 day
    size_t	          time, save_counter(1);

    // -----------------------
    // 9.0 Transient loop
    // -----------------------
    while ( model_time < max_time )
      {
         // monitor the regions
         monitor.ScalarPropertyRanges( model, model_time );
         monitor.ScalarPropertyIntegrals( model, model_time );

         // output variables
         if ( save_counter == save_frequency ) {
              time = static_cast<long>(model_time/day);
              // to Matlab files
              matlab.Write2DMatlabFile(  model, "saturation_oil", "saturation oil", time );
              matlab.Write2DMatlabFile(  model, "fluid_pressure", "fluid pressure", time );
              // to VTU files
              vtu.OutputDataToVTU( "fluid_pressure", "fluid pressure",    "Model", time );
              vtu.OutputDataToVTU( "saturation oil", "saturation oil",    "Model", time );
              save_counter = 0;
          }
         save_counter++;

         // compute advection of phases
         transport->TransportPhase( relperm_model, time_increment );

         // increment time
         model_time += time_increment;

         // compute 2phase properties
         computeTotalMobility( model, relperm_model );

         // compute pressure
         fluid_pressure.ComputeSteadyState( model.Region("Model") );

         // echo variables to screen
         printRangeOfVariable( model, "fluid pressure" );
         printRangeOfVariable( model, "velocity" );
         printRangeOfVariable( model, "saturation oil" );

         // runtime info
         cout <<"\n\nmain: RUNTIME (DAYS): "<< model_time/day << endl << endl;

      }

    // final output
    time = static_cast<uint32_t>(model_time/day);
    matlab.Write2DMatlabFile(  model, "saturation_oil", "saturation oil", time );
    matlab.Write2DMatlabFile(  model, "fluid_pressure", "fluid pressure", time );
    vtu.OutputDataToVTU( "fluid_pressure", "fluid pressure",    "Model", time );
    vtu.OutputDataToVTU( "saturation oil", "saturation oil",    "Model", time );

    // output monitored region properties
    monitor.Out( ranges.c_str() );

    // clocking the runtime
    clock_t end(clock());
    cerr << "\nmain: CPU time was " << static_cast<double>((end-start)/CLOCKS_PER_SEC) << " seconds " << endl;

    // terminate
    cerr << "\nmain: That's it..."<< endl;

    filesystem::current_path("../../example_inputs/");

} // Run()


void  Tutorial3_Example::computeTotalMobility( Model<2U>& mdl, TwoPhaseModel<2U>& relperm )
 {
    // keys to properties
    const Index  mobt_key(mdl.Database().StorageKey("total mobility"));
    const Index  satw_key(mdl.Database().StorageKey("saturation water"));
    const Index  sato_key(mdl.Database().StorageKey("saturation oil"));
    const double  one(1.);
    double        sw;
    ScalarVariable  mob_t;

    // 1. Computing the saturation of water = 1 - So
    //    loop over the FE nodes
    Region<2U>& mref = mdl.Region("Model");
    for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
      {
         // read in So, compute Sw and store back to nodes along with the flag of So
         sw = one - (*nit)->Read( sato_key );
         (*nit)->Store( satw_key, makeScalar((*nit)->Status(sato_key),sw));
      }


    // 2. Computing the multiphase flow properties
    //    loop over finite elements
    for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ )
      {
         // 1. setting up the relative permeability model
         // ---------------------------------------------
         relperm.Initialize( *(*eit) );
         relperm.InitializeForBaryCenter( *(*eit) );
         relperm.EffectiveSaturation();

         // 2. total mobility
         // -----------------
         mob_t = relperm.TotalMobility();
         (*eit)->Store( mobt_key, mob_t );
      }


} // end compute2PhaseFlowProperties

} // csmp
