#include "DESTwoPhaseFlow3D_Example.h"

// the CSMP model
#include "ANSYS_Model3D.h"

// FE algorithm
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// FV algorithms
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "TwoPhaseExplicitNodeCenteredFVTransport.h"
#include "TwoPhaseDESTransport.h"
#include "TwoPhaseTwoComponentDESTransport.h" 
#include "ExplicitStencilProcessor.h"

// relative permeability calculations
#include "BrooksCorey.h"
#include "CO2H2O_FunctionsModule1.h"

// monitoring individual regions
#include "RegionMonitor.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "MatlabInterface.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp{

void DESTwoPhaseFlow3D_Example::Specifications()
{
  SetTitle( "3D two phase flow simulation using DES (discrete event simulation)");
  SetDifficulty( 3 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "Qi Shao" );
  AddDescription( "3D two phase flow simulation with discrete event simulation (DES) or time-driven simulation (TDS).");
  AddDescription( "Output is written to to VTK files");
  AddRequirement( "source files: 'DESTwoPhaseFlow3D_Example.cpp' and '*.h'" );
  AddRequirement( "fracs4: files .dat, .asc, -regions.txt, -configuration.txt.");  
  AddRequirement( "variables(DES_2phase_variables.txt)" );
} 

/** 
    3D two phase flow simulation via CSMP's DES transport method 
    combining finite elements (for pressure) with finite volumes (for advection of non-wetting phase)

    Use models 'fracs4' (.dat, .asc, -regions.txt, -configuration.txt) as input file suites.
*/

void DESTwoPhaseFlow3D_Example::Run()
{
    // -------------------------------------
    // 0.0 Set variables used throughout the simulation
    // -------------------------------------

    double64 model_time =  0.0; // time

    // ---------------------------------------------------
    // 1.0 Create Model directly from ANSYS-ICEM mesh
    // ---------------------------------------------------
    string input_file;

    cerr<< "\nPlease enter the name of input mesh (default: fracs4):"<< endl;
    cin >> input_file;

    Standard_IO_Handler  stdio;
    bool  DES = stdio.YesNo("Do you want to solve the transport equation with DES? (y=DES, n=TDS)"); 
    bool  multi_component = stdio.YesNo("Do you want to perform multi-component transport?");
    double64 Courant_multiplier, PEP_parameter;
    cerr <<"\nEnter CFL multiplier (suggested value: 0.2) and PEP parameter (suggested value: 0.1)" << endl;
    cin >> Courant_multiplier >> PEP_parameter;  
     
    bool  with_capillary_spreading = stdio.YesNo("Do you want to include capillary effect (y/n)?"); 
    bool  with_gravity_forces = stdio.YesNo("Do you want to include gravity effect (y/n)?"); 
    
    ANSYS_Model3D                model( input_file.c_str(), "DES_2phase_variables.txt" );
    const PropertyDatabase<3>&   p_ref = model.Database();  

    // give the model dimensions
    printModelDimensions( model, true );

    // --------------------------------------------
    // 3.0 Configure the simulation from a file
    // --------------------------------------------
    InputDataManager<3U>  model_configuration;
    model_configuration.Configure_ANSYS_ModelFromFile( model, input_file.c_str() );

    // ---------------------------------------------------------------------
    // 4.0 Use the flow functions to compute the relative permeabilities
    // ---------------------------------------------------------------------
    CO2H2O_FunctionsModule1<3U> flowfunctions(model.Database());
    TwoPhaseDESTransport<3U,CO2H2O_FunctionsModule1>* DEStransport;
    if (!multi_component)
      DEStransport = new TwoPhaseDESTransport<3U,CO2H2O_FunctionsModule1>(model, "Model", with_capillary_spreading, with_gravity_forces, PEP_parameter, Courant_multiplier);
    else 
      DEStransport = new TwoPhaseTwoComponentDESTransport<3U,CO2H2O_FunctionsModule1>(model, "Model", with_capillary_spreading, with_gravity_forces, PEP_parameter, Courant_multiplier);

    computeTotalMobility( model, flowfunctions );

    // output the range of the result variable
    printRangeOfVariable( model, "total mobility permeability product" );

    // ------------------------------------------------------------------------------------------
    // 5.0 Setting up an FE algorithm to solve the diffusion equation 0 = div(lambda_t grad p)
    //
    //     lambda_t = total mobility
    //     p = fluid pressure
    // ------------------------------------------------------------------------------------------

    // create a steady-state CSMP FE Algorithm using a high-level class
    SteadyStateDiffusor<3U,Region> fluid_pressure( model,
                                                          "total mobility permeability product",
                                                          "fluid pressure",
                                                          "fluid volume source" );

    // operation to compute velocity
    VelocityAndVolumeFlux<3U,Element<3U> >  velo( model,
                               "total mobility permeability product",
                               "porosity",
                               "fluid pressure", true, "total velocity" );

    // add velocity calculation as post process
    fluid_pressure.AddPostProcess( &velo );

    // calculate initial pressure distribution
    fluid_pressure.ComputeSteadyState( model.Region("Model") );

    // show results
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "total velocity" );

    // -----------------------------
    // 6.0 Output initial conditions
    // -----------------------------
    VTU_Interface<3U>  vtu(model);  // binary VTK output, creates much smaller files than VTK

    vtu.OutputDataToVTU( "volume_flux",    "nodal volume flux", "Model", 0 );
    vtu.OutputDataToVTU( "fluid_pressure", "fluid pressure",    "Model", 0 );
    vtu.OutputDataToVTU( "saturation oil", "saturation carbonic phase",    "Model", 0 );
    vtu.OutputDataToVTU( "fluid_velocity", "total velocity",    "Model", 0 );
    if (multi_component) { 
      vtu.OutputDataToVTU( "dissolved CO2", "dissolved CO2",    "Model", 0 ); 
      vtu.OutputDataToVTU( "evaporated water", "evaporated water",    "Model", 0 );   
    }

    // -------------------------------------------------------------
    // 7.0 Construct the finite volume grid and DES transport algorithms
    // -------------------------------------------------------------

    model.InputPropertyValue( "nodal fluid volume source", makeScalar( PLAIN,0.0) ); // no FV sources/sinks
    //model.InputPropertyValue( "diffusivity coefficient carbonic phase", makeScalar(PLAIN,1.0e-25) ); // very small value
    //model.InputPropertyValue( "diffusivity coefficient aqueous phase", makeScalar(PLAIN,1.0e-25) ); // very small value

    // -----------------------
    // 8.0 Time Loop Variables
    // -----------------------

    // define some constant variables
    const double64    day(86400.0);
    const double64    max_time (60.0*day);     // run for 60 days
    double64          time_increment(0.6*day);      // timestep 0.6 day
    const long        save_frequency(3);       // write results to file every 3 days
    size_t	      time, save_counter(1);    
    
    // -----------------------
    // 9.0 Transient loop
    // -----------------------
    double64 solving_time = 0.;
    clock_t  T_begin;    
    while ( model_time < max_time )
      {
         // compute advection of phases
         T_begin = clock();
         if (DES) DEStransport->AdvectVariable_DES( model_time+time_increment, 1);
         else DEStransport->AdvectVariable_TDS( time_increment );
         solving_time += clock() - T_begin;
         
         // increment time
         model_time += time_increment;

         // compute 2phase properties
         computeTotalMobility( model, flowfunctions );

         // compute pressure
         fluid_pressure.ComputeSteadyState( model.Region("Model") );

         // echo variables to screen
         printRangeOfVariable( model, "fluid pressure" );
         printRangeOfVariable( model, "total velocity" );
         printRangeOfVariable( model, "saturation carbonic phase" );

         // output variables
         if ( save_counter == save_frequency ) {
              time = static_cast<long>(model_time/day);
              // to VTU files
              if(DES) {
                  vtu.OutputDataToVTU( "DES_fluid_pressure", "fluid pressure",    "Model", time );
                  vtu.OutputDataToVTU( "DES_saturation_oil", "saturation carbonic phase",    "Model", time );
                  vtu.OutputDataToVTU( "DES_volume_flux",    "nodal volume flux", "Model", time );
                  vtu.OutputDataToVTU( "DES_fluid_velocity", "total velocity",    "Model", time );
                  vtu.OutputDataToVTU( "DES_Update_count", "update count", "Model",  time );
                  vtu.OutputDataToVTU( "DES_CFL_multiplier", "cfl multiplier", "Model",  time );
                  //vtu.OutputDataToVTU( "DES_sn_shock", "shock saturation carbonic phase", "Model",  time );
                  vtu.OutputDataToVTU( "DES_sw_shock", "shock saturation aqueous phase", "Model",  time );
                  if (multi_component) { 
                    vtu.OutputDataToVTU( "DES_dissolved CO2", "dissolved CO2",    "Model", time );
                    vtu.OutputDataToVTU( "DES_evaporated water", "evaporated water",    "Model", time );
                  }                   
              } else {
                  vtu.OutputDataToVTU( "TDS_fluid_pressure", "fluid pressure",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_saturation_oil", "saturation carbonic phase",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_volume_flux",    "nodal volume flux", "Model", time );
                  vtu.OutputDataToVTU( "TDS_fluid_velocity", "total velocity",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_Update_count", "update count", "Model",  time );   
                  if (multi_component) {   
                    vtu.OutputDataToVTU( "TDS_dissolved CO2", "dissolved CO2",    "Model", time );
                    vtu.OutputDataToVTU( "TDS_evaporated water", "evaporated water",    "Model", time );
                  }                       
              }   
              save_counter = 0;
         }
         save_counter++;

         // runtime info
         cout <<"\n\nmain: RUNTIME (DAYS): "<< model_time/day << endl << endl;

      }

    // clocking the runtime
    if(DES) cerr << "\nmain: DES transport uses " << static_cast<double64>(solving_time/CLOCKS_PER_SEC) << " seconds " << endl;
    else    cerr << "\nmain: TDS transport uses " << static_cast<double64>(solving_time/CLOCKS_PER_SEC) << " seconds " << endl;
    
    // terminate
    cerr << "\nmain: That's it..."<< endl;

} // Run()


void DESTwoPhaseFlow3D_Example::computeTotalMobility( Model<3U>& mdl, CO2H2O_FunctionsModule1<3U>& flowfunctions )
 {
     
    static const Region<3U>& mref = mdl.Region("Model"); 
    // keys to properties
    static Index  mobt_key(mdl.Database().StorageKey("total mobility permeability product"));
    static Index  sw_key(mdl.Database().StorageKey("saturation aqueous phase"));
    static Index  snw_key(mdl.Database().StorageKey("saturation carbonic phase"));
    //csmp::INDEX<TENSOR,ELEMENT>  k_key(mdl.Database().StorageKey("permeability"));   
    //assert( k_key.type == TENSOR );  
    static Index  k_key(mdl.Database().StorageKey("permeability")); 

    // 1. Computing the saturation of water = 1 - So
    //    loop over the FE nodes
    vector<Node<3U>* >::const_iterator nit;
    for ( nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
    {
        double64 sw = 1. - (*nit)->Read(snw_key);
        (*nit)->Store( sw_key, makeScalar((*nit)->Status(snw_key),sw));
    }

    
    // 2. Computing the multiphase flow properties
    //    loop over finite elements
    vector<Element<3U>* >::const_iterator eit;
    for ( eit = mref.ElementsBegin(); eit!= mref.ElementsEnd(); eit++ )
    {
        //total mobility
        double64 mob_t = flowfunctions.TotalMobility(*eit);
        /*
        TensorVariable<3U> K;
        e.Read( k_key, K );
        double k = K.Trace() / 3.;
        */
        double64 k = (*eit)->Read(k_key);
        mob_t *= k;  
        (*eit)->Store( mobt_key, makeScalar(PLAIN, mob_t) );
    } 

} // end computeTotalMobility

} // csmp
