#include "DES2PhaseSlightlyCompressibleFlow_Example.h"

// FV algorithms
#include "DES2PhaseTransport.h"
#include "DES2PhaseSlightlyCompressibleTransport.h"

// monitoring individual regions
#include "RegionMonitor.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"

// solvers
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

// PDE
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"


using namespace std;

namespace csmp{

void DES2PhaseSlightlyCompressibleFlow_Example::Specifications()
{
  SetTitle( "Two phase slightly compressible flow simulation using DES (discrete event simulation)");
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "Qi Shao" );
  AddDescription( "Two phase slightly compressible flow simulation with discrete event simulation (DES) or time-driven simulation (TDS).");
  AddDescription( "Output is written to to VTK files");
  AddRequirement( "source files: 'DES2PhaseSlightlyCompressibleFlow_Example.cpp' and '*.h'" );
  AddRequirement( "2D: box2d_fault (CSMP native binary files), box2d_fault-configuration.txt, DES_2phase_variables.txt");
  AddRequirement( "3D: fracs4 (CSMP native binary files), frac4-configuration.txt, DES_2phase_variables.txt");
} 

/** 
    Two phase slightly compressible flow simulation via CSMP's DES transport method
    combining finite elements (for pressure) with finite volumes (for advection of non-wetting phase)

    2D case: use models 'box2d_fault' (CSMP native binary files, -configuration.txt) as input file suites.
    3D case: use models 'fracs4' (CSMP native binary files, -configuration.txt) as input file suites.
*/

void DES2PhaseSlightlyCompressibleFlow_Example::Run() {
    // reading in a csmp native format model, first determining whether it will be 2 or 3 dimensional.
    uint32_t dimension;
    cout << "\nPlease enter the dimension of the model (2 for 2D, 3 for 3D):" << endl;
    cin >> dimension;

    string model_name;
    //2D case
    if (dimension == 2U) {
      cout << "\nPlease enter the name of input model, or press ENTER to use the default model 'box2d_fault':" << endl;
      cin.ignore();
      getline(cin, model_name);
      if (model_name.length() == 0) model_name = "box2d_fault";
      //3D case
    } else if (dimension == 3U) {
      cout << "\nPlease enter the name of input model, or press ENTER to use the default model 'fracs4':" << endl;
      cin.ignore();
      getline(cin, model_name);
      if (model_name.length() == 0) model_name = "fracs4";
    } else {
      throw csmp::Exception(ERROR, "input dimension of model", "must be 2 or 3");
    }

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    file_name += ("_" + to_string(dimension) + "D");
    string variable_file = "DES_2phase_variables.txt";
    string config_file = model_name + "(DES_2phase)";
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);

    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    if (dimension == 2U) {
      Model<2U> model(model_name, variable_file);
      RunSimulation(model);
    } else if (dimension == 3U) {
      Model<3U> model(model_name, variable_file);
      RunSimulation(model);
    }
}


/**
    Performs either discrete-event simulation (DES) or time-driven simulation (TDS) of 2-phase flow through a porous medium.
*/
template<uint32_t dim>
void DES2PhaseSlightlyCompressibleFlow_Example::RunSimulation( Model<dim>& model )
{

    // -------------------------------------
    // 0.0 Set variables used throughout the simulation
    // -------------------------------------
    double model_time =  0.0; // time

    // -------------------------------------
    // 1.0 ive the model dimensions
    // -------------------------------------
    printModelDimensions( model, true );

    // ---------------------------------------------------
    // 2.0 Simulation settings
    // ---------------------------------------------------
    Standard_IO_Handler  stdio;
    bool  DES = stdio.YesNo("Do you want to solve the transport equation with DES? (y=DES, n=TDS)"); 
    double Courant_multiplier, PEP_parameter;
    cout <<"\nEnter CFL multiplier (suggested value: 0.3) and PEP parameter (suggested value: 0.1)" << endl;
    cin >> Courant_multiplier >> PEP_parameter;  
     
    bool  with_gravity_forces = stdio.YesNo("Do you want to include gravity effect (y/n)?"); 
    bool  with_capillary_spreading = stdio.YesNo("Do you want to include capillary effect (y/n)?");


    // --------------------------------------------
    // 3.0 Configure the simulation from a file
    // --------------------------------------------
    InputDataManager<dim>  model_configuration;
    string config_file = model.Name();
    config_file += "(DES_2phase)";
    model_configuration.ConfigureFromFile( model,
                                           config_file.c_str(),
                                           false,
                                           true,   // 2) default prop.values
                                           true,   // 3) regional prop.values
                                           true,   // 4) essential box-boundary conditions
                                           false,  // 5) essential flags
                                           false   // 6) boundary conditions
                                          );
    // check whether tensor k is used
    bool with_tensor_permeability(false);
    if(model.Database().IsDefined("tensor permeability")) {
      TensorVariable<dim> kk;
      static Index  key_kk(model.Database().StorageKey("tensor permeability"));
      model.Region("Model").E(0U)->Read(key_kk, kk);
      if(!isnan(kk(0U,0U))) with_tensor_permeability = true;
    }
    if(with_tensor_permeability) cout<<"\ntensor permeability is in use"<<endl;
    else  cout<<"\nscalar permeability is in use"<<endl;

    // ---------------------------------------------------------------------
    // 4.0 Use the flow functions to compute the relative permeabilities
    // ---------------------------------------------------------------------
    FlowFunctionsModule1<dim> flowfunctions(model.Database(), model.Read( model.Database().StorageKey("acceleration gravity") ));
    Compute2PhaseFlowProperties(model, flowfunctions, with_gravity_forces, with_tensor_permeability);

    // ------------------------------------------------------------------------------------------
    // 5.0 Setting up an FE algorithm to solve the steady state pressure
    // ------------------------------------------------------------------------------------------
    ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);
    printRangeOfVariable( model, stdio, "fluid pressure" );

    // -----------------------------
    // 6.0 Output initial conditions
    // -----------------------------
    VTU_Interface<dim>  vtu(model);  // binary VTK output, creates much smaller files than VTK
    vtu.OutputDataToVTU( "fluid_pressure", "fluid pressure",    "Model", 0 );
    vtu.OutputDataToVTU( "saturation carbonic phase", "saturation carbonic phase",    "Model", 0 );

    // -------------------------------------------------------------
    // 7.0 Construct the finite volume grid and DES transport algorithms
    // -------------------------------------------------------------
    DES2PhaseTransport<dim,FlowFunctionsModule1>* DEStransport{nullptr};
    DEStransport = new DES2PhaseSlightlyCompressibleTransport<dim,FlowFunctionsModule1> (model, "Model",
                                                                                         with_gravity_forces,
                                                                                         with_capillary_spreading,
                                                                                         Courant_multiplier,
                                                                                         PEP_parameter,
                                                                                         1., //relaxing factor
                                                                                         false, //tensor k
                                                                                         false, //2nd order in space
                                                                                         flowfunctions);

    // -----------------------
    // 8.0 Time Loop Variables
    // -----------------------
    // define some constant variables
    const double    day(86400.0);
    double          max_time (60.0*day);     // run for 60 days
    double          time_increment(0.5*day);      // timestep 0.5 day
    long            save_frequency(6);       // write results to file every 3 days
    size_t	        time, save_counter(1);
    size_t	        n_threads(1);


    // -----------------------
    // 9.0 Transient loop
    // -----------------------
    double solving_time = 0.;
    clock_t  T_begin;    
    while ( model_time < max_time )
      {
         // update pressure field
         Compute2PhaseFlowProperties( model, flowfunctions, with_gravity_forces, with_tensor_permeability );
         ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);

         // compute advection of phases
         T_begin = clock();
         if (DES) DEStransport->AdvectVariable_DES( time_increment, model_time+time_increment, n_threads );
         else DEStransport->AdvectVariable_TDS( time_increment, n_threads );
         solving_time += clock() - T_begin;
         
         // increment time
         model_time += time_increment;

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
                  vtu.OutputDataToVTU( "DES_saturation_carbonic_phase", "saturation carbonic phase",    "Model", time );
                  vtu.OutputDataToVTU( "DES_fluid_velocity", "total velocity",    "Model", time );
                  vtu.OutputDataToVTU( "DES_Update_count", "update count", "Model",  time );
                  vtu.OutputDataToVTU( "DES_CFL_multiplier", "cfl multiplier", "Model",  time );
                  vtu.OutputDataToVTU( "DES_sw_shock", "shock saturation aqueous phase", "Model",  time );                 
              } else {
                  vtu.OutputDataToVTU( "TDS_fluid_pressure", "fluid pressure",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_saturation_carbonic_phase", "saturation carbonic phase",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_fluid_velocity", "total velocity",    "Model", time );
                  vtu.OutputDataToVTU( "TDS_Update_count", "update count", "Model",  time );  
              }   
              save_counter = 0;
         }
         save_counter++;

         // runtime info
         cout <<"\n\nmain: RUNTIME (DAYS): "<< model_time/day << endl << endl;

      }

    // clocking the runtime
    if(DES) cerr << "\nmain: DES transport uses " << static_cast<double>(solving_time/CLOCKS_PER_SEC) << " seconds " << endl;
    else    cerr << "\nmain: TDS transport uses " << static_cast<double>(solving_time/CLOCKS_PER_SEC) << " seconds " << endl;
    
    // terminate
    delete DEStransport;
    cerr << "\nmain: That's it..."<< endl;

} // Run()




template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
void DES2PhaseSlightlyCompressibleFlow_Example::Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k )
{
    const size_t v( (dim==1u) ? 0u : 1u );

    //static const Region<dim>& mref = mdl.Region("Model");
    auto mref = mdl.Region("Model");
    // keys to properties
    static Index  sw_key(mdl.Database().StorageKey("saturation aqueous phase"));
    static Index  snw_key(mdl.Database().StorageKey("saturation carbonic phase"));
    static Index  thi_key(mdl.Database().StorageKey("thickness"));
    static Index  rhow_key(mdl.Database().StorageKey("density aqueous phase"));
    static Index  rhon_key(mdl.Database().StorageKey("density carbonic phase"));
    static Index  gt_key(mdl.Database().StorageKey("gravity term"));
    static Index  dip_key(mdl.Database().StorageKey("dip vector"));
    static Index  cw_key(mdl.Database().StorageKey("compressibility aqueous phase"));
    static Index  cn_key(mdl.Database().StorageKey("compressibility carbonic phase"));
    static Index  cR_key(mdl.Database().StorageKey("compressibility rock"));
    static Index  ct_key(mdl.Database().StorageKey("total system compressibility"));
    static Index  phi_key(mdl.Database().StorageKey("porosity"));
    static Index  swr_key(mdl.Database().StorageKey("residual saturation aqueous phase"));
    static Index  snr_key(mdl.Database().StorageKey("residual saturation carbonic phase"));
    static Index  muw_key(mdl.Database().StorageKey("viscosity aqueous phase"));
    static Index  mun_key(mdl.Database().StorageKey("viscosity carbonic phase"));

    const double gravity_acceleration = mdl.Read( mdl.Database().StorageKey("acceleration gravity") );

    // 1. Computing the saturation of water = 1 - So
    for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
    {
      if( (*nit)->Status(sw_key) != DIRICH ) {
        double sw = 1. - (*nit)->Read(snw_key);
        (*nit)->Store( sw_key, makeScalar((*nit)->Status(sw_key),sw));
      }
    }

    // 2. Computing the multiphase flow properties
    if(!with_tensor_k) { //scalar k
      static Index  k_key(mdl.Database().StorageKey("permeability"));
      static Index  mobt_key(mdl.Database().StorageKey("total mobility permeability product"));
      //vector<Element<dim>* >::const_iterator eit;
      for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ ) {
        Element<dim>* eptr = *eit;
        const size_t nodes(eptr->Nodes());
        vector<double> IPOL;
        eptr->N_AtBaryCenter( IPOL );
        double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.), e_cw(0.), e_cn(0.);
        for ( auto i{0U}; i<nodes; ++i ) {
          e_sw += IPOL[i] * eptr->N(i)->Read( sw_key );
          e_muw += IPOL[i] * eptr->N(i)->Read( muw_key );
          e_mun += IPOL[i] * eptr->N(i)->Read( mun_key );
          e_rhow += IPOL[i] * eptr->N(i)->Read( rhow_key );
          e_rhon += IPOL[i] * eptr->N(i)->Read( rhon_key );
          e_cw += IPOL[i] * eptr->N(i)->Read( cw_key );
          e_cn += IPOL[i] * eptr->N(i)->Read( cn_key );
          ipol_sum += IPOL[i];
        }
        if(ipol_sum > 0.) {
          e_sw *= 1. / ipol_sum;
          e_muw *= 1. / ipol_sum;
          e_mun *= 1. / ipol_sum;
          e_rhow *= 1. / ipol_sum;
          e_rhon *= 1. / ipol_sum;
          e_cw *= 1. / ipol_sum;
          e_cn *= 1. / ipol_sum;
        }
        e_sn = 1.0 - e_sw;

        const double e_swr = eptr->Read(swr_key);
        const double e_snr = eptr->Read(snr_key);
        double thickness = eptr->Read(thi_key);
        if(isnan(thickness)) thickness = 1.;
        double e_k = eptr->Read(k_key);

        //compute element mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = e_k * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = e_k * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun; //krn/mun
        e_lt = e_lw + e_ln;
        eptr->Store( mobt_key, makeScalar(PLAIN, e_lt) );

        //gravity term
        if ( with_gravity ) {
          VectorVariable<dim> gravity;
          eptr->Read( dip_key, gravity );
          if(isnan(gravity(1))) { //dip vector has not been initialised
            gravity(0u) = 0.; gravity(1u) = -1.;
            if(dim==3u) gravity(2u) = 0.;
            eptr->Store( dip_key, gravity );
          }
          double gravity_w(0.), gravity_n(0.);
          if(e_sw > e_swr) gravity_w = gravity_acceleration * e_k * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw * (e_rhow);
          if(e_sn > e_snr) gravity_n = gravity_acceleration * e_k * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun * (e_rhon);
          gravity_w *= fabs(gravity[v]);
          gravity_n *= fabs(gravity[v]);
          gravity *= (gravity_w + gravity_n);
          (*eit)->Store( gt_key, gravity );
        }

        /*
        // total system compressibility (weighted average approach)
        const double e_cR = eptr->Read( cR_key );
        double e_phi = eptr->Read(phi_key); //porosity
        double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        e_ct *= thickness;
        eptr->Store( ct_key, makeScalar(PLAIN, e_ct)  );
        */
      }

    } else { //tensor k
      static Index  kk_key(mdl.Database().StorageKey("tensor permeability"));
      static Index  LT_key(mdl.Database().StorageKey("tensor total mobility permeability product"));
      TensorVariable<dim> kk;
      VectorVariable<dim> kV;
      //vector<Element<dim>* >::const_iterator eit;
      for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ ) {
        Element<dim>* eptr = *eit;
        const size_t nodes(eptr->Nodes());
        vector<double> IPOL;
        eptr->N_AtBaryCenter( IPOL );
        double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.), e_cw(0.), e_cn(0.);
        for ( auto i{0U}; i<nodes; ++i ) {
          e_sw += IPOL[i] * eptr->N(i)->Read( sw_key );
          e_muw += IPOL[i] * eptr->N(i)->Read( muw_key );
          e_mun += IPOL[i] * eptr->N(i)->Read( mun_key );
          e_rhow += IPOL[i] * eptr->N(i)->Read( rhow_key );
          e_rhon += IPOL[i] * eptr->N(i)->Read( rhon_key );
          e_cw += IPOL[i] * eptr->N(i)->Read( cw_key );
          e_cn += IPOL[i] * eptr->N(i)->Read( cn_key );
          ipol_sum += IPOL[i];
        }
        if(ipol_sum > 0.) {
          e_sw *= 1. / ipol_sum;
          e_muw *= 1. / ipol_sum;
          e_mun *= 1. / ipol_sum;
          e_rhow *= 1. / ipol_sum;
          e_rhon *= 1. / ipol_sum;
          e_cw *= 1. / ipol_sum;
          e_cn *= 1. / ipol_sum;
        }
        e_sn = 1.0 - e_sw;

        const double e_swr = eptr->Read(swr_key);
        const double e_snr = eptr->Read(snr_key);
        double thickness = eptr->Read(thi_key); //thickness
        if(isnan(thickness)) thickness = 1.;

        //compute element mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun; //krn/mun
        e_lt = e_lw + e_ln;
        eptr->Read( kk_key, kk );
        kk *= e_lt;
        eptr->Store( LT_key, kk );

        //gravity term
        if ( with_gravity ) {
          VectorVariable<dim> gravity;
          eptr->Read( dip_key, gravity );
          if(isnan(gravity(1))) { //dip vector has not been initialised
            gravity(0u) = 0.; gravity(1u) = -1.;
            if(dim==3u) gravity(2u) = 0.;
            eptr->Store( dip_key, gravity );
          }
          // projection of tensor on the dip vector
          eptr->Read( kk_key, kk );
          kV = kk * gravity;
          double kV_magnitude = kV.Length();

          double gravity_w(0.), gravity_n(0.);
          if(e_sw > e_swr) gravity_w = gravity_acceleration * kV_magnitude * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw * (e_rhow);
          if(e_sn > e_snr) gravity_n = gravity_acceleration * kV_magnitude * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun * (e_rhon);
          gravity_w *= fabs(gravity[v]);
          gravity_n *= fabs(gravity[v]);
          gravity *= (gravity_w + gravity_n);
          (*eit)->Store( gt_key, gravity );
        }

        /*
        // total system compressibility (weighted average approach)
        const double e_cR = eptr->Read( cR_key );
        double e_phi = eptr->Read(phi_key); //porosity
        double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        e_ct *= thickness;
        eptr->Store( ct_key, makeScalar(PLAIN, e_ct)  );
        */
      }
    }
}

template void DES2PhaseSlightlyCompressibleFlow_Example::Compute2PhaseFlowProperties( Model<2U>&, FlowFunctionsModule1<2U>&, bool, bool );
template void DES2PhaseSlightlyCompressibleFlow_Example::Compute2PhaseFlowProperties( Model<3U>&, FlowFunctionsModule1<3U>&, bool, bool );




template<uint32_t dim>
void DES2PhaseSlightlyCompressibleFlow_Example::ComputeSteadyStatePressure( Model<dim>& mdl, bool with_gravity, bool with_tensor_k )
{
    bool verbose(false);

    if(verbose) {
      cout << "\nDES2PhaseSlightlyCompressibleFlow_Example::ComputeSteadyStatePressure: Input parameters: " << endl;
      printRangeOfVariable( mdl, "fluid pressure" );
      printRangeOfVariable( mdl, "total mobility permeability product" );
      printRangeOfVariable( mdl, "permeability" );
      printRangeOfVariable( mdl, "fluid volume source");
      printRangeOfVariable( mdl, "nodal fluid volume source" );
      if(with_gravity) printRangeOfVariable( mdl, "gravity term" );
    }

    std::string	conductance_operator;
    if(!with_tensor_k) conductance_operator = "total mobility permeability product";
    else conductance_operator = "tensor total mobility permeability product";

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                 samg_solver( &settings );
    PDE_Integrator<dim,Element>  steady_pressure(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
    PDE_Integrator<dim,Element>  steady_pressure(linear_solver);
#endif
    NumIntegral_dNT_lhsop_dN_dV<dim>  conductance( mdl.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<dim> elmt_volume_source( mdl.Database(), "fluid volume source", "fluid pressure" );
    //PointSource_rhsop<dim>          nodal_volume_source( mdl.Database(), "nodal fluid volume source", "fluid pressure" );
    steady_pressure.Add( &conductance );
    steady_pressure.Add( &elmt_volume_source );
    //steady_pressure.AddBoundaryIntegral( &influx ); //added
    steady_pressure.Verbose(verbose);

    NumIntegral_dNT_op_dV<dim>* gravity(nullptr);
    if ( with_gravity ) {
      gravity = new NumIntegral_dNT_op_dV<dim>( mdl.Database(), "gravity term", "fluid pressure" );
      steady_pressure.Add( gravity );
    }

/*
    // iout
    //settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
*/
    if(verbose) {
      cout << "\n\n\nDES2PhaseSlightlyCompressibleFlow_Example::ComputeSteadyStatePressure: ";
      cout << " Computing '" << "steady state fluid pressure" << "'" << endl;
    }

    mdl.Apply( steady_pressure );

    delete gravity;

} // end SteadyStatePressure

template void DES2PhaseSlightlyCompressibleFlow_Example::ComputeSteadyStatePressure(Model<2U>&, bool, bool);
template void DES2PhaseSlightlyCompressibleFlow_Example::ComputeSteadyStatePressure(Model<3U>&, bool, bool);


} // csmp
