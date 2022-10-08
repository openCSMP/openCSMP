#include "DES2PhaseFlowWithSplitBoundary_Example.h"

// the CSMP model
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"

// FV algorithms
#include "DES2PhaseSlightlyCompressibleTransport.h"
#include "SandPropertiesFor_VE_Model.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#endif

// monitoring individual regions
#include "RegionMonitor.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"

#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "FlowFunctionsModule.h"
#include "NumIntegral_NT_op_N_dS.h"

#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp{


void DES2PhaseFlowWithSplitBoundary_Example::Specifications()
{
  SetTitle( "Two phase slightly compressible flow simulation with splitboundaries, using DES (discrete event simulation)");
  SetDifficulty( 3 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "Qi Shao" );
  AddDescription( "Two phase slightly compressible flow simulation with splitboundaries, using discrete event simulation (DES) or time-driven simulation (TDS).");
  AddDescription( "Output is written to to VTK files");
  AddRequirement( "source files: 'DES2PhaseFlowWithSplitBoundary_Example.cpp' and '*.h'" );
  AddRequirement( "box2d_fault (CSMP native binary files), box2d_fault-configuration.txt, DES_2phase_variables.txt");
} 


/** 
    Two phase slightly compressible flow simulation with splitboundaries via CSMP's DES transport method
    combining finite elements (for pressure) with finite volumes (for advection of non-wetting phase)

    Use models 'box2d_fault' (CSMP native binary files, -configuration.txt) as input file suites.
*/

void DES2PhaseFlowWithSplitBoundary_Example::Run()
  {
    // reading in a csmp native format model, first determining whether it will be 2 or 3 dimensional.
    uint32_t dimension;
    cout << "\nPlease enter the dimension of the model (2 for 2D (default), 3 for 3D):" << endl;
    cin >> dimension;
    if (dimension != 2U and dimension != 3U)
      throw csmp::Exception(ERROR, "input dimension of model", "must be 2 or 3");

    string model_name;
    cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'box2d_fault':"<<endl;
    cin.ignore();
    getline(cin, model_name);
    if (model_name.length() == 0) model_name = "box2d_fault";

    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string variable_file = "DES_2phase_variables.txt";
    string config_file = model_name + "(DES_2phase_splitbdy)";
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);
    //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
    if (dimension == 2U) {
      Model<2U>  model(model_name, variable_file);
      RunSimulation(model);
    } else if (dimension == 3U) {
      Model<3U>  model(model_name, variable_file);
      RunSimulation(model);
    }

    fs::current_path("../../example_inputs/");
  } // end run
  
  

/**
    Performs either discrete-event simulation (DES) or time-driven simulation (TDS) of 2-phase flow through a porous medium.
*/
template<uint32_t dim>
void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation( Model<dim>& model )
 {
    // 1. MODEL CONFIGURATION
    // ----------------------
    // simulation settings
    Standard_IO_Handler stdio;
    VTU_Interface<dim>  vtu(model);

    bool  DES = stdio.YesNo("Do you want to solve the transport equation with DES? (y=DES, n=TDS)"); 
    double Courant_multiplier, PEP_parameter;
    cout <<"\nEnter CFL multiplier (suggested value: 0.3 for TDS, 0.1 for DES) and PEP parameter (suggested value: 0.1)" << endl;
    cin >> Courant_multiplier >> PEP_parameter;  
     
    const bool  with_gravity_forces = stdio.YesNo("Do you want to include gravity effect (y/n)?");
    const bool  with_capillary_spreading = stdio.YesNo("Do you want to include capillary effect (y/n)?");

    const bool  with_split_boundaries = stdio.YesNo("Do you want to create split boundaries (y/n)?");

    // give the model dimensions
    printModelDimensions( model, true );

    // Configure the simulation from a file
    InputDataManager<dim>  model_configuration;
    string config_file = model.Name();
    config_file += "(DES_2phase_splitbdy)";
    model_configuration.ConfigureFromFile( model, config_file.c_str(),
                                           false,  ///< region name from parameter range
                                           true,   ///< default property values
                                           true,   ///< regional property values
                                           false,  ///< boundary conditions for box-shaped model
                                           true,   ///< essential conditions for regions
                                           true );  ///< boundary conditions for arbitrary-shaped model

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

    /*
    // 2. PROVISIONS FOR LOWER-DIMENSIONAL REPRESENTATION OF SAND HORIZONS
    // -------------------------------------------------------------------
      {
         computeGravityDipVectors( model );
         if ( model.ContainsRegion("SAND") ) {
              const string dim_minus1_region("SAND");
              // not needed in current approach: 27/8/22
              //computeFV_Diameter_Normal_VerticalExtent( model, dim_minus1_region );
              //computeSpillPointSaturation( model, dim_minus1_region );
              //lowerDimensionalLayerDiagnostics( model, dim_minus1_region );
           }
         else cout <<"\nRunSimulation: no provisions made for lower-dimensional sandbodies."<< endl;
      }
    const double bcp{2.}, swr{0.15}, snr{0.};
    SandPropertiesFor_VE_Model  sandProperties( model, "SAND", bcp, swr, snr );
    */

    // 3. RELATIVE PERMEABILITY & CAPILLARY PRESSURE MODEL
    // ---------------------------------------------------
    // flow functions (Brooks Corey)
    FlowFunctionsModule1<dim> flowfunctions(model.Database(), model.Read( model.Database().StorageKey("acceleration gravity") ));

    // solve static pressure before split boundaries are created
    Compute2PhaseFlowProperties(model, flowfunctions, with_gravity_forces, with_tensor_permeability);
    /*
    if constexpr (dim == 2U )
      if ( sandProperties.HasLineElementRepresentation() )
         sandProperties.Compute2PhaseFlowPropertiesForSandLayer(model);
    */
    ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);
    printRangeOfVariable( model, stdio, "fluid pressure" );
    vtu.OutputDataToVTU( "steady_state_pressure", "fluid pressure", "Model", 0 );

    // create split boundaries between regions
    if( with_split_boundaries ) {
      //create splitboundaries for all interfaces
      size_t n_split_boundaries = model.SeparateUniqueRegionsBySplitBoundaries();
      cout<<"number of splitboundafies = "<<n_split_boundaries<<endl;
      model.SplitBoundariesOut();

      // using 'nodal variable' to visualise which nodes are manifolds
      Region<dim>  model_domain = model.Region("Model");
      if(!model.Database().IsDefined("number of collocated nodes"))
        model.CreateProperty( "number of collocated nodes", "n_nd", "uint", SCALAR, NODE, 1, 0, 100);;

      model.Region("Model").InputPropertyValue( "number of collocated nodes", makeScalar(PLAIN,1), COMPLETE);
      const csmp::Index nodes_key = model.Database().StorageKey("number of collocated nodes");
      for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
        if ( (*nit)->IsManifold() )
          // there should be as many branches as materials come together
          (*nit)->Store( nodes_key, makeScalar(PLAIN,(*nit)->Manifold()->Branches()));

      vtu.OmitZeroInFileName( true );
      vtu.OutputDataToVTU( "Collocated_nodes", "number of collocated nodes", "Model", static_cast<int>(0) );

    } //end create split boundaries

    //sort manifold nodes by entry pressure
    if ( with_split_boundaries ) {
      const csmp::Index pd_key = model.Database().StorageKey("entry pressure");
      for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++)
        (*mit).SortByVariableValue(pd_key);

      //check parent elements of duplicated nodes
      if(!model.Database().IsDefined("parent element id")) {
        model.CreateProperty( "parent element id", "eID", "uint", SCALAR, ELEMENT, 1, 0 ,1.00E+08);
        model.Region("Model").InputPropertyValue( "parent element id", makeScalar(PLAIN,0), COMPLETE);
      }
      const csmp::INDEX<SCALAR,ELEMENT> key_parent = csmp::INDEX<SCALAR,ELEMENT>( model.Database().StorageKey("parent element id") );
      size_t parent_id(1);
      for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++) {
        parent_id = 1;
        auto md = (*mit);
        auto master_node = md.N(0);
        for( auto e{0U}; e < master_node->Parents(); e++)
          master_node->Parent(e)->Store(key_parent, makeScalar(PLAIN, parent_id));

        for(size_t n=1;n<md.Branches();n++) {
          auto slave_node = md.N(n);
          parent_id++;
          for( auto e{0U}; e < slave_node->Parents(); e++)
            slave_node->Parent(e)->Store(key_parent, makeScalar(PLAIN, parent_id));
        }
      }
      vtu.OutputDataToVTU( "parent_elment_id", "parent element id",  "Model", 0 );
    }


    // construct DES transport
    auto DEStransport = new DES2PhaseSlightlyCompressibleTransport<dim,FlowFunctionsModule1> (model,
                                                                                             "Model",
                                                                                             with_gravity_forces,
                                                                                             with_capillary_spreading,
                                                                                             Courant_multiplier,
                                                                                             PEP_parameter,
                                                                                             1., //relaxation factor
                                                                                             with_tensor_permeability, //tensor k
                                                                                             false, //2nd order in space
                                                                                             flowfunctions);

    // solve static pressure after split boundaries are created
    if( with_split_boundaries ) {
      ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);
      printRangeOfVariable(model, stdio, "fluid pressure");
      vtu.OutputDataToVTU("steady_state_pressure_with_splitboundaries", "fluid pressure", "Model", 0);
      //vtu.OutputDataToVTU( "gravity_term", "gravity term", "Model", 0 );
    }

    // defining input properties
    std::list<string> input_properties;
    input_properties.emplace_back( "porosity" );
    input_properties.emplace_back( "fluid pressure" );
    input_properties.emplace_back( "permeability" );
    input_properties.emplace_back( "vertical permeability" );
    input_properties.emplace_back( "rocktype" );
    input_properties.emplace_back( "fluid volume source");
    input_properties.emplace_back( "nodal fluid volume source" );
    input_properties.emplace_back( "saturation carbonic phase" );
    input_properties.emplace_back( "saturation aqueous phase" );
    input_properties.emplace_back( "thickness" );
    input_properties.emplace_back( "entry pressure" );
    if ( with_split_boundaries ) {
      input_properties.emplace_back("pressure continuity status");
      input_properties.emplace_back("breakthrough status");
      /*
      // for model with lower dimensional representation of the sand horizons inside the split boundaries
      input_properties.emplace_back("finite volume diameter");
      input_properties.emplace_back("finite volume normal");
      input_properties.emplace_back("finite volume vertical extent");
      input_properties.emplace_back("finite volume diagnostics");
       */
    }

    //defining output properties
    std::list<string> output_properties;
    output_properties.emplace_back( "rocktype" );
    output_properties.emplace_back( "entry pressure" );
    output_properties.emplace_back( "fluid pressure" );
    output_properties.emplace_back( "saturation aqueous phase" );
    output_properties.emplace_back( "saturation carbonic phase" );
    output_properties.emplace_back( "update count" );
    if(!model.Database().IsDefined("out range value count"))
      model.CreateProperty( "out range value count", "orvc", "uint", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    model.Region("Model").InputPropertyValue( "out range value count", makeScalar(PLAIN,0), COMPLETE);
    output_properties.emplace_back( "out range value count" );
    if(with_split_boundaries) {
      output_properties.emplace_back("pressure continuity status");
      output_properties.emplace_back("breakthrough status");
    }

    //write input properties
    std::string model_name ( model.Name() );
    std::string file_name = model_name + "-input_data";
    vtu.OutputDataToVTU( file_name, input_properties, "Model", 0 );

    // global time for simulated runtime
    double	model_time(0.);

    // setting up region monitor
    RegionMonitor<dim>  monitor( model, "saturation carbonic phase", "saturation carbonic phase",  true );
    monitor.ScalarPropertyRanges( model, model_time );
    monitor.ScalarPropertyIntegrals( model, model_time );
    monitor.Out( (model_name + "-monitored_CO2_saturation").c_str() );

    RegionMonitor<dim>  monitor2( model, "fluid pressure", "fluid pressure",  true );
    monitor2.ScalarPropertyRanges( model, model_time );
    monitor2.ScalarPropertyIntegrals( model, model_time );
    monitor2.Out( (model_name + "-monitored_pressure").c_str() );

    // setting up time parameters (HARDWIRED PRESSURE STEPS!)
    const double day{86400.}, year{ 86400. * 365. };
    double max_time (5. * year);      // run for # years
    double time_increment(3600.); // timestep 
    double save_interval = 5. * day;  // save every # days

    size_t time;
    double end_time = model_time + max_time;
    double next_save_time = model_time + save_interval;

    // number of threads for use in transport
    size_t n_threads=1;

    // starting simulation
    if (DES) cerr <<"\n\nmain: Starting DES simulation: "<<endl;
    else cerr <<"\n\nmain: Starting TDS simulation: "<<endl;

    long output_count(1);
    while ( model_time < end_time )
    {
      Compute2PhaseFlowProperties( model, flowfunctions, with_gravity_forces, with_tensor_permeability );
      /*
      // setting multiphase flow properties for dim-1 elements representing sand
      if constexpr ( dim == 2U )
        if ( sandProperties.HasLineElementRepresentation() )
          sandProperties.Compute2PhaseFlowPropertiesForSandLayer( model );
      */
      // solve pressure
      ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);

      //transport steps
      cout << std::fixed << std::setprecision(3) <<"\nTransporting from "<<model_time<<" sec to "<<model_time+time_increment<<" sec, dt = "<<time_increment<<" sec\n"<<endl;
      if (DES) DEStransport->AdvectVariable_DES( time_increment, model_time+time_increment, n_threads );
      else DEStransport->AdvectVariable_TDS( time_increment, n_threads );

      //increment time
      model_time += time_increment;

      //output variables
      if ( model_time >= next_save_time ) {
        if (with_split_boundaries) {
          const csmp::Index p_key = model.Database().StorageKey("fluid pressure");
          const csmp::Index p_status_key = model.Database().StorageKey("pressure continuity status");
          for (auto mit = model.Mesh().NodeManifoldsBegin(); mit != model.Mesh().NodeManifoldsEnd(); mit++) {
            const auto n_branches{(*mit).Branches()};
            for (auto n{1U}; n < n_branches; n++) {
              auto node = (*mit).N(n);
              if (node->Status(p_key) == ROBIN) node->Store(p_status_key, makeScalar(PLAIN, 1));
              else node->Store(p_status_key, makeScalar(PLAIN, 0));
            }
          }
        }

        time = static_cast<long>(model_time/day);
        if (DES) {
          std::string runtime_file_name(model_name + "-DES_runtime_output");
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "Model", time);
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "SAND", time);
        } else {
          std::string runtime_file_name(model_name + "-TDS_runtime_output");
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "Model", output_count);
        }

        monitor.ScalarPropertyRanges(model, model_time);
        monitor.ScalarPropertyIntegrals(model, model_time);
        monitor.Out((model_name + "-monitored_CO2_saturation").c_str());
        monitor2.ScalarPropertyRanges(model, model_time);
        monitor2.ScalarPropertyIntegrals(model, model_time);
        monitor2.Out((model_name + "-monitored_pressure").c_str());

        next_save_time += save_interval;
        output_count++;

        model.Region("Model").InputPropertyValue("update count", makeScalar(PLAIN, 0), COMPLETE);

      }

      // runtime info
      cout <<"\n\nmain: RUNTIME (SEC): "<< model_time << endl << endl;

    }

    delete DEStransport;

    // terminate
    cerr << "\nmain: That's it..."<< endl;

  } // RunSimulation()

  template void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation(Model<2U>& model);
  template void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation(Model<3U>& model);





  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void DES2PhaseFlowWithSplitBoundary_Example::Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k )
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

template void DES2PhaseFlowWithSplitBoundary_Example::Compute2PhaseFlowProperties( Model<2U>&, FlowFunctionsModule1<2U>&, bool, bool );
template void DES2PhaseFlowWithSplitBoundary_Example::Compute2PhaseFlowProperties( Model<3U>&, FlowFunctionsModule1<3U>&, bool, bool );






  template<uint32_t dim>
  void DES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure( Model<dim>& mdl, bool with_gravity, bool with_tensor_k )
  {
    bool verbose(false);

    if(verbose) {
      cout << "\nDES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure: Input parameters: " << endl;
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

#ifdef USE_SAMG_SOLVER
    SAMG_Settings settings;
    // iout
    settings.ExplicitSecondary(true);
    settings.Set_iout1( -1 );
    settings.Set_iout2( -1 );
    settings.Set_idmp( -1 );
    settings.Set_mode_mess( -2 );
    SAMG_Solver                 samg_solver( &settings );
    PDE_Integrator<dim,Region>  steady_pressure(samg_solver);
#else
    EigenSolver linear_solver;
    PDE_Integrator<dim,Region>  steady_pressure(linear_solver);
#endif
    NumIntegral_dNT_op_dN_dV<dim>  conductance( mdl.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim>    elmt_volume_source( mdl.Database(), "fluid volume source", "fluid pressure" );
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
      cout << "\n\n\nDES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure: ";
      cout << " Computing '" << "steady state fluid pressure" << "'" << endl;
    }

    mdl.Apply( steady_pressure );

    delete gravity;

  } // end SteadyStatePressure

  template void DES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure(Model<2U>&, bool, bool);
  template void DES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure(Model<3U>&, bool, bool);




/**
      loops over the finite volumes of the supplied lower-dimensional supplied region, computing the diameters of the finite volumes and their normals
      @todo the "dip vector" variables and thickness attributes need to be computed elsewhere
*/
template<uint32_t dim>
void computeFV_Diameter_Normal_VerticalExtent( Model<dim>& model, const std::string& region_name )
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion( region_name ) )
      csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                      "target region does not exist");
                      
    Region<dim>& subdomain = model.Region( region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 )
      csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                      "target region is not lower dimensional (dim region != dim-1)");
                      
    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    const csmp::Index diam_key = model.Database().StorageKey("finite volume diameter");
    const csmp::Index nrml_key = model.Database().StorageKey("finite volume normal");
    const csmp::Index vext_key = model.Database().StorageKey("finite volume vertical extent");
    
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         // FV diameter/vertical extent
         const auto FV_props = diameterAndVerticalExtentOfLowerDimensional_FV( (*nit) );
         assert( FV_props.first > 0. ); // diameter must be greater than zero
         (*nit)->Store( diam_key, makeScalar(PLAIN,FV_props.first) );
         (*nit)->Store( vext_key, makeScalar(PLAIN,FV_props.second) );
         // FV normal
         Point<dim> nrml;
         bool was_able_to_compute_normal = (*nit)->UnitNormal( nrml );
         assert( was_able_to_compute_normal );
         // if this normal is not upward pointing, it is flipped
         if constexpr ( dim == 3U ) if ( dotProduct( nrml, Point<3U>(0.,1.,0.) ) < 0. ) nrml *= -1.;
         if constexpr ( dim == 2U ) if ( dotProduct( nrml, Point<2U>(0.,1.) ) < 0. )    nrml *= -1.;
         (*nit)->Store( nrml_key, move( VectorVariable<dim>(nrml) ) );
      }
 
 } // end computeFV_Diameter_Normal_VerticalExtent


template void computeFV_Diameter_Normal_VerticalExtent( Model<3U>&, const string& );
template void computeFV_Diameter_Normal_VerticalExtent( Model<2U>&, const string& );




/**
    compute gravity dip vectors borrowed from ACGS
*/
template<uint32_t dim>
void computeGravityDipVectors( Model<dim>& model )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. volumetric model region
  // --------------------------
  csmp::Region<dim>&   sgref = model.Region( "Model" );
  VectorVariable<dim>  gvt;
  Point<dim>           gravity_direction;
  if constexpr ( dim == 3U ) gravity_direction = { 0., -1., 0. };
  else                       gravity_direction = { 0., -1. };

  const csmp::Index key_dip  = model.Database().StorageKey("dip vector");
  const csmp::Index key_dipf = model.Database().StorageKey("face dip vector");

  for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); it++ )
  {
    // computing gravity direction vectors for lower dimensional elements
    // line elements
    if ( (*it)->IsLine() ) {
      Point<dim> line_vector( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
      line_vector.NormalizeLengthTo( 1. );
      if ( line_vector[1] > 0. ) line_vector *= -1.;
      gvt( 0 ) = line_vector[0];
      gvt( 1 ) = line_vector[1];
      if constexpr (dim == 3U ) gvt( 2 ) = line_vector[2];
    }
    // surface elements
    else if ( (*it)->IsSurface() ) {
      Point<dim> line_vector1( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
      Point<dim> line_vector2( (*it)->N( 1 )->Coordinate() - (*it)->N( 2 )->Coordinate() );

      // finding the normal to the surface element
      Point<dim> surface_normal( crossProduct( line_vector1, line_vector2 ) );

      // finding the gravity direction from the cross-product of the horizontal surface vector and the gravity direction
      Point<dim> gravity_projection_on_surface( crossProduct( surface_normal, crossProduct( surface_normal, gravity_direction ) ) );
      gravity_projection_on_surface.NormalizeLengthTo( 1. );

      if ( gravity_projection_on_surface[1] > 0. ) gravity_projection_on_surface *= -1.;

      gvt( 0 ) = gravity_projection_on_surface[0];
      gvt( 1 ) = gravity_projection_on_surface[1];
      if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_surface[2];
    }
    // volume elements
    else {
      gvt( 0 ) = 0.;
      gvt( 1 ) = -1.;
      if constexpr (dim == 3U ) gvt( 2 ) = 0.;
    }
    if(isnan(gvt( 0 )) || isnan(gvt( 1 )) || isnan(gvt( 2 ))) {
        gvt( 0 ) = 0.;
        gvt( 1 ) = -1.;
        if constexpr (dim == 3U ) gvt( 2 ) = 0.;
        cerr<<"gvt reset to [0, -1, 0]"<<endl;
    }
    (*it)->Store( key_dip, gvt );
  }

  // 2. for all model boundaries
  // ---------------------------
  for ( auto git = model.BoundariesBegin(); git != model.BoundariesEnd(); git++ )
  {
    //          cout <<"\n\tBoundary: "<< (*git).first.first << (*git).first.second;
    //          cout.flush();
    assert( (*git).second.Cells() > 0 );
    for ( auto it = (*git).second.CellsBegin(); it != (*git).second.CellsEnd(); it++ )
    {
      if ( (*it) != nullptr and (*it)->IsSurface() ) {
        // finding the normal to the surface element
        Point<dim> line_vector1( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
        Point<dim> line_vector2( (*it)->N( 1 )->Coordinate() - (*it)->N( 2 )->Coordinate() );
        Point<dim> surface_normal( crossProduct( line_vector1, line_vector2 ) );

        // finding the gravity direction from the cross-product of the hor izontal surface vector and the gravity direction
        Point<dim> gravity_projection_on_surface( crossProduct( surface_normal, crossProduct( surface_normal, gravity_direction ) ) );
        gravity_projection_on_surface.NormalizeLengthTo( 1. );

        if ( gravity_projection_on_surface[1] > 0. ) gravity_projection_on_surface *= -1.;

        gvt( 0 ) = gravity_projection_on_surface[0];
        gvt( 1 ) = gravity_projection_on_surface[1];
        if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_surface[2];

        (*it)->Store( key_dipf, gvt );
      }
      // line elements
      else if ( (*it) != nullptr and (*it)->IsLine() ) {
        Point<dim> gravity_projection_on_line( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
        // flip gravity vector if it is upward pointing
        if ( (*it)->N( 0 )->y() > (*it)->N( 1 )->y() ) gravity_projection_on_line *= -1.;
        gravity_projection_on_line.NormalizeLengthTo( 1. );

        gvt( 0 ) = gravity_projection_on_line[0];
        gvt( 1 ) = gravity_projection_on_line[1];
        if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_line[2];

        (*it)->Store( key_dipf, gvt );
      }
      else {
        if ( (*it) != nullptr ) (*it)->Out();
        csmp_error.Note( WARNING, "CO2_GeoSequestrationSimulator::ComputeGravityDipVectors:",
                        "could not compute boundary normal for (?line?) element." );
      }
    }
  }
  
 } // end gravityDipVectors

template void computeGravityDipVectors( Model<3U>& );
template void computeGravityDipVectors( Model<2U>& );



/**
    Computes the CO2 saturation value above which there a wedge of CO2  reaches from the lowest point of the FV to the highest point.
    Method uses 'vertical extent' and 'finite volume diameter' in this calculation.
    
    The scalar node variable 'spill-point saturation' is computed from the ratio of the (approximate) volume of CO2 required to create CO2 pool
    that spans its diameter in the dip direction and its pre-computed (exact total) 'finite volume'.
*/
template<uint32_t dim>
void computeSpillPointSaturation( Model<dim>& model, const std::string& region_name )
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion( region_name ) )
      csmp_error.Note( ERROR, "computeSpillPointSaturation",
                      "target region does not exist");
                      
    Region<dim>& subdomain = model.Region( region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 )
      csmp_error.Note( ERROR, "computeSpillPointSaturation",
                      "target region is not lower dimensional (dim region != dim-1)");
                      
    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    // input variables
    const csmp::Index fvol_key = model.Database().StorageKey("finite volume");
    const csmp::Index thi_key  = model.Database().StorageKey("thickness");
    const csmp::Index diam_key = model.Database().StorageKey("finite volume diameter");
    const csmp::Index vext_key = model.Database().StorageKey("finite volume vertical extent");
    const csmp::Index dip_key  = model.Database().StorageKey("dip vector");
    // output variables
    const csmp::Index spill_key = model.Database().StorageKey("spill-point saturation");
    
    // computation: finding the ratio between the approximate spill-point volume and the actual finite volume
    // (porosity is ignored in this geometric analysis)
    // (the approximate diameter of the finite volume is used 2R), see notes in MS Word.
    // ---------------------------------------------------------------------------------
    // V_wedge = 1/2 Pi r^2 (h1 + h2) -> h1=0, h2 = 'vertical extent'
    VectorVariable<dim> dip_vec, flat_vec;
    
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         const double h2 = (*nit)->Read( vext_key );
         // if the FV lies in the horizontal plane
         if ( fabs(h2) <= numeric_limits<double>::epsilon() * 100. ) {
              (*nit)->Store( spill_key, makeScalar(ANY,0.) );
              continue;
           }
         // if the FV is tilted
         const double r          = (*nit)->Read( diam_key ) / 2.;
         const double V_wedge    = 0.5 * PI * (r*r) * h2;
         const double sCO2_spill = V_wedge / (*nit)->Read( fvol_key );
         if ( sCO2_spill > 1. ) {
              cerr <<"\nNode "<< (*nit)->Idx() <<": sCO2_spill: "<< sCO2_spill <<", location: "<< (*nit)->Coordinate();
              cerr <<", dip of FV: ";
              // computing average dip of the lower-dim elmts making up the FV sectors
              double   avg_dip{ 0. };
              uint32_t n_surf_elmts{ 0u };
              for ( auto i{0U}; i<(*nit)->Parents(); i++ ) {
                   if constexpr ( dim == 3U ) {
                        if ( (*nit)->Parent(i)->IsSurface() ) {
                             (*nit)->Parent(i)->Read( dip_key, dip_vec );
                             // calculating dip
                             flat_vec    = dip_vec;
                             flat_vec(1) = 0.;
                             double dip = dip_vec.AngleTo( flat_vec );
                             avg_dip += dip;
                             n_surf_elmts++;
                          }
                     }
                   else if constexpr ( dim == 2U ) {
                        if ( (*nit)->Parent(i)->IsLine() ) {
                             (*nit)->Parent(i)->Read( dip_key, dip_vec );
                             // calculating dip
                             flat_vec    = dip_vec;
                             flat_vec(1) = 0.;
                             double dip = dip_vec.AngleTo( flat_vec );
                             avg_dip += dip;
                             n_surf_elmts++;
                          }
                     }
                   static_assert( dim != 1U, "computeSpillPointSaturation: method does not work in 1D." );
                   avg_dip /= static_cast<double>(n_surf_elmts);
                }
              cerr << avg_dip;
              csmp_error.Note( WARNING, "computeSpillPointSaturation:", "computed value greater than 1.");
           }
         else (*nit)->Store( spill_key, makeScalar(ANY,sCO2_spill) );
      }
 
  } // end computeSpillPointSaturation



/**
      Checks whether thickness is large enough to permit spillage.
      Checks whether:
      - the 'thickness' of the FV large enough to allow for spillage to occur
      - the FV diameter less than the thickness calling into question a lower-dimensional representation  (else a lower-dim representation may not be warranted
      
      The results are written to the "finite volume diagnostics" parameter.
*/
template<uint32_t dim>
void lowerDimensionalLayerDiagnostics( Model<dim>& model, const string& region_name )
  {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion( region_name ) )
      csmp_error.Note( ERROR, "lowerDimensionalLayerDiagnostics",
                      "target region does not exist");
                      
    Region<dim>& subdomain = model.Region( region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 )
      csmp_error.Note( ERROR, "lowerDimensionalLayerDiagnostics",
                      "target region is not lower dimensional (dim region != dim-1)");
                      
    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    // input variables
    const csmp::Index thi_key  = model.Database().StorageKey("thickness");
    const csmp::Index diam_key = model.Database().StorageKey("finite volume diameter");
    const csmp::Index nrml_key = model.Database().StorageKey("finite volume normal");
    const csmp::Index vext_key = model.Database().StorageKey("finite volume vertical extent");
    // output variables
    const csmp::Index fvdi_key = model.Database().StorageKey("finite volume diagnostics");
    
    // obtaining the diagnostics (options are mutually exclusive)
    enum FV_DIAGNOSTICS { OK, NOT_THICK_ENOUGH, THICKNESS_GREATER_THAN_DIAMETER };
    VectorVariable<dim> nrml;

    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         (*nit)->Read( nrml_key );
         const double FV_diameter        = (*nit)->Read( diam_key );
         const double FV_vertical_extent = (*nit)->Read( vext_key );
         (*nit)->Read( nrml_key );

         // Test 1: Is thickness' of the FV large enough to allow spillage to occur
         // -----------------------------------------------------------------------
         // computing the mininum thickness of the surface elements connected to FV
         double minimum_thickness{ 1.0e+30 }; // crazy high value
         for ( auto i{0U}; i<(*nit)->Parents(); i++ ) {
             if constexpr ( dim == 3U ) {
                 if ( (*nit)->Parent(i)->IsSurface() )
                   minimum_thickness = min( minimum_thickness, (*nit)->Parent(i)->Read( thi_key ) );
               }
             if constexpr ( dim == 2U ) {
                 if ( (*nit)->Parent(i)->IsLine() )
                   minimum_thickness = min( minimum_thickness, (*nit)->Parent(i)->Read( thi_key ) );
               }
             static_assert( dim != 1U, "lowerDimensionalLayerDiagnostics: method cannot be applied in 1D" );
           }
         // is 'minimum_thickness' >= vertical extent
         if ( minimum_thickness < FV_vertical_extent )
           (*nit)->Store( fvdi_key, makeScalar(ANY,NOT_THICK_ENOUGH) );

         // Test 2: Is thickness greater than the diameter of the FV
         // --------------------------------------------------------
         if ( minimum_thickness > FV_diameter )
           (*nit)->Store( fvdi_key, makeScalar(ANY,THICKNESS_GREATER_THAN_DIAMETER) );
      }

  } // end lowerDimensionalLayerDiagnostics

template void lowerDimensionalLayerDiagnostics( Model<2U>&, const string& );
template void lowerDimensionalLayerDiagnostics( Model<3U>&, const string& );



} // csmp
