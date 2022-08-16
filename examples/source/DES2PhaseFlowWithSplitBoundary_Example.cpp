#include "DES2PhaseFlowWithSplitBoundary_Example.h"

// the CSMP model
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"

// FV algorithms
#include "DES2PhaseSlightlyCompressibleTransport.h"

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
  AddRequirement( "box2d_fault: files .dat, .asc, -regions.txt, -configuration.txt.");  
  AddRequirement( "variables(DES_2phase_variables.txt)" );
} 

/** 
    Two phase slightly compressible flow simulation with splitboundaries via CSMP's DES transport method
    combining finite elements (for pressure) with finite volumes (for advection of non-wetting phase)

    Use models 'box2d_fault' (.dat, .asc, -regions.txt, -configuration.txt) as input file suites.
*/

  void DES2PhaseFlowWithSplitBoundary_Example::Run()
  {

    // model dimension
    uint32_t dimension;
    cerr << "\nPlease enter the dimension of the model (2 for 2D, 3 for 3D):" << endl;
    cin >> dimension;
    if (dimension != 2U and dimension != 3U)
      throw csmp::Exception(ERROR, "input dimension of model", "must be 2 or 3");

    // create model directly from ANSYS-ICEM mesh
    string input_file;
    cerr << "\nPlease enter the name of input mesh (default: box2d_fault):" << endl;
    cin >> input_file;

    string variables_file("DES_2phase_variables.txt");

    if (dimension == 2U) {
      ANSYS_Model2D model(input_file.c_str(), variables_file.c_str());
      RunSimulation(model);
    } else if (dimension == 3U) {
      ANSYS_Model3D model(input_file.c_str(), variables_file.c_str());
      RunSimulation(model);
    }

  }


  template<uint32_t dim>
  void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation(Model<dim>& model) {
  //void RunSimulation(Model<dim>& model) {

    // simulation settings
    Standard_IO_Handler stdio;
    VTU_Interface<dim>  vtu(model);

    bool  DES = stdio.YesNo("Do you want to solve the transport equation with DES? (y=DES, n=TDS)"); 
    double Courant_multiplier, PEP_parameter;
    cerr <<"\nEnter CFL multiplier (suggested value: 0.3 for TDS, 0.1 for DES) and PEP parameter (suggested value: 0.1)" << endl;
    cin >> Courant_multiplier >> PEP_parameter;  
     
    bool  with_gravity_forces = stdio.YesNo("Do you want to include gravity effect (y/n)?"); 
    bool  with_capillary_spreading = stdio.YesNo("Do you want to include capillary effect (y/n)?");

    bool  with_split_boundaries = stdio.YesNo("Do you want to create split boundaries (y/n)?");

    // give the model dimensions
    printModelDimensions( model, true );

    // Configure the simulation from a file
    InputDataManager<dim>  model_configuration;
    model_configuration.ConfigureFromFile( model, model.Name(),
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

    // flow functions (Brooks Corey)
    FlowFunctionsModule1<dim> flowfunctions(model.Database(), model.Read( model.Database().StorageKey("acceleration gravity") ));

    // solve static pressure before split boundaries are created
    Compute2PhaseFlowProperties(model, flowfunctions, with_gravity_forces, with_tensor_permeability);
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
        model.CreateProperty( "number of collocated nodes", "none", SCALAR, NODE, 1, 0, 100);;

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
    if( with_split_boundaries ) {
      const csmp::Index pd_key = model.Database().StorageKey("entry pressure");
      for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++)
        (*mit).SortByVariableValue(pd_key);

      //check parent elements of duplicated nodes
      if(!model.Database().IsDefined("parent element id")) {
        model.CreateProperty( "parent element id", "none", SCALAR, ELEMENT, 1, 0 ,1.00E+08);
        model.Region("Model").InputPropertyValue( "parent element id", makeScalar(PLAIN,0), COMPLETE);
      }
      csmp::INDEX<SCALAR,ELEMENT> key_parent = csmp::INDEX<SCALAR,ELEMENT>( model.Database().StorageKey("parent element id") );
      size_t parent_id(1);
      for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++) {
        parent_id = 1;
        auto md = (*mit);
        auto master_node = md.N(0);
        for(size_t e = 0; e < master_node->Parents(); e++)
          master_node->Parent(e)->Store(key_parent, makeScalar(PLAIN, parent_id));

        for(size_t n=1;n<md.Branches();n++) {
          auto slave_node = md.N(n);
          parent_id++;
          for(size_t e = 0; e < slave_node->Parents(); e++)
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
                                                                                             1., //relaxing factor
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
    if(with_split_boundaries) {
      input_properties.emplace_back("pressure continuity status");
      input_properties.emplace_back("breakthrough status");
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
      model.CreateProperty( "out range value count", "none", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
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

    // setting up time parameters
    const double    day(86400.0);
    double max_time (50.0*day); //run for 50 days
    double time_increment(0.5*day); //timestep 0.5 day
    double save_interval = 2.0*day; //save every 2 days

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
      // solve pressure
      Compute2PhaseFlowProperties( model, flowfunctions, with_gravity_forces, with_tensor_permeability );
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
  //template void RunSimulation(Model<2U>& model);
  //template void RunSimulation(Model<3U>& model);



  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void DES2PhaseFlowWithSplitBoundary_Example::Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k )
  //void Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k )
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
    //vector<Node<dim>* >::const_iterator nit;
    //cout<<"number of elements = "<<mref.Cells()<<endl;
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
        for ( size_t i=0U; i<nodes; ++i ) {
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
          //gravity_w = gravity_acceleration * e_k * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw * (e_rhow);
          //gravity_n = gravity_acceleration * e_k * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun * (e_rhon);
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
        for ( size_t i=0U; i<nodes; ++i ) {
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
  //template void Compute2PhaseFlowProperties( Model<2U>&, FlowFunctionsModule1<2U>&, bool, bool );
  //template void Compute2PhaseFlowProperties( Model<3U>&, FlowFunctionsModule1<3U>&, bool, bool );



  template<uint32_t dim>
  void DES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure(Model<dim>& mdl, bool with_gravity, bool with_tensor_k)
  //void ComputeSteadyStatePressure(Model<dim>& mdl, bool with_gravity, bool with_tensor_k)
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

    // SAMG Solver
    //SAMG_Settings settings;
    //SAMG_Solver                    solver(&settings);
    EigenSolver  solver;
    PDE_Integrator<dim,Region>      steady_pressure(solver);
    NumIntegral_dNT_op_dN_dV<dim>   conductance( mdl.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<dim>     elmt_volume_source( mdl.Database(), "fluid volume source", "fluid pressure" );
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
  //template void ComputeSteadyStatePressure(Model<2U>&, bool, bool);
  //template void ComputeSteadyStatePressure(Model<3U>&, bool, bool);



} // csmp
