#include "ThermalConvectionETHZ_Example.h"

// ------------------------------------------
// Created by James Patterson on 7.5.2015
// Modeified by Stephan Matthai, James Patterson 
// and Thomas Driesner, July 2015
// ------------------------------------------

// inclusion of class declarations necessary for this simulation
#include "CSMP_definitions.h"

#include "Region.h"
#include "ErrorHandler.h"
#include "VTU_Interface.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

#include "ModelTime.h"

// finite volumes 
#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"
#include "NodeCenteredFiniteVolumeTransport.h"

// finite elements
#include "PDE_Integrator.h"
// left-handside
#include "NumIntegral_dNT_op_dN_dV.h"  // mobility matrix
#include "NumIntegral_NT_lhsop_N_dV.h"
// righthand side
#include "NumIntegral_NT_op_N_dV.h"    // mass matrix
#include "NumIntegral_dNT_op_dV.h"     // gravity term
#include "NumIntegral_NT_op_N_dS.h"
#include "PointSource_rhsop.h"

#include "H2OLookup.h"

// utilities and monitoring
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"
//#include "RegionMonitor.h"


using namespace std;
namespace csmp {

  
  void ThermalConvectionETHZ_Example::Specifications()
  {
    SetTitle( "Thermal convection with H2O equation of state and advanced spatial discretization" );
    SetDifficulty( 3 );
    SetCategory( "Simulation of Physical Processes" );
    AddAuthor( "JWP+SKM+TD" );
    AddDescription( "Detailed description will follow" );
    AddDescription( "source file in: ThermalConvectionETHZ_Example" );
    AddRequirement( "mesh(3D_box), configuration(compressible_flow-variables.txt)" );
  }


  /*
    SKM QUESTIONS AND COMMENTS
 
    - cleanup of names in variable file may be necessary: in my mind, i associate 'content' with an absolute
    = integrated property, but this does not seem to be the case
      
    - 'nodal fluid volume source'  should perhaps be 'nodal fluid mass source' ?
  */

  //=================================================================================================
  void ThermalConvectionETHZ_Example::Run()
  {
    const  uint32_t DIM(3U);
    bool restart = false;
    std::string restart_name ("3D_box");    // 25 x 25 x 5-m small test module
    //  std::string restart_name ("3D_box2"); // as above but more refined
    //  std::string restart_name ("bu615"); // quadratic FEM test model
    //  std::string restart_name ("bwf");
    std::string geometry_name (restart_name);
    std::string regions_name (restart_name);
    const string config_name (restart_name); 
    std::string vars_name("compressible_flow-variables.txt");
  
    Model<DIM>* model(0);
    double& model_time( ModelTime::Instance().modelTime );

    if (!restart) {
      model = new ANSYS_Model3D(geometry_name.c_str(),
                                regions_name.c_str(),
                                vars_name.c_str(), true );
    }
    else model = new Model<DIM>(restart_name);
  
    long time_step(0);
    if (restart) {
      time_step = model->Read(model->Database().StorageKey("time step"));
      time_step++;
      cout<<"\n restarting computations at time =  "<<model_time;
      cin.get();
    }
  
    printModelDimensions<DIM>(*model, true );

    VTK_Interface<DIM> vtk_output;		//enable VTK output
  

    // -----------------------------------------------------------------------
    // 1. Assigning material properties, initial conditions, and boundary
    //    conditions from *-configuration.txt files.
    // -----------------------------------------------------------------------
    InputDataManager<DIM>  model_configuration;
    ComputationalSettings  run_settings;

    // NEEDS TO BE INITIALIZED HERE SO THAT finite-volume stencils get created and connected,
    // else FV-related variables will not be read from configuration file !
    model->InstantiateFiniteVolumes();
  
    if (!restart)
      model_configuration.ConfigureFromFile( *model, config_name.c_str(), false,    // groupname from parameter range
                                             true,    // default property values
                                             true,    // regional property values
                                             false,   // boundary conditions for box-shaped model
                                             true,    // essential conditions for groups
                                             true,    // csmp::Boundary properties
                                             run_settings );
  
    Region<DIM>& computation_domain = model->Region("Model");

 
    //! -------------------------------------------------------------
    //! 2. Temperature and pressure initialisation
    //! -------------------------------------------------------------
    H2OLookup  water;
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    #endif
    {
      #ifdef CSMP_WITH_SAMG_SOLVER
      PDE_Integrator<DIM,Element>  heat_conductor( *(new SAMG_Solver(&settings)));
      #else
      CSMP_DEFAULT_LINEAR_SOLVER   solver;
      PDE_Integrator<DIM,Element>  heat_conductor(solver);
      #endif
    
      NumIntegral_dNT_op_dN_dV<DIM> conductance( model->Database(), "thermal conductivity", "temperature", "temperature");
      NumIntegral_NT_op_N_dV<DIM>  heat_source( model->Database(), "heat source", "temperature" );
      NumIntegral_NT_op_N_dS<DIM>  basal_hfu( model->Database(), "basal heat flow", "temperature" );
       
      heat_conductor.Add( &conductance );
      heat_conductor.Add( &heat_source );
      heat_conductor.AddBoundaryIntegral(&basal_hfu);

      #ifdef CSMP_WITH_SAMG_SOLVER
      // solve down to round-off
      settings.Set_eps(0.); // absolute criterion
      #endif
     
      // Compute initial temperature profile
      heat_conductor.IntegrateOver( *model, computation_domain );
        
      printRangeOfVariable( *model, "temperature");

      //! computation of initial steady-state fluid pressure
      #ifdef CSMP_WITH_SAMG_SOLVER
      PDE_Integrator<DIM,Element>     steady_state_pressure( *(new SAMG_Solver(&settings) ));
      #else
      PDE_Integrator<DIM,Element>     steady_state_pressure( solver );
      #endif
      NumIntegral_dNT_op_dN_dV<DIM>  p_conductance( model->Database(), "mass conductivity", "fluid pressure", "fluid pressure" );
      NumIntegral_dNT_op_dV<DIM>     gravity( model->Database(), "gravity term times density", "fluid pressure" );

      steady_state_pressure.Add( &p_conductance );                           
      steady_state_pressure.Add( &gravity );

      #ifdef CSMP_WITH_SAMG_SOLVER
      // iout is reduced to bare minimum
      settings.Set_iout1( 0 );
      settings.Set_iout2( 0 );
      #endif
     
      // initialize gravity vector
      VectorVariable<DIM>  grav_vec(DIRICH,0.);
      grav_vec(1) = -9.80665;
      // TODO: this assignment is only correct if there are no lower dimensional elements; add l.d. in plane or line vectors
      model->InputPropertyValue( "gravity vector", grav_vec );
      
      // Iterate to find an initial vertical pressure profile
      // --------------------------------------------------------
      // initial pressure must be initialized in configuration file
      for ( int i=0; i<3; i++ ) {
        cout <<"\n\n\nmain: pressure initialisation loop; iteration: "<< i+1;
        steadyStatePressureOperands( *model, computation_domain, water );
        steady_state_pressure.IntegrateOver( *model, computation_domain );
        monitorVariableRange( *model );
        //printRangeOfVariable( *model, "fluid pressure" );
      }
      massFlux( *model, computation_domain );
    }
  
    outputToVTK( *model, vtk_output, ++time_step );
    //! -------------------------------------------------------------
    //! 3. Algorithms for transient pressure and temperature
    //! -------------------------------------------------------------
    // fluid pressure
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings pressure_settings;
    //pressure_settings.SetSolverInstance(1);
    PDE_Integrator<DIM,Element>  transient_pressure( *(new SAMG_Solver(&settings)) );
    #else
    CSMP_DEFAULT_LINEAR_SOLVER   solver;
    PDE_Integrator<DIM,Element>  transient_pressure(solver);
    #endif
    NumIntegral_dNT_op_dN_dV<DIM>  p_conductance( model->Database(), "mass conductivity", "fluid pressure", "fluid pressure" );

    NumIntegral_dNT_op_dV<DIM>     gravity( model->Database(), "gravity term times density", "fluid pressure" );
    gravity.AddAccumulateLater();
    gravity.MultiplyWithTimeIncrement(true);
    NumIntegral_NT_op_N_dV<DIM>    element_fluid_source( model->Database(), "fluid volume source", "fluid pressure" );
    element_fluid_source.AddAccumulateLater();
    // to account for absolute nodal fluid contributions due to PVT property effects
    PointSource_rhsop<DIM>         node_total_source( model->Database(), "nodal fluid volume source", "fluid pressure" );
    node_total_source.MultiplyWithTimeIncrement(false);
    node_total_source.LumpedFormulation( true );
    node_total_source.AddAccumulateLater();

    NumIntegral_NT_op_N_dV<DIM>    storage_rhs( model->Database(), "total compressibility", "fluid pressure" );
    storage_rhs.LumpedFormulation(true);
    storage_rhs.MultiplyWithTimeIncrement(true);

    NumIntegral_NT_lhsop_N_dV<DIM> storage_lhs( model->Database(), "total compressibility", "fluid pressure", "fluid pressure" );
    storage_lhs.LumpedFormulation(true);
    storage_lhs.MultiplyWithTimeIncrement(true);

    transient_pressure.Add( &p_conductance );
    transient_pressure.Add( &gravity);
    transient_pressure.Add( &element_fluid_source);
    transient_pressure.Add( &node_total_source);
    transient_pressure.Add( &storage_lhs );
    transient_pressure.Add( &storage_rhs );

    // temperature
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings temperature_settings;
    //temperature_settings.SetSolverInstance(2);
    PDE_Integrator<DIM,Element>  transient_temperature( *(new SAMG_Solver(&temperature_settings)) );
    #else
    CSMP_DEFAULT_LINEAR_SOLVER   solver2;
    PDE_Integrator<DIM,Element>  transient_temperature(solver2);
    #endif
    NumIntegral_dNT_op_dN_dV<DIM>  T_conductance( model->Database(), "thermal conductivity", "temperature", "temperature" );
    NumIntegral_NT_op_N_dV<DIM>    heat_source( model->Database(), "heat source", "temperature" );
    heat_source.AddAccumulateLater();

    PointSource_rhsop<DIM>         node_heat_source( model->Database(), "nodal heat source", "temperature" );
    node_heat_source.AddAccumulateLater();
    node_heat_source.MultiplyWithTimeIncrement(false);
    node_heat_source.LumpedFormulation( true );
    
    NumIntegral_NT_op_N_dV<DIM>    thermal_capacitance_rhs( model->Database(), "total heat capacity", "temperature" );
    thermal_capacitance_rhs.LumpedFormulation(true);
    thermal_capacitance_rhs.MultiplyWithTimeIncrement(true);

    NumIntegral_NT_lhsop_N_dV<DIM> thermal_capacitance_lhs( model->Database(), "total heat capacity", "temperature", "temperature" );
    thermal_capacitance_lhs.LumpedFormulation(true);
    thermal_capacitance_lhs.MultiplyWithTimeIncrement(true);

    NumIntegral_NT_op_N_dS<DIM>  basal_hfu( model->Database(), "basal heat flow", "temperature" );
    basal_hfu.AddAccumulateLater();

    transient_temperature.Add( &T_conductance );
    transient_temperature.Add( &heat_source);
    transient_temperature.Add( &node_heat_source);
    transient_temperature.Add( &thermal_capacitance_rhs );
    transient_temperature.Add( &thermal_capacitance_lhs );
    transient_temperature.AddBoundaryIntegral( &basal_hfu );


    // set-up explicit enthalpy transport scheme
    const bool second_order_accuracy(false);
    ExplicitMassBasedTransport<DIM,MassBasedStencilProcessor>
      mass_advector( "Model", //group
                     *model,  //model
                     "porosity", //porosity,
                     "fluid density", //advected_prop_lhs
                     "density liquid", //advected_prop_rhs
                     "Darcy velocity", //transp_velocity
                     "zero node", //nodal_source
                     second_order_accuracy ); //second_order_accuracy

    ExplicitMassBasedTransport<DIM,MassBasedStencilProcessor>
      heat_advector( "Model", *model,
                     "porosity",
                     "enthalpy content liquid", "volumetric enthalpy liquid",
                     "Darcy velocity",
                     "zero node",
                     second_order_accuracy );

    // initialising "finite volume"  and "effective FV porosity" and variables for transient calculation
    initialiseFiniteVolumeProperties( *model );
    printRangeOfVariable( *model, "finite volume" );
    printRangeOfVariable( *model, "FV pore volume" );
    printRangeOfVariable( *model, "sector volume" );
    printRangeOfVariable( *model, "sector weight" );
    // verifying model volume
    computation_domain.InputPropertyValue( "test variable", makeScalar(ANY,1.));
    cerr <<"\nmain: model volume: "<< computation_domain.VolumeIntegral( "test variable", false, false ) <<" m3.\n";

    // nodal fluid density and enthalpy content and 'previous sector enthalpy content'
    initialFluidProperties( *model, computation_domain, water );
    transientOperands( *model, computation_domain, water );
    
    //! -------------------------------------------------------------
    //! 4. Time evolution loop
    //! -------------------------------------------------------------
    monitorVariableRange( *model );
    double simulated_time(0.);
  
    while ( simulated_time < run_settings.Duration() )
      {
        cout <<"\nmain: time-stepping loop; step: "<< time_step << endl;
        double time_increment = mass_advector.AnisotropicCourantIncrement();
        // 1. transient temperature (starting with an initialized model with updated fluid properties)
        transient_temperature.TimeIncrement( 1. / time_increment );
        transient_temperature.IntegrateOver( *model, computation_domain );
        //printRangeOfVariable( *model, "temperature" );
        //printRangeOfVariable( *model, "fluid pressure" );
        
        // 2. advect mass and enthalpy
        cout <<"\ninput variables for transport scheme:";
        mass_advector.AdvectVariable( time_increment );
        heat_advector.AdvectVariable( time_increment );
        
        // 3. thermal equilibration
        thermalEquilibrationAtSectorIntegrationPoints( *model, computation_domain, water );
        monitorVariableRange( *model );
        // 4. calculation of source/sink coupling term with the pressure equation
        nodalSourceOperands( *model, computation_domain, water, time_increment );
        monitorVariableRange( *model );        
        // 5. solve transient pressure equation
        transient_pressure.TimeIncrement( 1. / time_increment );
        transient_pressure.IntegrateOver( *model, computation_domain );
        monitorVariableRange( *model );
        transientFluidProperties( *model, computation_domain, water );
        monitorVariableRange( *model );
        transientOperands( *model, computation_domain, water );
        monitorVariableRange( *model );
        massFlux( *model, computation_domain );
        cerr << "\nafter massFlux:\n";
        printRangeOfVariable( *model, "fluid pressure" );
        printRangeOfVariable( *model, "nodal fluid volume source" );
        printRangeOfVariable( *model, "density difference" );
        printRangeOfVariable( *model, "mass conductivity");
        printRangeOfVariable( *model, "gravity term times density");
        printRangeOfVariable( *model, "Darcy velocity" );
        printRangeOfVariable( *model, "total compressibility" );
        printRangeOfVariable( *model, "density liquid" );
        printRangeOfVariable( *model, "fluid density" );
        cerr << endl;
        // 6. output results at selected time intervals
        outputToVTK( *model, vtk_output, ++time_step );
      }
  
    delete model;
    cerr << "\n\nall done.\n\n";
    //return 0;
  
    //=================================================================================================
  } // end example
  //=================================================================================================






  // various tests used in the creation of this model

  /*
  // CHECKING BOUNDARY functionality
  // ----------------------------------------------------------------------
  cout <<"\nmain: model volume: "<< model->Region("Model").Volume() <<"\n";
  // @test testing bottom boundary
  Boundary<DIM>& bref(model->Boundary("BOTTOM"));
  const csmp::Index hfu_key = model->Database().StorageKey("basal heat flow");
  cout <<"\nmain: bottom boundary area: "<< bref.Area() <<", faces: "<< bref.Elements() <<"\n";
  for ( vector<Face<DIM>*>::const_iterator
  it=bref.ElementsBegin(); it!=bref.ElementsEnd(); it++ )
  if ( (*it)->Status( hfu_key ) != NEUMANN )
  cerr <<"\nmain: one of the b-faces was not flagged NEUMANN\n";
  printRangeOfVariable( *model, "BOTTOM", "basal heat flow" );
  // ----------------------------------------------------------------------

  // TESTING READING AND WRITING OF INTEGRATION POINT VARIABLES
  cerr <<"\nchecking pore volume computed from sector volumes and porosity\n";
  const csmp::Index phi_key = model->Database().StorageKey("porosity");
  const csmp::Index fv_key    = model->Database().StorageKey("finite volume");
  const csmp::Index fvphi_key = model->Database().StorageKey("FV pore volume");
  const csmp::Index sv_key    = model->Database().StorageKey("sector volume");
  const csmp::Index swt_key   = model->Database().StorageKey("sector weight");
  double total_volume(0.), total_pore_volume(0.);
  for ( vector<Element<3U>*>::iterator
  it=computation_domain.ElementsBegin(); it!=computation_domain.ElementsEnd(); ++it )
  {
  const double phi = (*it)->Read( phi_key );
  double evolume(0.);
  for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
  cerr <<"\n\t\t"<< (*it)->Read(i,0U,sv_key) <<" "<< (*it)->SectorVolume(i);
  evolume           += (*it)->Read( i, 0U, sv_key ); // (*it)->SectorVolume(i);
  total_volume      += (*it)->Read( i, 0U, sv_key );
  total_pore_volume += phi * (*it)->Read( i, 0U, sv_key );
  }
  if ( fabs(evolume - (*it)->Volume()) > numeric_limits<double>::epsilon() )
  cerr <<"\n\te"<< (*it)->Idx() <<": fv vs. e "<< evolume <<" "<< (*it)->Volume();
  }
  cerr <<"\ncomputed volume and  pore volume: "<< total_volume <<" "<< total_pore_volume << endl;

  */



  /**
     We store the volume of the finite volumes and their pore volumes
     
     - finite volume
     - finite volume (effective) pore volume
     - sector volume (stored at sector integration point)
     - sector weight: sector pore volume / finite volume pore volume = weighting factor
  */
  void ThermalConvectionETHZ_Example::initialiseFiniteVolumeProperties( Model<3U>& model )
  {
    const csmp::Index phi_key   = model.Database().StorageKey("porosity");

    const csmp::Index fv_key    = model.Database().StorageKey("finite volume");
    const csmp::Index fvphi_key = model.Database().StorageKey("FV pore volume");
    const csmp::Index sv_key    = model.Database().StorageKey("sector volume");
    const csmp::Index swt_key   = model.Database().StorageKey("sector weight");
   
    Region<3U>& ref = model.Region("Model");
   
    ref.InputPropertyValue( "finite volume", makeScalar( ANY, 0. ) );
    ref.InputPropertyValue( "FV pore volume", makeScalar( ANY, 0. ) );

    // initializing sector und FV volumes including pore volumes
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        const double phi = (*it)->Read( phi_key );
        
        for ( uint32_t i=0U; i<(*it)->Sectors(); ++i ) {
          // sector volume
          const double sec_vol = (*it)->SectorVolume(i);
          (*it)->Store( i, 0U, sv_key, makeScalar(PLAIN,sec_vol) );
          // accumulating the FV volume
          double FV_vol = (*it)->N(i)->Read( fv_key );
          FV_vol += sec_vol;
          (*it)->N(i)->Store( fv_key, makeScalar(PLAIN,FV_vol) );
          // accumulating FV pore volume
          double FVpore_vol = (*it)->N(i)->Read( fvphi_key );
          FVpore_vol += sec_vol * phi;
          (*it)->N(i)->Store( fvphi_key, makeScalar(PLAIN,FVpore_vol) );
        }
      }

    // initializing sector weights: (sector_V * phi) / FV pore volume
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        const double pore_volume( (*nit)->Read( fvphi_key ) );
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
          Element<3U>* const eptr((*nit)->Parent(i));
          const uint32_t sector = (*nit)->ParentNodeNumber(i);
          const double weight = (eptr->Read( sector, 0U, sv_key ) * eptr->Read( phi_key )) / pore_volume;
          //              const double weight = (eptr->SectorVolume(sector) * eptr->Read( phi_key )) / pore_volume;
          eptr->Store( sector, 0U, swt_key, makeScalar(PLAIN,weight) );
        }
      }

    cout <<"\ninitializeFiniteVolumeProperties: checking 'finite volume' and pore volume computed from sector volumes and porosity\n";
    double total_volume(0.), total_pore_volume(0.);
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        const double phi = (*it)->Read( phi_key );
        double evolume(0.);
        for ( uint32_t i=0U; i<(*it)->Sectors(); ++i ) {
          evolume           += (*it)->Read( i, 0U, sv_key ); // (*it)->SectorVolume(i);
          total_volume      += (*it)->Read( i, 0U, sv_key );
          total_pore_volume += phi * (*it)->Read( i, 0U, sv_key );
        }
        if ( fabs(evolume - (*it)->Volume()) > numeric_limits<double>::epsilon() )
          cerr <<"\n\te"<< (*it)->Idx() <<": fv vs. e "<< evolume <<" "<< (*it)->Volume();
      }
    cout <<"\ncomputed volume and  pore volume: "<< total_volume <<" "<< total_pore_volume << endl;


    cout <<"\ninitializeFiniteVolumeProperties: checking sector weights (should add up to 1).\n";
    // verifying that the weights add up to one
    double min_weight(1e30), max_weight(-1e30);
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        double total_weight(0.);
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
          const Element<3U>* eptr((*nit)->Parent(i));
          const uint32_t sector = (*nit)->ParentNodeNumber(i);
          const double weight = eptr->Read( sector, 0U, swt_key );
          total_weight += weight;
        }
        min_weight = std::min( min_weight, total_weight );
        max_weight = std::max( max_weight, total_weight );
      }
    cout <<"\n\tmin vs. max of the sector weights summed up over the finite volumes: "<< min_weight <<" vs. "<< max_weight << endl;
   
    cout <<"\nchecking writing of a sector variable (range should be 1).\n";
    const csmp::Index tv_key   = model.Database().StorageKey("test variable");
    ref.InputPropertyValue( "test variable", makeScalar( ANY, 0. ) );
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      for ( uint32_t i=0U; i<(*it)->Sectors(); ++i ) {
        double test = (*it)->N(i)->Read( tv_key );
        test += (*it)->Read( i, 0U, swt_key );
        (*it)->N(i)->Store( tv_key, makeScalar(PLAIN,test) );
      }
    printRangeOfVariable( model, "test variable" );
 
  } // end initializeFiniteVolumeProperties




  /* TESTING SECTOR INTEGRATION POINT STORAGE

  // sector storage: writing global node numbers to sector IP's and reading them out
  for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
  for ( size_t i=0U; i<(*it)->Sectors(); ++i )
  (*it)->Store( i, 0U, swt_key, makeScalar(PLAIN,(*it)->N(i)->Idx()) );
      
  // reading out node numbers and their double equivalents stored at the sector integration points
  for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it ) {
  cerr <<"\nelement: "<< (*it)->Idx() << endl;
  for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
  cerr << (*it)->N(i)->Idx() <<":";
  cerr << (*it)->Read( i, 0U, swt_key ) <<" ";
  }
  }
  */



  /*  TESTING FV volume calculations
      double volume(0.);
      for ( vector<Node<3U>*>::iterator it=ref.NodesBegin(); it!=ref.NodesEnd(); ++it )
      volume += (*it)->Read( fv_key );
      cerr <<"\nFV total volume: "<< volume;

      volume = 0.;
      for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      for ( size_t i=0U; i<(*it)->Sectors(); ++i )
      volume += (*it)->Read( i, 0U, sv_key );
      cerr <<"\nFV total volume: "<< volume;

      double pvolume(0.);
      for ( vector<Node<3U>*>::iterator it=ref.NodesBegin(); it!=ref.NodesEnd(); ++it )
      pvolume += (*it)->Read( fvphi_key );
      cerr <<"\nFV total volume: "<< pvolume;
  */





  /** 
      Uses the initial guesses of p, T to calculate fluid mass and total fluid energy
      content in the finite volume sectors.
  */
  void ThermalConvectionETHZ_Example::initialFluidProperties( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
  {
    // at the finite-element nodes
    const csmp::Index p_key   = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key   = model.Database().StorageKey("temperature");
    const csmp::Index FVT_key = model.Database().StorageKey("FV temperature");
    
    const csmp::Index mt_key  = model.Database().StorageKey("fluid density");
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index mu_key  = model.Database().StorageKey("fluid viscosity");
    const csmp::Index bf_key  = model.Database().StorageKey("fluid compressibility");
    const csmp::Index cpf_key = model.Database().StorageKey("fluid heat capacity");
    const csmp::Index hl_key  = model.Database().StorageKey("enthalpy liquid");
    const csmp::Index hVl_key = model.Database().StorageKey("volumetric enthalpy liquid");
    const csmp::Index hCl_key = model.Database().StorageKey("enthalpy content liquid");
    // at the finite element FV sector integration points
    const csmp::Index hCSl_key  = model.Database().StorageKey("sector enthalpy content liquid");
    const csmp::Index hCSlp_key = model.Database().StorageKey("previous sector enthalpy content liquid");
   
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        const double p  = (*nit)->Read( p_key );
        const double T  = (*nit)->Read( T_key );
        (*nit)->Store( FVT_key, makeScalar(PLAIN, T ) );
        (*nit)->Store( mt_key,  makeScalar(PLAIN, water.Density( T, p )) );
        (*nit)->Store( mu_key,  makeScalar(PLAIN, water.Viscosity( T, p ) ) );
        (*nit)->Store( bf_key,  makeScalar(PLAIN, water.Compressibility( T, p ) ) );
        (*nit)->Store( cpf_key, makeScalar(PLAIN, water.HeatCapacity( T, p ) ) );
        (*nit)->Store( hl_key,  makeScalar(PLAIN, water.Enthalpy( T, p ) ) );
        (*nit)->Store( hVl_key, makeScalar(PLAIN, water.Density( T, p ) * water.Enthalpy( T, p )) );
        (*nit)->Store( hCl_key, makeScalar(PLAIN, water.Density( T, p ) * water.Enthalpy( T, p )) );
      }

    // FV sector (specific=not integrated) enthalpy content computed at the sector integration points
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      for ( uint32_t i=0U; i<(*it)->Sectors(); ++i ) {
        const double p  = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, p_key );
        const double T  = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, T_key );
        const double h_content = water.Density( T, p ) * water.Enthalpy( T, p );
        (*it)->Store( i, 0U, hCSl_key,  makeScalar(PLAIN,h_content) );
        (*it)->Store( i, 0U, hCSlp_key, makeScalar(PLAIN,h_content) );
      }
   
  } // end initialFluidProperties



  /**
     A distinction is made between the actual fluid density (density liquid)
     and the current mass of fluid per pore-volume (in the transient calculation with advection).
  */
  
  void ThermalConvectionETHZ_Example::transientFluidProperties( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
  {
    const csmp::Index p_key   = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key   = model.Database().StorageKey("temperature");

    // density at given p, T
    const csmp::Index mu_key  = model.Database().StorageKey("fluid viscosity");
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index cf_key  = model.Database().StorageKey("fluid compressibility");
    const csmp::Index cpf_key = model.Database().StorageKey("fluid heat capacity");
       
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        const double p   = (*nit)->Read( p_key );
        const double T   = (*nit)->Read( T_key );
        //const double rho(1000.); // = water.Density( T, p );
        //const double mu(1e-3); //  = water.Viscosity( T, p );
        const double mu  = water.Viscosity( T, p );
        const double rho = water.Density( T, p ); 
        const double cf  = water.Compressibility( T, p );
        const double cpf = water.HeatCapacity( T, p );
                
        // storing the variables
        (*nit)->Store( mu_key,  makeScalar(PLAIN,mu) );
        (*nit)->Store( rho_key, makeScalar(PLAIN,rho) );
        (*nit)->Store( cf_key,  makeScalar(PLAIN,cf) );
        (*nit)->Store( cpf_key, makeScalar(PLAIN,cpf) );
      }
   
  } // end transientFluidProperties


  /**
     mobility, gravity term at element integration points
     TODO: update for lower-dimensional elements
  */
  void ThermalConvectionETHZ_Example::steadyStatePressureOperands( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
  {
    const csmp::Index p_key     = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key     = model.Database().StorageKey("temperature");
    
    const csmp::Index k_key     = model.Database().StorageKey("permeability");
    const csmp::Index rho_key   = model.Database().StorageKey("density liquid");
    const csmp::Index mu_key    = model.Database().StorageKey("fluid viscosity");
    const csmp::Index gv_key    = model.Database().StorageKey("gravity vector");

    const csmp::Index mc_key    = model.Database().StorageKey("mass conductivity");           // (rho k) / mu
    const csmp::Index rhogt_key = model.Database().StorageKey("gravity term times density");  // (rho^2 k g) / mu
    const csmp::Index gt_key    = model.Database().StorageKey("gravity term");
    
   
    VectorVariable<3U>  gv; // gravity vector
   
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        const double k = (*it)->Read( k_key );
                
        for ( uint32_t i=0U; i<(*it)->IntegrationPoints(); ++i ) {
          (*it)->Read( gv_key, gv ); // read here to get original value, not overwritten one          
          // compute fluid properties at integration point from interpolated T and p
          const double T = (*it)->PropertyValueAtIntegrationPoint( T_key, i );
          const double p = (*it)->PropertyValueAtIntegrationPoint( p_key, i );
         
          double rho, mu;
          rho  = water.Density( T, p );
          mu   = water.Viscosity( T,p );
          
          // computing and storing mass mobility value at the integration point
          (*it)->Store( i, mc_key, makeScalar(PLAIN, (k * rho) / mu ) );
           
          // computing the gravity term (mass-based p-equation, i.e. rho^2) for the righthand side
          gv *= (k / mu) * rho;
          (*it)->Store( i, gt_key, gv );
          gv *= rho;
          (*it)->Store( i, rhogt_key, gv );
        }
      }

  } // end steadyStateOperands



  /**
     Computing mass conductivity, the gravity term, total systems compressibility, and total heat capacity
  */
  void ThermalConvectionETHZ_Example::transientOperands( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
  {
    // rock properties
    const csmp::Index k_key      = model.Database().StorageKey("permeability");
    const csmp::Index phi_key    = model.Database().StorageKey("porosity");
    const csmp::Index betar_key  = model.Database().StorageKey("rock compressibility");
    const csmp::Index rhr_key    = model.Database().StorageKey("rock density");
    const csmp::Index cpr_key    = model.Database().StorageKey("rock heat capacity");
    
    // fluid mass per pore volume, to not "lose" or "gain" energy in total heat capacity
    const csmp::Index mt_key     = model.Database().StorageKey("fluid density");
    
    // transient operands that shall be computed
    const csmp::Index gv_key     = model.Database().StorageKey("gravity vector");
    const csmp::Index gt_key     = model.Database().StorageKey("gravity term");
    const csmp::Index rhogt_key  = model.Database().StorageKey("gravity term times density");
    const csmp::Index mc_key     = model.Database().StorageKey("mass conductivity"); // (rho k) / mu
    const csmp::Index betat_key  = model.Database().StorageKey("total compressibility"); // at integration point
    const csmp::Index cpt_key    = model.Database().StorageKey("total heat capacity");
    
    // needed to get the "correct" T and p at integration points
    const csmp::Index T_key      = model.Database().StorageKey("temperature");
    const csmp::Index p_key      = model.Database().StorageKey("fluid pressure");
    
    VectorVariable<3U>  gv;
   
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        const double k     = (*it)->Read( k_key );   // permeability
        const double phi   = (*it)->Read( phi_key ); // porosity
        const double betar = (*it)->Read( betar_key );  // rock compressibility
        const double rhr   = (*it)->Read( rhr_key ); // rock density
        const double cpr   = (*it)->Read( cpr_key ); // rock heat capacity
        
        for ( uint32_t i=0U; i<(*it)->IntegrationPoints(); ++i ) {
          (*it)->Read( gv_key, gv );
                    
          // Compute fluid properties from T and p at integration point
          const double T     = (*it)->PropertyValueAtIntegrationPoint( T_key, i );
          const double p     = (*it)->PropertyValueAtIntegrationPoint( p_key, i );
          
          const double rho   = water.Density( T, p );
          const double mu    = water.Viscosity( T, p );
                              
          // computing mass conductivity at element integration point
          (*it)->Store( i, mc_key, makeScalar(PLAIN, (k * rho) / mu ) );
           
          // computing the gravity term for the righthand side
          gv *= (k / mu) * rho;
          (*it)->Store( i, gt_key, gv );
          (*it)->Store( i, rhogt_key, gv*rho );
           
          // computing the total systems compressibility term
          const double betaf = water.Compressibility( T, p );
          const double betat = rho*( phi * betaf + (1. - phi) * betar ); //TD: rho or mt?
          (*it)->Store( i, betat_key, makeScalar(PLAIN, betat) );
           
          // computing the total heat capacity term
          ScalarVariable mt; //TD, it can be debated if rho or mt is to be taken, if taking rho instead -> compute from T&P
          (*it)->PropertyValueAtIntegrationPoint( mt_key, i, mt );
          const double cpf   = water.HeatCapacity( T, p );
          const double cpt = phi * mt() * cpf + (1. - phi) * cpr *  rhr;
          (*it)->Store( i, cpt_key, makeScalar(PLAIN, cpt) );
        }
      }
  } // end transientOperands








  /**
     Computes nodal heat source/sink by equilibrating the transported fluid mass/enthalpy
     by means of fluid-rock equilibration.
    
     1. the transported heat is distributed over the FV sectors
    
     2. the fluid properties are recomputed at each sector integration point
    
     3. the heating / cooling is computed sector by sector
    
     4.
    
     X. Also computes new finite volume "enthalpy content liquid" from sector contents
 
     Following a transport step, the total enthalpy on the node
     needs to be distributed between rock and fluid
     such that these are in thermal equilibrium.
    
     @attention  this method is concerned only with the computation of specific
     rather than volume-integrated quantities.
  */
  void ThermalConvectionETHZ_Example::thermalEquilibrationAtSectorIntegrationPoints( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
  {
    const csmp::Index pv_key         = model.Database().StorageKey("FV pore volume");
    const csmp::Index T_key          = model.Database().StorageKey("temperature");
    const csmp::Index FVT_key        = model.Database().StorageKey("FV temperature");
    // const csmp::Index ST_key         = model.Database().StorageKey("sector temperature");
    const csmp::Index pf_key         = model.Database().StorageKey("fluid pressure");
    const csmp::Index cp_fluid_key   = model.Database().StorageKey("fluid heat capacity");
    const csmp::Index cp_rock_key    = model.Database().StorageKey("rock heat capacity");
    const csmp::Index rho_rock_key   = model.Database().StorageKey("rock density");
    const csmp::Index hCl_key        = model.Database().StorageKey("enthalpy content liquid");
    const csmp::Index hVl_key        = model.Database().StorageKey("volumetric enthalpy liquid");
    const csmp::Index hVSl_key       = model.Database().StorageKey("sector volumetric enthalpy liquid");
    const csmp::Index hCSl_key       = model.Database().StorageKey("sector enthalpy content liquid");
    const csmp::Index hCSlp_key      = model.Database().StorageKey("previous sector enthalpy content liquid");
    const csmp::Index CPTS_key       = model.Database().StorageKey("total sector heat capacity");
    const csmp::Index sv_key         = model.Database().StorageKey("sector volume");
    const csmp::Index sW_key         = model.Database().StorageKey("sector weight");
    const csmp::Index smW_key        = model.Database().StorageKey("sector mass weight");
    const csmp::Index sm_key         = model.Database().StorageKey("sector mass");
    // mt = fluid mass / pore volume
    const csmp::Index mt_key         = model.Database().StorageKey("fluid density"); // at node
    const csmp::Index rln_key        = model.Database().StorageKey("density liquid");
    const csmp::Index rl_key         = model.Database().StorageKey("sector density liquid");
    const csmp::Index phi_key        = model.Database().StorageKey("porosity");
    const csmp::Index dhs_key        = model.Database().StorageKey("sector enthalpy change");
   
    // looping over FV sectors determining the amount by which their heat content has changed
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        // element-related properties
        const double phi = (*it)->Read( phi_key );
        const double cpr = (*it)->Read( cp_rock_key );
        const double rhr = (*it)->Read( rho_rock_key );
        
        // sector by sector calculations
        for ( uint32_t i=0U; i<(*it)->Sectors(); ++i )
          {
            // 1. interpolating old variable values p, T to sector integration points
            // ----------------------------------------------------------------------
            const double p = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, pf_key );
            double       T = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, T_key );
             
            // 2. recalculating heat capacitiy at the sector integration point
            // ---------------------------------------------------------------
            // NB: water heat capacitity is not exact, but iteration is neglected
            const double cp_fluid   = water.HeatCapacity( T, p );
            const double mt         = (*it)->N(i)->Read( mt_key );
            const double CPT_sector = (1. - phi) * cpr * rhr + phi * cp_fluid * mt; //TD check consistency
            //                                                                    ^^
            (*it)->Store( i, 0U, CPTS_key, makeScalar(PLAIN,CPT_sector) ); // for later computation of FV temperature
            
            // 3. computing (specific) fluid enthalpy and sector temperature change and new sector temperature
            // -----------------------------------------------------------------------------------------------
            double hCSl  = (*it)->N(i)->Read( hCl_key );    //TD: hCSl is never being updated during transport, rather need to read hCl from FV!
            double hCSlp = (*it)->Read( i, 0U, hCSlp_key ); // previous enthalpy content from sector
            
            // specific 'sector enthalpy change'
            const double dhCl = hCSl - hCSlp;
            (*it)->Store( i, 0U, dhs_key, makeScalar(PLAIN,dhCl) );
            // simple version of T-change without phase changes and assuming linearity
            const double dT = phi*dhCl / CPT_sector;
            T += dT;
            // (*it)->Store( i, 0U, ST_key, makeScalar(PLAIN,T) ); //TD: unnecessary?
            // 4. recomputing fluid properties at old p, new T
            // -----------------------------------------------
            // new 'sector density liquid'
            const double rl = water.Density( T, p );
            (*it)->Store( i, 0U, rl_key, makeScalar( PLAIN, rl ) );
            // new 'sector enthalpy content liquid'
            const double hl = water.Enthalpy( T, p );
            // new sector enthalpy content liquid - accumulated lhs
            hCSl  = hl * mt;
            // 'previous sector enthalpy content liquid' previous heat content at sector ip
            (*it)->Store( i, 0U, hCSlp_key, makeScalar(PLAIN,hCSl) );
            (*it)->Store( i, 0U, hVSl_key,  makeScalar(PLAIN,rl*hl) );
            
            const double sector_mass = (phi*mt + (1.0-phi)*rhr);
            (*it)->Store( i, 0U, sm_key,  makeScalar(PLAIN,sector_mass) );
          }
      }
 
    // Looping over sectors to determine their 'sector mass weight'
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        double  FV_mass(0.);
          // determing FV mass
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
            const uint32_t sector = (*nit)->ParentNodeNumber( i );
            const Element<3U>* const eptr((*nit)->Parent(i));
            const double sector_mass = eptr->Read( sector, 0U, sm_key )*eptr->Read( sector, 0U, sv_key );
            FV_mass += sector_mass;
          }
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
            const uint32_t sector = (*nit)->ParentNodeNumber( i );
            Element<3U>* const eptr((*nit)->Parent(i));
            const double sector_mass = eptr->Read( sector, 0U, sm_key )*eptr->Read( sector, 0U, sv_key );
            eptr->Store( sector, 0U, smW_key, makeScalar( PLAIN, sector_mass/FV_mass ) );
          }
      }
    
    // 5. recalculating the effective FV 'enthalpy content liquid' (J m-3) after temperature change for next advection step
    // --------------------------------------------------------------------------------------------------------------------
    // Notice:
    // - hCl is the 'enthalpy content liquid [J m-3]' of the FV, stored at the node, the LHS property for heat advection
    // - here, we integrate how it has changed after individual thermal equilibration in the sector; because water gained 
    //   or lost enthalpy by thermally equilibrating with rock
    
    /* Forget this for the moment; in the  logic that we currently have, we need to update T at the nodes, in order
     to get as accurate interpolations as possible in transientOperands; I therefore introduce the variable "FV temperature"
    // - T at the node is NOT being updated here; like in the pressure equation, a nodal source term will be added to the
    //   transient_temperature FE algorithm to compute a new temperature, similar to how the p-equation is solved; 
    //   I wonder if this is the best way, we could also thermally equilibrate the FV after summing to get the "exact" T ... 
    */
     
     for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        if((*nit)->Status(T_key) != DIRICH )
        {
          (*nit)->Store( hCl_key, makeScalar(PLAIN,0.) );
          (*nit)->Store( hVl_key, makeScalar(PLAIN,0.) );
          
          double  hCl(0.);
          double  hVl(0.);
          
          // drawing together the contributions of the sectors that make up the finite volume
          // the sector weight is sector pore volume divided by finite volume pore volume
          for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
            const uint32_t sector = (*nit)->ParentNodeNumber( i );
            const Element<3U>* const eptr((*nit)->Parent(i));
            const double weight = eptr->Read( sector, 0U, sW_key );
            hCl += weight * eptr->Read( sector, 0U, hCSlp_key );
            hVl += weight * eptr->Read( sector, 0U, hVSl_key );
            //T   += weight * eptr->Read( sector, 0U, ST_key );
          }
        // 'enthalpy content liquid' = specific rather than volume-integrated quantity
          (*nit)->Store( hCl_key, makeScalar(PLAIN,hCl) );
          (*nit)->Store( hVl_key, makeScalar(PLAIN,hVl) );
        }
      }
    // 6. computing, rhs (effective'density liquid' in the FV), from integration point values of 'sector density liquid'
    // ----------------------------------------------------------------------------------------------------------------
    // (= RHS 'density liquid' FV equation)
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        (*nit)->Store( rln_key, makeScalar(PLAIN,0.) );
        double  rhol_node(0.);
        // looping over the sectors of the current FV
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
          const uint32_t sector = (*nit)->ParentNodeNumber(i);
          const Element<3U>* const eptr((*nit)->Parent(i));
          // weighting the sector contributions taking their porosity into account
          const double weight = eptr->Read( sector, 0U, sW_key );
          rhol_node += weight * eptr->Read( sector, 0U, rl_key );
        }
        (*nit)->Store( rln_key, makeScalar(PLAIN,rhol_node) );
      }
    
    // 6. computing, fake "FV temperature"
    // ----------------------------------------------------------------------------------------------------------------
    // (= RHS 'density liquid' FV equation)
    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        (*nit)->Store( FVT_key, makeScalar(PLAIN,0.) );
        double FVT =(*nit)->Read( T_key );
        
        double  FV_CPT(0.);
        double  FV_dH(0.);
        // looping over the sectors of the current FV
        for ( uint32_t i=0U; i<(*nit)->Parents(); ++i ) {
          const uint32_t sector = (*nit)->ParentNodeNumber(i);
          const Element<3U>* const eptr((*nit)->Parent(i));
          // weighting the sector contributions taking their porosity into account
          const double mweight = eptr->Read( sector, 0U, smW_key );
          const double vweight = eptr->Read( sector, 0U, sW_key );
          FV_CPT += mweight * eptr->Read( sector, 0U, CPTS_key );
          FV_dH  += vweight * eptr->Read( sector, 0U, hCSlp_key ) * eptr->Read( phi_key );
        }
        (*nit)->Store( FVT_key, makeScalar(PLAIN, FVT + FV_CPT/FV_dH ) );
      }
    
  } // end thermalEquilibrationAtSectorIntegrationPoints







  /**
     computes nodal source-sink terms due to changes in the  density and energy content of the fluid 
     following (as are established during the thermal equilibration step)
    
     Effective terms at the nodes were already computed in thermal equilibrate so that these do not 
     need to be repeated here
 
     The 'nodal fluid volume source' is a total (integrated) point source, which is calculated from
     the difference in the density liquid (actual fluid density at onset of transport step)
     and the 'fluid density' calculated after the enthalpy transport step.
  */
  void ThermalConvectionETHZ_Example::nodalSourceOperands( Model<3U>& model, Region<3U>& ref, H2OLookup& water, const double& time_increment )
  {
    // for the nodal source term in the pressure equation
    // mt = fluid mass / pore volume
    const csmp::Index mt_key     = model.Database().StorageKey("fluid density");
    // actual thermodynamic density at current T,p
    const csmp::Index FVphi_key  = model.Database().StorageKey("FV pore volume");
    const csmp::Index sV_key     = model.Database().StorageKey("sector volume");
    const csmp::Index rho_key    = model.Database().StorageKey("density liquid");
    const csmp::Index dH_key     = model.Database().StorageKey("sector enthalpy change");
    const csmp::Index sW_key     = model.Database().StorageKey("sector weight");
    
    // output variables
    const csmp::Index QVn_key    = model.Database().StorageKey("nodal fluid volume source");
    const csmp::Index QHn_key    = model.Database().StorageKey("nodal heat source");
    const csmp::Index drho_key   = model.Database().StorageKey("density difference");

    for ( auto nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
        // 'nodal fluid volume source' term for pressure equation
        // ------------------------------------------------------
        const double mt    = (*nit)->Read( mt_key );
        // the 'density liquid' at the node, rl, was obtained from that at the sector integration points by
        // weighting the sector values by their contribution to the total pore volume
        const double rl    = (*nit)->Read( rho_key );
        const double FVphi = (*nit)->Read( FVphi_key );
        (*nit)->Store( QVn_key, makeScalar( PLAIN, FVphi * (mt - rl) ) );
        (*nit)->Store( drho_key, makeScalar( PLAIN, (mt - rl) ) );
        // zeroing the 'nodal heat source' variable so that it can be accumulated in next loop
        (*nit)->Store( QHn_key, makeScalar(PLAIN,0.) );
      }

    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      // 'nodal heat source' term for temperature equation
      // -------------------------------------------------
      for ( uint32_t i=0U; i<(*it)->Sectors(); ++i ) {
        // integration of specific term by multiplication with sector weight
        const double sQHn = (*it)->Read( i, 0U, dH_key ) * (*it)->Read( i, 0U, sW_key );//(*it)->Read( i, 0U, sV_key );
        double hcontribution = (*it)->N(i)->Read( QHn_key );
        hcontribution += sQHn * (*it)->N(i)->Read( FVphi_key );
        (*it)->N(i)->Store( QHn_key, makeScalar( PLAIN, hcontribution/*time_increment*/ ) );
      }
   
  } // end nodalSourceSinkOperand





  /**
     Computing the mass flux qm at the element integration points.
  */
  void ThermalConvectionETHZ_Example::massFlux( Model<3U>& model, Region<3U>& ref )
  {
    // required input variables
    const csmp::Index p_key   = model.Database().StorageKey("fluid pressure");
    const csmp::Index gt_key  = model.Database().StorageKey("gravity term");
    const csmp::Index mc_key  = model.Database().StorageKey("mass conductivity"); // (rho k) / mu
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    // output variables of this function
    const csmp::Index qm_key  = model.Database().StorageKey("mass flux"); // (rho k) / mu (grad p + rho g)
    const csmp::Index qme_key = model.Database().StorageKey("average mass flux"); // element
    const csmp::Index vD_key  = model.Database().StorageKey("Darcy velocity");
   
    VectorVariable<3U>  vflux, gflux, avg_flux;
    DenseMatrix<DM_MIN> DN;
   
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); ++it )
      {
        avg_flux = 0.;
        for ( uint32_t i=0U; i<(*it)->IntegrationPoints(); ++i )
          {
            // computing pressure gradient term
            // --------------------------------
            // collecting mass conductivity for calculation
            const double K = (*it)->Read( i, mc_key );
           
            vflux = 0.;
            (*it)->dN_AtIntegrationPoint( DN, i );
            for ( uint32_t j=0U; j<(*it)->Nodes(); ++j )
              for ( uint32_t k=0U; k<3U; ++k )
                vflux(k) += K * -DN(k,j) * (*it)->N(j)->Read( p_key );
         
            // gravity-related component of the flux
            // -------------------------------------
            (*it)->Read( i, gt_key, gflux );

            // adding pressure components together obtaining the mass flux term
            vflux += gflux;
         
            // storing the mass flux at the integration point
            (*it)->Store( i, qm_key, vflux );
           
            avg_flux += vflux;
          }
        
        // storing averaged value on the element
        avg_flux /= static_cast<double>((*it)->IntegrationPoints());
        (*it)->Store( qme_key, avg_flux );
        
        // computing the Darcy velocity using an average density at the element barycenter
        ScalarVariable rho_l;
        (*it)->PropertyValueAtBaryCenter( rho_key, rho_l );
        (*it)->Store( vD_key, avg_flux / rho_l() );
      }

  } // end massFlux









  void ThermalConvectionETHZ_Example::outputToVTK( const Model<3U>& model, VTK_Interface<3U>& vtk_output, long time_step )
  {
    vtk_output.OutputDataToVTK( model,
                                "temperature",
                                "temperature",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "fluid_pressure",
                                "fluid pressure",
                                time_step);
     
    vtk_output.OutputDataToVTK( model,
                                "mass_flux",
                                "mass flux",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "fluid_density",
                                "fluid density",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "density_liquid",
                                "density liquid",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "permeability",
                                "permeability",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "fluid_compressibility",
                                "fluid compressibility",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "nodal_fluid_volume_source",
                                "nodal fluid volume source",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "nodal_heat_source",
                                "nodal heat source",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "gravity_term",
                                "gravity term",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "enthalpy_liquid",
                                "enthalpy liquid",
                                time_step);

    vtk_output.OutputDataToVTK( model,
                                "test_variable",
                                "test variable",
                                time_step);
  
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "temperature" );
    printRangeOfVariable( model, "density liquid" );
    printRangeOfVariable( model, "fluid density" );
    printRangeOfVariable( model, "fluid viscosity" );
    printRangeOfVariable( model, "thermal conductivity");
    printRangeOfVariable( model, "total heat capacity");
    printRangeOfVariable( model, "total compressibility");
    printRangeOfVariable( model, "mass conductivity");
  }

  
void ThermalConvectionETHZ_Example::monitorVariableRange( Model<3U>& model )
  {
    printRangeOfVariable( model, "temperature" );
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "density liquid" );
    printRangeOfVariable( model, "fluid density" );
    printRangeOfVariable( model, "fluid viscosity" );
    printRangeOfVariable( model, "fluid compressibility" );
    printRangeOfVariable( model, "mass conductivity");
    printRangeOfVariable( model, "gravity term times density");    
    printRangeOfVariable( model, "FV temperature" );
    printRangeOfVariable( model, "FV pore volume" );
    printRangeOfVariable( model, "finite volume" );
    printRangeOfVariable( model, "nodal fluid volume source" );
    printRangeOfVariable( model, "density difference" );
    printRangeOfVariable( model, "Darcy velocity" );
    printRangeOfVariable( model, "total compressibility" );
    printRangeOfVariable( model, "sector enthalpy change" );
  }
} // end csmp
  

  
  
  
  
  

