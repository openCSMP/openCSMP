// ------------------------------------------
// Created by James Patterson on 7.5.2015
// ------------------------------------------

#include "Model.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "VTU_Interface.h"
#include "VTK_Interface.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"

#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "ModelTime.h"

// finite volumes 
#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"
#include "NodeCenteredFiniteVolumeTransport.h"

// finite elements
#include "PDE_Integrator.h"
#include "PDE_IntegratorExperimental.h"
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

#include "SAMG_Settings.h"

namespace csmp {

/// sector volume, finite volume, FV pore volume
void initializeFiniteVolumeProperties( Model<3U>& );

void initialMassAndEnthalpyContent( Model<3U>&, Region<3U>&, H2OLookup& );

void fluidProperties( Model<3U>&, Region<3U>&, H2OLookup& );

/// mobility, gravity term
void steadyStateOperands( Model<3U>&, Region<3U>& );

/// system compressibility etc.
void transientOperands( Model<3U>&, Region<3U>& );

/// calculation of new density liquid and enthalpy in FV sectors
void thermalEquilibrationAtSectorIntegrationPoints( Model<3U>&, Region<3U>&, H2OLookup& );

/// nodal source-sink terms due to density amd energy content changes upon thermal equilibration
void nodalSourceOperands( Model<3U>&, Region<3U>&, H2OLookup& );

/// qm at integration point
void massFlux( Model<3U>&, Region<3U>& );

void outputToVTK( const Model<3U>&, VTK_Interface<3U>&, long time_step );

}

using namespace std;
using namespace csmp;

/*
    SKM QUESTIONS AND COMMENTS
 
    - cleanup of names in variable file may be necessary: in my mind, i associate 'content' with an absolute
      = integrated property, but this does not seem to be the case
      
    - 'nodal fluid volume source'  should perhaps be 'nodal fluid mass source' ?
*/

//=================================================================================================
int main()
//=================================================================================================
{
  const  size_t DIM(3U);
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
  double64& model_time( ModelTime::Instance().modelTime );

  if (!restart) {
        model = new ANSYS_Model3D(geometry_name.c_str(), 
                                  regions_name.c_str(), 
                                  vars_name.c_str() , 
                                  true, true, true);
    }
  else model = new ANSYS_Model3D(restart_name);
  
  long64 time_step(0);
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
    SAMG_Settings settings;
     {
        PDE_IntegratorExperimental<DIM,Region>  heat_conductor(new SAMG_Solver(&settings));
    
        NumIntegral_dNT_op_dN_dV<DIM> conductance( model->Database(), "thermal conductivity", "temperature", "temperature");
        NumIntegral_NT_op_N_dV<DIM>  heat_source( model->Database(), "heat source", "temperature" );
        NumIntegral_NT_op_N_dS<DIM,Face>  basal_hfu( model->Database(), "basal heat flow", "temperature" );
       
        heat_conductor.Add( &conductance );
        heat_conductor.Add( &heat_source );
        heat_conductor.AddBoundaryIntegrals( &basal_hfu );
       
        // solve down to round-off
        settings.Set_eps(0.); // absolute criterion
        
        heat_conductor.IntegrateOver( *model, computation_domain );
        
        printRangeOfVariable( *model, "temperature");

        //! computation of initial steady-state fluid pressure
        PDE_Integrator<DIM,Region>  steady_state_pressure( new SAMG_Solver(&settings) );
        NumIntegral_dNT_op_dN_dV<DIM>  p_conductance( model->Database(), "mass conductivity", "fluid pressure", "fluid pressure" );
        NumIntegral_dNT_op_dV<DIM>     gravity( model->Database(), "gravity term", "fluid pressure" );

        steady_state_pressure.Add( &p_conductance );                           
        steady_state_pressure.Add(&gravity);
     
         // iout is reduced to bare minimum
        settings.Set_iout1( 0 );
        settings.Set_iout2( 0 );
     
        // initialize gravity vector
        VectorVariable<DIM>  grav_vec(DIRICH,0.);
        grav_vec(1) = -9.80665;

        // TODO: this assignment is only correct if there are no lower dimensional elements; add l.d. in plane or line vectors
        model->InputPropertyValue("gravity vector", grav_vec );
      
        // Iterate to find an initial vertical pressure profile
        // --------------------------------------------------------
        // initial pressure must be initialized in configuration file
        for ( int i=0; i<3; i++ ) {
             cout <<"\n\n\nmain: pressure initialisation loop; iteration: "<< i+1;
             fluidProperties( *model, computation_domain, water );
             steadyStateOperands( *model, computation_domain );
             steady_state_pressure.IntegrateOver( computation_domain );
             printRangeOfVariable( *model, "fluid pressure" );
          }
        massFlux( *model, computation_domain );
        printRangeOfVariable( *model, "mass flux" );
        printRangeOfVariable( *model, "average mass flux" );
     }
  
  
    //! -------------------------------------------------------------
    //! 3. Algorithms for transient pressure and temperature
    //! -------------------------------------------------------------
    // fluid pressure
    SAMG_Settings pressure_settings;
    pressure_settings.SetSolverInstance(1);
    PDE_Integrator<DIM,Region>     transient_pressure( new SAMG_Solver(&pressure_settings) );
    NumIntegral_dNT_op_dN_dV<DIM>  p_conductance( model->Database(), "mass conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_dNT_op_dV<DIM>     gravity( model->Database(), "gravity term", "fluid pressure" );
    NumIntegral_NT_op_N_dV<DIM>    element_fluid_source( model->Database(), "fluid volume source", "fluid pressure" );
                                   element_fluid_source.AddAccumulateLater();
    // to account for absolute nodal fluid contributions due to PVT property effects
    PointSource_rhsop<DIM>         node_total_source( model->Database(), "nodal fluid volume source", "fluid pressure" );
                                   node_total_source.AddAccumulateLater();
    NumIntegral_NT_op_N_dV<DIM>    storage_rhs( model->Database(), "total compressibility", "fluid pressure" );
                                   storage_rhs.LumpedFormulation();
                                   storage_rhs.MultiplyWithTimeIncrement();
    NumIntegral_NT_lhsop_N_dV<DIM> storage_lhs( model->Database(), "total compressibility", "fluid pressure", "fluid pressure" );
                                   storage_lhs.LumpedFormulation();
                                   storage_lhs.MultiplyWithTimeIncrement();

    transient_pressure.Add( &p_conductance );
    transient_pressure.Add( &gravity);
    transient_pressure.Add( &element_fluid_source);
    transient_pressure.Add( &node_total_source);
    transient_pressure.Add( &storage_lhs );
    transient_pressure.Add( &storage_rhs );

    // temperature
    SAMG_Settings temperature_settings;
    temperature_settings.SetSolverInstance(2);
    PDE_IntegratorExperimental<DIM,Region>  transient_temperature( new SAMG_Solver(&temperature_settings) );
    NumIntegral_dNT_op_dN_dV<DIM>  T_conductance( model->Database(), "thermal conductivity", "temperature", "temperature" );
    NumIntegral_NT_op_N_dV<DIM>    heat_source( model->Database(), "heat source", "temperature" );
                                   heat_source.AddAccumulateLater();

    PointSource_rhsop<DIM>         node_heat_source( model->Database(), "nodal heat source", "temperature" );
                                   node_heat_source.AddAccumulateLater();

    NumIntegral_NT_op_N_dV<DIM>    thermal_capacitance_rhs( model->Database(), "total heat capacity", "temperature" );
                                   thermal_capacitance_rhs.LumpedFormulation();
                                   thermal_capacitance_rhs.MultiplyWithTimeIncrement();
    NumIntegral_NT_lhsop_N_dV<DIM> thermal_capacitance_lhs( model->Database(), "total heat capacity", "temperature", "temperature" );
                                   thermal_capacitance_lhs.LumpedFormulation();
                                   thermal_capacitance_lhs.MultiplyWithTimeIncrement();

    NumIntegral_NT_op_N_dS<DIM,Face>  basal_hfu( model->Database(), "basal heat flow", "temperature" );
                                      basal_hfu.AddAccumulateLater();

    transient_temperature.Add( &T_conductance );
    transient_temperature.Add( &heat_source);
    transient_temperature.Add( &node_heat_source);
    transient_temperature.Add( &thermal_capacitance_rhs );
    transient_temperature.Add( &thermal_capacitance_lhs );
    transient_temperature.AddBoundaryIntegrals( &basal_hfu );


    // set-up explicit enthalpy transport scheme
    const bool second_order_accuracy(false);
    ExplicitMassBasedTransport<DIM,MassBasedStencilProcessor>
                                     mass_advector( "Model", *model,
                                               "porosity",
                                               "fluid density", "density liquid",
                                               "Darcy velocity",
                                               "zero node",
                                               second_order_accuracy );

    ExplicitMassBasedTransport<DIM,MassBasedStencilProcessor>
                                     heat_advector( "Model", *model,
                                               "porosity",
                                               "enthalpy content liquid", "volumetric enthalpy liquid",
                                               "Darcy velocity",
                                               "zero node",
                                               second_order_accuracy );

    // initialising "finite volume"  and "effective FV porosity" and variables for transient calculation
    initializeFiniteVolumeProperties( *model );
    printRangeOfVariable( *model, "finite volume" );
    printRangeOfVariable( *model, "FV pore volume" );
    printRangeOfVariable( *model, "sector volume" );
    printRangeOfVariable( *model, "sector weight" );
    // verifying model volume
    computation_domain.InputPropertyValue( "test variable", makeScalar(ANY,1.));
    cerr <<"\nmain: model volume: "<< computation_domain.VolumeIntegral( "test variable", false, false ) <<" m3.\n";

    // nodal fluid density and enthalpy content and 'previous sector enthalpy content'
    initialMassAndEnthalpyContent( *model, computation_domain, water );
    transientOperands( *model, computation_domain );


    //! -------------------------------------------------------------
    //! 4. Time evolution loop
    //! -------------------------------------------------------------
    printRangeOfVariable( *model, "Darcy velocity" );
    double64 simulated_time(0.), time_increment(mass_advector.AnisotropicCourantIncrement());
  
    while ( simulated_time < run_settings.Duration() )
      {
         cout <<"\nmain: time-stepping loop; step: "<< time_step << endl;
         // 1. transient temperature (starting with an initialized model with updated fluid properties)
         transient_temperature.TimeIncrement( 1. / time_increment );
         transient_temperature.IntegrateOver( computation_domain );
         printRangeOfVariable( *model, "temperature" );
        
         // 2. advect mass and enthalpy
         cout <<"\ninput variables for transport scheme:";
         printRangeOfVariable( *model, "Darcy velocity" );
         printRangeOfVariable( *model, "fluid density" );
         mass_advector.AdvectVariable( time_increment );
         heat_advector.AdvectVariable( time_increment );
        
         // 3. thermal equilibration
         thermalEquilibrationAtSectorIntegrationPoints( *model, computation_domain, water );
         printRangeOfVariable( *model, "sector enthalpy change" );
         printRangeOfVariable( *model, "temperature" );
        
         // 4. calculation of source/sink coupling term with the pressure equation
         nodalSourceOperands( *model, computation_domain, water );
         printRangeOfVariable( *model, "nodal fluid volume source" );
         printRangeOfVariable( *model, "nodal heat source" );
        
         // 5. solve transient pressure equation
         cout <<"\ninput variables for transient pressure computation:";
         printRangeOfVariable( *model, "nodal fluid volume source" );
         printRangeOfVariable( *model, "total compressibility" );
         printRangeOfVariable( *model, "density liquid" );
         transient_pressure.TimeIncrement( 1. / time_increment );
         transient_pressure.IntegrateOver( computation_domain );
         printRangeOfVariable( *model, "fluid pressure" );
         fluidProperties( *model, computation_domain, water );
         transientOperands( *model, computation_domain );
         massFlux( *model, computation_domain );
        
          // 6. output results at selected time intervals
          outputToVTK( *model, vtk_output, ++time_step );
      }
  
  delete model;
  cerr << "\n\nall done.\n\n";
  return 0;
  
//=================================================================================================
} // end main
//=================================================================================================





namespace csmp {


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
         const double64 phi = (*it)->Read( phi_key );
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
void initializeFiniteVolumeProperties( Model<3U>& model, Region<dim> ref )
 {
    const csmp::Index phi_key = model.Database().StorageKey("porosity");

    const csmp::Index fv_key    = model.Database().StorageKey("finite volume");
    const csmp::Index fvphi_key = model.Database().StorageKey("FV pore volume");
    const csmp::Index sv_key    = model.Database().StorageKey("sector volume");
    const csmp::Index swt_key   = model.Database().StorageKey("sector weight");
   
    ref.InputPropertyValue( "finite volume", makeScalar( ANY, 0. ) );
    ref.InputPropertyValue( "FV pore volume", makeScalar( ANY, 0. ) );

    // initializing sector und FV volumes including pore volumes
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         const double64 phi = (*it)->Read( phi_key );
        
         for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
              // sector volume
              const double64 sec_vol = (*it)->SectorVolume(i);
              (*it)->Store( i, 0U, sv_key, makeScalar(PLAIN,sec_vol) );
              // accumulating the FV volume
              double64 FV_vol = (*it)->N(i)->Read( fv_key );
              FV_vol += sec_vol;
              (*it)->N(i)->Store( fv_key, makeScalar(PLAIN,FV_vol) );
              // accumulating FV pore volume
              double64 FVpore_vol = (*it)->N(i)->Read( fvphi_key );
              FVpore_vol += sec_vol * phi;
              (*it)->N(i)->Store( fvphi_key, makeScalar(PLAIN,FVpore_vol) );
           }
      }

    // initializing sector weights: (sector_V * phi) / FV pore volume
    for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         const double64 pore_volume( (*nit)->Read( fvphi_key ) );
         for ( size_t i=0U; i<(*nit)->Parents(); ++i ) {
              Element<3U>* const eptr((*nit)->Parent(i));
              const size_t sector = (*nit)->ParentNodeNumber(i);
              const double64 weight = (eptr->Read( sector, 0U, sv_key ) * eptr->Read( phi_key )) / pore_volume;
//              const double64 weight = (eptr->SectorVolume(sector) * eptr->Read( phi_key )) / pore_volume;
              eptr->Store( sector, 0U, swt_key, makeScalar(PLAIN,weight) );
           }
      }

    cout <<"\ninitializeFiniteVolumeProperties: checking 'finite volume' and pore volume computed from sector volumes and porosity\n";
    double total_volume(0.), total_pore_volume(0.);
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         const double64 phi = (*it)->Read( phi_key );
         double evolume(0.);
         for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
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
     double64 min_weight(1e30), max_weight(-1e30);
     for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         double64 total_weight(0.);
         for ( size_t i=0U; i<(*nit)->Parents(); ++i ) {
              const Element<3U>* eptr((*nit)->Parent(i));
              const size_t sector = (*nit)->ParentNodeNumber(i);
              const double64 weight = eptr->Read( sector, 0U, swt_key );
              total_weight += weight;
           }
         min_weight = std::min( min_weight, total_weight );
         max_weight = std::max( max_weight, total_weight );
      }
    cout <<"\n\tmin vs. max of the sector weights summed up over the finite volumes: "<< min_weight <<" vs. "<< max_weight << endl;
   
    cout <<"\nchecking writing of a sector variable (range should be 1).\n";
    const csmp::Index tv_key   = model.Database().StorageKey("test variable");
    ref.InputPropertyValue( "test variable", makeScalar( ANY, 0. ) );
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
           double64 test = (*it)->N(i)->Read( tv_key );
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
void initialMassAndEnthalpyContent( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
 {
    // at the finite-element nodes
    const csmp::Index p_key   = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key   = model.Database().StorageKey("temperature");
    const csmp::Index mt_key  = model.Database().StorageKey("fluid density");
    const csmp::Index hCl_key = model.Database().StorageKey("enthalpy content liquid");
    // at the finite element FV sector integration points
    const csmp::Index hCSlp_key = model.Database().StorageKey("previous sector enthalpy content liquid");
   
    for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         const double64 p  = (*nit)->Read( p_key );
         const double64 T  = (*nit)->Read( T_key );
         (*nit)->Store( mt_key,  makeScalar(PLAIN,water.Density( T, p )) );
         (*nit)->Store( hCl_key, makeScalar(PLAIN,water.Density( T, p ) * water.Enthalpy( T, p )) );
      }

    // FV sector (specific=not integrated) enthalpy content computed at the sector integration points
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
            const double64 p  = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, p_key );
            const double64 T  = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, T_key );
            const double64 h_content = water.Density( T, p ) * water.Enthalpy( T, p );
            (*it)->Store( i, 0U, hCSlp_key, makeScalar(PLAIN,h_content) );
         }
   
 } // end initialMassAndEnthalpyContent




/**
     A distinction is made between the actual fluid density (density liquid)
     and the current mass of fluid per pore-volume (in the transient calculation with advection).
*/
void fluidProperties( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
 {
    const csmp::Index p_key = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key = model.Database().StorageKey("temperature");

    // density at given p, T
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index mu_key  = model.Database().StorageKey("fluid viscosity");
    const csmp::Index cf_key  = model.Database().StorageKey("fluid compressibility");
    const csmp::Index cpf_key = model.Database().StorageKey("fluid heat capacity");
    const csmp::Index hl_key  = model.Database().StorageKey("enthalpy liquid");
    const csmp::Index hvl_key = model.Database().StorageKey("volumetric enthalpy liquid");
   
    for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         const double64 p   = (*nit)->Read( p_key );
         const double64 T   = (*nit)->Read( T_key );
//         const double64 rho(1000.); // = water.Density( T, p );
//         const double64 mu(1e-3); //  = water.Viscosity( T, p );
         const double64 rho = water.Density( T, p );
         const double64 mu  = water.Viscosity( T, p );
         const double64 cf  = water.Compressibility( T, p );
         const double64 hl  = water.Enthalpy( T, p );
         const double64 cpf = water.HeatCapacity( T, p );
         const double64 hvl = hl * rho;
        
         // storing the variables
         (*nit)->Store( rho_key, makeScalar(PLAIN,rho) );
         (*nit)->Store( mu_key, makeScalar(PLAIN,mu) );
         (*nit)->Store( cf_key, makeScalar(PLAIN,cf) );
         (*nit)->Store( hl_key, makeScalar(PLAIN,hl) );
         (*nit)->Store( cpf_key, makeScalar(PLAIN,cpf) );
         (*nit)->Store( hvl_key, makeScalar(PLAIN,hvl) );
      }
   
 } // end fluidProperties





/**
     mobility, gravity term at element integration points
     TODO: update for lower-dimensional elements
*/
void steadyStateOperands( Model<3U>& model, Region<3U>& ref )
 {
    const csmp::Index k_key = model.Database().StorageKey("permeability");
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index mu_key = model.Database().StorageKey("fluid viscosity");
    const csmp::Index gv_key = model.Database().StorageKey("gravity vector");

    const csmp::Index gt_key = model.Database().StorageKey("gravity term");
    const csmp::Index mc_key = model.Database().StorageKey("mass conductivity"); // (rho k) / mu
   
    VectorVariable<3U>  gv;
   
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         const double64 k = (*it)->Read( k_key );
        
         for ( size_t i=0; i<(*it)->IntegrationPoints(); ++i ) {
              (*it)->Read( gv_key, gv );
              ScalarVariable rho, mu;
              (*it)->PropertyValueAtIntegrationPoint( rho_key, i, rho );
              (*it)->PropertyValueAtIntegrationPoint( mu_key, i, mu );

              // computing mobility value at the integration point
              (*it)->Store( i, mc_key, makeScalar(PLAIN, (k * rho()) / mu() ) );
           
              // computing the gravity term for the righthand side
              gv *= (k / mu()) * rho() * rho();
              (*it)->Store( i, gt_key, gv );
           }
      }

 } // end steadyStateOperands



/**
     Computing mass conductivity, the gravity term, total systems compressibility, and total heat capacity
*/
void transientOperands( Model<3U>& model, Region<3U>& ref )
 {
    const csmp::Index k_key   = model.Database().StorageKey("permeability");
    const csmp::Index phi_key = model.Database().StorageKey("porosity");
    const csmp::Index mt_key  = model.Database().StorageKey("fluid density");
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index mu_key  = model.Database().StorageKey("fluid viscosity");
    const csmp::Index cf_key  = model.Database().StorageKey("fluid compressibility");
    const csmp::Index cr_key  = model.Database().StorageKey("rock compressibility");
    const csmp::Index cpf_key  = model.Database().StorageKey("fluid heat capacity");
    const csmp::Index rhr_key  = model.Database().StorageKey("rock density");
    const csmp::Index cpr_key  = model.Database().StorageKey("rock heat capacity");
    const csmp::Index gv_key  = model.Database().StorageKey("gravity vector");

    const csmp::Index gt_key  = model.Database().StorageKey("gravity term");
    const csmp::Index mc_key  = model.Database().StorageKey("mass conductivity"); // (rho k) / mu
    const csmp::Index ct_key  = model.Database().StorageKey("total compressibility"); // at integration point
    const csmp::Index tcp_key = model.Database().StorageKey("total heat capacity");
   
    VectorVariable<3U>  gv;
   
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         const double64 k   = (*it)->Read( k_key );   // permeability
         const double64 phi = (*it)->Read( phi_key ); // porosity
         const double64 cr  = (*it)->Read( cr_key );  // rock compressibility
         const double64 rhr = (*it)->Read( rhr_key ); // rock density
         const double64 cpr = (*it)->Read( cpr_key ); // rock heat capacity
        
         for ( size_t i=0; i<(*it)->IntegrationPoints(); ++i ) {
             (*it)->Read( gv_key, gv );
              ScalarVariable rho, mu, cf;
              (*it)->PropertyValueAtIntegrationPoint( rho_key, i, rho );
              (*it)->PropertyValueAtIntegrationPoint( mu_key, i, mu );
              (*it)->PropertyValueAtIntegrationPoint( cf_key, i, cf );
           
              // computing mass conductivity at element integration point
              (*it)->Store( i, mc_key, makeScalar(PLAIN, (k * rho()) / mu() ) );
           
              // computing the gravity term for the righthand side
              gv *= (k / mu()) * rho() * rho();
              (*it)->Store( i, gt_key, gv );
           
              // computing the total systems compressibility term
              const double64 ct = phi * cf() + (1. - phi) * cr;
              (*it)->Store( i, ct_key, makeScalar(PLAIN,ct) );
           
              // computing the total heat capacity term
              ScalarVariable mt, cpf;
              (*it)->PropertyValueAtIntegrationPoint( mt_key, i, mt );
              (*it)->PropertyValueAtIntegrationPoint( cpf_key, i, cpf );
              const double64 tcp = phi * mt() * cpf() + (1. - phi) * cpr *  rhr;
              (*it)->Store( i, tcp_key, makeScalar(PLAIN,tcp) );
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
void thermalEquilibrationAtSectorIntegrationPoints( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
 {
    const csmp::Index pv_key         = model.Database().StorageKey("FV pore volume");
    const csmp::Index T_key          = model.Database().StorageKey("temperature");
    const csmp::Index pf_key         = model.Database().StorageKey("fluid pressure");
    const csmp::Index cp_fluid_key   = model.Database().StorageKey("fluid heat capacity");
    const csmp::Index cp_rock_key    = model.Database().StorageKey("rock heat capacity");
    const csmp::Index rho_rock_key   = model.Database().StorageKey("rock density");
    const csmp::Index hCl_key        = model.Database().StorageKey("enthalpy content liquid");
    const csmp::Index hCSl_key       = model.Database().StorageKey("sector enthalpy content liquid");
    const csmp::Index hCSlp_key      = model.Database().StorageKey("previous sector enthalpy content liquid");
    const csmp::Index sW_key         = model.Database().StorageKey("sector weight");
    // mt = fluid mass / pore volume
    const csmp::Index mt_key         = model.Database().StorageKey("fluid density"); // at node
    const csmp::Index rln_key        = model.Database().StorageKey("density liquid");
    const csmp::Index rl_key         = model.Database().StorageKey("sector density liquid");
    const csmp::Index phi_key        = model.Database().StorageKey("porosity");
    const csmp::Index dhs_key        = model.Database().StorageKey("sector enthalpy change");
   
    // looping over FV sectors determining the amount by which their heat content has changed
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         // element-related properties
         const double64 phi = (*it)->Read( phi_key );
         const double64 cpr = (*it)->Read( cp_rock_key );
         const double64 rhr = (*it)->Read( rho_rock_key );
        
         // sector by sector calculations
         for ( size_t i=0U; i<(*it)->Sectors(); ++i )
           {
              // 1. interpolating old variable values p, T to sector integration points
              // ----------------------------------------------------------------------
              const double64 p = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, pf_key );
              double64       T = (*it)->PropertyValueAtSectorIntegrationPoint( i, 0U, T_key );
             
              // 2. recalculating heat capacitiy at the sector integration point
              // ---------------------------------------------------------------
              // NB: water heat capacitity is not exact, but iteration is neglected
              const double64 cp_fluid = water.HeatCapacity( T, p );
              const double64 mt = (*it)->N(i)->Read( mt_key );
              const double64 CPT_sector = (1. - phi) * cpr * rhr + phi * cp_fluid * mt;
              //                                                                    ^^
            
              // 3. computing (specific) fluid enthalpy and sector temperature change and new sector temperature
              // -----------------------------------------------------------------------------------------------
              double64 hCSl  = (*it)->Read( i, 0U, hCSl_key );  // from FV sector
              double64 hCSlp = (*it)->Read( i, 0U, hCSlp_key ); // previous enthalpy content from sector
              // (computing specific quantity)
              const double64 dhCl = phi * (hCSl - hCSlp);
              // specific 'sector enthalpy change'
              (*it)->Store( i, 0U, dhs_key, makeScalar(PLAIN,dhCl) );
              // simple version of T-change without phase changes and assuming linearity
              const double64 dT = dhCl / CPT_sector;
              T += dT;

              // 4. recomputing fluid properties at old p, new T
              // -----------------------------------------------
              // new 'sector density liquid'
              const double64 rl = water.Density(T,p);
              (*it)->Store( i, 0U, rl_key, makeScalar(PLAIN,rl) );
              // new 'sector enthalpy content liquid'
              const double64 hl = water.Enthalpy(T,p);
              // new enthalpy content liquid - accumulated lhs
              hCSl  = hl * mt;
              // 'previous sector enthalpy content liquid' previous heat content at sector ip
              (*it)->Store( i, 0U, hCSlp_key, makeScalar(PLAIN,hCSl) );
          }
      }
 
    // 5. recalculating the effective FV 'enthalpy content liquid' (J m-3) after temperature change for next advection step
    // --------------------------------------------------------------------------------------------------------------------
    // (=LHS enthalpy property for FV equation)
    for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         (*nit)->Store( hCl_key, makeScalar(PLAIN,0.) );
         double64  hCl(0.);
         // drawing together the contributions of the sectors that make up the finite volume
         for ( size_t i=0U; i<(*nit)->Parents(); ++i ) {
              const size_t sector = (*nit)->ParentNodeNumber( i );
              const Element<3U>* const eptr((*nit)->Parent(i));
              const double64 weight = eptr->Read( sector, 0U, sW_key );
              hCl += weight * eptr->Read( sector, 0U, hCSlp_key );
           }
         // 'enthalpy content liquid' = specific rather than volume-integrated quantity
         (*nit)->Store( hCl_key, makeScalar(PLAIN,hCl) );
      }
 
     // 6. computing, rh (effective'density liquid' in the FV), from integration point values of 'sector density liquid'
     // ----------------------------------------------------------------------------------------------------------------
     // (= RHS 'density liquid' FV equation)
     for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         (*nit)->Store( rln_key, makeScalar(PLAIN,0.) );
         double64  rhol_node(0.);
         // looping over the sectors of the current FV
         for ( size_t i=0U; i<(*nit)->Parents(); ++i ) {
              const size_t sector = (*nit)->ParentNodeNumber(i);
              const Element<3U>* const eptr((*nit)->Parent(i));
              // weighting the sector contributions taking their porosity into account
              const double64 weight = eptr->Read( sector, 0U, sW_key );
              rhol_node += weight * eptr->Read( sector, 0U, rl_key );
           }
         (*nit)->Store( rln_key, makeScalar(PLAIN,rhol_node) );
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
void nodalSourceOperands( Model<3U>& model, Region<3U>& ref, H2OLookup& water )
 {
    // for the nodal source term in the pressure equation
    // mt = fluid mass / pore volume
    const csmp::Index mt_key = model.Database().StorageKey("fluid density");
    // actual thermodynamic density at current T,p
    const csmp::Index FV_key  = model.Database().StorageKey("finite volume");
    const csmp::Index sV_key  = model.Database().StorageKey("sector volume");
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    const csmp::Index dH_key  = model.Database().StorageKey("sector enthalpy change");
    // output variables
    const csmp::Index QVn_key = model.Database().StorageKey("nodal fluid volume source");
    const csmp::Index QHn_key = model.Database().StorageKey("nodal heat source");

    for ( vector<Node<3U>*>::iterator nit=ref.NodesBegin(); nit!=ref.NodesEnd(); ++nit )
      {
         // 'nodal fluid volume source' term for pressure equation
         // ------------------------------------------------------
         const double64 mt = (*nit)->Read( mt_key );
         // the 'density liquid' at the node, rl, was obtained from that at the sector integration points by
         // weighting the sector values by their contribution to the total pore volume
         const double64 rl = (*nit)->Read( rho_key );
         const double64 FV = (*nit)->Read( FV_key );
         (*nit)->Store( QVn_key, makeScalar( PLAIN, FV * (mt - rl) ) );

         // zeroing the 'nodal heat source' variable so that it can be accumulated in next loop
         (*nit)->Store( QHn_key, makeScalar(PLAIN,0.) );
      }

    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      // 'nodal heat source' term for temperature equation
      // -------------------------------------------------
      for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
           // integration of specific term by multiplication with sector volume
           const double64 sQHn = (*it)->Read( i, 0U, dH_key ) * (*it)->Read( i, 0U, sV_key );
           double64 hcontribution = (*it)->N(i)->Read( QHn_key );
           hcontribution += sQHn;
           (*it)->N(i)->Store( QHn_key, makeScalar(PLAIN,hcontribution) );
        }
   
 } // end nodalSourceSinkOperand





/**
    Computing the mass flux qm at the element integration points.
*/
void massFlux( Model<3U>& model, Region<3U>& ref )
 {
    // required input variables
    const csmp::Index p_key = model.Database().StorageKey("fluid pressure");
    const csmp::Index gt_key = model.Database().StorageKey("gravity term");
    const csmp::Index mc_key = model.Database().StorageKey("mass conductivity"); // (rho k) / mu
    const csmp::Index rho_key = model.Database().StorageKey("density liquid");
    // output variables of this function
    const csmp::Index qm_key  = model.Database().StorageKey("mass flux"); // (rho k) / mu (grad p + rho g)
    const csmp::Index qme_key = model.Database().StorageKey("average mass flux"); // element
    const csmp::Index vD_key  = model.Database().StorageKey("Darcy velocity");
   
    VectorVariable<3U>  vflux, gflux, avg_flux;
    DenseMatrix<DM_MIN> DN;
   
    for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
      {
         avg_flux = 0.;
         for ( size_t i=0; i<(*it)->IntegrationPoints(); ++i )
           {
              // computing pressure gradient term
              // --------------------------------
              // collecting mass conductivity for calculation
              const double64 K = (*it)->Read( i, mc_key );
           
              vflux = 0.;
              (*it)->dN_AtIntegrationPoint( DN, i );
              for ( size_t j=0U; j<(*it)->Nodes(); ++j )
                for ( size_t k=0U; k<3U; ++k )
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
         avg_flux /= static_cast<double64>((*it)->IntegrationPoints());
         (*it)->Store( qme_key, avg_flux );
        
         // computing the Darcy velocity using an average density at the element barycenter
         ScalarVariable rho_l;
         (*it)->PropertyValueAtBaryCenter( rho_key, rho_l );
         (*it)->Store( vD_key, avg_flux / rho_l() );
      }

 } // end massFlux









void outputToVTK( const Model<3U>& model, VTK_Interface<3U>& vtk_output, long time_step )
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
                          "fluid_volume_source",
                          "fluid volume source",
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

} // end csmp
  

  
  
  
  
  

