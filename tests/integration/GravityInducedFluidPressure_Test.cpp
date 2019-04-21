//
//  GravityInducedFluidPressure_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 18/2/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#include "GravityInducedFluidPressure_Test.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "LineElementMesher.h"
#include "TextInterface.h"
#include "VTU_Interface.h"

#include "EOS_CO2H2ONaCl_Spycher2004.h"

#include "PDE_Integrator.h"
#include "PDE_Integrator_UoM.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "GaussJordan_Solver.h"

#include "ANSYS_Model3D.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"

using namespace std;

namespace csmp {

GravityInducedFluidPressure_Test::GravityInducedFluidPressure_Test()
 : top_(3000.), verbose_(true)
 {
    InitialiseModel1D( 3000. );
 }


GravityInducedFluidPressure_Test::~GravityInducedFluidPressure_Test()
 {
    delete model1D_;
 }





/**
    Builds a 2km-tall 1D (100-element) model
*/
void GravityInducedFluidPressure_Test::InitialiseModel1D( double64 model_height )
 {
    assert( model_height > 0. );
    VSet<1U>                mesh_container;
    const uint32            N_ELEMENTS(100);
    LineElementMesher<1U>   mesher;
    top_ = model_height; 
    mesher.BuildUniformMesh( mesh_container, model_height, N_ELEMENTS );
    delete model1D_;
    model1D_ = new Model<1U>( mesh_container, "GravityInducedFluidPressure_Test-variables.txt", true, false );
    Region<1U>& model_domain(model1D_->Region("Model"));
    model_domain.UpdateMemberIndexes();
    printModelDimensions( *model1D_ );
    // NOTE: Y points upward, else gravity will act in opposite direction
    const size_t topNode   = findNode( *model1D_, top_, 1.0e-7 );
    const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
    assert( topNode < N_ELEMENTS+1 );
    assert( bottomNode0m < 101 );
    cout <<"\nGravityInducedFluidPressure_Test: model end points:\n";
    cout <<"\t"<< model_domain.N(bottomNode0m)->Coordinate() << endl;
    cout <<"\t"<< model_domain.N(topNode)->Coordinate() << endl;

    // Input of material properties and initial conditions
    model1D_->InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-10) );
    model1D_->InputPropertyValue( "total mobility permeability product", makeScalar(PLAIN,1.0e-6) );
    model1D_->InputPropertyValue( "fluid pressure",        makeScalar(PLAIN,1.0e7) ); // 1 bar
    model1D_->InputPropertyValue( "thermal conductivity",  makeScalar(PLAIN,2.2) );
    model1D_->InputPropertyValue( "density carbonic phase",  makeScalar(PLAIN,200.) );
    model1D_->InputPropertyValue( "viscosity carbonic phase",  makeScalar(PLAIN,1.0e-5) );
    model1D_->InputPropertyValue( "acceleration gravity",  makeScalar(PLAIN,9.8061) );
   
 } // end InitialiseModel1D





/**
    Builds a 2km-tall 1D (100-element) model
*/
void GravityInducedFluidPressure_Test::InitialiseModel3D( const char* model )
 {
    // checking whether a binary version of the model already exists by looking for the prop database file
    ifstream ifs( string(model) + "_variables.dat" );
    const bool restart = (ifs.good()) ? true : false;
 
    if (!restart)
      model3D_ = new ANSYS_Model3D( model, model, "CO2-geo-sequestration-variables.txt", true, true, true );
    else
      model3D_ = new Model<3U>( string(model) );

    Region<1U>& model_domain(model1D_->Region("Model"));
    model_domain.UpdateMemberIndexes();
    printModelDimensions( *model1D_ );

    InputDataManager<3U>   model_configuration;
    ComputationalSettings  run_settings;

    if (!restart)
      model_configuration.ConfigureFromFile( *model3D_, model, false,    // groupname from parameter range
                                             true,    // default property values
                                             true,    // regional property values
                                             false,   // boundary conditions for box-shaped model
                                             true,    // essential conditions for groups
                                             true,    // csmp::Boundary properties
                                             run_settings );
   
 } // end InitialiseModel3D





/**
    From top temperature and a positive temperature gradient in oC/m.
*/
void GravityInducedFluidPressure_Test::InitialiseTemperatureProfile1D( double64 T_top, double64 grad_T )
 {
    assert( T_top  > 8. );
    assert( grad_T > 0. );
    const csmp::Index key_T(model1D_->Database().StorageKey("temperature"));
    Region<1U>& model_domain = model1D_->Region("Model");

    for ( vector<Node<1U>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 ToC = T_top + (top_ - (*nit)->x()) * grad_T;
         (*nit)->Store( key_T, makeScalar(PLAIN,ToC) );
      }
    printRangeOfVariable( *model1D_, "temperature" );

 } // end InitialiseTemperatureProfile1D



/**
    Finds the shallowest point in the model and uses that as a reference for the initialisation.
    Since the thermal gradient is linear and constant, the node coordinate Y values are sufficient to make the assignment.
*/
void GravityInducedFluidPressure_Test::InitialiseTemperatureProfile3D( double64 T_top, double64 grad_T )
 {
    assert( T_top  > 8. );
    assert( grad_T > 0. );
   
    // finds the shallowest (largest-Y point in the model)
    Point<3U> xyz_min, xyz_max;
    model3D_->MinMaxCoordinates( xyz_min, xyz_max );
    top_ = xyz_max[1];
   
    const csmp::Index key_T(model3D_->Database().StorageKey("temperature"));
    Region<3U>& model_domain = model3D_->Region("Model");

    for ( vector<Node<3U>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 ToC = T_top + (top_ - (*nit)->y()) * grad_T;
         (*nit)->Store( key_T, makeScalar(PLAIN,ToC) );
      }
    printRangeOfVariable( *model3D_, "temperature" );

 } // end InitialiseTemperatureProfile3D





/**
   As above, but for fluid pressure
*/
void GravityInducedFluidPressure_Test::InitialisePressure3D( double64 pf_top, double64 fluid_density )
 {
    assert( pf_top  >= patm_ );
    assert( fluid_density >= 1. );
   
    // finds the shallowest (largest-Y point in the model)
    Point<3U> xyz_min, xyz_max;
    model3D_->MinMaxCoordinates( xyz_min, xyz_max );
    top_ = xyz_max[1];
   
    // recovering acceleration of gravity
    const csmp::Index  g_key(model3D_->Database().StorageKey("acceleration gravity"));
    const double64     acc_gravity(model3D_->Read(g_key));

    const csmp::Index key_pf(model3D_->Database().StorageKey("fluid pressure"));
    Region<3U>& model_domain = model3D_->Region("Model");

    for ( vector<Node<3U>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 pf = pf_top + (top_ - (*nit)->y()) * fluid_density * acc_gravity;
         (*nit)->Store( key_pf, makeScalar(PLAIN,pf) );
      }
    printRangeOfVariable( *model3D_, "fluid pressure" );

 } // end InitialisePressure3D




/**
    TESTING IS DONE HERE
*/
void GravityInducedFluidPressure_Test::run()
 {
    // 1. testing pressure computation for a constant fluid density
    // ------------------------------------------------------------
    const double64 ref_density(1000.);
    ReferencePressureForFixedDensity( ref_density );
    ReferencePressureByTopDownIntegration( ref_density );
    _test( TestComputedWithReferencePressure() );
    if ( verbose_ ) OutputResultsToText( "GravityInducedFluidPressure_Test0" );


    // 2. testing computation for compressible CO2 (Spycher EOS)
    // ------------------------------------------------------------
    const double64 T_model_top(25.);
    const double64 grad_T_K_per_m(0.02); // 20oC/km
    InitialiseTemperatureProfile1D( T_model_top, grad_T_K_per_m );
    const double pf_model_top( 4.0e6 ); // hydrostatic at 400-m depth
    ReferencePressureByTopDownIntegrationCO2( pf_model_top );
    if ( verbose_ ) OutputResultsToText( "GravityInducedFluidPressure_Test1" );
    // testing
    ComputeCO2Pressure_PDE_Integrator( pf_model_top );
    ComputeCO2Pressure_PDE_Integrator2( pf_model_top );
    ComputeCO2PressureFromReducedPressure_PDE_Integrator2( model1D_, pf_model_top );
//    _test( TestComputedWithReferencePressure() );
    if ( verbose_ ) OutputResultsToText( "GravityInducedFluidPressure_Test2" );
   
   
    // 3. testing the same computation for a 3D model with surface topography
    //    and versions with difference element types
    // ---------------------------------------------
    const double64 fluid_density(1000.);
    const string   modelName("Greenshank_prism_mini"); // 2800 nodes
    InitialiseModel3D( modelName.c_str() );
    InitialiseTemperatureProfile3D( T_model_top, grad_T_K_per_m );
    InitialisePressure3D( pf_model_top,  fluid_density );
    model3D_->InputPropertyValue( "density carbonic phase", makeScalar(PLAIN,200.) );
    model3D_->InputPropertyValue( "fluid mixture density", makeScalar(PLAIN,200.) );
    ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan( model3D_, pf_model_top );
    // creating test dataset for comparison with other solver
    model3D_->CopyReplace( "fluid pressure", "previous fluid pressure" );
    if ( verbose_ ) {
         list<string> propertyNames({"temperature","fluid pressure", "reduced fluid pressure"});
         VTU_Interface<3U> vtu_output( *model3D_, "GravityInducedFluidPressure_Test" );
         vtu_output.OutputDataToVTU( modelName, propertyNames, "Model", 0 );
      }
   // compariing the Gauss-Jordan solver solution with that of SAMG
   ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG( model3D_, pf_model_top );
    _test( TestResultsByComparison( model3D_, "fluid pressure", "previous fluid pressure", ptol_ ) );

 } // end run
  
  
  
  
  
/**
    Integrating from the surface downward to the bottom of the model, the pressure is accumulated.
*/
void GravityInducedFluidPressure_Test::ReferencePressureByTopDownIntegration( double64 fluid_density )
  {
     assert( fluid_density > 1. );
     const csmp::Index key_pf(model1D_->Database().StorageKey("fluid pressure"));
     const csmp::Index g_key(model1D_->Database().StorageKey("acceleration gravity"));
     Region<1U>& model_domain(model1D_->Region("Model"));
     model_domain.UpdateMemberIndexes();
     const size_t topNode   = findNode( *model1D_, top_, 1.0e-7 );
     const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
     const double64 acc_gravity(model1D_->Read(g_key));

     // finding the top node
     Node<1U>* nptr = model_domain.N(topNode);
     Node<1U>* const bottom_node_ptr = model_domain.N(bottomNode0m);
     assert( nptr != nullptr );
     Element<1U>* eptr = nptr->Parent(0);
     assert( eptr != nullptr );
     // initialising first node
     nptr->Store( key_pf, makeScalar(DIRICH,patm_) );
 
     while ( nptr != bottom_node_ptr )
       {
          // computing the fluid properties at the next deeper node
          // ------------------------------------------------------
          Node<1U>* nptr2 = (eptr->N(0) == nptr) ? eptr->N(1) : eptr->N(0);
          assert( nptr2 != nullptr );
          assert( nptr2->x() < nptr->x() );
          // estimating pressure at deeper node using fluid density from current node
          double64 elength = eptr->Volume(); // =length in 1D
          // pressure calculation and storage
          double64 pf = nptr->Read(key_pf) + elength * fluid_density * acc_gravity;
         
          // recording the calculated fluid pressure value on the lower node
          // ---------------------------------------------------------------
          nptr2->Store( key_pf, makeScalar(DIRICH,pf) );
         
          // continue to next element
          if ( nptr2->Parents() > 1U ) {
               eptr = (nptr2->Parent(0) == eptr) ? nptr2->Parent(1) : nptr2->Parent(0);
            }
          else break;
          nptr = nptr2;
       }
 
     // get first idea about the results
     if ( verbose_ )
       assert( printRangeOfVariable( *model1D_, "fluid pressure" ) >= patm_ );

 } // end ReferencePressureByTopDownIntegration







/**
    Integrating from the surface downward to the bottom of the model, the pressure is accumulated.
    CO2 density is taken from the Spycher et al equation of state.
*/
void GravityInducedFluidPressure_Test::ReferencePressureByTopDownIntegrationCO2( double64 pf_top )
  {
     const csmp::Index key_T(model1D_->Database().StorageKey("temperature"));
     const csmp::Index key_pf(model1D_->Database().StorageKey("reference pressure"));
     const csmp::Index g_key(model1D_->Database().StorageKey("acceleration gravity"));
     const double64    acc_gravity(model1D_->Read(g_key));
     const csmp::Index key_rho(model1D_->Database().StorageKey("density carbonic phase"));
     const csmp::Index key_mu(model1D_->Database().StorageKey("viscosity carbonic phase"));
     const csmp::Index key_rhom(model1D_->Database().StorageKey("fluid mixture density"));

     Region<1U>& model_domain(model1D_->Region("Model"));
     model_domain.UpdateMemberIndexes();
     const size_t topNode   = findNode( *model1D_, top_, 1.0e-7 );
     const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
     // finding the top node
     Node<1U>* nptr = model_domain.N(topNode);
     Node<1U>* const bottom_node_ptr = model_domain.N(bottomNode0m);
     assert( nptr != nullptr );
     Element<1U>* eptr = nptr->Parent(0);
     assert( eptr != nullptr );
     // initialising EOS on first node
     nptr->Store( key_pf, makeScalar(DIRICH,pf_top) );
    
     // Spycher et al 2003 EOS
     EOS_CO2H2ONaCl_Spycher04  eos;
     nptr->Store( key_rho, makeScalar( nptr->Status(key_rho), eos.Rho_CarbonicPhase( pf_top, nptr->Read(key_T) ) ) );
     nptr->Store( key_mu, makeScalar( nptr->Status(key_mu), eos.mu_CarbonicPhase( pf_top, nptr->Read(key_T) ) ) );
 
     while ( nptr != bottom_node_ptr )
       {
          // computing the fluid properties at the deeper node:
          Node<1U>* nptr2 = (eptr->N(0) == nptr) ? eptr->N(1) : eptr->N(0);
          assert( nptr2 != nullptr );
          assert( nptr2->x() < nptr->x() );
          // estimating pfluid pressure at deeper node using fluid density at barycentre
          double64 elength  = eptr->Volume(); // =length in 1D
          double64 ToC      = eptr->PropertyValueAtBaryCenter( key_T );
          double64 pf_above = nptr->Read(key_pf);
          // guessing fluid density
          double64 pf_guess = pf_above + (elength/2.) * acc_gravity * nptr->Read(key_rho);
          double64 rho_mix = eos.Rho_CarbonicPhase( pf_guess, ToC );
          // refining the guess
          pf_guess = pf_above + (elength/2.) * acc_gravity * rho_mix;
          rho_mix = eos.Rho_CarbonicPhase( pf_guess, ToC );
          pf_guess = pf_above + (elength/2.) * acc_gravity * rho_mix;
          rho_mix = eos.Rho_CarbonicPhase( pf_guess, ToC );
          // storing the estimates
          double64 pf_below = pf_above + elength * acc_gravity * rho_mix;
          nptr2->Store( key_pf, makeScalar( nptr2->Status(key_pf), pf_below ) );
          nptr2->Store( key_rho, makeScalar( nptr2->Status(key_rho), eos.Rho_CarbonicPhase( pf_below, ToC ) ) );
          nptr2->Store( key_mu, makeScalar( nptr2->Status(key_mu), eos.mu_CarbonicPhase( pf_below, ToC ) ) );
          eptr->Store( key_rhom, makeScalar(eptr->Status(key_rhom),rho_mix) );
          // continue to next element
          if ( nptr2->Parents() > 1U ) {
               eptr = (nptr2->Parent(0) == eptr) ? nptr2->Parent(1) : nptr2->Parent(0);
            }
          else break;
          nptr = nptr2;
       }
 
     // for comparison, you may want to introduce a new variable
     if ( verbose_ ) {
          printRangeOfVariable( *model1D_, "density carbonic phase" );
          printRangeOfVariable( *model1D_, "viscosity carbonic phase" );
          printRangeOfVariable( *model1D_, "fluid mixture density" );
          printRangeOfVariable( *model1D_, "reference pressure" );
       }

 } // end ReferencePressureByTopDownIntegration







 
/**
    Assumes datum of zero as reference level at atmospheric pressure
*/
void GravityInducedFluidPressure_Test::ReferencePressureForFixedDensity( double64 ref_density )
 {
    const csmp::Index p_key(model1D_->Database().StorageKey("reference pressure"));
    const csmp::Index g_key(model1D_->Database().StorageKey("acceleration gravity"));
    Region<1U>&       model_domain = model1D_->Region("Model");
    const double64    datum(top_); // model top with Y-axis pointing upwards
    const double64    acc_gravity(model1D_->Read(g_key));
   
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ )
      (*nit)->Store( p_key, makeScalar( (*nit)->Status(p_key), ((*nit)->y() - datum) * acc_gravity * ref_density + patm_ ) );
   
    if ( verbose_ ) printRangeOfVariable( *model1D_, "fluid pressure" );
 }





/**
    Checks whether the difference between 'fluid pressure'
    and 'reference pressure' is smaller than the allowed tolerance.
*/
bool GravityInducedFluidPressure_Test::TestComputedWithReferencePressure()
 {
    const csmp::Index pr_key(model1D_->Database().StorageKey("reference pressure"));
    const csmp::Index pf_key(model1D_->Database().StorageKey("fluid pressure"));
    const Region<1U>& model_domain = model1D_->Region("Model");
    size_t out_of_range_vals(0U);

    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 abs_val_pf_diff = fabs( (*nit)->Read(pr_key) - (*nit)->Read(pf_key ) );
         _test( abs_val_pf_diff <= ptol_ );
         if ( abs_val_pf_diff > ptol_ )
           out_of_range_vals++;
      }
   
    if ( out_of_range_vals > 0U ) return false;
    return true;
   
 } // end TestComputedWithReferencePressure







/**
    Compares the node variable values with one-another reporting false if they differ by more than the user-specified tolerance
*/
template<size_t dim>
bool GravityInducedFluidPressure_Test::TestResultsByComparison( const Model<dim>* const model, const char* test_variable, const char* reference_variable, double64 tolerance )
 {
    const csmp::Index test_key = model->Database().StorageKey(test_variable);
    const csmp::Index refv_key = model->Database().StorageKey(reference_variable);
    assert( test_key.type == SCALAR );
    assert( refv_key.type == SCALAR );
    assert( test_key.place == NODE );
    assert( refv_key.place == NODE );

    const Region<dim>& model_domain = model->Region("Model");
    bool values_within_tolerance(true);
    for ( typename vector<Node<dim>*>::const_iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 difference = fabs( (*nit)->Read(refv_key) - (*nit)->Read(test_key) );
         _test( difference <= tolerance );
         if ( difference > tolerance )
           values_within_tolerance = false;
      }
   
    return values_within_tolerance;
   
 } // end TestResultsByComparison
 
template bool GravityInducedFluidPressure_Test::TestResultsByComparison( const Model<1U>* const, const char*, const char*, double64 );
template bool GravityInducedFluidPressure_Test::TestResultsByComparison( const Model<2U>* const, const char*, const char*, double64 );
template bool GravityInducedFluidPressure_Test::TestResultsByComparison( const Model<3U>* const, const char*, const char*, double64 );






/**
    Uses SAMG solver and gravity term to compute vertical fluid pressure profile
 
    This computation needs to read 'fluid mixture density' as input values.
*/
void GravityInducedFluidPressure_Test::ComputeCO2Pressure_PDE_Integrator( double64 pf_top )
 {
    cout << "\n\n\nComputeCO2Pressure_PDE_Integrator: Computing CO2-static pressure..." << endl;

    const csmp::Index rho_key = model1D_->Database().StorageKey("fluid mixture density");
    // the last node (id=n_nodes) is located at the top of the model (by anology with the Y-axis)
    Region<1U>& model_domain(model1D_->Region("Model"));
    model_domain.UpdateMemberIndexes();
    const size_t topNode = findNode( *model1D_, top_, 1.0e-7 );

    model_domain.ChangePropertyStatus( "fluid pressure", ANY );
    model_domain.N(topNode)->Store( model1D_->Database().StorageKey("fluid pressure"), makeScalar(DIRICH,pf_top) );

    // 7. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
    GaussJordan_Solver         GJ_solver;
    PDE_Integrator<1U,Region>  hydrostatic_pressure(GJ_solver);
  
    NumIntegral_dNT_op_dN_dV<1U,Element<1U> >  hydrostatic_conductance( model1D_->Database(), "total mobility permeability product",  "fluid pressure", "fluid pressure" );

    const csmp::Index g_key(model1D_->Database().StorageKey("acceleration gravity"));
    const double64    acc_gravity(model1D_->Read(g_key));
    // unless specified otherwise, in a 1D model, gravity will automatically act in the x-direction
    NumIntegral_NT_op_dNi_dV<1U,Element<1U> >  hydrostatic_gravity( model1D_->Database(), "fluid mixture density",
                                                                   "total mobility permeability product", "fluid pressure", acc_gravity );
    hydrostatic_pressure.Add( &hydrostatic_conductance );
    hydrostatic_pressure.Add( &hydrostatic_gravity );


    // 8. Compute the initial hydrostatic pressure. Since fluid properties will change after the pressure was computed
    //     recompute the fluid properties and fluid pressure 3 times such that fluid pressure and fluid properties converge
    // --------------------------------------------------------------------------------------------------------------------
//    model1D_->InputPropertyValue( "total mobility permeability product", makeScalar(PLAIN,1.0e-6) );
    printRangeOfVariable( *model1D_, "total mobility permeability product" );
    printRangeOfVariable( *model1D_, "fluid mixture density" );

    hydrostatic_pressure.IntegrateOver( model_domain );
   
    printRangeOfVariable( *model1D_, "fluid pressure" );

 } // end ComputeCO2Pressure_PDE_Integrator






/**
    Version that computes the fluid-static pressure exactly like the ACGSS simulator
*/
void GravityInducedFluidPressure_Test::ComputeCO2Pressure_PDE_Integrator2( double64 pf_top )
 {
    cout << "\n\n\nComputeCO2Pressure_PDE_Integrator2: Computing CO2-static fluid pressure..." << endl;

    Region<1U>& model_domain(model1D_->Region("Model"));
    model_domain.UpdateMemberIndexes();
    const size_t topNode = findNode( *model1D_, top_, 1.0e-7 );

    model_domain.ChangePropertyStatus( "fluid pressure", ANY );
    model_domain.N(topNode)->Store( model1D_->Database().StorageKey("fluid pressure"), makeScalar(DIRICH,pf_top) );

    // 1. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
    const size_t                   dim(1U);
    GaussJordan_Solver             GJ_solver;
    PDE_Integrator<dim,Region>     hydrostatic_pressure(GJ_solver);
    NumIntegral_dNT_op_dN_dV<dim>  p_conductance( model1D_->Database(), "total mobility permeability product", "fluid pressure", "fluid pressure" );
    NumIntegral_dNT_op_dV<dim>     gravity( model1D_->Database(), "gravity term", "fluid pressure" );

    hydrostatic_pressure.Add( &p_conductance );
    hydrostatic_pressure.Add( &gravity );
   
    // 1.1 initialising 'total mobility permeability product' including the density term
    const csmp::Index mobt_key = model1D_->Database().StorageKey("total mobility permeability product");
    const csmp::Index k_key    = model1D_->Database().StorageKey("permeability");
    const csmp::Index rho_key  = model1D_->Database().StorageKey("fluid mixture density");
    const csmp::Index mu_key   = model1D_->Database().StorageKey("viscosity carbonic phase");
    // rho k/mu
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         double64 mobt = (*it)->Read( rho_key ) * (*it)->Read( k_key ) / (*it)->PropertyValueAtBaryCenter( mu_key );
         (*it)->Store( mobt_key, makeScalar(PLAIN,mobt) );
      }
    printRangeOfVariable( *model1D_, "total mobility permeability product" );

    // 1.2 initialize gravity term
    VectorVariable<dim>  grav_vec(DIRICH,0.);
    grav_vec(1) = -1.;
    model1D_->InputPropertyValue( "dip vector", grav_vec );
    // recovering acceleration of gravity
    const csmp::Index  g_key(model1D_->Database().StorageKey("acceleration gravity"));
    const double64     acc_gravity(model1D_->Read(g_key));
    // computing the gravity term: k g (rho_CO2 * lambda_CO2 * rho_CO2 + rho_water * lambda water * rho water) for the PDE operator
    const csmp::Index  gv_key(model1D_->Database().StorageKey("dip vector"));
    const csmp::Index  gt_key(model1D_->Database().StorageKey("gravity term"));
    // single phase version: rho^2 k/mu g
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        (*it)->Read( gv_key, grav_vec );
        grav_vec(1) = -acc_gravity * (*it)->Read(rho_key) * (*it)->Read(rho_key) * ((*it)->Read(k_key) / (*it)->PropertyValueAtBaryCenter(mu_key));
        (*it)->Store( gt_key, grav_vec );
     }
    printRangeOfVariable( *model1D_, "gravity term" );


    // 2. Compute hydrostatic pressure
    // -------------------------------------------------------------
    printRangeOfVariable( *model1D_, "fluid mixture density" );

    hydrostatic_pressure.IntegrateOver( model_domain );
   
    printRangeOfVariable( *model1D_, "fluid pressure" );

 } // end ComputeCO2Pressure_PDE_Integrator2




void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_Integrator2( Model<1U>* const model, double64 pf_top )
 {
    cout << "\n\n\nComputeCO2PressureFromReducedPressure_PDE_Integrator2: Computing reduced CO2-static pressure..." << endl;
    const size_t dim(1U);
    Region<dim>& model_domain(model->Region("Model"));
    model_domain.UpdateMemberIndexes();
    const size_t topNode = findNode( *model, top_, 1.0e-7 );

    model_domain.ChangePropertyStatus( "reduced fluid pressure", ANY );
    model_domain.N(topNode)->Store( model->Database().StorageKey("reduced fluid pressure"), makeScalar(DIRICH,pf_top) );

    // 1. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
    GaussJordan_Solver             GJ_solver;
    PDE_Integrator<dim,Region>     hydrostatic_pressure(GJ_solver);
    NumIntegral_dNT_op_dN_dV<dim>  p_conductance( model->Database(), "total mobility permeability product", "reduced fluid pressure", "reduced fluid pressure" );
    NumIntegral_dNT_op_dV<dim>     gravity( model->Database(), "gravity term", "reduced fluid pressure" );

    hydrostatic_pressure.Add( &p_conductance );
    hydrostatic_pressure.Add( &gravity );
   
    // 1.1 initialising 'total mobility permeability product' including the density term
    const csmp::Index mobt_key = model->Database().StorageKey("total mobility permeability product");
    const csmp::Index k_key    = model->Database().StorageKey("permeability");
    const csmp::Index rho_key  = model->Database().StorageKey("fluid mixture density");
    const csmp::Index mu_key   = model->Database().StorageKey("viscosity carbonic phase");
    // rho k/mu
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         double64 mobt = (*it)->Read( rho_key ) * (*it)->Read( k_key ) / (*it)->PropertyValueAtBaryCenter( mu_key );
         (*it)->Store( mobt_key, makeScalar(PLAIN,mobt) );
      }
    printRangeOfVariable( *model1D_, "total mobility permeability product" );
   
    // 1.2 choosing a reference fluid density
    const double64 reference_density(100.);

    // 1.3 initialize gravity term
    VectorVariable<dim>  grav_vec(DIRICH,0.);
    grav_vec(1) = -1.;
    model1D_->InputPropertyValue( "dip vector", grav_vec );
    // recovering acceleration of gravity
    const csmp::Index  g_key(model->Database().StorageKey("acceleration gravity"));
    const double64     acc_gravity(model->Read(g_key));
    // computing the gravity term: k g (rho_CO2 * lambda_CO2 * rho_CO2 + rho_water * lambda water * rho water) for the PDE operator
    const csmp::Index  gv_key(model->Database().StorageKey("dip vector"));
    const csmp::Index  gt_key(model->Database().StorageKey("gravity term"));
    // single phase version: rho^2 k/mu g
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        (*it)->Read( gv_key, grav_vec );
        grav_vec(1) = -acc_gravity * (*it)->Read(rho_key) * ((*it)->Read(rho_key) - reference_density) * ((*it)->Read(k_key) / (*it)->PropertyValueAtBaryCenter(mu_key));
        (*it)->Store( gt_key, grav_vec );
     }
    printRangeOfVariable( *model, "gravity term" );


    // 2. Compute hydrostatic pressure
    // -------------------------------------------------------------
    printRangeOfVariable( *model, "fluid mixture density" );

    hydrostatic_pressure.IntegrateOver( model_domain );
   
    printRangeOfVariable( *model, "reduced fluid pressure" );
   
    // post-processing the absolute fluid pressure
    const csmp::Index pr_key   = model->Database().StorageKey("fluid pressure");
    const csmp::Index pf_key   = model->Database().StorageKey("reduced fluid pressure");
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         // adding ref-density-static pressure gradient to 'reduced pressure'
         double64 pf = (*nit)->Read( pr_key ) + (top_ - (*nit)->x()) * reference_density * acc_gravity;
         (*nit)->Store( pf_key, makeScalar( (*nit)->Status(pf_key),pf) );
      }

    printRangeOfVariable( *model, "fluid pressure" );

 } // end ComputeCO2PressureFromReducedPressure_PDE_Integrator2






/**
    Version that computes the fluid-static pressure exactly like the ACGSS simulator
*/
template<size_t dim>
void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan( Model<dim>* const model, double64 pf_top )
 {
    cout << "\n\n\nComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan<dim>: Computing reduced CO2-static pressure..." << endl;
   
    Region<dim>& model_domain(model->Region("Model"));
    model_domain.ChangePropertyStatus( "reduced fluid pressure", ANY );
   
    // 0. looping over the top boundary to make topography adjustments
    // ---------------------------------------------------------------
    Boundary<dim>& top_boundary = model->Boundary("TOP");
    const csmp::Index  g_key(model3D_->Database().StorageKey("acceleration gravity"));
    const double64     acc_gravity(model3D_->Read(g_key));
    const csmp::Index key_rf(model3D_->Database().StorageKey("reduced fluid pressure"));
    const csmp::Index key_rho(model3D_->Database().StorageKey("density carbonic phase"));
    Point<3U> xyz_min, xyz_max;
    model3D_->MinMaxCoordinates( xyz_min, xyz_max );
    top_ = xyz_max[1];
    for ( typename vector<Node<dim>*>::iterator nit=top_boundary.NodesBegin(); nit!=top_boundary.NodesEnd(); nit++ ) {
         double64 pf = pf_top + (top_ - (*nit)->y()) * (*nit)->Read(key_rho) * acc_gravity;
         assert( !isnan(pf) );
         (*nit)->Store( key_rf, makeScalar(DIRICH,pf) );
      }
    printRangeOfVariable( *model, "TOP", "reduced fluid pressure" );

    // 1. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
    GaussJordan_Solver             GJ_solver;
    PDE_Integrator<dim,Region>     hydrostatic_pressure(GJ_solver);
    NumIntegral_dNT_op_dN_dV<dim>  p_conductance( model->Database(), "total mobility permeability product", "reduced fluid pressure", "reduced fluid pressure" );
    NumIntegral_dNT_op_dV<dim>     gravity( model->Database(), "gravity term", "reduced fluid pressure" );

    hydrostatic_pressure.Add( &p_conductance );
    hydrostatic_pressure.Add( &gravity );
   
    // 1.1 initialising 'total mobility permeability product' including the density term
    const csmp::Index mobt_key = model->Database().StorageKey("total mobility permeability product");
    const csmp::Index k_key    = model->Database().StorageKey("permeability");
    const csmp::Index rho_key  = model->Database().StorageKey("fluid mixture density");
    const csmp::Index mu_key   = model->Database().StorageKey("viscosity carbonic phase");
    // rho k/mu
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         double64 mobt = (*it)->Read( rho_key ) * (*it)->Read( k_key ) / (*it)->PropertyValueAtBaryCenter( mu_key );
         (*it)->Store( mobt_key, makeScalar(PLAIN,mobt) );
      }
    printRangeOfVariable( *model, "total mobility permeability product" );
   
    // 1.2 choosing a reference fluid density
    const double64 reference_density(100.);

    // 1.3 initialize gravity term
    VectorVariable<dim>  grav_vec(DIRICH,0.);
    grav_vec(1) = -1.;
    model->InputPropertyValue( "dip vector", grav_vec );
    // computing the gravity term: k g (rho_CO2 * lambda_CO2 * rho_CO2 + rho_water * lambda water * rho water) for the PDE operator
    const csmp::Index  gv_key(model->Database().StorageKey("dip vector"));
    const csmp::Index  gt_key(model->Database().StorageKey("gravity term"));
    // single phase version: rho^2 k/mu g
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        (*it)->Read( gv_key, grav_vec );
        grav_vec(1) = -acc_gravity * (*it)->Read(rho_key) * ((*it)->Read(rho_key) - reference_density) * ((*it)->Read(k_key) / (*it)->PropertyValueAtBaryCenter(mu_key));
        (*it)->Store( gt_key, grav_vec );
     }
    printRangeOfVariable( *model, "gravity term" );


    // 2. Compute hydrostatic pressure iteratively
    // -------------------------------------------------------------
    printRangeOfVariable( *model, "fluid mixture density" );

    hydrostatic_pressure.IntegrateOver( model_domain );
   
/* FOR ITERATIVE APPROXIMATION OF PROPERTIES

    cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid properties... " << endl;
    for ( uint32 i=0; i<=5U; i++ ) {
         cout <<"\n\titeration "<< i+1U <<":"<< endl;
         model1D_->Apply( hydrostatic_pressure );
         model1D_->Accept( properties_visitor );
         model1D_->InterpolateNodeToElementProperty( "fluid density", "element fluid density" );
         rhof += total_dissolved_solids;
         printRangeOfVariable( *model1D_, "fluid pressure" );
         printRangeOfVariable( *model1D_, "element fluid density" );
      }

*/

    printRangeOfVariable( *model, "reduced fluid pressure" );
   
    // post-processing the absolute fluid pressure
    const csmp::Index pr_key   = model->Database().StorageKey("fluid pressure");
    const csmp::Index pf_key   = model->Database().StorageKey("reduced fluid pressure");
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         // adding ref-density-static pressure gradient to 'reduced pressure'
         double64 pf = (*nit)->Read( pr_key ) + (top_ - (*nit)->y()) * reference_density * acc_gravity;
         (*nit)->Store( pf_key, makeScalar( (*nit)->Status(pf_key),pf) );
      }

    printRangeOfVariable( *model, "fluid pressure" );

 } // end ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan

template void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan( Model<2U>* const, double64 );
template void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan( Model<3U>* const, double64 );










/**
    Version that computes the fluid-static pressure exactly like the ACGSS simulator
*/
template<size_t dim>
void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG( Model<dim>* const model, double64 pf_top )
 {
    cout << "\n\n\nComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG<"<< dim <<">: Computing reduced CO2-static pressure..." << endl;
   
    Region<dim>& model_domain(model->Region("Model"));
    model_domain.ChangePropertyStatus( "reduced fluid pressure", ANY );
   
    // 0. looping over the top boundary to make topography adjustments
    // ---------------------------------------------------------------
    Boundary<dim>& top_boundary = model->Boundary("TOP");
    const csmp::Index  g_key(model3D_->Database().StorageKey("acceleration gravity"));
    const double64     acc_gravity(model3D_->Read(g_key));
    const csmp::Index key_rf(model3D_->Database().StorageKey("reduced fluid pressure"));
    const csmp::Index key_rho(model3D_->Database().StorageKey("density carbonic phase"));
    Point<3U> xyz_min, xyz_max;
    model3D_->MinMaxCoordinates( xyz_min, xyz_max );
    top_ = xyz_max[1];
    for ( typename vector<Node<dim>*>::iterator nit=top_boundary.NodesBegin(); nit!=top_boundary.NodesEnd(); nit++ ) {
         double64 pf = pf_top + (top_ - (*nit)->y()) * (*nit)->Read(key_rho) * acc_gravity;
         assert( !isnan(pf) );
         (*nit)->Store( key_rf, makeScalar(DIRICH,pf) );
      }
    printRangeOfVariable( *model, "TOP", "reduced fluid pressure" );

    // 1. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
    SAMG_Settings                  settings;
    SAMG_Solver                    samg_solver( &settings );
    PDE_Integrator_UoM<dim,Region> hydrostatic_pressure(samg_solver);
    NumIntegral_dNT_op_dN_dV<dim>  p_conductance( model->Database(), "total mobility permeability product", "reduced fluid pressure", "reduced fluid pressure" );
    NumIntegral_dNT_op_dV<dim>     gravity( model->Database(), "gravity term", "reduced fluid pressure" );

    hydrostatic_pressure.Add( &p_conductance );
    hydrostatic_pressure.Add( &gravity );
   
    // specify the settings for SAMG here
   
    // 1.1 initialising 'total mobility permeability product' including the density term
    const csmp::Index mobt_key = model->Database().StorageKey("total mobility permeability product");
    const csmp::Index k_key    = model->Database().StorageKey("permeability");
    const csmp::Index rho_key  = model->Database().StorageKey("fluid mixture density");
    const csmp::Index mu_key   = model->Database().StorageKey("viscosity carbonic phase");
    // rho k/mu
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         double64 mobt = (*it)->Read( rho_key ) * (*it)->Read( k_key ) / (*it)->PropertyValueAtBaryCenter( mu_key );
         (*it)->Store( mobt_key, makeScalar(PLAIN,mobt) );
      }
    printRangeOfVariable( *model, "total mobility permeability product" );
   
    // 1.2 choosing a reference fluid density
    const double64 reference_density(100.);

    // 1.3 initialize gravity term
    VectorVariable<dim>  grav_vec(DIRICH,0.);
    grav_vec(1) = -1.;
    model->InputPropertyValue( "dip vector", grav_vec );
    // computing the gravity term: k g (rho_CO2 * lambda_CO2 * rho_CO2 + rho_water * lambda water * rho water) for the PDE operator
    const csmp::Index  gv_key(model->Database().StorageKey("dip vector"));
    const csmp::Index  gt_key(model->Database().StorageKey("gravity term"));
    // single phase version: rho^2 k/mu g
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        (*it)->Read( gv_key, grav_vec );
        grav_vec(1) = -acc_gravity * (*it)->Read(rho_key) * ((*it)->Read(rho_key) - reference_density) * ((*it)->Read(k_key) / (*it)->PropertyValueAtBaryCenter(mu_key));
        (*it)->Store( gt_key, grav_vec );
     }
    printRangeOfVariable( *model, "gravity term" );


    // 2. Compute hydrostatic pressure iteratively
    // -------------------------------------------------------------
    printRangeOfVariable( *model, "fluid mixture density" );

    hydrostatic_pressure.IntegrateOver( model_domain );
   
/* FOR ITERATIVE APPROXIMATION OF PROPERTIES

    cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid properties... " << endl;
    for ( uint32 i=0; i<=5U; i++ ) {
         cout <<"\n\titeration "<< i+1U <<":"<< endl;
         model1D_->Apply( hydrostatic_pressure );
         model1D_->Accept( properties_visitor );
         model1D_->InterpolateNodeToElementProperty( "fluid density", "element fluid density" );
         rhof += total_dissolved_solids;
         printRangeOfVariable( *model1D_, "fluid pressure" );
         printRangeOfVariable( *model1D_, "element fluid density" );
      }

*/

    printRangeOfVariable( *model, "reduced fluid pressure" );
   
    // post-processing the absolute fluid pressure
    const csmp::Index pr_key   = model->Database().StorageKey("fluid pressure");
    const csmp::Index pf_key   = model->Database().StorageKey("reduced fluid pressure");
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         // adding ref-density-static pressure gradient to 'reduced pressure'
         double64 pf = (*nit)->Read( pr_key ) + (top_ - (*nit)->y()) * reference_density * acc_gravity;
         (*nit)->Store( pf_key, makeScalar( (*nit)->Status(pf_key),pf) );
      }

    printRangeOfVariable( *model, "fluid pressure" );

 } // end ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG

template void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG( Model<2U>* const, double64 );
template void GravityInducedFluidPressure_Test::ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG( Model<3U>* const, double64 );





/**
    Outputs the variables discretised on the model:
 
    - fluid pressure
    - temperature
    - fluid densities and viscosities
 
    appends "-results.txt" to file name.
*/
void GravityInducedFluidPressure_Test::OutputResultsToText( const char* file_name ) const
 {
    TextInterface  text_output;
    list<string>   output_variables_node({"fluid pressure", "temperature", "reference pressure",
                                          "density carbonic phase", "viscosity carbonic phase" });

    list<string>   output_variables_elmt({"fluid mixture density"});

    // plot these results with MS Excel or similar
    CoordinateTransformer<1U>  no_transformations;
    string filename(file_name);
    filename += "-results";
    text_output.OutputDataAsTextColumns( *model1D_, filename.c_str(), "Model",
                                         no_transformations,
                                         output_variables_node );

    text_output.OutputDataAsTextColumns( *model1D_, filename.c_str(), "Model",
                                         no_transformations,
                                         output_variables_elmt );
 } // end OutputResultsToText



} // end csmp
