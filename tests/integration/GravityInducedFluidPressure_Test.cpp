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
#include "EOS_CO2H2ONaCl_Spycher2004.h"

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
    // mesh_container.Out(); // OK - neighbor connectivity etc.
 
    model1D_ = new Model<1U>( mesh_container, "GravityInducedFluidPressure_Test-variables.txt", true, false );
    Region<1U>& model_domain(model1D_->Region("Model"));
    model_domain.UpdateMemberIndexes();
    printModelDimensions( *model1D_ );
    // NOTE: Y points upward, else gravity will act in opposite direction
    const size_t topNode2km   = findNode( *model1D_, top_, 1.0e-7 );
    const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
    assert( topNode2km < N_ELEMENTS+1 );
    assert( bottomNode0m < 101 );
    cout <<"\nGravityInducedFluidPressure_Test: model end points:\n";
    cout <<"\t"<< model_domain.N(bottomNode0m)->Coordinate() << endl;
    cout <<"\t"<< model_domain.N(topNode2km)->Coordinate() << endl;

    // Input of material properties and initial conditions
    model1D_->InputPropertyValue( "total mobility permeability product", makeScalar(PLAIN,1.0e-13) );
    model1D_->InputPropertyValue( "fluid pressure",        makeScalar(PLAIN,1.0e7) ); // 1 bar
    model1D_->InputPropertyValue( "thermal conductivity",  makeScalar(PLAIN,2.2) );
//    model1D_->InputPropertyValue( "saturation aqueous phase",  makeScalar(PLAIN,0.) );
//    model1D_->InputPropertyValue( "saturation carbonic phase",  makeScalar(PLAIN,1.) );
    model1D_->InputPropertyValue( "density carbonic phase",  makeScalar(PLAIN,200.) );
    model1D_->InputPropertyValue( "viscosity carbonic phase",  makeScalar(PLAIN,1.0e-5) );
    model1D_->InputPropertyValue( "acceleration gravity",  makeScalar(PLAIN,9.8061) );
   
 } // end InitialiseModel1D





/**
    From top temperature and a positive temperature gradient in oC/m.
*/
void GravityInducedFluidPressure_Test::InitialiseTemperatureProfile( double64 T_top, double64 grad_T )
 {
    assert( T_top  > 8. );
    assert( grad_T > 0. );
    const csmp::Index key_T(model1D_->Database().StorageKey("temperature"));
    Region<1U>& model_domain = model1D_->Region("Model");
    //             oC           K/m                W/m3
    double64 heat_source(0.);
    model1D_->InputPropertyValue( "energy source", makeScalar(PLAIN,heat_source) );
    // default: Y = 15oC @ 2000m, increasing downward on a gradient of 45C/km
    for ( vector<Node<1U>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
         double64 ToC = T_top + (top_ - (*nit)->x()) * grad_T;
         (*nit)->Store( key_T, makeScalar(PLAIN,ToC) );
      }
    printRangeOfVariable( *model1D_, "temperature" );

 } // end InitialiseTemperatureProfile




/**
    TESTING IS DONE HERE
*/
void GravityInducedFluidPressure_Test::run()
 {
    // testing pressure computation for a constant fluid density
    const double64 ref_density(1000.);
    ReferencePressureForFixedDensity( ref_density );
    ReferencePressureByTopDownIntegration( ref_density );
    _test( TestComputedWithReferencePressure() );
    if ( verbose_ ) OutputResultsToText( "GravityInducedFluidPressure_Test0" );

    // testing computation for compressible CO2 (Spycher EOS)
    const double64 grad_T_K_per_m(0.02); // 20oC/km
    InitialiseTemperatureProfile( 12., grad_T_K_per_m );
    ReferencePressureByTopDownIntegrationCO2();
    if ( verbose_ ) OutputResultsToText( "GravityInducedFluidPressure_Test1" );
 }
  
  
  
  
  
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
     const size_t topNode2km   = findNode( *model1D_, top_, 1.0e-7 );
     const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
     const double64 acc_gravity(model1D_->Read(g_key));

     // finding the top node
     Node<1U>* nptr = model_domain.N(topNode2km);
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
          double64 rho_mix = fluid_density;
          double64 pf      = nptr->Read(key_pf) + elength * rho_mix * acc_gravity;
         
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
     assert( printRangeOfVariable( *model1D_, "fluid pressure" ) >= patm_ );

 } // end ReferencePressureByTopDownIntegration







/**
    Integrating from the surface downward to the bottom of the model, the pressure is accumulated.
    CO2 density is taken from the Spycher et al equation of state.
*/
void GravityInducedFluidPressure_Test::ReferencePressureByTopDownIntegrationCO2( double64 pf_top )
  {
     const csmp::Index key_T(model1D_->Database().StorageKey("temperature"));
     const csmp::Index key_pf(model1D_->Database().StorageKey("fluid pressure"));
     const csmp::Index g_key(model1D_->Database().StorageKey("acceleration gravity"));
     const double64    acc_gravity(model1D_->Read(g_key));
     const csmp::Index key_rho(model1D_->Database().StorageKey("density carbonic phase"));
     const csmp::Index key_mu(model1D_->Database().StorageKey("viscosity carbonic phase"));
     const csmp::Index key_rhom(model1D_->Database().StorageKey("fluid mixture density"));

     Region<1U>& model_domain(model1D_->Region("Model"));
     model_domain.UpdateMemberIndexes();
     const size_t topNode2km   = findNode( *model1D_, top_, 1.0e-7 );
     const size_t bottomNode0m = findNode( *model1D_, 0., 1.0e-7 );
     // finding the top node
     Node<1U>* nptr = model_domain.N(topNode2km);
     Node<1U>* const bottom_node_ptr = model_domain.N(bottomNode0m);
     assert( nptr != nullptr );
     Element<1U>* eptr = nptr->Parent(0);
     assert( eptr != nullptr );
     // initialising EOS on first node
     nptr->Store( key_pf, makeScalar(DIRICH,pf_top) );
    
     // Spycher et al 2003 EOS
     EOS_CO2H2ONaCl_Spycher04  eos;
 
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
     printRangeOfVariable( *model1D_, "density carbonic phase" );
     printRangeOfVariable( *model1D_, "viscosity carbonic phase" );
     printRangeOfVariable( *model1D_, "fluid mixture density" );
     printRangeOfVariable( *model1D_, "fluid pressure" );

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
 }





bool GravityInducedFluidPressure_Test::TestComputedWithReferencePressure()
 {
    const csmp::Index pr_key(model1D_->Database().StorageKey("reference pressure"));
    const csmp::Index pf_key(model1D_->Database().StorageKey("reference pressure"));
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
    Outputs the variables discretised on the model:
 
    - fluid pressure
    - temperature
    - fluid densities and viscosities
 
    appends "-results.txt" to file name.
*/
void GravityInducedFluidPressure_Test::OutputResultsToText( const char* file_name ) const
 {
    TextInterface  text_output;
    list<string>   output_variables_node({"fluid pressure", "temperature",
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
