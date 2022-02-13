//
//  SlopeMechanics_Example
//  CSMP_API_library
//
//  Created by Stephan Matthai on 1/27/17.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "SlopeMechanics_Example.h"
#include "compareFloats.h"
#include "ANSYS_Model2D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "PDE_Integrator_UoM.h"
#include "PDE_IntegratorExperimental.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#endif

#include "InputDataManager.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"

// include any header files that you need here...
#include "SteadyStateDiffusor.h"
#include "CSMP_highLevelUtilities.h"
#include "ComputationalSettings.h"
#include "IAPWS_H2OPropertiesVisitor.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "StressesAndStrains.h"
#include "StressesAndStrains2.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "SinglePhaseVelocityVisitor.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;

namespace csmp {

void SlopeMechanics_Example::Specifications()
  {
     SetTitle( "Fluid pressure and stress distribution in hill slope" );
     SetDifficulty( 2 );
     SetCategory( "Simulation of Physical Processes" );
     AddAuthor( "Stephan Matthai" );
     AddDescription( "source in: SlopeMechanics_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "Input file suite: Jura-slope1" );
//     AddRequirement( "no predefined model or variables file" );
  }





/**
     Put CSMP code that you would like to test here and run it as part of the
     example suite.
 
     Systems Modelling and Design model that calculates the effective stress in a dam,
     using the gravitational loading and the plane stress assumption.
 
     input models:
     - slope_model1
*/
void SlopeMechanics_Example::Run()
{
  // ----------------------------------------------------------
  // 1. building and configering model from ANSYS - csp dataset
  // ----------------------------------------------------------
   string input_file("Jura-slope1");
   ANSYS_Model2D   model( input_file.c_str(), "CSMP_field_scale_mechanics_variables.txt" );
   Region<DIM>& model_domain(model.Region("Model"));

   printModelDimensions( model, true );

   InputDataManager<DIM>  model_configuration;
   ComputationalSettings settings;

   model_configuration.ConfigureFromFile( model, input_file.c_str(),
                                          false,           ///< region name from parameter range
                                          true,            ///< default property values
                                          true,            ///< regional property values
                                          true,            ///< boundary conditions for box-shaped model
                                          true,           ///< essential conditions for regions
                                          true,           ///< boundary conditions for arbitrary-shaped model
                                          settings );

  
  // ----------------------------------------------------------
  // 2. computing initial steady state temperature distribution
  // ----------------------------------------------------------
    {
      SteadyStateDiffusor<DIM,Region>  temperature( model, "thermal conductivity", "temperature", "energy source");
      printRangeOfVariable( model, "thermal conductivity" );
      printRangeOfVariable( model, "energy source");
      // TODO: consider potential stress changes due to insolation of slope etc.
      model.InputBoundaryValue( TOP, "temperature", makeScalar(DIRICH,15.));
      //model.InputBoundaryValue( BOTTOM, "temperature", makeScalar(DIRICH,80.));
      temperature.ComputeSteadyState( model.Region("Model") );
      printRangeOfVariable( model, "temperature");
    }
   VTK_Interface<DIM>  vtk_output;
   vtk_output.OutputDataToVTK( model, "temperature", "temperature", 0, true );


  
  // -------------------------------------------------------------------------
  // 3. computing initial hydrostatic pressure (no loading via grain skeleton)
  // -------------------------------------------------------------------------
  //  Visitor computes thermodynamic properties of water; density is interpolated to
  //     element barycentre for later vertical integration
  // --------------------------------------------------------------------------------------
  IAPWS_H2OPropertiesVisitor<DIM>  properties_visitor( model, "fluid pressure" ,"fluid density", "fluid viscosity");
  // computing fluid density and viscosity
  model.Accept( properties_visitor );
  printRangeOfVariable( model, "fluid density");
  printRangeOfVariable( model, "fluid viscosity");

  // hydraulic conductivity K = k/mu
  const csmp::Index perm_key(model.Database().StorageKey("permeability"));
  const csmp::Index visc_key(model.Database().StorageKey("fluid viscosity"));
  const csmp::Index cond_key(model.Database().StorageKey("conductivity"));
  ScalarVariable  visc;
  for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
       (*it)->PropertyValueAtBaryCenter( visc_key, visc );
       (*it)->Store( cond_key, makeScalar( (*it)->Status(cond_key), (*it)->Read(perm_key) / visc() ) );
    }
  printRangeOfVariable( model, "conductivity" );

  model.InterpolateNodeToElementProperty( "fluid density", "element fluid density" );
  
  // Setting up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
  SAMG_Settings  samg_settings;
  SAMG_Solver    samg_solver(&samg_settings);
  PDE_Integrator_UoM<DIM,Region>  hydrostatic_pressure(samg_solver);
  // minimizing screen output
  samg_settings.Set_iout1( 0 );
  samg_settings.Set_iout2( 0 );
  samg_settings.Set_idmp( -1 );
  samg_settings.SetSolverInstance(2);

  NumIntegral_dNT_op_dN_dV<DIM>  hydrostatic_conductance( model.Database(),
                                                         "conductivity",
                                                         "fluid pressure", "fluid pressure" );

  // cin >> "\nmain: Enter the acceleration of gravity (kg/m.s2): in the area of interest: ";
  double acc_gravity(9.81);
  // unless specified otherwise, in a 1D model, gravity will automatically act in the x-direction
  NumIntegral_NT_op_dNi_dV<DIM>  hydrostatic_gravity( model.Database(), "element fluid density",
                                                                 "conductivity", "fluid pressure", acc_gravity );
  hydrostatic_pressure.Add( &hydrostatic_conductance );
  hydrostatic_pressure.Add( &hydrostatic_gravity );

  // "\nmain: Enter the amount of total dissolved solids (ppm = g/tonne; normal seawater=12000 g/t): ";
  double total_dissolved_solids(0.); // g->kg (157500-ppm = 157kg salt)
  total_dissolved_solids /= 1000.; // gets kg/m3
  printRangeOfVariable( model, "element fluid density" );
  PropertyHandle<DIM>  rhof( model, "element fluid density", SCALAR, ELEMENT );
  rhof.OutputCondition(ANY);
  rhof += total_dissolved_solids/1000.;
  printRangeOfVariable( model, "element fluid density" );
  vtk_output.OutputDataToVTK( model, "element-fluid-density", "element fluid density", 0, true );
  
  cout <<"\nmain: initial guess of fluid pressure.";
  printRangeOfVariable( model, "fluid pressure" );
/*
  cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid density and viscosity... " << endl;
  for ( uint32_t i=0; i<=5U; i++ ) {
       cout <<"\n\titeration "<< i+1U <<":"<< endl;
       model.Apply( hydrostatic_pressure );
//printRangeOfVariable( model, "fluid pressure" );
//vtk_output.OutputDataToVTK( model, "fluid pressure", "fluid pressure", 0, true );
       model.Accept( properties_visitor );
       for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
            (*it)->PropertyValueAtBaryCenter( visc_key, visc );
            (*it)->Store( cond_key, makeScalar( (*it)->Status(cond_key), (*it)->Read(perm_key) / visc() ) );
         }
       model.InterpolateNodeToElementProperty( "fluid density", "element fluid density" );
       rhof += total_dissolved_solids;
       printRangeOfVariable( model, "fluid pressure" );
       printRangeOfVariable( model, "element fluid density" );
    }
   vtk_output.OutputDataToVTK( model, "fluid pressure", "fluid pressure", 0, true );
*/



  // -----------------------------------------------------------
  // 4. initial rock effective stress due to gravitional loading
  // -----------------------------------------------------------
  // mechanical properties are placed on the element integration points
  // computes the increase in pore-pressure due to the gravitational loading
    SAMG_Settings  samg_mechanics_settings;
    samg_mechanics_settings.Set_napproach(2); // sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
    samg_mechanics_settings.SetSolverInstance(2);
                               // which is needed for deformation simulations
    SAMG_Solver                 samg_solver2(&samg_mechanics_settings);
    PDE_Integrator<DIM,Region>  deformation(samg_solver2);
    deformation.ScaleEssentialConditions(1.0e15);

    // this will also include boundary stresses translated into nodal forces
    const bool principal_vectors(true), plane_strain(false);
    PT_op<DIM>                  bforces( model.Database(), "force", "displacement" );
    NumIntegral_BT_D_B_dV<DIM>  stiffness( model.Database(),
                                          "Youngs modulus", "Poissons ratio", "displacement", "displacement", plane_strain);
    NumIntegral_PT_op_dV<DIM>   bodyforce( model.Database(), "gravity force", "displacement");
    NumIntegral_BT_D_op_dV<DIM> volstrain( model.Database(),
                                          "fluid volume source", "Youngs modulus", "Poissons ratio", "displacement");
    NumIntegral_BT_op_dV<DIM>   porepressure( model.Database(), "Biot term", "displacement");

    // computes the 'Biot term'
    porePressureBiotAlphaProduct( model, "Model" );
    cout <<"\nmain: taking into account (in the displacement equation), the effects of pore pressure.";
    printRangeOfVariable( model, "Biot alpha" );
    printRangeOfVariable( model, "Biot term" );
  
    // computing the gravity force
    gravityForce( model, acc_gravity );
    printRangeOfVariable( model, "gravity force" );
  
    deformation.Add( &porepressure );
    deformation.Add( &stiffness );
    deformation.Add( &bforces );  // force vector must always be there so that Dirichlet conditions are accumulated
    deformation.Add( &bodyforce );
    deformation.Add( &volstrain );
    deformation.Add( &porepressure );

    //StressesAndStrains<DIM>  postpro( model, "Youngs modulus", "Poissons ratio", "displacement", principal_vectors );
    StressesAndStrains<DIM>  postpro( model, "Youngs modulus", "Poissons ratio", "displacement", plane_strain, principal_vectors );
    deformation.AddPostProcess( &postpro );
 
    // -------------------------------------------
    // displacement, strain and stress computation
    // -------------------------------------------
    printRangeOfVariable( model, "Youngs modulus" );
    printRangeOfVariable( model, "Poissons ratio" );
    printRangeOfVariable( model, "gravity force" );
    deformation.IntegrateOver( model_domain );
    // among other things the stresses and strains operator computes the (FE-based) (elastic) dilatation
    printRangeOfVariable( model, "displacement" );
  
  
   // ---------------------------------------------------------------------------------------
   // extracting the vertical stress component
   // ---------------------------------------------------------------------------------------
    if ( DIM == 2 ) {
         ExtractTensorVariableComponent<DIM>  ystress(   model.Database(), "stress", "stress-y",  1,1 );
         ExtractTensorVariableComponent<DIM>  stress_xy( model.Database(), "stress", "stress-xy", 0,1 );
         model.Apply( ystress );
         model.Apply( stress_xy );
      }
    model.MoveNodeCoordinatesBy("displacement");
  
    vtk_output.OutputDataToVTK( model, "displacement", "displacement", 1, true );
    vtk_output.OutputDataToVTK( model, "strain",       "strain",       1, true );
    vtk_output.OutputDataToVTK( model, "stress",       "stress",       1, true );
    vtk_output.OutputDataToVTK( model, "mean-stress",  "mean stress",  1, true );
    vtk_output.OutputDataToVTK( model, "dilatation",   "dilatation",   1, true );
    vtk_output.OutputDataToVTK( model, "sigma1_",      "sigma1",       1, true );
    vtk_output.OutputDataToVTK( model, "sigma2_",      "sigma2",       1, true );
    vtk_output.OutputDataToVTK( model, "strain1_",     "strain1",      1, true );
    vtk_output.OutputDataToVTK( model, "strain2_",     "strain2",      1, true );
    vtk_output.OutputDataToVTK( model, "stress-y",     "stress-y",     1, true );
    vtk_output.OutputDataToVTK( model, "stress-xy",    "stress-xy",    1, true );
  
   // get the divergence
    dilatationInducedChangeInPorePressure( model );
    printRangeOfVariable( model, "fluid pressure" );
  
    SinglePhaseVelocityVisitor<DIM>  velo( model, "porosity", "conductivity",
                                          "fluid density", "fluid pressure", "total velocity" );
    model.Accept( velo );
    printRangeOfVariable( model, "total velocity" );

  // divergence of fluid flow to get volume strains (none yet)
  // ---------------------------------------------------------
  // loop over the finite volumes and measure the divergence of the facet fluxes
  // ignoring the FV at the model boundaries
   const bool second_order_in_space(false);
   const bool second_order_in_time(false);
   NodeCenteredFiniteVolumeTransport<DIM>  ncfvt( "Model", model, "porosity", "concentration", "total velocity",
                                                  "nodal fluid volume source",
                                                   second_order_in_space, second_order_in_time );
   // the divergence is assigned to the fvs
   ncfvt.Divergence( "total velocity", "nodal fluid volume source");

   model.MoveNodeCoordinatesBy("displacement");
  
   cout <<"\nSlopeMechanics_Example: That's it!\n";
  
} // end Run




// OTHER FUNCTION RELATIONSHIPS

/**
    Since fluid pressure is a node property we are looking
    for the dilatation at the node which is computed using the 
    finite volume framework.
 
    @note Notes

       1. Loading the rock skeleton results in an increase in pore pressure.
       2. the corresponding velocity field is no longer divergence free
       3. computing ensuing fluxes over the time-interval of interest, reveals
          how much fluid is leaving the rock
       4. the rock volume must be reduced accordingly, but without creating a stress response
*/
void dilatationInducedChangeInPorePressure( Model<DIM>& model )
 {
    const csmp::Index dil_key(model.Database().StorageKey("dilatation"));
    const csmp::Index pf_key(model.Database().StorageKey("fluid pressure"));
    const csmp::Index bf_key(model.Database().StorageKey("fluid compressibility"));
 
    Region<DIM>& model_domain(model.Region("Model"));
    ScalarVariable fluid_compressibility;

    // TODO: we need to get the nodal dilatation using the FVM discretisation
    // volume-weighted average of the dilatation of the FV sectors surrounding the node
   
    for ( auto it=model_domain.NodesBegin(); it!=model_domain.NodesEnd(); ++it )
      {
         // get the inputs
         double pf         = (*it)->Read( pf_key );
         double fluid_compressibility = (*it)->Read(bf_key);
         // finite volume average dilatation
         double fv_dilatation(0.), fv_volume(0.);
         for ( auto i{0}; i<(*it)->Parents(); ++i ) {
             // needed: sector volumes and elemental dilatation values
             Element<DIM>*  eptr((*it)->Parent(i));
             double sector_dilatation = eptr->Read( dil_key );
             size_t esector = (*it)->ParentNodeNumber(i);
             double sector_volume  = eptr->SectorVolume(esector);
             fv_dilatation += sector_dilatation * sector_volume;
             fv_volume     += sector_volume;
           }
         fv_dilatation /= fv_volume;
        
         // compute pressure change and store the new pore pressure
         (*it)->Store( pf_key, makeScalar((*it)->Status(pf_key), pf + fv_dilatation / fluid_compressibility) );
      }
 
 } // end dilatationInducedChangeInPorePressure




void permeabilityPorosityCorrelation( Model<DIM>& model )
 {
    const csmp::Index phi_key(model.Database().StorageKey("porosity"));
    const csmp::Index k_key(model.Database().StorageKey("permeability"));
 
    Region<DIM>& model_domain(model.Region("Model"));
   
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
      {
         double porosity = (*it)->Read( phi_key );
         // TODO: introduce proper relationship here
         double permeability = porosity * 1.0e-12;
         (*it)->Store( k_key, makeScalar((*it)->Status(k_key), permeability) );
      }
 
 } // end permeabilityPorosityCorrelation





/**
    Computes the product of the Biot coefficient alpha with the pore pressure 
    for assignment to the pre-stress term in the governing equation
    
    The result is returned into the scalar element variable 'Biot term'
    
    @attention the fluid pressure effect is averaged on the element barycenter
*/
void porePressureBiotAlphaProduct( Model<DIM>& model, const char* target_region )
 {
    const csmp::Index pf_key    = model.Database().StorageKey("fluid pressure");
    const csmp::Index alpha_key = model.Database().StorageKey("Biot alpha");
    assert( alpha_key.place == ELEMENT );
    assert( alpha_key.type  == SCALAR );
    // result
    const csmp::Index res_key = model.Database().StorageKey("Biot term");
    assert( res_key.place == ELEMENT );
    assert( res_key.type  == SCALAR );

    csmp::Region<DIM>& ref = model.Region(target_region);
    ScalarVariable sc;
    
    for ( auto it=ref.ElementsBegin(); it!=ref.ElementsEnd(); it++ )
      {
         const double alpha((*it)->Read( alpha_key ));
         if ( alpha < 0. or alpha > 1. ) {
              cerr <<"\n\tBiot coefficient alpha (1 - K_dry/K_grain): "<< alpha;
              throw csmp::Exception( ERROR, "porePressureBiotAlphaProduct:", "Biot coefficient alpha is out of range." );
           }
         (*it)->PropertyValueAtBaryCenter( pf_key, sc );
         sc *= alpha;
         (*it)->Store( res_key, makeScalar( (*it)->Status(res_key), sc() ) );
      }
      
 } // end porePressureBiotAlphaProduct





/**
    Computing the specific gravity force that acts on the rock skeleton from the "dry rock density"
*/
void gravityForce( Model<DIM>& model, double acc_gravity )
 {
    csmp::Region<DIM>& ref = model.Region("Model");

    const csmp::Index drd_key = model.Database().StorageKey("dry rock density");
    const csmp::Index gf_key  = model.Database().StorageKey("gravity force");
   
    VectorVariable<DIM> gforce(ANY,0. );
    model.InputPropertyValue( "gravity force", gforce );
   
    for ( auto it=ref.ElementsBegin(); it!=ref.ElementsEnd(); it++ )
      {
         (*it)->Read( gf_key, gforce );
         gforce(1) += (*it)->Read( drd_key ) * acc_gravity;
         (*it)->Store( gf_key, gforce );
      }
   
 } // end gravity force



} // csmp
