//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "Experimental_Example.h"
#include "CompareFloats.h"
#include "ANSYS_Model3D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "PDE_IntegratorExperimental.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include <iostream>
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
#include "NodeCenteredFiniteVolumeTransport.h"
#include "SinglePhaseVelocityVisitor.h"

using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
{
   SetTitle( "Experimental_Example" );
   SetDifficulty( 1 );
   SetCategory( "Software Functionality" );
   AddAuthor( "Irina Sin" );
   AddDescription( "source in: Experimental_Example.cpp" );
   AddDescription( "Empty example for the user to experiment with" );
   AddRequirement( "none" );
   AddRequirement( "no predefined model or variables file" );
}

// OTHER RELATIONSHIPS

/**
    Since fluid pressure is a node property we are looking
    for the dilatation at the node which is computed using the 
    finite volume framework.
*/
void dilatationInducedChangeInPorePressure( Model<3U>& model )
 {
    const csmp::Index dil_key(model.Database().StorageKey("dilatation"));
    const csmp::Index pf_key(model.Database().StorageKey("fluid pressure"));
    const csmp::Index bf_key(model.Database().StorageKey("fluid compressibility"));
 
    Region<3U>& model_domain(model.Region("Model"));
    ScalarVariable fluid_compressibility;

    // TODO: we need to get the nodal dilatation using the FVM discretisation
    // volume-weighted average of the dilatation of the FV sectors surrounding the node
   
    for ( auto it=model_domain.NodesBegin(); it!=model_domain.NodesEnd(); ++it )
      {
         // get the inputs
         double64 pf         = (*it)->Read( pf_key );
         double64 fluid_compressibility = (*it)->Read(bf_key);
         // finite volume average dilatation
         double64 fv_dilatation(0.), fv_volume(0.);
         for ( size_t i=0U; i<(*it)->Parents(); ++i ) {
             // needed: sector volumes and elemental dilatation values
             Element<3U>*  eptr((*it)->Parent(i));
             double64 sector_dilatation = eptr->Read( dil_key );
             size_t esector = (*it)->ParentNodeNumber(i);
             double64 sector_volume  = eptr->SectorVolume(esector);
             fv_dilatation += sector_dilatation * sector_volume;
             fv_volume     += sector_volume;
           }
         fv_dilatation /= fv_volume;
        
         // compute pressure change and store the new pore pressure
         (*it)->Store( pf_key, makeScalar((*it)->Status(pf_key), pf + fv_dilatation / fluid_compressibility) );
      }
 
 } // end dilatationInducedChangeInPorePressure




void permeabilityPorosityCorrelation( Model<3U>& model )
 {
    const csmp::Index phi_key(model.Database().StorageKey("porosity"));
    const csmp::Index k_key(model.Database().StorageKey("permeability"));
 
    Region<3U>& model_domain(model.Region("Model"));
   
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
      {
         double64 porosity = (*it)->Read( phi_key );
         // TODO: introduce proper relationship here
         double64 permeability = porosity * 1.0e-12;
         (*it)->Store( k_key, makeScalar((*it)->Status(k_key), permeability) );
      }
 
 } // end permeabilityPorosityCorrelation





/**
    Computes the product of the Biot coefficient alpha with the pore pressure 
    for assignment to the pre-stress term in the governing equation
    
    The result is returned into the scalar element variable 'Biot term'
    
    @attention the fluid pressure effect is averaged on the element barycenter
*/
template<size_t dim>
void porePressureBiotAlphaProduct( Model<dim>& model, const char* target_region )
 {
    const csmp::Index pf_key    = model.Database().StorageKey("fluid pressure");
    const csmp::Index alpha_key = model.Database().StorageKey("Biot alpha");
    assert( alpha_key.place == ELEMENT );
    assert( alpha_key.type  == SCALAR );
    // result
    const csmp::Index res_key = model.Database().StorageKey("Biot term");
    assert( res_key.place == ELEMENT );
    assert( res_key.type  == SCALAR );

    csmp::Region<dim>& ref = model.Region(target_region);
    ScalarVariable sc;
    
    for ( typename vector<Element<dim>*>::iterator
          it=ref.ElementsBegin(); it!=ref.ElementsEnd(); it++ )
      {
         const double64 alpha((*it)->Read( alpha_key ));
         if ( alpha < 0. or alpha > 1. ) {
              cerr <<"\n\tBiot coefficient alpha (1 - K_dry/K_grain): "<< alpha;
              throw csmp::Exception( ERROR, "porePressureBiotAlphaProduct:", "Biot coefficient alpha is out of range." );
           }
         (*it)->PropertyValueAtBaryCenter( pf_key, sc );
         sc *= alpha;
         (*it)->Store( res_key, makeScalar( (*it)->Status(res_key), sc() ) );
      }
      
 } // end porePressureBiotAlphaProduct

template void porePressureBiotAlphaProduct( Model<3U>&, const char* );



/**
    Computing the specific gravity force that acts on the rock skeleton from the "dry rock density"
*/
void gravityForce( Model<3U>& model, double64 acc_gravity )
 {
    csmp::Region<3U>& ref = model.Region("Model");

    const csmp::Index drd_key = model.Database().StorageKey("dry rock density");
    const csmp::Index gf_key  = model.Database().StorageKey("gravity force");
   
    VectorVariable<3U> gforce(ANY,ANY,ANY, 0., 0., 0. );
    model.InputPropertyValue( "gravity force", gforce );
   
    for ( typename vector<Element<3U>*>::iterator
          it=ref.ElementsBegin(); it!=ref.ElementsEnd(); it++ )
      {
         (*it)->Read( gf_key, gforce );
         gforce(1) += (*it)->Read( drd_key ) * acc_gravity;
         (*it)->Store( gf_key, gforce );
      }
   
 } // end gravity force







/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.
*/
void Experimental_Example::Run()
{
  // compaction model, draft 1
  
  // 1. build model
  const size_t   DIM(3U);
  const string   model_name("3D_box2");
  ANSYS_Model3D  model( model_name.c_str(), "CSMP_field_scale_mechanics_variables.txt");
  Region<3U>&    model_domain(model.Region("Model"));

  printModelDimensions( model, true );

  // 2. configure it from file
  InputDataManager<3U>  model_configuration;
  ComputationalSettings settings;

  model_configuration.ConfigureFromFile( model, "3D_box2",
                                         false,           ///< region name from parameter range
                                         true,            ///< default property values
                                         true,            ///< regional property values
                                         true,            ///< boundary conditions for box-shaped model
                                         false,           ///< essential conditions for regions
                                         false,           ///< boundary conditions for arbitrary-shaped model
                                         settings );

  VTK_Interface<3U>  vtk_output;
  
  // --------------------------------------------------------
  // 3. compute initial steady state temperature distribution
  // --------------------------------------------------------
  {
    SteadyStateDiffusor<3U,Region>  temperature( model, "thermal conductivity", "temperature", "heat source");
    printRangeOfVariable( model, "thermal conductivity" );
    printRangeOfVariable( model, "heat source");
    //model.InputBoundaryValue( TOP, "temperature", makeScalar(DIRICH,20.));
    //model.InputBoundaryValue( BOTTOM, "temperature", makeScalar(DIRICH,80.));
    temperature.ComputeSteadyState( model );
    printRangeOfVariable( model, "temperature");
  }
  
  // -------------------------------------------------------
  // 4. initial hydrostatic pressure (no loading via grain skeleton)
  // --------------------------------------------------------
  //  Visitor computes thermodynamic properties of water; density is interpolated to
  //     element barycentre for later vertical integration
  // --------------------------------------------------------------------------------------
  IAPWS_H2OPropertiesVisitor<DIM>  properties_visitor( model, "fluid pressure" ,"fluid density", "fluid viscosity");
  // computing fluid density and viscosity
  model.Accept( properties_visitor );
  printRangeOfVariable( model, "fluid density");
  printRangeOfVariable( model, "fluid viscosity");

  model.InterpolateNodeToElementProperty( "fluid density", "element fluid density" );

  // Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
  SAMG_Settings  samg_settings;
  SAMG_Solver    samg_solver(&samg_settings);
  PDE_Integrator<DIM,Region>  hydrostatic_pressure(samg_solver);
  // minimizing screen output
  samg_settings.Set_iout1( 0 );
  samg_settings.Set_iout2( 0 );
  samg_settings.Set_idmp( -1 );

  NumIntegral_dNT_op_dN_dV<DIM>  hydrostatic_conductance( model.Database(),
                                                         "conductivity",
                                                         "fluid pressure", "fluid pressure" );

  // acceleration of gravity (kg/m.s2): in the area of interest: ";
  double64 acc_gravity(9.81);
  // unless specified otherwise, in a 1D model, gravity will automatically act in the x-direction
  NumIntegral_NT_op_dNi_dV<DIM>  hydrostatic_gravity( model.Database(), "element fluid density",
                                                                 "conductivity", "fluid pressure", acc_gravity );
  hydrostatic_pressure.Add( &hydrostatic_conductance );
  hydrostatic_pressure.Add( &hydrostatic_gravity );

  // "\nmain: Enter the amount of total dissolved solids (ppm = g/tonne; normal seawater=12000 g/t): ";
  double64 total_dissolved_solids(12000./1000.); // g->kg (157500-ppm = 157kg salt)
  total_dissolved_solids /= 1000.; // gets kg/m3
  printRangeOfVariable( model, "element fluid density" );
  PropertyHandle<DIM>  rhof( model, "element fluid density", SCALAR, ELEMENT );
  rhof.OutputCondition(ANY);
  rhof += total_dissolved_solids/1000.;
  printRangeOfVariable( model, "element fluid density" );

  cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid properties... " << endl;
  for ( uint32 i=0; i<=5U; i++ ) {
       cout <<"\n\titeration "<< i+1U <<":"<< endl;
       model.Apply( hydrostatic_pressure );
       model.Accept( properties_visitor );
       model.InterpolateNodeToElementProperty( "fluid density", "element fluid density" );
       rhof += total_dissolved_solids;
       printRangeOfVariable( model, "fluid pressure" );
       printRangeOfVariable( model, "element fluid density" );
    }

  // -----------------------------------------------------------
  // 5. initial rock effective stress due to gravitional loading
  // -----------------------------------------------------------
  // computes the increase in pore-pressure due to the gravitational loading
    SAMG_Settings  samg_mechanics_settings;
    samg_mechanics_settings.Set_napproach(2); // sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
    samg_mechanics_settings.SetSolverInstance(2);
                               // which is needed for deformation simulations
    SAMG_Solver                 samg_solver2(&samg_mechanics_settings);
    PDE_Integrator<DIM,Region>  deformation(samg_solver2);
    deformation.ScaleEssentialConditions(1.0e15);

    // this will also include boundary stresses translated into nodal forces
    PT_op<DIM>                  bforces( model.Database(), "force", "displacement" );
    NumIntegral_BT_D_B_dV<DIM>  stiffness( model.Database(),
                                          "Youngs modulus", "Poissons ratio", "displacement", "displacement");
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

    const bool principal_vectors(true);
    StressesAndStrains<DIM>  postpro( model, "Youngs modulus", "Poissons ratio", "displacement", principal_vectors );
    deformation.AddPostProcess( &postpro );
 
    // computation
    deformation.IntegrateOver( model_domain );
    // among other things the stresses and strains operator computes the (FE-based) (elastic) dilatation
    printRangeOfVariable( model, "displacement" );
  
    // from the dilatation (compression) we can compute the change in pore pressure, fluid flow and then
    model.InstantiateFiniteVolumes();
    // get the divergence
    dilatationInducedChangeInPorePressure( model );
    printRangeOfVariable( model, "fluid pressure" );
  
    SinglePhaseVelocityVisitor<DIM>  velo( model, "porosity", "conductivity",
                                          "fluid density", "fluid pressure", "total velocity" );
    model.Accept( velo );
    printRangeOfVariable( model, "total velocity" );

  // divergence of fluid flow to get volume strains (none yet)
  // --------------------------------------------------------
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
  
   // NOTE:
   /*
       1. Loading the rock skeleton results in an increase in pore pressure.
       2. the corresponding velocity field is no longer divergence free
       3. computing ensuing fluxes over the time-interval of interest, reveals
          how much fluid is leaving the rock
       4. the rock volume must be reduced accordingly, but without creating a stress response
   */



  // -------------------------------------------------------
  // 6. UNDRAINED SPLIT COMPACTION COMPUTATION
  // --------------------------------------------------------




  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run



} // csmp
