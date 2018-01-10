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
#include "ANSYS_Model2D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"
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

#include "Gauss_Solver.h"

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
#include "StressesAndStrains2.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "SinglePhaseVelocityVisitor.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;

namespace csmp {

const size_t DIM(2U);

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

template void porePressureBiotAlphaProduct( Model<DIM>&, const char* );



/**
    Computing the specific gravity force that acts on the rock skeleton from the "dry rock density"
*/
void gravityForce( Model<DIM>& model, double64 acc_gravity )
 {
    csmp::Region<DIM>& ref = model.Region("Model");

    const csmp::Index drd_key = model.Database().StorageKey("dry rock density");
    const csmp::Index gf_key  = model.Database().StorageKey("gravity force");
   
    VectorVariable<DIM> gforce(ANY,0. );
    model.InputPropertyValue( "gravity force", gforce );
   
    for ( typename vector<Element<DIM>*>::iterator
          it=ref.ElementsBegin(); it!=ref.ElementsEnd(); it++ )
      {
         (*it)->Read( gf_key, gforce );
         // minus (-) because gravity is acting in the opposite direction of the Y axis
         gforce(1) += (*it)->Read( drd_key ) * -acc_gravity;
         (*it)->Store( gf_key, gforce );
      }
   
 } // end gravity force




/** 
    Computation of gravity-equilibrated fixed pressure at vertical boundaries of 2D models,
    by top down piecewise integration. The new pressure values are flagged Dirichlet.
*/
void computeAndFreezeBoundaryPressure( Model<DIM>& model, const char* edge_boundary, double64 pressure_at_top, double64 acc_gravity )
 {
    const Boundary<DIM>& bref(model.Boundary(edge_boundary));
   
    // 1. making a map of elements on the basis of the Y coordinates
    map<double64,Node<DIM>*,greater<double64> > nodes_map;
    set<double64> pressure_range;

    // ordering the boundary nodes vertically, largest Y-value first
    for ( vector<Node<DIM>*>::const_iterator nit=bref.NodesBegin(); nit!=bref.NodesEnd(); ++nit )
      nodes_map.insert( make_pair( (*nit)->y(), (*nit) ) );

    // 2. performing the vertical integration approximating fluid density as piecewise averages
    //    on the elements (trapezoidal rule), top down
    const csmp::Index rho_key(model.Database().StorageKey("fluid density"));
    const csmp::Index pf_key(model.Database().StorageKey("fluid pressure"));
    map<double64,Node<DIM>*,greater<double64> >::iterator nit(nodes_map.begin());
    double64 density1((*nit).second->Read(rho_key)), density2;
    double64 pressure1(pressure_at_top), elevation((*nit).first);
    pressure_range.insert( pressure1 );
    nit++;
   
    while ( nit != nodes_map.end() ) {
         density2 = (*nit).second->Read(rho_key);
         double64 rho_avg = (density1 + density2) / 2.;
         double64 pressure2 = pressure1 + acc_gravity * rho_avg * (elevation - (*nit).first);
         (*nit).second->Store(pf_key, makeScalar(DIRICH,pressure2) );
         pressure_range.insert( pressure2 );
         pressure1 = pressure2;
         density1  = density2;
         elevation = (*nit).first;
         nit++;
      }
 
    cout <<"\ncomputeAndFreezeBoundaryPressure: new pressure range at boundary '"<< edge_boundary <<"' ";
    cout << (*pressure_range.begin()) <<" to "<< (*pressure_range.rbegin()) <<" Pa.\n";
   
 } // ComputeAndFixLateralBoundaryPressure













/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.
     
     Systems Modelling and Design model that calculates the effective stress in a dam, 
     using the gravitational loading and the plane stress assumption.
     
     input models:
     - slope_model1
*/
void Experimental_Example::Run()
{
  // ----------------------------------------------------------
  // 1. building and configering model from ANSYS - csp dataset
  // ----------------------------------------------------------
   string input_file("Jura-slope1");
   ANSYS_Model2D   model( input_file.c_str(), "CSMP_field_scale_mechanics_variables.txt" );
  // ANSYS_Model3D  model( input_file.c_str(), "CSMP_field_scale_mechanics_variables.txt");
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
       SAMG_Settings  samg_settings;
       SAMG_Solver    samg_solver(&samg_settings);
       //Gauss_Solver  linear_solver;
       PDE_Integrator<2U,Region>  temperature(samg_solver);
       NumIntegral_dNT_op_dN_dV<2U,Element<2U> >  conductance( model.Database(), "thermal conductivity", "temperature",  "temperature" );
       NumIntegral_NT_op_N_dV<2U,Element<2U> >    source( model.Database(),  "energy source", "temperature" );
       temperature.Add( &conductance );
       temperature.Add( &source );
       printRangeOfVariable( model, "thermal conductivity" );
       printRangeOfVariable( model, "energy source");
       model.Apply( temperature );
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
  vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0, true );
  
  // Setting up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
  SAMG_Settings  samg_settings;
  
  SAMG_Solver    samg_solver(&samg_settings);
  PDE_Integrator<DIM,Region>  hydrostatic_pressure(samg_solver);
  
  // minimizing screen output
  samg_settings.Set_iout1( 0 );
  samg_settings.Set_iout2( 0 );
  samg_settings.Set_idmp( -1 );

  samg_settings.Set_iswit(5);
  // SAMG solution criterion
  samg_settings.Set_eps(0.);
  // Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
//    bpressure.GetSolverSettings().Set_nred(1);
  // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
  samg_settings.Set_a_cmplx(2);
  // Pre-adjust SAMG mesh complexity, based on solver output
  samg_settings.Set_g_cmplx(1.5);
  samg_settings.Set_w_avrge(2);

  samg_settings.SetSolverInstance(2);
  
//  Gauss_Solver  linear_solver;
//  PDE_Integrator<2U,Region>  hydrostatic_pressure(linear_solver);

  NumIntegral_dNT_op_dN_dV<DIM>  hydrostatic_conductance( model.Database(),
                                                         "conductivity",
                                                         "fluid pressure", "fluid pressure" );
//  NumIntegral_dNT_dN_dV<DIM>  hydrostatic_conductance( model.Database(), "fluid pressure", "fluid pressure" );

  // cin >> "\nmain: Enter the acceleration of gravity (kg/m.s2): in the area of interest: ";
  // reading this from the model
  const csmp::Index g_key(model.Database().StorageKey("acceleration gravity"));
  const double64 acc_gravity = model.Read( g_key );
  // unless specified otherwise, gravity will automatically act in the opposite direction of the Y axis
//  NumIntegral_NT_op_dNi_dV<DIM>  hydrostatic_gravity( model.Database(), "fluid density", "one", "fluid pressure", acc_gravity );
  NumIntegral_NT_op_dNi_dV<DIM>  hydrostatic_gravity( model.Database(), "fluid density", "conductivity", "fluid pressure", acc_gravity );
  hydrostatic_pressure.Add( &hydrostatic_conductance );
  hydrostatic_pressure.Add( &hydrostatic_gravity );

  // "\nmain: Enter the amount of total dissolved solids (ppm = g/tonne; normal seawater=12000 g/t): ";
  double64 total_dissolved_solids(0.); // g->kg (157500-ppm = 157kg salt)
  total_dissolved_solids /= 1000.; // gets kg/m3
  printRangeOfVariable( model, "fluid density" );
  PropertyHandle<DIM>  rhof( model, "fluid density", SCALAR, NODE );
  rhof.OutputCondition(ANY);
  rhof += total_dissolved_solids/1000.;
  printRangeOfVariable( model, "fluid density" );
  vtk_output.OutputDataToVTK( model, "fluid-density", "fluid density", 0, true );
  
  cout <<"\nmain: initial guess of fluid pressure.";
  const double64 p_atm(100325.);
  computeAndFreezeBoundaryPressure( model, "RIGHT", p_atm, acc_gravity );
  printRangeOfVariable( model, "fluid pressure" );

  cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid density and viscosity... " << endl;
  for ( uint32 i=0; i<=3U; i++ ) {
       cout <<"\n\titeration "<< i+1U <<":"<< endl;
       model.Apply( hydrostatic_pressure );
       printRangeOfVariable( model, "fluid pressure" );
//vtk_output.OutputDataToVTK( model, "fluid pressure", "fluid pressure", i, true );
       model.Accept( properties_visitor );
       computeAndFreezeBoundaryPressure( model, "RIGHT", p_atm, acc_gravity );
       for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
            (*it)->PropertyValueAtBaryCenter( visc_key, visc );
            (*it)->Store( cond_key, makeScalar( (*it)->Status(cond_key), (*it)->Read(perm_key) / visc() ) );
         }
       rhof += total_dissolved_solids;
       printRangeOfVariable( model, "fluid density" );
    }
   vtk_output.OutputDataToVTK( model, "fluid pressure", "fluid pressure", 0, true );
  
  // -----------------------------------------------------------
  // 4. computing the Darcy velocity
  // -----------------------------------------------------------
   assignNodeCoordinatesTo( model, 'Y', "elevation" );
   vtk_output.OutputDataToVTK( model, "elevation", "elevation", 0, true );
  // 4.1 computing the hydraulic head, h = pf/rho*g + z   from the fluid pressure distribution
    const csmp::Index pf_key(model.Database().StorageKey("fluid pressure"));
    const csmp::Index rhof_key(model.Database().StorageKey("fluid density"));
    const csmp::Index E_key(model.Database().StorageKey("elevation"));
    const csmp::Index h_key(model.Database().StorageKey("hydraulic head"));
    const csmp::Index K_key(model.Database().StorageKey("conductivity"));
    const csmp::Index vD_key(model.Database().StorageKey("Darcy velocity"));

    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      {
          double64 head = (*nit)->Read(pf_key) / ((*nit)->Read(rhof_key) * acc_gravity) + (*nit)->Read(E_key);
         (*nit)->Store( h_key, makeScalar( (*nit)->Status(h_key), head ) );
      }
    printRangeOfVariable( model, "hydraulic head" );
    vtk_output.OutputDataToVTK( model, "hydraulic head", "hydraulic head", 0, true );

    // 4.2 computing hydraulic head, its gradients and the ensuing Darcy velocity
    DenseMatrix<DM_MIN> DERIV;
    VectorVariable<DIM> velo;
    // 4.3 computing ensuing fluid flow
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
      {
         const double64 K = (*it)->Read( K_key );
         velo = 0.;
         (*it)->dN_AtBaryCenter( DERIV, 1U );
         for ( size_t i=0U; i<(*it)->Nodes(); i++ ) {
              double64 h = (*it)->N(i)->Read( h_key );
              velo(0)  += h  * -DERIV(0,i) * K;
              velo(1)  += h  * -DERIV(1,i) * K;
              velo.Flag(0) = (*it)->Status(vD_key,0);
              velo.Flag(1) = (*it)->Status(vD_key,1);
           }
        
         (*it)->Store( vD_key, velo );
      }
    printRangeOfVariable( model, "Darcy velocity" );
    vtk_output.OutputDataToVTK( model, "Darcy-velocity", "Darcy velocity", 0, true );
    model.ExtrapolateElementToNodeProperty( "Darcy velocity", "Darcy velocity node" );
    printRangeOfVariable( model, "Darcy velocity node" );
    vtk_output.OutputDataToVTK( model, "Darcy-velocity-node", "Darcy velocity node", 0, true );


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
                                          "Youngs modulus", "Poissons ratio", "displacement", "displacement", plane_strain );
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
    deformation.IntegrateOver( model_domain );
    // among other things the stresses and strains operator computes the (FE-based) (elastic) dilatation
    printRangeOfVariable( model, "displacement" );
  
  
   // ---------------------------------------------------------------------------------------
   // extracting the vertical stress component
   // ---------------------------------------------------------------------------------------
    ExtractTensorVariableComponent<DIM>  ystress(   model.Database(), "stress", "stress-y",  1,1 );
    ExtractTensorVariableComponent<DIM>  stress_xy( model.Database(), "stress", "stress-xy", 0,1 );
    model.Apply( ystress );
    model.Apply( stress_xy );
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

  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run

} // csmp
