//
//  EffectiveStressDilatation2D_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/24/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "EffectiveStressDilatation2D_Example.h"

#include "Model.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "ModelTime.h"
#include "PDE_Integrator_CRM.h"
#include "TransientDiffusor.h"
#include "ConstantFactor.h"

// interfaces
#include "ANSYS_Interface.h"
#include "VSetConverter.h"
#include "VTK_Interface.h"

// algebraic multigrid solvers
#include "LinearSolver.h"
// fluid pressure algorithm and velocity computation
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

// transport scheme
#include "NodeCenteredFiniteVolumeTransport.h"

// utilities
#include "CSMP_highLevelUtilities.h"
#include "parallelPlatePermeabilityFromChannelWidth.h"
#include "RegionBoundaryFluxVisitor.h"

using namespace std;

namespace csmp {

// setting the dimensionality of the model
const size_t DIM(2U);

bool createLowerDimensionalRegion( Model<DIM>& model, const char* name_of_new_region );
void computeGasFlowProperties( Model<DIM>& model, const char* target_region );


EffectiveStressDilatation2D_Example::EffectiveStressDilatation2D_Example()
 : fluid_pressure_(0), verbose_(true)
{
}


void EffectiveStressDilatation2D_Example::Specifications()
{
  SetTitle( "Pressure diffusion through fracture network in low permeability brittle shale." );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "Two regions: 2 fracture sets SET1 SET and rock  MATRIX, steady state pf and velo output to VTK" );
  AddDescription( "Quadratic finite element interpolation, Dirichlet essential conditions" );
  AddDescription( "Region saved to a binary vset, then a second Region is built from it" );
  AddDescription( "source in: EffectiveStressDilatation2D_Example.cpp" );
  AddRequirement( "file set: 'fractured_slate'");
}


/** *****************************************************************************************
 
   Diffusion out of ultra-low permeability shale through a fracture network.

   Exercise: on the basis of the code presented herein and the linear elasticity example, implement:

   1. model of over-pressured fractured layer (is there a gas generation rate?)
      - no flow boundary conditions, fixed initial pressure

   2. computation of total system compressibility, ct taking into account gas compressibility
   
   3. simulation of the evolution of the bhp and reservoir pressure over time
   
   4. visualisation with Paraview

   5. well: flow rate monitoring for constant bh pressure
   
   6. introduce split boundaries to do mechanics

   suppporting functionality to illustrate:
   - implement: gas EOS
   - Knudson diffusion in rock matrix
   - implement: conversion of linear bar elements into quadratic bar elements
   - split boundaries
   - mechanics (change property placement to element integration points)
 
******************************************************************************************* */
void EffectiveStressDilatation2D_Example::Run()
{
  // ___________________________________________________________________________
  //
  // MODEL CREATION & PROPERTY ASSIGNMENT
  // ___________________________________________________________________________
  
  // name of input file set
  const char* bin_file="fractured_slate";
  
  // 1.0 create an input interface for ANSYS 2D meshes
  // ---------------------------------------------------------------------------
  ANSYS_Interface    mesh_interface(true); // true = isoparametric elements
  ModelTopology      mesh_topology(true);
  VSet<DIM>          mesh_container;

  // reading mesh
  const bool binary_file( true );
  const bool irregular_mesh( false );
  mesh_interface.Read_ANSYS_Mesh( bin_file, mesh_container, mesh_topology, binary_file, irregular_mesh );
  // mesh_topology.Out();
  // eliminating potentially unwanted regions
  mesh_topology.ReduceToRegions( bin_file );
  map<size_t,size_t>  old_and_new_elmtids;
  mesh_topology.CreateNewElementNumbers( old_and_new_elmtids );
  mesh_container.ReduceTo( old_and_new_elmtids );
  old_and_new_elmtids.clear();

// Optional: convert elements with linear- basis functions into ones with quadratic interpolation ones
//    VSetConverter<DIM>  mesh_converter;
//    mesh_converter.ConvertLinearToQuadraticTriangles( mesh_container );


  // 2.0 Building the "Model" with associated property storage (see variables file: example2.txt
  // -------------------------------------------------------------------------------------------
  Model<DIM>  model( mesh_topology, mesh_container, "EffectiveStressDilatation2D_Example_variables.txt" );
  mesh_container.Erase();
  mesh_topology.Erase();
  
  // 2.1 split Boundaries
  model.InsertSplitBoundary("SET1");
  model.InsertSplitBoundary("SET2");
  bool delete_underlying_region(true);
  model.InsertSplitBoundary("WELL_FRACTURE");
  
  cout <<"\nmain: created split boundaries: ";
  for ( Model<2>::splitBoundaryConstIterator it=model.SplitBoundariesBegin(); it!=model.SplitBoundariesEnd(); it++ )
    cout << (*it).first <<" ";
  cout <<"\n";


  // 3.0 Assigning material properties and initial valiues to the model and its subregions
  // -------------------------------------------------------------  300 nD  --------------
  // the names of the subregions are listed in *.asc file
  model.InputPropertyValue( "permeability",        makeScalar(PLAIN,3.0e-19) );
  model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.0) );
  model.InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN,0.0) );
  model.InputPropertyValue( "porosity",            makeScalar(PLAIN,0.1) );
  model.InputPropertyValue( "viscosity gas",       makeScalar(PLAIN,1.0e-5) );
  model.InputPropertyValue( "rock compressibility", makeScalar(PLAIN,1.0e-11) );
  model.InputPropertyValue( "gas compressibility",  makeScalar(PLAIN,1.0e-6) );
  // getting access to the regions via references
  Region<DIM>& set1 = model.Region("SET1");
  Region<DIM>& set2 = model.Region("SET2");
  SplitBoundary<DIM>& well = model.SplitBoundary("WELL_FRACTURE");

  // parallel plate law, ~1-mm aperture
  const double frac_perm((pow(1.0e-3,2)/12.) * 1.0e-3);
  set1.InputPropertyValue( "permeability",  makeScalar(PLAIN,frac_perm*1.) );
  set2.InputPropertyValue( "permeability",  makeScalar(PLAIN,frac_perm*2.) );

  // 'ct' must be computed from water, gas, and rock compressibility
  computeGasFlowProperties( model, "Model" );
  printRangeOfVariable( model, "conductivity" );


  // 4. Initial uniform fluid pressure in reservoir layer
  // ---------------------------------------------------------------
  // hydrostatic pressure at 2300 m depth
  const double64 pf(2300. * 9.81 * 1000. + 100325.);
  model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,pf) );
  // gas static pressure in the well
  const double64 pg(2300. * 9.81 * 100. + 100325.);
  well.InputNodePropertyValue( "fluid pressure", makeScalar(PLAIN,pg), COMPLETE, INSIDE );
  well.ChangePropertyStatus( "fluid pressure", DIRICH );
  well.InputPropertyValue( "fracture permeability", makeScalar(PLAIN,1.0e-7) );
  well.InputPropertyValue( "fracture porosity",     makeScalar(PLAIN,1.0) );


  // 5. Form a Region named "fractures" from all lower-dimensional elements
  // ----------------------------------------------------------------------
  createLowerDimensionalRegion( model, "fractures" );

  printRangeOfVariable( model, "permeability", true );

  VTK_Interface<DIM>  vtk_output;
  vtk_output.OutputDataToVTK( model, "fractures", "permeability", "permeability", 0, true );

  // output of splitboundary properties
  vtk_output.OutputDataToVTK( model, "WELL", "fracperm", "fracture permeability", 0, true );
  

    // 7. Output the results to VTK files
    // ---------------------------------------------------------------------------
    vtk_output.OutputDataToVTK( model, "permeability",   "permeability",      0, true );
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure",    0, true );
    vtk_output.OutputDataToVTK( model, "velocity",       "velocity",          0, true );
    vtk_output.OutputDataToVTK( model, "nvelocity",      "nodal velocity",    0, true );
    vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",       0, true );
    vtk_output.OutputDataToVTK( model, "nvolume-flux",   "nodal volume flux", 0, true );


  // ___________________________________________________________________________
  //
  //   COMPUTATION OF TRANSIENT FLUID PRESSURE
  // ___________________________________________________________________________

    // watch out with the time unit because spurious oscillations may result...
    const double64 well_life(86400000.), time_unit(600.); // 1000 days vs. 10-min
    double64       time_increment(60.);       // 60 sec to start with
    double64&      model_time( ModelTime::Instance().modelTime );

    // to invoke the finite-volume framework a transport algorithm is constructed
    NodeCenteredFiniteVolumeTransport<DIM>  solute_transport( "fractures", model, "porosity", "concentration",
                                                              "velocity", "nodal fluid volume source", false, false );

    RegionBoundaryFluxVisitor<DIM,Region> well_influx( model, "WELL_FRACTURE", "velocity" );
    cout <<"\nRun: length of the well (m): "<< well.Area() <<"\n";
    double64 cumulative_production(0.); // m3 at well pressure
  
    while ( model_time <= well_life )
      {
         model_time += time_increment;
         ComputeTransientFluidPressure( model, "WELL_FRACTURE", time_increment ); // timestep report to file in hours
         well_influx.TimeIncrement(time_increment);
         well.Accept(well_influx);
         cout <<"\nRun: Cumulative flux into the well (t="<< model_time <<"): "<< well_influx.InFlux() <<" (m3).\n";
         cumulative_production = well_influx.InFlux() * time_increment;
         vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", static_cast<size_t>(rint(model_time/time_unit)) );
         time_increment *= 1.2;
      }
  
    cout <<"\nRun: a total of "<< cumulative_production <<" (m3), was produced after "<< model_time/86400. <<" days\n";
  
} // end Run







/**
    Slightly compressible formulation for pressure equation, using total system compressibility-based hydraulic diffusivity.
*/
void EffectiveStressDilatation2D_Example::ComputeTransientFluidPressure( Model<DIM>& model, double64 time_increment )
 {
    static bool first_call(true);
    
     if ( first_call ) {
         fluid_pressure_ = new TransientDiffusor<DIM,Region>( model,
                                                           "conductivity",
                                                           "fluid pressure", "total system compressibility",
                                                           "fluid volume source" );
       
         // post-processing operation to compute flow velocities
         static VelocityAndVolumeFlux<DIM,Element<DIM> >  velocity( model, "conductivity", "porosity", "fluid pressure", false );
         fluid_pressure_->AddPostProcess( &velocity );

         #ifdef CSMP_WITH_SAMG_SOLVER
         // targeting SAMG DLL 2 for this pressure solver (same as for steady state solver, but now reusing solution from previous timestep if available)
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().SetSolverInstance(2);
         // io
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iout1( 0 );
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iout2( 0 );
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_idmp( -1 );
         // SAMG solution criteria
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_eps(0.);
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_rel_eps(1.E-10);
         
         // Re-use solver setup from previous timestep, internal checks force setup when required
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iswit(7);

         // Re-use solver setup from previous timestep, internal checks force setup when required
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iswit(7);
         // Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
      //    dynamic_cast<SteadyStateDiffusionSolver<3U,Region>*>(fluid_pressure_)->GetSolverSettings().GetSolverSettings().Set_nred(1);
         // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_a_cmplx(2);
         // Pre-adjust SAMG mesh complexity, based on solver output
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_g_cmplx(1.5);
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_w_avrge(2);
         #else
         /// add extra functionality for alternative solver if needed
         #endif
       }

    double64& model_time( ModelTime::Instance().modelTime );

#ifdef SAMG_OUTPUT_TO_FILE
    // trigger output to file
    dynamic_cast<TransientDiffusionSolver<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_idmp( 8 );
    // // define SAMG file output format for reduced file size, idmp > 1 is required
    dynamic_cast<TransientDiffusionSolver<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_ioform( "f" );
    // set filename for SAMG file output other than default "level", idmp > 1 is required
    string currentDumpFileName( "transient_pf" );
    currentDumpFileName.append( numberToString( static_cast<size_t>( model_time ) ) );
    dynamic_cast<TransientDiffusionSolver<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_filnam_dump( currentDumpFileName );
#endif

    if ( verbose_ ) {
         cout <<"\nComputeTransientFluidPressure: Input parameter ranges: "<< endl;
         printRangeOfVariable( model, "conductivity", true );
         printRangeOfVariable( model, "fluid volume source", true );
         printRangeOfVariable( model, "total system compressibility", true );
      }
    
    cout <<"\n\n\n\nComputeTransientFluidPressure: ";
    cout <<"Computing transient '"<< "fluid pressure" <<"' at t="<< fixed << setprecision(0) << model_time <<" secs.";
    cout <<"\n\t\t(FD Backward-Euler time increment, delta_t="<< time_increment <<" secs)."<< endl;
    fluid_pressure_->TimeIncrement( 1. / time_increment );

    model.Apply( *fluid_pressure_ );

    // setting relative solution tolerance after first solve
    if ( first_call )  {
        #ifdef CSMP_WITH_SAMG_SOLVER
         // suppress verbose output
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iout1(-1);
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iout2(-1);
         // absolute tolerance set to 1e-14
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_eps( 1.E-14 );
         // use solution from last step as an initial guess
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_itypu(0);
         // reuse solver setup from last call
         dynamic_cast<TransientDiffusor<DIM,Region>*>(fluid_pressure_)->GetSolverSettings().Set_iswit(3);
        #else
        /// add extra functionality for alternative solver if needed
        #endif
      }
  
//    if ( verbose_ ) printRangeOfVariable( model, stdio_, "fluid pressure", true );
   
    first_call = false;     
   
 } // ComputeTransientFluidPressure







/**
    For SplitBoundary computations.
*/
void EffectiveStressDilatation2D_Example::ComputeTransientFluidPressure( Model<2U>& model, const char* split_boundary_name, double64 time_increment )
 {
    static bool first_call(true);


    // ([C] + dt[K]){p}t+dt = [C]{p}t + dt {Q}t+dt
    PDE_Integrator_CRM<2U,SplitBoundary>  transient_pressure;
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    SAMG_Solver    samg_solver(&settings);
    transient_pressure.SetSolver(&samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    /// add extra functionality for alternative solver if needed
    transient_pressure.SetSolver(&linear_solver);
#endif

    NumIntegral_dNT_op_dN_dV<2U,InterFace<2U> >  conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                               conductance.MultiplyWithTimeIncrement(true);

    NumIntegral_NT_lhsop_N_dV<2U,InterFace<2U> > capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                               capacitance_lhs.LumpedFormulation(true);

    NumIntegral_NT_op_N_dV<2U,InterFace<2U> >  capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                               capacitance_rhs.LumpedFormulation(true);

    NumIntegral_NT_op_N_dV<2U,InterFace<2U> >  source( model.Database(), "fluid volume source",  "fluid pressure" );
                                               source.MultiplyWithTimeIncrement(true);
                                               source.AddAccumulateLater();
                                               source.LumpedFormulation(true);
   
    PointSource_rhsop<2U>                      nsource( model.Database(), "nodal fluid volume source",  "fluid pressure" );
                                               source.MultiplyWithTimeIncrement(true);
                                               source.AddAccumulateLater();
                                               source.LumpedFormulation(true);

    VelocityAndVolumeFlux<2U,Element<2U> >    velocity( model,  "conductivity", "porosity", "fluid pressure", false );

    transient_pressure.Add( &conductance );
    transient_pressure.Add( &capacitance_lhs );
    transient_pressure.Add( &capacitance_rhs );
    transient_pressure.Add( &source );
    transient_pressure.Add( &nsource );
    transient_pressure.AddPostProcess( &velocity );
    transient_pressure.TimeIncrement( time_increment );
   
     if ( first_call ) {
#ifdef CSMP_WITH_SAMG_SOLVER
         // targeting SAMG DLL 2 for this pressure solver (same as for steady state solver, but now reusing solution from previous timestep if available)
         settings.SetSolverInstance(2);
         // io
         settings.Set_iout1( 0 );
         settings.Set_iout2( 0 );
         settings.Set_idmp( -1 );
         // SAMG solution criteria
         settings.Set_eps(0.);
         settings.Set_rel_eps(1.E-10);
         
         // Re-use solver setup from previous timestep, internal checks force setup when required
         settings.Set_iswit(7);

         // Re-use solver setup from previous timestep, internal checks force setup when required
         settings.Set_iswit(7);
         // Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
      //    dynamic_cast<SteadyStateDiffusionSolver<3U,Region>*>(fluid_pressure_)->GetSolverSettings().GetSolverSettings().Set_nred(1);
         // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
         settings.Set_a_cmplx(2);
         // Pre-adjust SAMG mesh complexity, based on solver output
         settings.Set_g_cmplx(1.5);
         settings.Set_w_avrge(2);
#else
         /// add extra functionality for alternative solver if needed
#endif
       }

    double64& model_time( ModelTime::Instance().modelTime );

    if ( verbose_ ) {
         cout <<"\nComputeTransientFluidPressure: Input parameter ranges: "<< endl;
         printRangeOfVariable( model, "fracture conductivity", true );
         printRangeOfVariable( model, "nodal fluid volume source", true );
         printRangeOfVariable( model, "total system compressibility", true );
      }
    
    cout <<"\n\n\n\nComputeTransientFluidPressure: ";
    cout <<"Computing transient '"<< "fluid pressure" <<"' at t="<< fixed << setprecision(0) << model_time <<" secs.";
    cout <<"\n\t\t(FD Backward-Euler time increment, delta_t="<< time_increment <<" secs)."<< endl;
    transient_pressure.TimeIncrement( 1. / time_increment );

    transient_pressure.IntegrateOver( model.SplitBoundary(split_boundary_name) );

    // setting relative solution tolerance after first solve
    if ( first_call )  {
        #ifdef CSMP_WITH_SAMG_SOLVER
         // suppress verbose output
         settings.Set_iout1(-1);
         settings.Set_iout2(-1);
         // absolute tolerance set to 1e-14
         settings.Set_eps( 1.E-14 );
         // use solution from last step as an initial guess
         settings.Set_itypu(0);
         // reuse solver setup from last call
         settings.Set_iswit(3);
        #else
         /// add extra functionality for alternative solver if needed
        #endif
      }
  
    if ( verbose_ ) printRangeOfVariable( model, split_boundary_name, "fluid pressure", true );
   
    first_call = false;     
   
 } // ComputeTransientFluidPressure(SplitBoundary)








/**
    creates a region from finite elements that have a lower dimension than the model dimension.
    The user has the option to define the name of this region.
    
    function returns true if it succeeded in forming the new region.
*/
bool createLowerDimensionalRegion( Model<DIM>& model, const char* name_of_new_region )
 {
    Region<DIM>& mref = model.Region("Model");
    set<size_t> element_idx;
    for ( vector<Element<DIM>*>::const_iterator it=mref.ElementsBegin(); it!=mref.ElementsEnd(); ++it )
      if ( (*it)->FE()->IsLineElement() ) element_idx.insert( (*it)->Idx() );
    vector<size_t> unique_idx( element_idx.begin(), element_idx.end() );
   
    if ( model.ContainsRegion( name_of_new_region ) ) return false;
   
    model.FormRegionFrom( name_of_new_region, unique_idx );
   
    return true;
 }





/**
    Computes the gas conductivity, Kg, total systems compressibility, ct.
    
       ct = phi * cf + (1 - phi) * cr

       multiphase for later: ct(sw) = phi * (cw * sw + (1 - sw) * cg)  + (1 - phi) * cr

    hydraulic diffusivity, kappa = k / (phi * ct * mu);
*/
void computeGasFlowProperties( Model<DIM>& model, const char* target_region )
 {
    csmp::Index ct_key = model.Database().StorageKey("total system compressibility");
    csmp::Index cr_key = model.Database().StorageKey("rock compressibility");
    csmp::Index cf_key = model.Database().StorageKey("gas compressibility");
    csmp::Index phi_key = model.Database().StorageKey("porosity");
    csmp::Index k_key = model.Database().StorageKey("permeability");
    csmp::Index mu_key = model.Database().StorageKey("viscosity gas");
    csmp::Index K_key = model.Database().StorageKey("conductivity");
 
    Region<DIM>& mref = model.Region(target_region);
    ScalarVariable  cf, mu;

    for ( vector<Element<DIM>*>::iterator
          it=mref.ElementsBegin(); it!=mref.ElementsEnd(); ++it )
      {
         // retrieve input properties (knowing where they are placed)
         double64 phi = (*it)->Read( phi_key );
         double64 k   = (*it)->Read( k_key );
         double64 cr   = (*it)->Read( cr_key );
         (*it)->PropertyValueAtBaryCenter( cf_key, cf );
         (*it)->PropertyValueAtBaryCenter( mu_key, mu );
        
         // calculation
         double64 ct = phi * cf() + (1 - phi) * cr;
         double64 K  = k / mu();
        
         // store results
         (*it)->Store( ct_key, makeScalar(INIT_GUESS,ct) );
         (*it)->Store( K_key, makeScalar(INIT_GUESS,K) );
      }

 } // end computeGasFlowProperties

} // csmp
