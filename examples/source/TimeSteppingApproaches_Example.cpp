#include "TimeSteppingApproaches_Example.h"

#include "Model.h"
#include "ModelTopology.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "ModelTime.h"

// fluid pressure algorithm and velocity computation
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "SAMG_Solver.h"

#include "ConstantFactor.h"

#include "TextInterface.h"
#include "MapleInterface.h"

#include "LineElementMesher.h"

#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

using namespace std;

namespace csmp{

void TimeSteppingApproaches_Example::Specifications()
{
  SetTitle( "Finite-difference time-stepping approaches" );
  SetDifficulty( 2 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "SKM" );
  AddDescription( "illustration of Backward-Euler-, Crank-Nicholson, and explicit time discretization of pressure equation in 1D" );
  AddDescription( "source in: TimeSteppingApproaches_Example.cpp" );
}



/** *****************************************************************************************

   (17) 1D fluid pressure diffusion comparing Crank-Nicholson time-stepping with other
        approaches.

        Use Maple or Excel to visualize the output.

 ***************************************************************************************** */
void TimeSteppingApproaches_Example::Run()
{
	// ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  enum{ DIM=1 };
  // 1. 1D mesh with 100 elements
   // --------------------------------------------
   VSet<DIM>  mesh_container;
   //         min_element_size, max_size, width_transition_zone, model_length
   LineElementMesher<DIM>  mesher;
   mesher.BuildUniformMesh( mesh_container, 10., 10 ); // 10 m sample

   // creating a model topological region
   set<string>    fem_types; fem_types.insert("LINEAR_BAR");
   ModelTopology  mesh_topology(true);  // isoparametric

   // making a rock and a fracture region of elements
   vector<size_t>  elms;
   for ( size_t i=0; i<mesh_container.Elements(); i++ ) elms.push_back(i);
   mesh_topology.AddDomain( "ROCK", fem_types, elms );
   elms.erase( elms.begin(), elms.end() );

   mesh_container.Out();


   // 2. Building the 1D Region named 'model'
   // --------------------------------------------
   Model<DIM>  model( mesh_topology, mesh_container, "example22.txt", false );


   // 3. Input of material properties and initial conditions
   // ------------------------------------------------------
   model.InputPropertyValue( "permeability",          		    makeScalar(PLAIN,1.0e-14) );
   model.InputPropertyValue( "porosity",              		    makeScalar(PLAIN,0.2) );
   model.InputPropertyValue( "storativity",           		    makeScalar(PLAIN,1.0e-9) );
   model.InputPropertyValue( "previous fluid volume source", makeScalar(PLAIN,0.) );
   model.InputPropertyValue( "fluid volume source",   		    makeScalar(PLAIN,0.) );


   // 4. Calculating hydraulic conductivity
   // --------------------------------------
   const double fluid_viscosity(1.0e-03);
   ConstantFactor<DIM,divides>  conductivity( model.Database(),
                                                      "conductivity", "permeability",
                                                       fluid_viscosity );
   model.Apply( conductivity );
   printRangeOfVariable( model, "conductivity" );


   // 5. Dirichlet boundary conditions for pressure and temperature
   // -------------------------------------------------------------
   cout <<"\nmain: Enter left and right fluid pressure: ";
   double pleft, pright;
   cin >> pleft >> pright;

   model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,pright) );
   model.InputBoundaryValue( CNR1, "fluid pressure", makeScalar(DIRICH, pleft) );
   model.InputBoundaryValue( CNR2, "fluid pressure", makeScalar(DIRICH, pright) );
   ScalarVariable  bpleft( PLAIN, pleft ), bpright( PLAIN, pright );
 //  model.AssignConstraint( "fluid pressure", mesh_container.Vertices(), bpright );
   model.Region("Model").N(mesh_container.Vertices()-1U)->Store( model.Database().StorageKey("fluid pressure"), bpright );

   // giving the first element a high permeability
   ScalarVariable  K_fracture( PLAIN, 1.0e-9 / 1.0e-3 );
 //  model.AssignConstraint( "conductivity", 1U, k_fracture );
   model.Region("Model").E(mesh_container.Elements()-1U)->Store( model.Database().StorageKey("conductivity"), K_fracture );

   mesh_container.Erase();


   // ------------------------------------------------------------------------------------
   //
   // Build a transient fluid pressure algorithm using Crank-Nicholson time-stepping:
   //
   //  Vs.1:   ([C] + dt/2 [K]){p}t+dt = ([C] - dt/2 [K]){p}t + dt/2 ({Q}t + {Q}t+dt)
   //
   // ------------------------------------------------------------------------------------
   PDE_Integrator<DIM,Region>  CrankNicholson1_pf;
   #ifdef CSMP_WITH_SAMG_SOLVER
   SAMG_Solver solver;
   CrankNicholson1_pf.SetSolver( solver );
   #else
   CSMP_DEFAULT_LINEAR_SOLVER solver;
   CrankNicholson1_pf.SetSolver( solver );
   #endif

   // conductance matrix dt/2 * [K] at pressure t+dt
   NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  CN1_conductance1( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                            CN1_conductance1.MultiplyWithTimeIncrement(true);
                                            CN1_conductance1.LumpedFormulation(true);

   // conductance matrix dt/2 * [K] at pressure t
   NumIntegral_DNT_rhsop_DN_dV<DIM,Element<DIM> >  CN1_conductance2( model.Database(), "conductivity", "fluid pressure" );
                                               CN1_conductance2.MultiplyWithTimeIncrement(true);
                                               CN1_conductance2.SubtractAccumulate();

   // left-hand side capacitance matrix [C]
   NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> > CN1_capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                            CN1_capacitance_lhs.LumpedFormulation(true);

   // right-hand side capacitance matrix [C]
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN1_capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                          CN1_capacitance_rhs.LumpedFormulation(true);

   // source vector at dt * {Q} at current pressure
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN1_source1( model.Database(), "previous fluid volume source",  "fluid pressure" );
                                          CN1_source1.MultiplyWithTimeIncrement(true);
                                          CN1_source1.AddAccumulateLater();
                                          CN1_source1.LumpedFormulation(true);

   // source vector at dt * {Q} at current pressure
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN1_source2( model.Database(), "fluid volume source",  "fluid pressure" );
                                          CN1_source2.MultiplyWithTimeIncrement(true);
                                          CN1_source2.AddAccumulateLater();
                                          CN1_source2.LumpedFormulation(true);

   VelocityAndVolumeFlux<DIM,Element<DIM> >   velocity( model, "conductivity", "porosity", "fluid pressure", false );

   // add PDE_Operators and post-processor to the FE Algorithm
   // and multiply with time-increment
   CrankNicholson1_pf.Add( &CN1_conductance1 );
   CrankNicholson1_pf.Add( &CN1_conductance2 );
   CrankNicholson1_pf.Add( &CN1_capacitance_lhs );
   CrankNicholson1_pf.Add( &CN1_capacitance_rhs );
   CrankNicholson1_pf.Add( &CN1_source1 );
   CrankNicholson1_pf.Add( &CN1_source2 );
   CrankNicholson1_pf.AddPostProcess( &velocity );






   // ------------------------------------------------------------------------------------
   // Build a transient fluid pressure algorithm using Crank-Nicholson time-stepping:
   //
   // Vs. 2    ([C]/dt + 1/2 [K]){p}t+dt = ([C]/dt - 1/2 [K]){p}t + 1/2 ({Q}t + {Q}t+dt)
   //
   // ------------------------------------------------------------------------------------
   PDE_Integrator<DIM,Region>  CrankNicholson2_pf;
#ifdef CSMP_WITH_SAMG_SOLVER
   SAMG_Solver  samg_solver;
   CrankNicholson2_pf.SetSolver(samg_solver);
#else
   CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
   CrankNicholson2_pf.SetSolver(linear_solver);
#endif

   // conductance matrix dt/2 * [K] at pressure t+dt
   NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  CN2_conductance1( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                            CN2_conductance1.MultiplyBy( 1./2. );

   // conductance matrix dt/2 * [K] at pressure t
   NumIntegral_DNT_rhsop_DN_dV<DIM,Element<DIM> >  CN2_conductance2( model.Database(), "conductivity", "fluid pressure" );
                                               CN2_conductance2.MultiplyBy( 1./2. );
                                               CN2_conductance2.SubtractAccumulate();

   // left-hand side capacitance matrix [C]
   NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> >  CN2_capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                             CN2_capacitance_lhs.MultiplyWithTimeIncrement(true);
                                             CN2_capacitance_lhs.LumpedFormulation(true);

   // right-hand side capacitance matrix [C]
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN2_capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                          CN2_capacitance_rhs.MultiplyWithTimeIncrement(true);
                                          CN2_capacitance_rhs.LumpedFormulation(true);

   // source vector at dt * {Q} at current pressure
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN2_source1( model.Database(), "previous fluid volume source",  "fluid pressure" );
                                          CN2_source1.MultiplyBy( 1./2. );
                                          CN2_source1.AddAccumulateLater();
                                          CN2_source1.LumpedFormulation(true);

   // source vector at dt * {Q} at current pressure
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  CN2_source2( model.Database(), "fluid volume source",  "fluid pressure" );
                                          CN2_source2.MultiplyBy( 1./2. );
                                          CN2_source2.AddAccumulateLater();
                                          CN2_source2.LumpedFormulation(true);

   // add PDE_Operators and post-processor to the FE Algorithm
   // and multiply with time-increment
   CrankNicholson2_pf.Add( &CN2_conductance1 );
   CrankNicholson2_pf.Add( &CN2_conductance2 );
   CrankNicholson2_pf.Add( &CN2_capacitance_lhs );
   CrankNicholson2_pf.Add( &CN2_capacitance_rhs );
   CrankNicholson2_pf.Add( &CN2_source1 );
   CrankNicholson2_pf.Add( &CN2_source2 );
   CrankNicholson2_pf.AddPostProcess( &velocity );



 // THIS IS THE ALGORITHM WHICH GIVES THE BEST ACCURACY BECAUSE MATRIX LUMPING IS RESTRICTED
 // TO THE CAPACITANCE MATRICES


   // ------------------------------------------------------------------------------------
   // 6. Backward Euler but with divided capacitance (WORKS BEST)
   //
   //     ([C]/dt + [K]){p}t+dt = ([C]/dt){p}t + ({Q}t+dt)
   //
   // NB: It is not a good idea to lump the conductance matrix
   //
   // ------------------------------------------------------------------------------------
   PDE_Integrator<DIM,Region>  BackwardEuler_pf;
   #ifdef CSMP_WITH_SAMG_SOLVER
   BackwardEuler_pf.SetSolver( solver );
   #else
   BackwardEuler_pf.SetSolver( solver );
   #endif

   // conductance matrix dt * [K] at pressure t+dt
   NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  BE_conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );

   // left-hand side capacitance matrix [C]
   NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> > BE_capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                            BE_capacitance_lhs.MultiplyWithTimeIncrement(true);
                                            BE_capacitance_lhs.LumpedFormulation(true);

   // right-hand side capacitance matrix [C]
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  BE_capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                          BE_capacitance_rhs.MultiplyWithTimeIncrement(true);
                                          BE_capacitance_rhs.LumpedFormulation(true);

   // source vector at dt * {Q} at current pressure
   NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  BE_source( model.Database(), "fluid volume source",  "fluid pressure" );
                                          BE_source.AddAccumulateLater();
                                          BE_source.LumpedFormulation(true);

   // add PDE_Operators and post-processor to the FE Algorithm
   // and multiply with time-increment
   BackwardEuler_pf.Add( &BE_conductance );
   BackwardEuler_pf.Add( &BE_capacitance_lhs );
   BackwardEuler_pf.Add( &BE_capacitance_rhs );
   BackwardEuler_pf.Add( &BE_source );
   BackwardEuler_pf.AddPostProcess( &velocity );




   // 7. Setting up the time-stepping loop for the transient calculation
   // ------------------------------------------------------------------
   double&        model_time( ModelTime::Instance().modelTime );
   double 	      maxtime, time_increment, ramp_factor;
   const double   year(31536000.0);  // 1 year in seconds
   uint32_t           timestep(1);

   cout << "\nEnter the run-time in years: ";
   cin  >> maxtime;
   maxtime *= year;
   cout << "\nEnter the time increment in seconds and ramp factor: ";
   cin  >> time_increment >> ramp_factor;


   // 7. Transient loop: Compute fluid pressure during each time-step and output the results for each time step
   // ---------------------------------------------------------------------------------------------------------
   cout <<"\nmain: Choose timestepping method Cranck-Nicholson vs. 1 (1), Cranck-Nicholson vs. 2 (2), Backward-Euler (3): ";
   int32_t  answer(0);
   cin >> answer;

   while ( model_time <= maxtime )
     {
       cout << "\n\nmain: COMPUTING TIMESTEP " << timestep << endl;

       switch( answer ) {
           case 1:
               CrankNicholson1_pf.TimeIncrement( time_increment/2. );
           model.Apply( CrankNicholson1_pf );
             break;

           case 2:
           CrankNicholson2_pf.TimeIncrement( 1. / time_increment );
           model.Apply( CrankNicholson2_pf );
             break;

           case 3:
                 BackwardEuler_pf.TimeIncrement( 1. / time_increment );
                 model.Apply( BackwardEuler_pf );
             break;
           default:
               cout <<"\nmain: error: choice of time stepping scheme not recognized (1-3), vs.: "<< answer;
               return;
          }

       // output variables to file and screen
       printRangeOfVariable( model, "fluid pressure" );
       printRangeOfVariable( model, "velocity" );

       writeVariableToMapleTextFile( model, "fluid pressure", timestep, model_time );

       // Preparing next Time Step
       model_time    += time_increment;
       time_increment *= ramp_factor;
       timestep++;

       cout << "\nmain: ELAPSED TIME " << model_time / year << " years " << endl;
     }

   cout <<"\nmain: That's it..."<< endl;
} // Run()

} // csmp
