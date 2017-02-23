//
//  GenericFiniteVolumeTransport_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/01/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

// generic transport scheme
#include "GenericFiniteVolumeTransport_Test.h"
#include "finiteVolumeFunctions.h"
#include "FacetFlux_TracerTransferExplicit.h"
#include "TimeStepEvaluator.h"
#include "ExplicitTransport.h"

// model
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "ComputationalSettings.h"
#include "InputDataManager.h"
#include "ANSYS_Model3D.h"
#include "ModelTopology.h"
#include "VTK_Interface.h"

// integration od PDEs and post-processing
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// interrelations
#include "ConstantFactor.h"

// legacy finite-volume transport scheme
#include "FiniteVolumeStencil.h"
#include "StencilProcessor.h"
#include "ExplicitStencilProcessor.h"
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"
#include "NodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {


// *************************************************************************************************
//
// definitions of auxiliary functions
//
// *************************************************************************************************

/**
    Tests performed:
    
    Speed comparison between projections made in parametric versus physical space while checking accuracy
    at same time.
*/
void GenericFiniteVolumeTransport_Test::run()
 {
     BenchmarkGlobalVersusParametricIntegration();
   
 } // end run



/**
    Testing finite volume projection calculations for the underpinning stencils.
    
    1. get model and compute pressure gradients
    
    2. verify (for the interior elements only) that the flux balance is indeed zero using established approach
    
    3. verify face by face computations
    
    4. verify flux balance cell-by-cell
    
    5. verify overall flux balances for hybrid element mesh
    
    6. compare computation times for flux balances
 
*/
void GenericFiniteVolumeTransport_Test::BenchmarkGlobalVersusParametricIntegration()
 {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
      string  model_name("prism_test");
      ANSYS_Model3D  model( model_name.c_str(), "example25.txt");
      printModelDimensions( model, true );

     // ------------------------------------------------------------
     // 2. configuring the model
     // ------------------------------------------------------------
      InputDataManager<3U>   model_configuration;
      ComputationalSettings  run_settings;
      model_configuration.ConfigureFromFile( model,
                                             model_name.c_str(),
                                             false, 
                                             true,   // 2) default prop.values
                                             true,   // 3) region prop.values
                                             true,   // 4) essential box-boundary conditions
                                             true,   // 5) essential flags
                                             true,   // 6) boundary conditions
                                             run_settings );
      Standard_IO_Handler  stdio;
      printRangeOfVariable( model, stdio, "permeability" );

      // optional visualization of the input permeability and boundary conditions
      VTK_Interface<3U>  vtk_output;
      vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );


     // -----------------------------------------------------------------------
     // 3. hydraulic conductivity computation
     // -----------------------------------------------------------------------
      const double64  fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model.Apply( conductivity );
      printRangeOfVariable( model, "conductivity" );
      vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0 );


     // -----------------------------------------------------------------------
     // 4. computing a steady-state fluid pressure distribution in the model
     // -----------------------------------------------------------------------
      SteadyStateDiffusor<3U,Region> steady_state_pressure( model,
                                                            "conductivity", "fluid pressure",
                                                            "fluid volume source" );
    
      // postprocessing of pressure gradients and flow velocities
      VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model, "conductivity", "porosity", "fluid pressure" );
      steady_state_pressure.AddPostProcess( &postpro0 );

      // the calculation of fluid pressure
      steady_state_pressure.ComputeSteadyState( model );

      // results: the pore velocity is the Darcy velocity divided by the porosity
      printRangeOfVariable( model, stdio, "fluid pressure" );
      printRangeOfVariable( model, stdio, "velocity" );
      printRangeOfVariable( model, stdio, "pore velocity" );

      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
      vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1 );

     // -----------------------------------------------------------------------
     // 4. stepping over the model comparing facet by facet flux calculations
     // -----------------------------------------------------------------------
     model.InstantiateFiniteVolumes();
     const Region<3>& model_domain(model.Region("Model"));
     const csmp::Index p_key(model.Database().StorageKey("fluid pressure")),
                       K_key(model.Database().StorageKey("conductivity")),
                       v_key(model.Database().StorageKey("velocity"));
   
     std::vector<double64> DNR, DNS, DNT;
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
       {
          const size_t nodes((*it)->Nodes());
          // 1. computing facet velocity in parametric space
          // -----------------------------------------------
          // 1.1 getting ipol-functions at barycentre and computing the pressure gradient in parametric space
          Point<3U> bctr = (*it)->FV_Stencil()->Barycenter();
          (*it)->FE()->dNr( bctr[0], bctr[1], bctr[2], DNR );
          (*it)->FE()->dNs( bctr[0], bctr[1], bctr[2], DNS );
          (*it)->FE()->dNt( bctr[0], bctr[1], bctr[2], DNT );
          // pressure gradients / velocities
          const double64 K((*it)->Read(K_key));
          double64  dpdr(0.), dpds(0.), dpdt(0.);
          Point<3U> vD(0.);
          for ( size_t i=0U; i<nodes; ++i ) {
               const double64 p_node((*it)->N(i)->Read(p_key));
               dpdr  += DNR[i] * p_node;
               dpds  += DNS[i] * p_node;
               dpdt  += DNT[i] * p_node;
               vD[0] += -K * dpdr;
               vD[1] += -K * dpds;
               vD[2] += -K * dpdt;
            }
         
          // 2. facet projections
          // --------------------
          const size_t facets=(*it)->Facets();
          for ( size_t i=0U; i<facets; ++i ) {
               // 2.1 classic way of calculating facet fluxes in physical space
               // ------------------------------------------------------------
               double64 flux = (*it)->ProjectionOnFacetNormal( i, v_key );
               // Darcy velocity computation
               double64 flux_physical = (*it)->FacetArea(i) * flux;
            
               // 2.2 parametric space computation
               // --------------------------------
               Point<3U> fnu = (*it)->ParametricFacetNormal( i );
               double64  local_facet_area = (*it)->ParametricFacetArea( i );
               // computing the facet projection
               double64 flux_from_parametric = dotProduct( fnu, vD );
               // performing local integration via multiplication with integration weight
               flux_from_parametric *= local_facet_area;
               // transforming the result to physical space
               assert( (*it)->IsVolumeElement() );
               // get Jacobian and its determinant at the facet integration point
               Point<3U> fip( (*it)->FV_Stencil()->FacetIntegrationPoint(i,0U) );
               (*it)->FE()->dNr( fip[0], fip[1], fip[2], DNR );
               (*it)->FE()->dNs( fip[0], fip[1], fip[2], DNS );
               (*it)->FE()->dNt( fip[0], fip[1], fip[2], DNT );
               (*it)->CoordinateMatrix();
               (*it)->FE()->Jacobian( DNR, DNS, DNT );
               double64 scale_factor = (*it)->FE()->JacobianDeterminant();
               flux_from_parametric *= scale_factor;
            
               // 3. testing that the fluxes are the same
               // ---------------------------------------
               _equal( flux_physical, flux_from_parametric, numeric_limits<double64>::epsilon() );
            }
       }

 } // end





/**
    Executing and testing all finite-volume related methods
*/
bool test_NCFVT_methods( Model<3U>& model3D, NodeCenteredFiniteVolumeTransport<3U>& advector3D, VTK_Interface<3U>& vtk_output )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  //----------------------------------------------------------------------
  //------------------------- TESTING METHODS ----------------------------
  //----------------------------------------------------------------------


          /*
          * -- TESTING OF NCFVT METHODS-----------------------
          *
          * -- 1.)AdvectVariable -----------------------------
          */
          const double64 timeInterval(1.e3);
          const double64 courantMultiplier(1.e5);
          double64 courantIncrement;
          cout << "\n\n\n\n\n";
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 1.)AdvectVariable -----------------------------\n";
          cout << "\nInputAguments - timeInterval: " << timeInterval <<
                  ", courantMultiplier: " << courantMultiplier;

          courantIncrement = advector3D.AdvectVariable(timeInterval,courantMultiplier);
          //                            ^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;



          /*
          * -----------   2.)TransportPhase    -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 2.)TransportPhase 1D ONLY ---------------------\n";
          //cout << "\nInputAguments - timeInterval: " << timeInterval;

          //advector1D.TransportPhase(relperms, timeInterval);
          //         ^^^^^^^^^^^^^^




          /*
          * -----------  3.)CourantIncrement() -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 3.)CourantIncrement() 1D ONLY -----------------\n";

          //courantIncrement = advector1D.CourantIncrement();
          //                            ^^^^^^^^^^^^^^^^

          //cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * -----------   4.)CourantIncrement    -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 4.)CourantIncrement 1D ONLY -------------------\n";

          //courantIncrement = advector1D.CourantIncrement(relperms);
          //                            ^^^^^^^^^^^^^^^^

          //cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * ----------- 5.)AnisotropicCourantIncrement() -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 5.)AnisotropicCourantIncrement() --------------\n";

          courantIncrement = advector3D.AnisotropicCourantIncrement();
          //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * ----------- 6.)AnisotropicCourantIncrement -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 6.)AnisotropicCourantIncrement ----------------\n";

          //courantIncrement = advector3D.AnisotropicCourantIncrement(relperms);
          //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * -----------    7.)CFL_Multiplier   -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 7.)CFL_Multiplier -----------------------------\n";

          double64 CFLcheck(1.e2);
          advector3D.CFL_Multiplier(CFLcheck);
          //         ^^^^^^^^^^^^^^
          assert(CFLcheck == advector3D.CFL_Multiplier());
          cout << "\nManual CFL input value: " << CFLcheck ;
          cout << "\nReturn function for CFL multiplier gives: " << advector3D.CFL_Multiplier();
          //                                                                   ^^^^^^^^^^^^^^


          /*
          * -----------    8.)Inflow/Outflow   -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 8.)Model Inflow/Outflow -----------------------\n";

          double64 outflow(advector3D.ModelOutflow());
          //                          ^^^^^^^^^^^^
          double64 inflow(advector3D.ModelInflow());
          //                         ^^^^^^^^^^^^

          cout << "\nadvector.ModelOutflow(): " << outflow ;
          cout << "\nadvector.ModelInflow(): " << inflow ;


          /*
          * -----------  9.)BoundaryFluxes  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 9.)BoundaryFluxes -----------------------------\n";

          double64 boundaryFluxes(advector3D.BoundaryFluxes(inflow, outflow));
          //                                 ^^^^^^^^^^^^^^
          cout << "\nBoundary Fluxes with previous as Input";
          cout << " arguments returns: " << boundaryFluxes ;


          /*
          * -----------  10.)FluxBalance  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 10.)FluxBalance -------------------------------\n";

          double64 fmin, fmax;
          advector3D.FluxBalance(fmin, fmax);
          //         ^^^^^^^^^^^
          cout << "\nFluxBalance returns " << fmin << " as minimum and ";
          cout << fmax << " as maximum FV flux balance.";


          /*
          * ----- 11.)MultiplyScalarNodePropertyByFiniteVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 11.)MultiplyScalarNodePropertyByFiniteVolume --\n";

          advector3D.MultiplyScalarNodePropertyByFiniteVolume("fluid pressure");
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


          /*
          * ----- 12.)VolumeIntegrateScalarFiniteElementVariable  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS------------------------";
          cout << "\n\n/** -- 12.)VolumeIntegrateScalarFiniteVolumeVariable --\n";

          double64 poreVolume(0.);
          poreVolume =
          advector3D.VolumeIntegrateScalarFiniteElementVariable("porosity", true);
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nIntegration of porosity  returns " << poreVolume << " for FE Integration.";



          /*
          * ----- 13.)FiniteVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 13.)FiniteVolumes -------------------------------\n";

          poreVolume =
          advector3D.FiniteVolume("concentration");
          //         ^^^^^^^^^^^^
          vtk_output.OutputDataToVTK( model3D, "finite-volume", "concentration", 0 );

          cout << "\nFinite Volumes of porosity  returns " << poreVolume;

          /*
          * ----- 14.)VolumeIntegrate  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 14.)VolumeIntegrate -----------------------------\n";

          advector3D.VolumeIntegrate("porosity", "poreVolume");
          //         ^^^^^^^^^^^^^^^
          advector3D.VolumeIntegrate("porosity", "saturation oil", "oilVolume");
          //         ^^^^^^^^^^^^^^^

          /*
          * ----- 15.)AssignScalarBoundaryValues  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 15.)AssignScalarBoundaryValues ------------------\n";
          double64 customPressure(2.e7);
          advector3D.AssignScalarBoundaryValues(LEFT, "fluid pressure", DIRICH, customPressure, true);
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^


          /*
          * ----- 16.)PoreVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 16.)PoreVolume -----------------------------\n";

          poreVolume =
          advector3D.PoreVolume(50);
          //         ^^^^^^^^^^
          cout << "\nPoreVolume returns " << poreVolume;


          /*
          * ----- 17.)Out  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 17.)Out -----------------------------------------\n";

          advector3D.Out();
          //         ^^^

  return true;

} // test_NCFVT_methods





void testNodeCenteredFiniteVolumeStencils( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  //Testing NodeCenteredFiniteVolumeTransport.h
  //method testFiniteVolumeStencil

  cout << "Testing NodeCenteredFiniteVolumeTransport.h method testFiniteVolumeStencil " << endl;
  testFiniteVolumeStencil( sg.Database(), sg.Region( "Model" ), advector );

  test_NCFVT_methods( sg, advector, vtkOut );
  cout << "End of Testing!" << endl;
  
} //end TestNodeCenteredFiniteVolumeStencil






/**
  Calculates flux mismatch for predefined velocity field and normalizes it
  by finite volume. The finite volume and absolute value of the flux mismatch
  are reported to the variables "finite volume" and "nodal flux mismatch", respectively.
  do not use when surface elements are also present in model!
*/
double64  testNodeCenteredFiniteVolumeTransport_PrescribedVelocity( Model<3U>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  VectorVariable<3U>  velo(PLAIN,PLAIN,PLAIN, 3., 7., 1. );
  velo /= velo.Length(); // unit length
  velo.Out();
  sg.InputPropertyValue( "velocity", velo );

  cout <<"\n\tMeasuring the time required to build basic transport algorithm."<< endl;
  clock_t ticks = clock();
  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

  PropertyHandle<3U>  fv( sg,"finite volume",SCALAR,NODE);
  advector.FiniteVolume( "finite volume" );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Measuring the divergence of fluxes."<< endl;
  advector.Divergence( "velocity", "nodal flux mismatch" );

  // identifying the Dirichlet boundaries (since they will have in or outflow)
  csmp::Index  pf_key = sg.Database().StorageKey("fluid pressure");

  // zapping result values at model boundaries and normalizing divergence by finite volume
  csmp::Index          fv_key   = sg.Database().StorageKey("finite volume");
  csmp::Index          prop_key = sg.Database().StorageKey("nodal flux mismatch");
  ScalarVariable       sc;
  double64             emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  // for all interior nodes we calculate the normalised flux balance
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.PerimeterNodesBegin(); it++ ) {
       sc = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       (*it)->Store( prop_key, sc );
       emax = std::max( emax, fabs(sc()) );
    }

  // for all boundary nodes we set the balance to zero because we cannot evaluate it
  for ( vector<Node<3U>*>::iterator it=gref.PerimeterNodesBegin(); it!=gref.NodesEnd(); it++ )
    (*it)->Store( prop_key, sc=0. );

  // finding the worst finite volume and analyzing it
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - fabs((*it)->Read( prop_key ))) <= numeric_limits<double64>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( size_t i=0U; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

  return emax;

 } // end testNodeCenteredFiniteVolumeTransport_PrescribedVelocity










void advectVariableExplicit( Model<3U>& sg, bool second_order )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  ExplicitNodeCenteredFiniteVolumeTransport<3U,ExplicitStencilProcessor>  explicit_advector(
                                                                            "Model", sg,
                                                                            "porosity",
                                                                            "concentration",
                                                                            "velocity",
                                                                            "nodal fluid volume source",
                                                                             second_order );

   cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
   if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
   else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
   cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement() << endl;
   cout <<"\nEnter advection time: ";
   double64 time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
   clock_t ticks = clock();
   explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
   ticks = clock() - ticks;
   cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit



void advectVariableExplicit( Model<3U>& sg, const char* region, bool second_order )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  ExplicitNodeCenteredFiniteVolumeTransport<3U,ExplicitStencilProcessor>  explicit_advector(
                                                                            region, sg,
                                                                            "porosity",
                                                                            "concentration",
                                                                            "velocity",
                                                                            "nodal fluid volume source",
                                                                             second_order );

   cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
   if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
   else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
   cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement() << endl;
   cout <<"\nEnter advection time: ";
   double64 time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock();
  explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit








void advectVariableFirstOrderImplicit( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
 {
    // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
    //ostream &cout = *GetStream();

    NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source",false, false );

    cout <<"\n\nadvectVariableFirstOrderImplicit: Configuring TRANSPORT simulation: IMPIMS"<< endl;
    cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement();
    cout.flush();
    cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
    double64 time_interval, Courant_multiplier;
    cin >> time_interval >> Courant_multiplier;
    cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
    clock_t ticks = clock();
    for ( int i=0; i<5; ++i ) {
        advector.AdvectVariable( time_interval/5, Courant_multiplier );
        vtkOut.OutputDataToVTK( sg, "concentration", "concentration", i+1 );
     }
    ticks = clock() - ticks;
    cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableFirstOrderImplicit





// restricted to a group
void advectVariableFirstOrderImplicit( Model<3U>& sg, const char* group )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( group, sg,
                                                  "porosity", "concentration", "velocity",
                                                  "nodal fluid volume source", false, false );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Configuring TRANSPORT simulation: IMPIMS for region'"<< group <<"'"<< endl;
  cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;

  cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
  double64 time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //         fluxbalancecorrection=true, updateporevols=false
  advector.AdvectVariable( time_interval, Courant_multiplier, false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableFirstOrderImplicit





void advectVariableSecondOrderImplicit( Model<3U>& sg, bool bijective_mapping )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source", true, false );

  cout <<"\n\nadvectVariableSecondOrderImplicit: Configuring TRANSPORT simulation: ";
  if ( bijective_mapping ) cout <<" IMPIMS with BIJECTIVE MAPPING."<< endl;
  else                     cout <<" IMPIMS without BIJECTIVE MAPPING."<< endl;
  cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;
  cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
  double64 time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //                  flux_balance_correction  update_pore_volumes
  advector.AdvectVariable( time_interval, Courant_multiplier, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableSecondOrderImplicit




 void advectVariableSecondOrderImplicitSecondOrderInTime( Model<3U>& sg, bool bijective_mapping )
   {
     // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
     //ostream &cout = *GetStream();

    NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source", true, true );

    cout <<"\n\nadvectVariableSecondOrderImplicitSecondOrderInTime: Configuring TRANSPORT simulation: ";
    if ( bijective_mapping ) cout <<" IMPIMS with BIJECTIVE MAPPING."<< endl;
    else                     cout <<" IMPIMS without BIJECTIVE MAPPING."<< endl;
    cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;
    cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
    double64 time_interval, Courant_multiplier;
    cin >> time_interval >> Courant_multiplier;

    cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
    clock_t ticks = clock();
    advector.AdvectVariable( time_interval, Courant_multiplier, bijective_mapping );
    ticks = clock() - ticks;
    cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

   } // end advectVariableSecondOrderImplicit









void testNodeCenteredFiniteVolumeTransport( Model<3U>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  cout <<"\n\tMeasuring the time required to build basic transport algorithm."<< endl;
  clock_t ticks = clock();
  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

  advector.FiniteVolume( "finite volume" );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Measuring the divergence of fluxes."<< endl;
  advector.Divergence( "velocity", "nodal flux mismatch" );

  // identifying the Dirichlet boundaries (since they will have in or outflow)
  csmp::Index  pf_key = sg.Database().StorageKey("fluid pressure");

  // zapping result values at model boundaries and normalizing divergence by finite volume
  csmp::Index     fv_key   = sg.Database().StorageKey("finite volume");
  csmp::Index     prop_key = sg.Database().StorageKey("nodal flux mismatch");
  ScalarVariable  sc;
  double64        emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ ) {
       sc = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       if ( (*it)->AtBoundary() != NOT and (*it)->Status( pf_key ) == DIRICH )
         (*it)->Store( prop_key, sc=0. );
       else
         (*it)->Store( prop_key, sc );
       emax = std::max( emax, sc() );
    }

  // finding the worst finite volume and analyzing it
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - (*it)->Read( prop_key )) <= numeric_limits<double64>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( size_t i=0U; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

 } // end TestNodeCenteredFiniteVolumeTransport












// *************************************************************************************************
//
// definitions of main test function:  run()
//
// *************************************************************************************************


void testSchemeAsComponent()
  {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
      string  model_name("prism_test");
      //cout <<"\nmain: Enter name of 'ANSYS TETRA' input file (binary): ";
      //cin >> model_name;

      ANSYS_Model3D  model3D( model_name.c_str(), "example25.txt");
      printModelDimensions( model3D, true );


     // ------------------------------------------------------------
     // 2. configuring the model
     // ------------------------------------------------------------
      InputDataManager<3U>  model_configuration;
      ComputationalSettings  run_settings;
      model_configuration.ConfigureFromFile( model3D,
                                             model_name.c_str(),
                                             false, 
                                             true,   // 2) default prop.values
                                             true,   // 3) group prop.values
                                             true,   // 4) essential box-boundary conditions
                                             true,   // 5) essential flags
                                             true,   // 6) boundary conditions
                                             run_settings );
      Standard_IO_Handler  stdio;
      printRangeOfVariable( model3D, stdio, "permeability" );

      // optional visualization of the input permeability and boundary conditions
      VTK_Interface<3U>  vtk_output;
      vtk_output.OutputDataToVTK( model3D, "permeability", "permeability", 0 );
      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );


     // -----------------------------------------------------------------------
     // 3. hydraulic conductivity computation
     // -----------------------------------------------------------------------
      const double64  fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model3D.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model3D.Apply( conductivity );
      printRangeOfVariable( model3D, "conductivity" );

      vtk_output.OutputDataToVTK( model3D, "conductivity", "conductivity", 0 );


     // -----------------------------------------------------------------------
     // 4. computing a steady-state fluid pressure distribution in the model
     // -----------------------------------------------------------------------
      SteadyStateDiffusor<3U,Region> steady_state_pressure( model3D,
                                                            "conductivity", "fluid pressure",
                                                            "fluid volume source" );
    
      // postprocessing of pressure gradients and flow velocities
      VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model3D, "conductivity", "porosity", "fluid pressure" );
      steady_state_pressure.AddPostProcess( &postpro0 );

      // the calculation of fluid pressure
      steady_state_pressure.ComputeSteadyState( model3D );

      // results: the pore velocity is the Darcy velocity divided by the porosity
      printRangeOfVariable( model3D, stdio, "fluid pressure" );
      printRangeOfVariable( model3D, stdio, "velocity" );
      printRangeOfVariable( model3D, stdio, "pore velocity" );

      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 1 );
      vtk_output.OutputDataToVTK( model3D, "velocity",       "velocity",       1 );


     // ---------------------------------------------------------------------------------
     // 5. Testing the generic transport scheme for tracer transport in single phase flow
     // ---------------------------------------------------------------------------------
     //   (here you can compare different schemes with one another and
     //    overstep CFL to see how this adds numerical diffusion to the solution)
     /*  
          Options
          -------
          A. prescribed divergence free velocity field
            1=explicit,
            2=explicit, O(2),
            3=implicit, 
            4=implicit O(2), 
            5=4+bijective mapping, 
            6,7=tests of volume integration,
            8=TestNodeCenteredFiniteVolumeStencil
      
          B. computed velocity field
     */
     // -----------------------------------------------------------------------

      int32     tmethod(1);
      double64  max_error(1.0e-5);
    
      // 5.1 setting up the transport scheme
      /*
           We store the volume of the finite volumes and their pore volumes
           
           - finite volume
           - finite volume (effective) pore volume
           - sector volume (stored at sector integration point)
           - sector weight: sector pore volume / finite volume pore volume = weighting factor
      */
      Region<3U>&  flow_domain(model3D.Region("Model"));
      initializeFiniteVolumeProperties( model3D, flow_domain );
    
// TESTING SECTOR INTEGRATION POINT STORAGE
      const csmp::Index swt_key(model3D.Database().StorageKey("node number"));
      // sector storage: writing global node numbers to sector IP's and reading them out
      for ( vector<Element<3U>*>::iterator it=flow_domain.ElementsBegin(); it!=flow_domain.ElementsEnd(); ++it )
        for ( size_t i=0U; i<(*it)->Sectors(); ++i )
          (*it)->Store( i, 0U, swt_key, makeScalar(PLAIN,(*it)->N(i)->Idx()) );
        
      // reading out node numbers and their double equivalents stored at the sector integration points
      for ( vector<Element<3U>*>::iterator it=flow_domain.ElementsBegin(); it!=flow_domain.ElementsEnd(); ++it ) {
           cerr <<"\nelement: "<< (*it)->Idx() << endl;
           for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
                cerr << (*it)->N(i)->Idx() <<":";
                cerr << (*it)->Read( i, 0U, swt_key ) <<" ";
             }
        }

//  TESTING FV VOLUME CALCULATIONS
      const csmp::Index fv_key(model3D.Database().StorageKey("finite element volume"));
      const csmp::Index sv_key(model3D.Database().StorageKey("fv sector volume"));
      const csmp::Index fvphi_key(model3D.Database().StorageKey("fv pore volume"));
      // accumulating matching sector volumes with finite element volumes
      double volume(0.);
      for ( vector<Node<3U>*>::iterator it=flow_domain.NodesBegin(); it!=flow_domain.NodesEnd(); ++it )
        volume += (*it)->Read( fv_key );
      cerr <<"\nFV total volume: "<< volume;

      volume = 0.;
      for ( vector<Element<3U>*>::iterator it=flow_domain.ElementsBegin(); it!=flow_domain.ElementsEnd(); ++it )
        for ( size_t i=0U; i<(*it)->Sectors(); ++i )
          volume += (*it)->Read( i, 0U, sv_key );
      cerr <<"\nFV total volume: "<< volume;

      double pvolume(0.);
      for ( vector<Node<3U>*>::iterator it=flow_domain.NodesBegin(); it!=flow_domain.NodesEnd(); ++it )
        pvolume += (*it)->Read( fvphi_key );
      cerr <<"\nFV total volume: "<< pvolume;



      // first order version only
      ExplicitTransport<3U>  explicit_advector( model3D, "Model" );
      printRangeOfVariable( model3D, "sector volume" );
      printRangeOfVariable( model3D, "sector pore volume" );
      printRangeOfVariable( model3D, "finite volume" );
      printRangeOfVariable( model3D, "FV pore volume" );

      cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
    //  if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
      cout <<" FIRST ORDER SCHEME."<< endl;
      cout <<"\nThe grid Courant number is "<< explicit_advector.TimeIncrement() << endl;
      cout <<"\nEnter advection time: ";
      double64 time_interval;
      cin >> time_interval;

      cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
      clock_t ticks = clock();
      explicit_advector.AdvectVariable( time_interval );
      ticks = clock() - ticks;
      cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;








      switch( tmethod ) {
           case 1:
              cout <<"\nmain: Would you like to restrict computation to group (yes=1, 0=no)? ";
              cin >> tmethod;
              if ( tmethod != 1 ) advectVariableExplicit( model3D, false /* second order=false */ );
              else {
                   string group_name;
                   cout <<"\nmain: Enter name of model region: ";
                   cin >> group_name;
                   advectVariableExplicit( model3D, group_name.c_str() );
                   vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "new-concentration", "new concentration", 1, true );
                }
             break;
           case 2:  advectVariableExplicit( model3D, true );
             break;
           case 3:
              cout <<"\nmain: Would you like to restrict computation to model region (yes=1, 0=no)? ";
              cin >> tmethod;
              if ( tmethod != 1 ) advectVariableFirstOrderImplicit( model3D, vtk_output );
              else {
                   string group_name;
                   cout <<"\nmain: Enter name of region: ";
                   cin >> group_name;
                   advectVariableFirstOrderImplicit( model3D, group_name.c_str() );
                   vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "concentration", "concentration", 1, true );
                }
             break;
           case 4:  advectVariableSecondOrderImplicit( model3D, false );
             break;
           case 5:  advectVariableSecondOrderImplicit( model3D, true );
             break;
           case 6:
               cout <<"\nmain: Calculated flux mismatch: ";
               cout << testNodeCenteredFiniteVolumeTransport_PrescribedVelocity( model3D ) << endl;
               printRangeOfVariable( model3D, stdio, "finite volume" );
               vtk_output.OutputDataToVTK( model3D, "finite-volume", "finite volume", 1 );

               max_error = printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
               vtk_output.OutputDataToVTK( model3D, "nodal-flux-mismatch", "nodal flux mismatch", 1 );

               if ( max_error > 1.0e-7 ) {
                    model3D.FormRegionFrom( "corrupted-flux", "nodal flux mismatch", 1e-7, 100. );
                    vtk_output.OutputDataToVTK( model3D, "corrupted-flux", "nodal-flux-mismatch", "nodal flux mismatch", 1, true );
                }
             break;
           case 7:
              testNodeCenteredFiniteVolumeTransport( model3D );
              printRangeOfVariable( model3D, stdio, "finite volume" );
              vtk_output.OutputDataToVTK( model3D, "finite-volume", "finite volume", 1 );

              max_error = printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
              vtk_output.OutputDataToVTK( model3D, "nodal-flux-mismatch", "nodal flux mismatch", 1 );

              if ( max_error > 1.0e-7 ) {
                   model3D.FormRegionFrom( "corrupted-flux", "nodal flux mismatch", 1e-7, 100. );
                   vtk_output.OutputDataToVTK( model3D, "corrupted-flux", "nodal-flux-mismatch", "nodal flux mismatch", 1, true );
               }
             break;

           case 8:
             break;

           default:
               cout <<"\nmain: Transport method not recognized."<< endl;
             return;
        }

      printRangeOfVariable( model3D, stdio, "concentration" );
      vtk_output.OutputDataToVTK( model3D, "concentration", "concentration", 99 );

      cout <<"\nmain: That's it."<< endl;

  } // end run



} // end csmp
