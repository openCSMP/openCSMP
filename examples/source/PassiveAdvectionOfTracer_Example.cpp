#include "PassiveAdvectionOfTracer_Example.h"

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

#include "MeshDiagnostics.h"

using namespace std;

namespace csmp {

void PassiveAdvectionOfTracer_Example::Specifications()
  {
    SetTitle( "FVM passive advection of concentration in FEM velocity field" );
    SetDifficulty( 3 );
    SetCategory( "Numerical Methods" );
    AddAuthor( "SKM" );
    AddDescription( "3D passive tracer advection, choice of different advection schemes" );
    AddRequirement( "source files: 'PassiveAdvectionOfTracer_Example.cpp' and '*.h'" );
    AddRequirement( "prism_test: files .dat, .asc, -regions.txt, -configuration.txt.");
    AddRequirement( "variables(example25.txt)" );
  }


/** 
    Illustration of the  NodeCenteredFiniteVolumeTransport  scheme based family of transport schemes.
    Explicit vs implicit, first- vs. second-order accurate in space and time.
    With and without fluid volume sources.
    
    3D passive tracer advection combining finite elements (for pressure) with finite volumes (for advection of concentration profile)

    User can test degree of CFL overstepping that the model can cope with
    and what the consequences are for the shape of the advection front.

    Use models 'hex2_3', 'prism_test' or 'fracs2000' (.dat, .asc, -regions.txt, -configuration.txt)
    as input file suites.
*/
void PassiveAdvectionOfTracer_Example::Run()
{
 // ------------------------------------------------------------
 // 1. building model from ANSYS data files
 // ------------------------------------------------------------
  string  model_name;
  cout <<"\nmain: Enter name of 'ANSYS TETRA' input file (binary): ";
  cin >> model_name;

  ANSYS_Model3D  model3D( model_name.c_str(), "example25.txt");


 // ------------------------------------------------------------
 // 2. checking model and mesh quality
 // ------------------------------------------------------------
  printModelDimensions( model3D, true );
  // different ways of printing what a model actually consists off
  model3D.Out();
  model3D.RegionsOut();
  model3D.BoundariesOut();
  model3D.SplitBoundariesOut();

  MeshDiagnostics<3U>  mesh_check;
  double             vol_min, vol_max;
  mesh_check.ElementVolumeRange( model3D, vol_min, vol_max );
  cout <<"\nmain: element volume range: "<< vol_min <<" to "<< vol_max << endl;
  cout <<"\tlarge element-volume range (>10^3) can pose problems for the linear solver.2\n";


 // ------------------------------------------------------------
 // 3. configuring the model
 // ------------------------------------------------------------
  InputDataManager<3U>  model_configuration;
//  model_configuration.ConfigureFromFile( model3D, model_name.c_str(),
//                                         false, true, true, true, false );

  // the boolean variables determine which blocks in the input file shall be read
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

  // visualizing the input permeability and boundary conditions
  VTK_Interface<3U>  vtk_output;
  vtk_output.OutputDataToVTK( model3D, "permeability", "permeability", 0 );
  vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );


 // -----------------------------------------------------------------------
 // 4. hydraulic conductivity and other interrelations
 // -----------------------------------------------------------------------
  const double  fluid_viscosity(1.0e-03);
  ConstantFactor<3U,divides>  conductivity( model3D.Database(),
                                           "conductivity", "permeability",
                                            fluid_viscosity );
  model3D.Apply( conductivity );
  printRangeOfVariable( model3D, "conductivity" );

  vtk_output.OutputDataToVTK( model3D, "conductivity", "conductivity", 0 );


 // -----------------------------------------------------------------------
 // 5. computing a steady-state fluid pressure distribution in the model
 // -----------------------------------------------------------------------
  SteadyStateDiffusor<3U,Region> steady_state_pressure( model3D,
                                                        "conductivity", "fluid pressure",
                                                        "fluid volume source" );
  // postprocessing of pressure gradients and flow velocities
  VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model3D, "conductivity", "porosity", "fluid pressure" );
  steady_state_pressure.AddPostProcess( &postpro0 );

  // the calculation of fluid pressure
  steady_state_pressure.ComputeSteadyState( model3D.Region("Model") );

  // results: the pore velocity is the Darcy velocity divided by the porosity
  printRangeOfVariable( model3D, stdio, "fluid pressure" );
  printRangeOfVariable( model3D, stdio, "velocity" );
  printRangeOfVariable( model3D, stdio, "pore velocity" );

  vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 1 );
  vtk_output.OutputDataToVTK( model3D, "velocity",       "velocity",       1 );


 // -----------------------------------------------------------------------
 // 6. Demonstration of the generic transport algorithm
 //   (here you can compare different schemes with one another and
 //    overstep CFL to see how this adds numerical diffusion to the solution)
 // -----------------------------------------------------------------------
  printModelDimensions( model3D, true );

  cout <<"\nmain: Choose method of transport (1=explicit, 2=explicit, O(2), ";
  cout <<"3=implicit, 4=implicit O(2), 5=4+bijective mapping, 6,7=tests of volume integration. ";
  cout <<"8=TestNodeCenteredFiniteVolumeStencil ";
  cout <<" Try (3)=implicit single-step computation to get an idea of transport distance: ";
  int32_t     tmethod;
  double  max_error(1.0e-5);
  cin >>    tmethod;

  switch( tmethod ) {
       case 1:
          cout <<"\nmain: Would you like to restrict computation to group (yes=1, 0=no)? ";
          cin >> tmethod;
          if ( tmethod != 1 ) AdvectVariableExplicit( model3D, false /* second order=false */ );
          else {
               string group_name;
               cout <<"\nmain: Enter name of model region: ";
               cin >> group_name;
               AdvectVariableExplicit( model3D, group_name.c_str() );
               vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "new-concentration", "new concentration", 1, true );
            }
         break;
       case 2:  AdvectVariableExplicit( model3D, true );
         break;
       case 3:
          cout <<"\nmain: Would you like to restrict computation to model region (yes=1, 0=no)? ";
          cin >> tmethod;
          if ( tmethod != 1 ) AdvectVariableFirstOrderImplicit( model3D, vtk_output );
          else {
               string group_name;
               cout <<"\nmain: Enter name of region: ";
               cin >> group_name;
               AdvectVariableFirstOrderImplicit( model3D, group_name.c_str() );
               vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "concentration", "concentration", 1, true );
            }
         break;
       case 4:  AdvectVariableSecondOrderImplicit( model3D, false );
         break;
       case 5:  AdvectVariableSecondOrderImplicit( model3D, true );
         break;
       case 6:
           cout <<"\nmain: Calculated flux mismatch: ";
           cout << TestNodeCenteredFiniteVolumeTransport_PrescribedVelocity( model3D ) << endl;
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
          TestNodeCenteredFiniteVolumeTransport( model3D );
          printRangeOfVariable( model3D, stdio, "finite volume" );
          vtk_output.OutputDataToVTK( model3D, "finite-volume", "finite volume", 1 );

          max_error = printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
          vtk_output.OutputDataToVTK( model3D, "nodal-flux-mismatch", "nodal flux mismatch", 1 );

          if ( max_error > 1.0e-7 ) {
               model3D.FormRegionFrom( "corrupted-flux", "nodal flux mismatch", 1e-7, 100. );
               vtk_output.OutputDataToVTK( model3D, "corrupted-flux", "nodal-flux-mismatch", "nodal flux mismatch", 1, true );
           }
         break;

       case 8:  TestNodeCenteredFiniteVolumeStencil( model3D, vtk_output );
         break;

       default:
           cout <<"\nmain: Transport method not recognized."<< endl;
         return;
    }

  printRangeOfVariable( model3D, stdio, "concentration" );
  vtk_output.OutputDataToVTK( model3D, "concentration", "concentration", 99 );

  cout <<"\nmain: That's it."<< endl;

} // end Run




// *************************************************************************************************
//
// definitions of auxiliary functions
//
// *************************************************************************************************

void PassiveAdvectionOfTracer_Example::AdvectVariableExplicit( Model<3U>& sg, bool second_order )
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
   double time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
   clock_t ticks = clock();
   explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
   ticks = clock() - ticks;
   cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit



void PassiveAdvectionOfTracer_Example::AdvectVariableExplicit( Model<3U>& sg, const char* region, bool second_order )
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
   double time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock();
  explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit








void PassiveAdvectionOfTracer_Example::AdvectVariableFirstOrderImplicit( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
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
    double time_interval, Courant_multiplier;
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
void PassiveAdvectionOfTracer_Example::AdvectVariableFirstOrderImplicit( Model<3U>& sg, const char* group )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( group, sg,
                                                  "porosity", "concentration", "velocity",
                                                  "nodal fluid volume source", false, false );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Configuring TRANSPORT simulation: IMPIMS for region'"<< group <<"'"<< endl;
  cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;

  cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
  double time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //         fluxbalancecorrection=true, updateporevols=false
  advector.AdvectVariable( time_interval, Courant_multiplier, false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableFirstOrderImplicit





void PassiveAdvectionOfTracer_Example::AdvectVariableSecondOrderImplicit( Model<3U>& sg, bool bijective_mapping )
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
  double time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //                  flux_balance_correction  update_pore_volumes
  advector.AdvectVariable( time_interval, Courant_multiplier, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableSecondOrderImplicit




 void PassiveAdvectionOfTracer_Example::AdvectVariableSecondOrderImplicitSecondOrderInTime( Model<3U>& sg, bool bijective_mapping )
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
    double time_interval, Courant_multiplier;
    cin >> time_interval >> Courant_multiplier;

    cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
    clock_t ticks = clock();
    advector.AdvectVariable( time_interval, Courant_multiplier, bijective_mapping );
    ticks = clock() - ticks;
    cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

   } // end advectVariableSecondOrderImplicit





/**
  Calculates flux mismatch for predefined velocity field and normalizes it
  by finite volume. The finite volume and absolute value of the flux mismatch
  are reported to the variables "finite volume" and "nodal flux mismatch", respectively.
  do not use when surface elements are also present in model!
*/
double  PassiveAdvectionOfTracer_Example::TestNodeCenteredFiniteVolumeTransport_PrescribedVelocity( Model<3U>& sg )
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
  double             emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  // for all interior nodes we calculate the normalised flux balance
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.PerimeterNodesBegin(); it++ ) {
       sc() = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       (*it)->Store( prop_key, sc );
       emax = std::max( emax, fabs(sc()) );
    }

  // for all boundary nodes we set the balance to zero because we cannot evaluate it
  for ( vector<Node<3U>*>::iterator it=gref.PerimeterNodesBegin(); it!=gref.NodesEnd(); it++ )
    (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), 0.) );

  // finding the worst finite volume and analyzing it
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - fabs((*it)->Read( prop_key ))) <= numeric_limits<double>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( size_t i=0U; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

  return emax;

 } // end TestNodeCenteredFiniteVolumeTransport_PrescribedVelocity






void PassiveAdvectionOfTracer_Example::TestNodeCenteredFiniteVolumeTransport( Model<3U>& sg )
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
  double        emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ ) {
       sc() = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       if ( (*it)->AtBoundary() != NOT and (*it)->Status( pf_key ) == DIRICH )
         (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key),0.) );
       else
         (*it)->Store( prop_key, sc );
       emax = std::max( emax, sc() );
    }

  // finding the worst finite volume and analyzing it
  for ( vector<Node<3U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - (*it)->Read( prop_key )) <= numeric_limits<double>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( size_t i=0U; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

 } // end TestNodeCenteredFiniteVolumeTransport





void PassiveAdvectionOfTracer_Example::TestNodeCenteredFiniteVolumeStencil( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
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

  NCFVT_methods( sg, advector, vtkOut );
  cout << "End of Testing!" << endl;
  
} //end TestNodeCenteredFiniteVolumeStencil




/**
    Executing and testing all finite-volume related methods
*/
bool PassiveAdvectionOfTracer_Example::NCFVT_methods( Model<3U>& model3D,
                                                      NodeCenteredFiniteVolumeTransport<3U>& advector3D,
                                                      VTK_Interface<3U>& vtk_output )
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
          const double timeInterval(1.e3);
          const double courantMultiplier(1.e5);
          double courantIncrement;
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

          double CFLcheck(1.e2);
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

          double outflow(advector3D.ModelOutflow());
          //                          ^^^^^^^^^^^^
          double inflow(advector3D.ModelInflow());
          //                         ^^^^^^^^^^^^

          cout << "\nadvector.ModelOutflow(): " << outflow ;
          cout << "\nadvector.ModelInflow(): " << inflow ;


          /*
          * -----------  9.)BoundaryFluxes  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 9.)BoundaryFluxes -----------------------------\n";

          double boundaryFluxes(advector3D.BoundaryFluxes(inflow, outflow));
          //                                 ^^^^^^^^^^^^^^
          cout << "\nBoundary Fluxes with previous as Input";
          cout << " arguments returns: " << boundaryFluxes ;


          /*
          * -----------  10.)FluxBalance  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 10.)FluxBalance -------------------------------\n";

          double fmin, fmax;
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

          double poreVolume(0.);
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
          double customPressure(2.e7);
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

} // NCFVT_methods




} // csmp


