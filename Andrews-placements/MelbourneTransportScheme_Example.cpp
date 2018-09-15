#include "MelbourneTransportScheme_Example.h"
#include "Variables_TracerTransfer.h"

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
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"
#include "ExplicitTransport.h"
#include "ImplicitTransport.h"

#include "MeshDiagnostics.h"

using namespace std;

namespace csmp {

void MelbourneTransportScheme_Example::Specifications()
  {
    SetTitle( "Melbourne transport scheme: FVM passive advection of concentration in FEM velocity field" );
    SetDifficulty( 2 );
    SetCategory( "Numerical Methods" );
    AddAuthor( "SKM and Andrew Bromage" );
    AddDescription( "3D passive tracer advection, choice of different advection schemes" );
    AddRequirement( "source files: 'MelbourneTransportScheme_Example.cpp' and '*.h'" );
    AddRequirement( "prism_test: files .dat, .asc, -regions.txt, -configuration.txt.");
    AddRequirement( "variables(example15.txt)" );
  }


/** 
    Illustration of refactored transport scheme taking into account the insights from the 2012 Colleoli CSMP workshop.

    Use models 'hex2_3', 'prism_test' or 'fracs2000' (.dat, .asc, -regions.txt, -configuration.txt)
    as input file suites.
*/
void MelbourneTransportScheme_Example::Run()
{
 // ------------------------------------------------------------
 // 1. building model from ANSYS data files
 // ------------------------------------------------------------
  string  model_name;
  cout <<"\nmain: Enter name of 'ANSYS TETRA' input file (binary): ";
  cin >> model_name;

  ANSYS_Model3D  model3D( model_name.c_str(), "MelbourneTransportScheme_Example-variables.txt");
  printModelDimensions( model3D, true );


 // ------------------------------------------------------------
 // 2. checking the quality of the mesh
 // ------------------------------------------------------------
  MeshDiagnostics<3U>  mesh_check;
  double64             vol_min, vol_max;
  mesh_check.ElementVolumeRange( model3D, vol_min, vol_max );
  cout <<"\nmain: element volume range: "<< vol_min <<" to "<< vol_max << endl;
  cout <<"\tlarge element-volume range (>10^3) can pose problems for the linear solver.2\n";


 // ------------------------------------------------------------
 // 3. configuring the model
 // ------------------------------------------------------------
  InputDataManager<3U>  model_configuration;
  
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
  const double64  fluid_viscosity(1.0e-03);
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
  
  /** Comments
  
      - philosophy is to remove the need for the reader to know whehere a variable is placed. 
      - the scheme only needs to know whether this is an Element stencil or a Finite Volume loop
        and what the placement of the result variable should be
        
        Obtain() instead of Interpolate() since these operations may or may not be interpolations
        
        use compile-time asserts to make compile-time error message more telling
  */

  // EXPLICIT TRANSPORT SCHEME
  // -------------------------
  bool second_order_in_space(false);
  
  
  // Variable placement
  // ------------------
  INDEX<SCALAR,NODE>    pf_key(model3D.Database().StorageKey("fluid pressure"));
  INDEX<SCALAR,ELEMENT> K_key(model3D.Database().StorageKey("conductivity"));
  Region<3U>&           model_domain(model3D.Region("Model"));
  
  for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
  
       // the old way
       double64 K = (*it)->Read(K_key);
       ScalarVariable pf;
       (*it)->PropertyValueAtBaryCenter( pf_key, pf );
    
       // using the finite element placement
       auto e = (*it)->AtBarycenter();
       double64 pf_val = e.Obtain(pf_key);
    
       // comparisons
       cerr <<"\n\tconductivity:   "<< K <<" vs. "<< e.Read(K_key);
       cerr <<"\n\tfluid pressure: "<< pf <<" vs. "<< pf_val <<"\n";
    }
  
  // 6.1 Most basic case: explicit, first-order, no fluid sources and sinks, prescribed velocity field
  // -------------------------------------------------------------------------------------------------
  
  ExplicitTransport<3U> advector( model3D, "Model", false);
  advector.StepSizeReductionFactor(0.1);
  
   cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
   if ( second_order_in_space ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
   else                         cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
   //cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;
   cout <<"\nEnter advection time: ";
   double64 time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
   clock_t ticks = clock();
   for (unsigned i = 1; i < 10; ++i) advector.AdvectVariable(time_interval/10.);
   ticks = clock() - ticks;
   cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;
  
  
  
  
  // 6.2 explicit, first-order, no fluid sources and sinks, heterogeneity and computed velocity field
  // -------------------------------------------------------------------------------------------------

  // 6.3 explicit, second-order, no fluid sources and sinks, heterogeneity and computed velocity field
  // -------------------------------------------------------------------------------------------------

  // 6.4 explicit, second-order, fluid sources and sinks, heterogeneity and computed velocity field
  // -------------------------------------------------------------------------------------------------

  // 6.5 explicit, second-order, transient pressure field, heterogeneity and computed velocity field
  // -------------------------------------------------------------------------------------------------

/*
  cout <<"\nmain: Choose method of transport (1=explicit, 2=explicit, O(2), ";
  cout <<"3=implicit, 4=implicit O(2), 5=4+bijective mapping, 6,7=tests of volume integration. ";
  cout <<"8=TestNodeCenteredFiniteVolumeStencil ";
  cout <<" Try (3)=implicit single-step computation to get an idea of transport distance: ";
  int32     tmethod;
  double64  max_error(1.0e-5);
  cin >>    tmethod;

  switch( tmethod ) {
       case 1:
          cout <<"\nmain: Would you like to restrict computation to group (yes=1, 0=no)? ";
          cin >> tmethod;
//          if ( tmethod != 1 ) AdvectVariableExplicit( model3D, false ); // second order=false
/*
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
*/
  printRangeOfVariable( model3D, stdio, "concentration" );
  vtk_output.OutputDataToVTK( model3D, "concentration", "concentration", 99 );

  cout <<"\nmain: That's it."<< endl;

} // end Run



} // csmp


