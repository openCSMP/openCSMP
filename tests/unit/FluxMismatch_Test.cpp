#include "FluxMismatch_Test.h"
#include "AP_testUtilities.h"
#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"
#include "ConstantFactor.h"
#include "VTK_Interface.h"
#include "FiniteVolumeStencil.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"

using namespace std;

namespace csmp {

// testing whether there is a flux mismatch when integrating velocity fields 
// with the node-centered FV scheme
int flux_mismatch( bool bPrescribedVelocity ) 
  {
    const bool verbose(false);
    
    VSet<3>  mesh_container;
   
   //(a)test_Create_One_Hexahedra_VSet(mesh_container); //0
   
   //(b)test_Create_One_Prism_VSet(mesh_container,true); //0

   //(c)
   test_Create_Pyramid_Hexa_VSet(mesh_container,true);

   //(d)
   //test_Create_Hexahedra_VSet(mesh_container, true); //10-17

   //(e)test_Create_Prism_VSet(mesh_container, true); //10-16

   //(f) 
   //test_Create_Pyramid_VSet(mesh_container, true); //10-16
   
   //(g)
  // test_Create_Prism_Hexa_VSet(mesh_container, true);
   
   Model<3U>  sg( mesh_container, "example15.txt" ); 
   mesh_container.Erase();

   //Initialize properties (depends on what you need...):
   //model configuration
   sg.InputPropertyValue ( "fluid pressure", makeScalar(PLAIN,0.) );
   sg.InputPropertyValue ( "permeability", makeScalar(PLAIN,1.0e-12) );
   sg.InputPropertyValue ( "porosity", makeScalar(PLAIN,0.25) );
   sg.InputPropertyValue ( "fluid volume source", makeScalar(PLAIN,0.) );
   sg.InputPropertyValue ( "nodal fluid volume source", makeScalar(PLAIN,0.) );
   sg.InputPropertyValue ( "concentration", makeScalar(PLAIN,0.) );
   
   //Pressure boundary conditions
   const double pressure (10*101325.);
   sg.InputBoundaryValue( LEFT,  "fluid pressure",   makeScalar(DIRICH, 0.) ); 
   sg.InputBoundaryValue( RIGHT, "fluid pressure",   makeScalar(DIRICH,pressure) ); // 1 bar

    // -----------------------------------------------------------------------
   // 1. hydraulic conductivity and other interrelations
   // -----------------------------------------------------------------------
    const double fluid_viscosity(1.0e-03);
    ConstantFactor<3U,divides>  conductivity( sg.Database(), 
                                             "conductivity", "permeability", 
                                              fluid_viscosity );
    sg.Apply( conductivity );
    printRangeOfVariable( sg, "conductivity" );

   // -----------------------------------------------------------------------
   // 2. steady-state fluid pressure 
   // -----------------------------------------------------------------------
    SteadyStateDiffusor<3U,Region> steady_state_pressure( sg, "conductivity", "fluid pressure", "fluid volume source" );

    VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( sg, "conductivity", "porosity", "fluid pressure" );
    
    steady_state_pressure.AddPostProcess( &postpro0 );
    steady_state_pressure.ComputeSteadyState( sg.Region("Model") );

   // -----------------------------------------------------------------------
   // 3. Measure speed of the NodeCenteredFiniteVolumeTransport constructor
   // -----------------------------------------------------------------------    
    if ( bPrescribedVelocity )
      {
        VectorVariable<3U> velo(PLAIN,PLAIN,PLAIN, 3., 7., 1. );
        velo /= velo.Length();
        sg.InputPropertyValue( "velocity", velo );
      }
  
    if ( verbose ) cout << "\nflux_mismatch: Measuring the time required to build basic transport algorithm: ";
	  clock_t ticks = clock();
    NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg, "porosity", "concentration", "velocity", 
                                                         "nodal fluid volume source",false, false );
    ticks = clock() - ticks;
	  if ( verbose ) cout << ticks << endl << endl;
      
    if ( verbose ) cout <<"\nflux_mismatch: total volume of the model: "<< advector.FiniteVolume( "finite volume" ) << endl;

    if ( verbose ) cout <<"\n\nflux_mismatch: Measuring the divergence of fluxes."<< endl;
    const size_t compareSpeedTIMES(50U);
    ticks = clock() - ticks;
    for(size_t times = 0; times < compareSpeedTIMES; times++)
     advector.Divergence( "velocity", "nodal flux mismatch" );
	  if ( verbose ) cout <<"\ncompleted 50 divergence computation in: "<<  clock() - ticks << endl << endl;
        
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( sg, "nodal flux mismatch", "nodal flux mismatch", 0 );
  
    // identifying the Dirichlet boundaries (since they will have in or outflow)
    csmp::Index  pf_key = sg.Database().StorageKey("fluid pressure");
    
    ofstream ofs("vals_by_areas.txt");
    
    // zapping result values at model boundaries and normalizing divergence by finite volume
    csmp::Index  fv_key   = sg.Database().StorageKey("finite volume");
    csmp::Index  prop_key = sg.Database().StorageKey("nodal flux mismatch");
    ScalarVariable  sc;
    double  emax(0.);
    Region<3>&  sgref(sg.Region("Model"));
    for ( auto it=sgref.NodesBegin(); it!=sgref.PerimeterNodesBegin(); it++ ) {
         sc = fabs( (*it)->Read( prop_key ));// / cross_section_fv[it->ID()-1U] );
         (*it)->Store( prop_key, sc );           
         emax = std::max( emax, sc() );
      }
    for ( auto it=sgref.PerimeterNodesBegin(); it!=sgref.NodesEnd(); it++ )
      (*it)->Store( prop_key, sc=0. );
    
    // finding the worst finite volume and analyzing it
    for ( auto it=sgref.NodesBegin(); it!=sgref.PerimeterNodesBegin(); it++ )
      if ( fabs(emax - (*it)->Read( prop_key )) <= 1e-15 )
        if ( verbose )
        {
           cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: " << (*it)->Read( prop_key ) << endl;
           (*it)->Out();
           cout <<"\ncomposed of the element types: "<< endl;
           for ( auto i{0}; i<(*it)->Parents(); i++ )
             {
               cout << (*it)->Parent(i) << " ";             
               cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;       
             }
           cout << endl << endl;
       }
    
    if ( verbose ) {
        cout << "\nModel Inflow: " << setprecision(15) << advector.ModelInflow();
        cout << "\nModel Outflow: " << setprecision(15) << advector.ModelOutflow();
        cout << "\nDiff: " << setprecision(15) << advector.ModelOutflow()-advector.ModelInflow() << endl;
      }
    
    Standard_IO_Handler  stdio;
    if ( verbose ) {
        printRangeOfVariable( sg, stdio, "finite volume" );
        printRangeOfVariable( sg, stdio, "nodal flux mismatch" );
        printRangeOfVariable( sg, stdio, "velocity" );
      }
    
    vtk_output.OutputDataToVTK( sg, "nodal-flux-mismatch", "nodal flux mismatch", 0 );
  
//    rhinoOutput(sg);
//    printFiniteVolumes( sg );
                                    
	if ( verbose ) cout << "\n Finished comparing nodal flux mismatch, and speed.";
    
  return 0;
  
} // end 

 
 
} // namespace csmp
























