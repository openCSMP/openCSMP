#include "RhinoMesh_Example.h"

#include "Exception.h"
#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_definitions.h"
#include "Standard_IO_Handler.h"

// PDE Operators
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

// Interrelations
#include "ConstantFactor.h"

// finite-volume based transport for triangular FE
#include "StencilProcessor.h"
#include "NodeCenteredFiniteVolumeTransport.h"

// output
#include "PropertyData.h"
#include "VTK_Interface.h"
#include "RhinoSurfaceReader.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

void RhinoMesh_Example::Specifications()
{
  SetTitle( "RhinoSurfaceReader: Input of triangulated surfaces from Rhinoceros (McNeel&Assocs.) as mesh." );
  SetDifficulty( 2 );
  SetCategory( "Software Functionality" );
  AddAuthor( "SKM" );
  AddDescription( "source in: RhinoMesh_Example.cpp" );
  AddDescription( "3D fracture-only flow & transport simulation using Rhino meshes as input" );
  AddDescription( "transport is computed with the higher-order theta-limited implicit scheme" );
  AddRequirement( "file set: Rhino output '.raw' file called 'example20.raw'; variables file: 'example20.txt'");
} 




void RhinoMesh_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();
      char  file_name[200];
      cout <<"\nmain: Enter name of input file (ASCII format): ";
      cin >> file_name;

     // -----------------------------------------------------------------------
     // 0. reading in '.raw' meshes that were written to file as labeled
     //    entities
     // -----------------------------------------------------------------------
      SKM_RhinoSurfaceReader  rhino_surface( file_name );

      VSet<3U>  mesh_container;
      rhino_surface.OutputObjectTo( "fault1", mesh_container );
      // setting to a numerically integrated element as this is later needed for the transport calculation
      mesh_container.SingleElementType(ISOPARAMETRIC_LINEAR_TRIANGLE);

      PropertyData  permdata( ELEMENT, SCALAR, 3U );
      permdata.Reserve( mesh_container.Elements() );
      for ( auto i{0}; i<mesh_container.Elements(); ++i )
        pushBack( permdata, makeScalar( PLAIN, 1.0e-12 ) );
      mesh_container.AddData( "permeability", permdata );

      const bool isoparametric(true);
      Model<3U>  model3D( mesh_container, "example20.txt", isoparametric );
      mesh_container.Erase();

      printModelDimensions( model3D, true );

      // provisions for single surface model
      model3D.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
      model3D.InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN,0.) );
      model3D.InputPropertyValue( "fluid pressure",      makeScalar(PLAIN,0.) );
      // parallel plate permeability of a 1-mm fracture
      const double frac_perm(std::pow(0.001,2.)/12.);
      model3D.InputPropertyValue( "permeability",        makeScalar(PLAIN,frac_perm) );
      model3D.InputPropertyValue( "porosity",            makeScalar(PLAIN,1.) );
      model3D.InputPropertyValue( "concentration",       makeScalar(PLAIN,0.) );

      // boundary conditions applied on opposite sides of the model
      SideBoundaryConditions( model3D );
      // rectangular region with a concentration of 3 on the Rhino surface
      ConcentrationRectangle( model3D, 3. );

      VTK_Interface<3U>  vtk_output;

      vtk_output.OutputDataToVTK( model3D, "concentration",  "concentration", 0 );
      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );

      Standard_IO_Handler  stdio;

      printRangeOfVariable( model3D, stdio, "concentration" );
      printRangeOfVariable( model3D, stdio, "permeability" );


     // -----------------------------------------------------------------------
     // 0. hydraulic conductivity and other interrelations
     // -----------------------------------------------------------------------
      const double fluid_viscosity(1.6e-3);
      ConstantFactor<3U,divides>  conductivity( model3D.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model3D.Apply( conductivity );
      printRangeOfVariable( model3D, "conductivity" );


     // -----------------------------------------------------------------------
     // 1. steady-state fluid pressure
     // -----------------------------------------------------------------------
      #ifdef CSMP_WITH_SAMG_SOLVER
      PDE_Integrator<3U,Region>  steady_state_pressure(new SAMG_Solver());
      #else
      PDE_Integrator<3U,Region>  steady_state_pressure(new CSMP_DEFAULT_LINEAR_SOLVER());
      #endif

      NumIntegral_dNT_op_dN_dV<3U,Element<3U> >  conductance0( model3D.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
      NumIntegral_NT_op_N_dV<3U,Element<3U> >    source0( model3D.Database(), "fluid volume source", "fluid pressure" );
      VelocityAndVolumeFlux<3U,Element<3U> >     postpro0( model3D, "conductivity", "porosity", "fluid pressure" );

      steady_state_pressure.Add( &conductance0 );
      steady_state_pressure.Add( &source0 );
      steady_state_pressure.AddPostProcess( &postpro0 );

      model3D.Apply( steady_state_pressure );

      // output of results
      printRangeOfVariable( model3D, stdio, "fluid pressure" );
      printRangeOfVariable( model3D, stdio, "velocity" );
      printRangeOfVariable( model3D, stdio, "volume flux" );

      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );
      vtk_output.OutputDataToVTK( model3D, "velocity",       "velocity",       0 );
      vtk_output.OutputDataToVTK( model3D, "volume-flux",    "volume flux", 0 );


     // -----------------------------------------------------------------------
     // 2. advection of concentration field on fracture surfaces
     // -----------------------------------------------------------------------
     // TODO: include construction enabling simultaneous solution of diffusion problem
      NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", model3D, "porosity", "concentration",
                                                       "velocity", "nodal fluid volume source", true, true );

      double duration, cfl_mult, time_increment = advector.AnisotropicCourantIncrement();
      cout <<"\nmain: Enter advection time and CFL overstepping multiplier (CFL="<< time_increment <<" s): ";
      cin >> duration >> cfl_mult;

      advector.AdvectVariable( cfl_mult, 3, true, false );

      // testing the validity of the FV cells
      printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
      vtk_output.OutputDataToVTK( model3D, "flux-mismatch", "nodal flux mismatch", 0 );

      vtk_output.OutputDataToVTK( model3D, "concentration",   "concentration", 1 );

      cout <<"\nmain: That's it..."<< endl;

} // Run()





/**
     For a simple Rhino model, this method assigns nodal boundary values for 
     fluid pressure and concentration.
 
     The boundary nodes are identified on the basis of their spatial position.
*/
void RhinoMesh_Example::SideBoundaryConditions( Model<3U>& sg )
 {
   const csmp::Index  Skey = sg.Database().StorageKey("concentration");
   const csmp::Index  pkey = sg.Database().StorageKey("fluid pressure");
   Region<3>&   gref(sg.Region("Model"));
   const double tol(5.0e-1); // 50-cm match of the position of the nodes

   for ( auto nit=gref.NodesBegin(); nit!=gref.NodesEnd(); nit++ )
    {
      // for all boundary nodes
      // if x,y=0 a boundary value of 0. is assigned to the fluid pressure
      if ( fabs((*nit)->x() - 0.) <= tol and
           fabs((*nit)->y() - 0.) <= tol ) {
           (*nit)->Store( pkey, makeScalar(DIRICH,0.) );
           cout <<".";
        }
      // if x,y=50 a boundary value of 100 is assigned to the fluid pressure,
      // and the concentration is set to 1.
      if ( fabs((*nit)->x() - 100.25) <= tol and
           fabs((*nit)->y() - 100.25) <= tol ) {
           (*nit)->Store( pkey, makeScalar(DIRICH,141.42136) );
           (*nit)->Store( Skey, makeScalar(DIRICH,1.) );
           cout <<"~";
        }
    }
   cout <<"\n";

 } // end



/**
     Creates a rectangular high concentration region for later
     advection on the Rhino surface object.
*/
void  RhinoMesh_Example::ConcentrationRectangle( Model<3U>& sg, double concentration )
 {
   csmp::Index  Skey = sg.Database().StorageKey("concentration");
   Region<3>&  gref(sg.Region("Model"));

   for ( auto nit=gref.NodesBegin(); nit!=gref.NodesEnd(); nit++ )
     {
        if ( (*nit)->x() >= 20. && (*nit)->x() <= 80. &&
             (*nit)->y() >= 20. && (*nit)->y() <= 80. &&
             (*nit)->z() >= 10. && (*nit)->z() <= 40. )
         (*nit)->Store( Skey, makeScalar(PLAIN,concentration) );
     }
 }

} // csmp
