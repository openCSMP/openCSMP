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
  SetCategory( "Meshing Interfaces" );
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
     string model_name;
     cout << "\nPlease enter the name of input model, or press ENTER to use the default model 'example20':" << endl;
     cin.ignore();
     getline(cin, model_name);
     if (model_name.length() == 0) model_name = "example20";

     std::string variables_file("example20.txt");

     // create a working directory with current example name, go into this directory, and copy input files into it.
     CopyInputFiles(model_name, variables_file);

     auto file_name = model_name.c_str();

     // -----------------------------------------------------------------------
     // 0. reading in '.raw' meshes that were written to file as labeled
     //    entities
     // -----------------------------------------------------------------------
      RhinoSurfaceReader  rhino_surface( file_name );

      VSet<3U>  mesh_container;
      rhino_surface.OutputObjectTo( "fault1", mesh_container );
      // setting to a numerically integrated element as this is later needed for the transport calculation
      mesh_container.SingleElementType(ISOPARAMETRIC_LINEAR_TRIANGLE);

      PropertyData  permdata( ELEMENT, SCALAR, 3U );
      permdata.Reserve( mesh_container.Elements() );
      for ( auto i{0}; i<mesh_container.Elements(); ++i )
        pushBack( permdata, makeScalar( PLAIN, 1.0e-12 ) );
      mesh_container.AddData( "permeability", permdata );

      Model<3U>  model3D( mesh_container, "example20.txt" );
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
      SAMG_Solver solver;
      PDE_Integrator<3U,Element>  steady_state_pressure(solver);
      #else
      CSMP_DEFAULT_LINEAR_SOLVER solver;
      PDE_Integrator<3U,Element>  steady_state_pressure(solver);
      #endif

      NumIntegral_dNT_op_dN_dV<3U>  conductance0( model3D.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
      NumIntegral_NT_op_N_dV<3U>    source0( model3D.Database(), "fluid volume source", "fluid pressure" );
      VelocityAndVolumeFlux<3U>     postpro0( model3D, "conductivity", "porosity", "fluid pressure" );

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

      filesystem::current_path("../../example_inputs/");

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


  void RhinoMesh_Example::CopyInputFiles(std::string& model_name, std::string& variable_file) {
    //find the name of current example source file
    string example_name = GetExampleFileName(__FILE__);
    //create a working directory with the name of this example and go into it
    filesystem::create_directory("../example_outputs");
    filesystem::current_path("../example_outputs");
    if (filesystem::is_directory(example_name)) filesystem::remove_all(example_name); //if directory already exists, delete it
    filesystem::create_directory(example_name);
    filesystem::current_path(example_name);

    //copy model files into working directory
    string input_directory = (filesystem::current_path().parent_path().parent_path()).string();
    input_directory += "/example_inputs/input_meshes/";

    string path = "../../example_inputs/input_meshes/";
    string name = model_name + ".raw";
    string file_name = path + name;
    if (filesystem::exists(file_name)) filesystem::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      error_message += (name + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }

    name = model_name + ".3dm";
    file_name = path + name;
    if (filesystem::exists(file_name)) filesystem::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      error_message += (name + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }

    //copy variable file into working directory
    input_directory = (filesystem::current_path().parent_path().parent_path()).string();
    input_directory += "/example_inputs/variables_and_configuration_files/";

    path = "../../example_inputs/variables_and_configuration_files/";
    file_name = path + variable_file;
    if (filesystem::exists(file_name)) filesystem::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      error_message += (variable_file + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }
  }

} // csmp
