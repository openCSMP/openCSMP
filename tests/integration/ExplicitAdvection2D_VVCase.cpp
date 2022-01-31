#include "ExplicitAdvection2D_VVCase.h"
#include "ModelTopology.h"

// CSMP Files
#include "Model.h"
#include "Region.h"
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
//#include "PropertyHandle.h"

// File I/O and Initialization
#include "ANSYS_Interface.h"
#include "MeshDiagnostics.h"

// PDE operators & solvers
#include "VelocityAndVolumeFlux.h"

// outputting
#include "VTK_Interface.h"

// finite-volume stuff
// explicit scheme
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {
  ExplicitAdvection2D_VVCase::ExplicitAdvection2D_VVCase()
  {
  }

  ExplicitAdvection2D_VVCase::ExplicitAdvection2D_VVCase(const char* prefix)
  {
    this->setName("ExplicitAdvection2D_VVCase");
    prefix_=prefix;
  }

  ExplicitAdvection2D_VVCase::~ExplicitAdvection2D_VVCase()
  {
  }


  // ****************************************************************************************
  // Testing Advection via CSMP's generic CVFE method
  // *****************************************************************************************
  void ExplicitAdvection2D_VVCase::run()
  {
    //tolerance
    double fTolerance(1.e-7);
    // isoparametric parameter
    const bool       isoparametric(true);
    ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
    VSet<2U>         mesh_container;
    ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements
    string input_file_name(prefix_);

    cout <<"\nRunning Test Case Simulation..."<< endl;
    // Building Region object from ANSYS data files
    cout <<"Reading mesh..."<<endl;
    const bool binary_file( true );
    mesh_interface.Read_ANSYS_Mesh( prefix_, mesh_container, mesh_topology, binary_file, true );
    cout <<"Finished reading mesh..."<<endl;
    cout <<"Building Model..."<<endl;

    Model<2U>  model2D( mesh_topology, mesh_container, (input_file_name+".txt").c_str() );
    cout <<"Finished building model."<<endl;

    mesh_container.Erase();
    mesh_topology.Erase();
    cout <<"Cleared mesh_container and mesh_topology..."<<endl;

    // Checks
    Standard_IO_Handler  stdio;
    printModelDimensions( model2D, true );
    MeshDiagnostics<2U>  mesh_check;
    double             vol_min, vol_max;
    mesh_check.ElementVolumeRange( model2D, vol_min, vol_max );

    // Configuring the built model (adding boundary condition info, etc, from input files:
    // *-configuration.txt
    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile( model2D, prefix_,
                                           false, true, true, true, false );

    printRangeOfVariable( model2D, stdio, "permeability" );
    printModelDimensions( model2D, true );

    // -----------------------------------------------------------------------
    // Build prescribed velocity fuild to be used in tests.
    // -----------------------------------------------------------------------
    VectorVariable<2U>  v_const(PLAIN,PLAIN,5.0e-6, 0.0e-6);
    model2D.InputPropertyValue("velocity", v_const );
    VTK_Interface<2U>  vtk_output;
    vtk_output.OutputDataToVTK( model2D, "velocity", "velocity", 1, true );

    // -----------------------------------------------------------------------
    // 1. Test of the Explicit transport algorithm (1rst order)
    // -----------------------------------------------------------------------

    ExplicitNodeCenteredFiniteVolumeTransport<2U,ExplicitStencilProcessor>  explicit_advector(
        "Model", model2D,
        "porosity",
        "concentration",
        "velocity",
        "nodal fluid volume source",
        true );
    printRangeOfVariable( model2D, stdio, "concentration" );
    cout <<"\nConfiguring Explicit TRANSPORT simulation (1rst Order): ";
    cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement() << endl;
    //cout <<"\nEnter advection time: ";
    double time_interval;
    //cin >> time_interval;
    time_interval=10000U;

    for(int i=0U; i<5; ++i)
    {
      explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
      vtk_output.OutputDataToVTK( model2D, "concentration", "concentration", static_cast<int>(i+1U), true );
    }
    // Difference diagnostics will go here.
    //call diagnostics function here..
    _equal( 0, 0., fTolerance );
    return;

  } // end


} // namespace csmp
























