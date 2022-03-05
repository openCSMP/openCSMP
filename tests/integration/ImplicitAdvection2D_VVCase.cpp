#include "ImplicitAdvection2D_VVCase.h"

// CSMP Files
#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"
#include "InputDataManager.h"
#include "PropertyHandle.h"

// File I/O and Initialization
#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "MeshDiagnostics.h"

// PDE operators & solvers
#include "VelocityAndVolumeFlux.h"

// Interrelations
#include "ConstantFactor.h"

// outputting
#include "VTK_Interface.h"

// finite-volume scheme
#include "NodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {

/** 
    *****************************************************************************************
      2D advection via CSMP's generic CVFE method
    ***************************************************************************************** 
*/
int test_implicit_advection_2D()
  {
    // isoparametric parameter
    const bool       isoparametric(true);
    ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
    VSet<2U>         mesh_container;
    ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements

    // Building Region object from ANSYS data files
    cerr <<"Reading mesh..."<<endl;
    const bool binary_file( true );
    mesh_interface.Read_ANSYS_Mesh( "square", mesh_container, mesh_topology, binary_file, true );
    cerr <<"Finished reading mesh..."<<endl;
    cerr <<"Building Model..."<<endl;
    Model<2U>  model2D( mesh_topology, mesh_container, "square.txt" );
    cerr <<"Finished building model."<<endl;

    mesh_container.Erase();
    mesh_topology.Erase();
    cerr <<"Cleared mesh_container and mesh_topology..."<<endl;

    // Checks
    Standard_IO_Handler  stdio;
    printModelDimensions( model2D, true );
    MeshDiagnostics<2U>  mesh_check;
    double               vol_min, vol_max;
    mesh_check.ElementVolumeRange( model2D, vol_min, vol_max );

    // Configuring the built model (adding boundary condition info, etc, from input files:
    // *-configuration.txt
    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile( model2D, "square",
                                           false, true, true, true, false );

    printRangeOfVariable( model2D, stdio, "permeability" );
    printModelDimensions( model2D, true );

    // -----------------------------------------------------------------------
    // Build prescribed velocity field to be used in tests.
    // -----------------------------------------------------------------------
    VectorVariable<2U>  v_const(PLAIN,PLAIN,5.0e-6, 0.0e-6);
    model2D.InputPropertyValue("velocity", v_const );
    VTK_Interface<2U>  vtk_output;
    vtk_output.OutputDataToVTK( model2D, "velocity", "velocity", 1 );

    // -----------------------------------------------------------------------
    // 1. Test of the Implicit transport algorithm (1rst order)
    // -----------------------------------------------------------------------

    //-------------------------------

    NodeCenteredFiniteVolumeTransport<2U>  advector( "Model", model2D,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source",false, false );

    cout <<"\nConfiguring TRANSPORT simulation... Prescribed velocity"<< endl;
    cout <<"\nThe grid Courant number is "<< advector .AnisotropicCourantIncrement() << endl;
    cout <<"\nEnter advection time and Courant multiplier: ";
    double Courant_multiplier, time_interval;
    cin >> time_interval >> Courant_multiplier;
    //for(int i=0; i<5; ++i)
    //{
    advector.AdvectVariable( time_interval, Courant_multiplier );
    //vtk_output.OutputDataToVTK( model2D, "concentration", "concentration", i+1 );
    //}




    //        FiniteVolumeTransport_PrescribedVelocity( model2D );

    //        vtk_output.OutputDataToVTK( model2D, "fluid-pressure", "fluid pressure",    1 );

    //        vtk_output.OutputDataToVTK( model2D, "concentration", "concentration", 1 );

    //        FiniteVolumeTransport_PrescribedVelocity( model2D );

    //        vtk_output.OutputDataToVTK( model2D, "fluid-pressure", "fluid pressure",    2 );
    //        vtk_output.OutputDataToVTK( model2D, "velocity",       "velocity",          2 );
    //        vtk_output.OutputDataToVTK( model2D, "concentration", "concentration", 2 );

    //        FiniteVolumeTransport_PrescribedVelocity( model2D );

    //        vtk_output.OutputDataToVTK( model2D, "fluid-pressure", "fluid pressure",    3 );
    //        vtk_output.OutputDataToVTK( model2D, "velocity",       "velocity",          3 );
    //        vtk_output.OutputDataToVTK( model2D, "concentration", "concentration", 3 );

    //        printRangeOfVariable( model2D, stdio, "fluid pressure" );
    //        printRangeOfVariable( model2D, stdio, "velocity" );
    //        printRangeOfVariable( model2D, stdio, "pore velocity" );
    //        printRangeOfVariable( model2D, stdio, "concentration" );


    cout <<"\nShort_Tests: ADvection 2D_test is done."<< endl;

    return 0;

  } // end



} // namespace csmp
























