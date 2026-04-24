#include "TemperatureDensityPressure_Example.h"

#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "ModelTopology.h"
#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"

#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "LinearSolver.h"

// H2O properties
#include "IAPWS_H2OPropertiesVisitor.h"
#include "CSMP_physical_constants.h"

#include "TextInterface.h"
#include "MapleInterface.h"

#include "LineElementMesher.h"
#include "VSetConverter.h"

#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

using namespace std;

namespace csmp {

void TemperatureDensityPressure_Example::Specifications()
{
  SetTitle( "1D vertical hydrostatic fluid-pressure profile taking into acount p,T-dependent H2O density" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "vertical (1D) temperature - fluid density - fluid pressure distribution" );
  AddDescription( "this example also illustrates the output of text files and Maple plots to paste into worksheets" );
  AddDescription( "source in: 'TemperatureDensityPressure_Example.cpp'" );
  AddRequirement( "variable file (example17.txt)");

}




/** *****************************************************************************************
 
   Computation of a hydrostatic (pore) pressure profile in a sedimentary basin.

   Illustration of how to use CSMP for one-dimensional models.
   Illustration of how to use the fluid equation of state (EOS) module.
 
    Steps: 
    
    1. First a vertical (1D) temperature profile is computed. This can be done in a single
       step because it is essentially independent on fluid density.
       
    2. Using an initial guess of the water density of 1000 kg/m3 an iterative process
       is started where rho(p,T) is adjusted incrementally using the new pressure estimates.
       If no further change occurs (within a given tolerance), 
       the (Picard) iteration is stopped.
        
   @attention Use Maple or Excel to visualize '*.mpl' and '*.txt' output files, respectively.
 
   *****************************************************************************************

*/
void TemperatureDensityPressure_Example::Run()
{
    // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
    //ostream &cout = *GetStream();

  //0. Create of directory with current example name, go into this directory, and copy input files into it.
    //find the name of current example source file
    string file_name = GetExampleFileName(__FILE__);
    string model_name; //empty - no copy required
    string variable_file = "example17.txt";

    CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file);

  // 1. builds 4km-tall 1D model
    // -------------------------------------------
    VSet<1U>       mesh_container;
    const size_t   N_ELEMENTS(4000);  // 4,000 meter tall model
    LineElementMesher<1U>   mesher;
    mesher.BuildUniformMesh( mesh_container, 4000., N_ELEMENTS );

    // creating a model topological region
    set<string>    fem_types; fem_types.insert("ISOPARAMETRIC_LINEAR_BAR");
    ModelTopology  mesh_topology(true);  // isoparametric

    // creating elements into the positive x direction (meaning upwards, see later)
    vector<size_t>  elms;
    for ( size_t i{0}; i<N_ELEMENTS; ++i ) elms.push_back(i);
    mesh_topology.AddDomain( "ROCK", fem_types, elms );
    elms.erase( elms.begin(), elms.end() );

    // to get isoparametric elements
    VSetConverter().ConvertElementTypesToOnesUsingLocalCoordinateSystem( mesh_container );
    Model<1U>  model( mesh_topology, mesh_container, "example17.txt", true );
    printModelDimensions( model );

    // 2. Input of material properties and initial conditions
    // ------------------------------------------------------
    model.InputPropertyValue( "conductivity",          makeScalar(PLAIN,1.0e-15) );
    model.InputPropertyValue( "fluid pressure",        makeScalar(PLAIN,101325.) ); // 1 bar
    model.InputPropertyValue( "thermal conductivity",  makeScalar(PLAIN,2.2) );


    // 3. Dirichlet boundary conditions for pressure and temperature
    // -------------------------------------------------------------
    cout <<"\nmain: Enter top and bottom temperature, and radiogenic heat source for the model: ";
    double T_top, T_bottom, heat_source;
    cin >> T_top >> T_bottom >> heat_source;

    model.InputBoundaryValue( CNR1, "temperature", makeScalar(DIRICH,T_bottom) );
    model.InputBoundaryValue( CNR2, "temperature", makeScalar(DIRICH,T_top) );
    model.InputPropertyValue( "energy source", makeScalar(PLAIN,heat_source) );


    // 4. compute vertical temperature profile
    // -------------------------------------------------------------
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver solver;
    PDE_Integrator<1U,Element>  temperature( solver );
    #else
    EigenSolver solver;
    PDE_Integrator<1U,Element>  temperature( solver );
    #endif

    NumIntegral_dNT_op_dN_dV<1U>  temperature_conductance( model.Database(), "thermal conductivity",  "temperature", "temperature" );
    NumIntegral_NT_op_N_dV<1U>    energy_source( model.Database(), "energy source", "temperature" );

    temperature.Add( &temperature_conductance );
    temperature.Add( &energy_source );

    model.Apply( temperature );
    printRangeOfVariable( model, "temperature" );


    // 5.  Visitor computes thermodynamic properties of water; density is interpolated to
    //     element barycentre for later vertical integration
    // --------------------------------------------------------------------------------------
    IAPWS_H2OPropertiesVisitor<1U>  properties_visitor( model, "fluid pressure" ,"fluid density","fluid viscosity");

    model.Accept( properties_visitor );
    model.InterpolateNodeToCellProperty( "fluid density", "element fluid density" );
  
  
    // 6. Setting up boundary conditions for the pressure iteration
    // --------------------------------------------------------------------------------------------
    csmp::Index     rhow_key = model.Database().StorageKey("element fluid density");
    csmp::Index     pfa_key  = model.Database().StorageKey("fluid pressure analytic");
    ScalarVariable  pf_top(DIRICH,101325.); // 1 bar
    // the last node (id=n_nodes) is located at the top of the model (by anology with the Y-axis)
    model.InputBoundaryValue( CNR1, "fluid pressure", pf_top );
    model.InputBoundaryValue( CNR2, "fluid pressure analytic", pf_top );
  

    // 7. Set up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
    // --------------------------------------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    SAMG_Solver    samg_solver(&settings);
    PDE_Integrator<1U,Element>  hydrostatic_pressure(samg_solver);
    // minimizing screen output
    settings.Set_iout1( 0 );
    settings.Set_iout2( 0 );
    settings.Set_idmp( -1 );
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<1U,Element>  hydrostatic_pressure(linear_solver);
#endif

    NumIntegral_dNT_op_dN_dV<1U>  hydrostatic_conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );

    cout <<"\nmain: enter the acceleration of gravity (kg/m.s2): in the area of interest: ";
    double acc_gravity(9.81);
    cin >> acc_gravity;
    // unless specified otherwise, in a 1D model, gravity will automatically act in the x-direction
    NumIntegral_NT_op_dNi_dV<1U>  hydrostatic_gravity( model.Database(), "element fluid density",
                                                      "conductivity", "fluid pressure", acc_gravity );
                                                      
    hydrostatic_pressure.Add( &hydrostatic_conductance );
    hydrostatic_pressure.Add( &hydrostatic_gravity );


    // 8. Compute the initial hydrostatic pressure. Since fluid properties will change after the pressure was computed
    //     recompute the fluid properties and fluid pressure 3 times such that fluid pressure and fluid properties converge
    // --------------------------------------------------------------------------------------------------------------------
    // fluid salinity is accounted for in the crudest fashion: the weight of the
    // total dissolved solids is simply added to the fluid density
    cout <<"\nmain: Enter the amount of total dissolved solids (ppm = g/tonne; normal seawater=35000 g/t): ";
    double total_dissolved_solids(12000./1000.); // g->kg (157500-ppm = 157kg salt)
    cin >> total_dissolved_solids;
    total_dissolved_solids /= 1000.; // gets kg/m3
    printRangeOfVariable( model, "element fluid density" );
    PropertyHandle<1U>  rhof( model, "element fluid density", SCALAR, ELEMENT );
    rhof.OutputCondition(ANY);
    rhof += total_dissolved_solids/1000.;
    printRangeOfVariable( model, "element fluid density" );

    cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid properties... " << endl;
    for ( uint32_t i=0; i<=5U; i++ ) {
         cout <<"\n\titeration "<< i+1U <<":"<< endl;
         model.Apply( hydrostatic_pressure );
         model.Accept( properties_visitor );
         model.InterpolateNodeToCellProperty( "fluid density", "element fluid density" );
         rhof += total_dissolved_solids;
         printRangeOfVariable( model, "fluid pressure" );
         printRangeOfVariable( model, "element fluid density" );
      }


    // 9. top down integration over density
    // -------------------------------------------
    Region<1U>& model_domain = model.Region("ROCK");
    Element<1U>* eptr = (*model_domain.CellVector().rbegin()); // top-most element
    const csmp::Index rhon_key  = model.Database().StorageKey("fluid density");
    ScalarVariable pf_analytic = pf_top;
    while(eptr->Neighbor(1))
      {
          pf_analytic += eptr->Volume()*eptr->Read(rhon_key) * 9.8;
          eptr->N(0)->Store(pfa_key, pf_analytic);
          eptr=eptr->Neighbor(1);
      }
    // first node of first element
    pf_analytic += eptr->Volume()*eptr->Read(rhon_key)*9.8;
    eptr->N(0)->Store(pfa_key, pf_analytic);


    // 10. Output of variables to screen and files
    // -------------------------------------------
    TextInterface  text_output;

    // plot these results with MS Excel or similar
    text_output.OutputDataAsTextColumns( model, "temperature",    "temperature", 0 );
    text_output.OutputDataAsTextColumns( model, "fluid-pressure", "fluid pressure", 0 );
    text_output.OutputDataAsTextColumns( model, "fluid-density",  "element fluid density", 0 );

    // plot these results with Maple
    writeVariableToMapleTextFile( model,  "fluid pressure analytic", 0, 0. );
    writeVariableToMapleTextFile( model,  "fluid pressure", 0, 0. );
    writeVariableToMapleTextFile( model, "temperature",     0, 0. );
    writeVariableToMapleTextFile( model, "fluid density",   0, 0. ); // without salt
    writeVariableToMapleTextFile( model, "element fluid density",   0, 0. ); // with salt
    writeVariablesToMapleTextFile( model, "temperature", "fluid density", 0, 0. );

    cout <<"\nmain: That's it..."<< endl;

    filesystem::current_path("../../example_inputs/");
  
} // end Run

} // csmp
