#include "RegionProperties_Example.h" 
#include "vsetMakers.h"
#include "PropertyHandle.h"
#include "Model.h"
#include "Region.h"

#include "ArithmeticMean.h"
#include "ConstantFactor.h"

#include "VSet.h"
#include "ModelTopology.h"

#include "CSMP_highLevelUtilities.h"
#include "TextInterface.h"
#include "VTK_Interface.h"



using namespace std;

namespace csmp {

void RegionProperties_Example::Specifications()
{
   SetTitle( "Properties associated with csmp::Region objects" );
   SetDifficulty( 2 );
   SetCategory( "Software Functionality" );
   AddAuthor( "SKM" );
   AddDescription( "source in: RegionProperties_Example.cpp" );
   AddDescription( "application of interrelations that use REGION variables etc." );
   AddRequirement( "example1.txt" );
} 


void RegionProperties_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

    enum {DIM=2};

    //find the name of current example source file
    string example_name = GetExampleFileName(__FILE__);
    //create a working directory with the name of this example and go into it
    fs::create_directory("../example_outputs");
    fs::current_path("../example_outputs");
    if (fs::is_directory(example_name)) fs::remove_all(example_name); //if directory already exists, delete it
    fs::create_directory(example_name);
    fs::current_path(example_name);

    //copy variable file into working directory
    string variable_file ("example1.txt" );
    string path = "../../example_inputs/variables_and_configuration_files/";
    string file_name = path + variable_file;
    if (fs::exists(file_name)) fs::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      string input_directory = fs::current_path().parent_path().parent_path();
      input_directory += "/example_inputs/variables_and_configuration_files/";
      error_message += (variable_file + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }

    // ---------------------------------------------------------------------------------------------
    // 1. making a Model of squares to examine interpolations to element barycentres more readily
    // ---------------------------------------------------------------------------------------------
    VSet<DIM>     mesh_container;
    const size_t  n_squares_on_side(4U);
    const bool    skewed(false);
    const bool    isoparametric(true); // quadrilateral exists only in isoparametric form
    test_Create_Square_VSet( mesh_container, n_squares_on_side, DIM, skewed );
    
    // making two extra regions
    set<string>  fem_type;
    fem_type.insert("ISOPARAMETRIC_LINEAR_QUADRILATERAL");
    ModelTopology   mesh_topology( "variable-access-test model", isoparametric );
    vector<size_t>  elms; 
    elms.reserve(n_squares_on_side);
    for ( auto i{0}; i<n_squares_on_side; i++ ) elms.push_back(i);
    mesh_topology.AddDomain( "region1", fem_type, elms );
    elms.erase( elms.begin(), elms.end() );
    elms.reserve(mesh_container.Elements()-n_squares_on_side);
    for ( size_t i=n_squares_on_side; i<mesh_container.Elements(); i++ ) elms.push_back(i);
    mesh_topology.AddDomain( "region2", fem_type, elms );
    elms.erase( elms.begin(), elms.end() );
    mesh_topology.Out();
 
    // constructing the model with the constructor for ANSYS meshes 
    Model<DIM>  model( mesh_topology, mesh_container, "example1.txt", true );


    // assigning material properties 
    // -----------------------------------------------------  
    model.InputPropertyValue( "porosity", makeScalar(PLAIN,0.) );
    // element property
    model.Region("region1").InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-12) );
    model.Region("region2").InputPropertyValue( "permeability", makeScalar(PLAIN,2.0e-12) );
    // node property
    model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,1e5) );
    model.Region("region1").InputPropertyValue( "fluid pressure", makeScalar(PLAIN,1.0e6) );
    
    cout <<"\nmain: creating and initialising the 'hydraulic conductivity' K as constraint point property:\n";
    PropertyHandle<DIM>  K( model, "hydraulic conductivity", SCALAR, ELEMENT_INTEGRATION_POINT );
    
    // step by step initialisation of the new 'hydraulic conductivity' variable that lives on the element integratioon points
    cout <<"\nmain: setting K to to permeability:\n";
    // assign element property 'permeability' to integration points
    ArithmeticMean<DIM,ScalarVariable>  avg( model.Database(), "hydraulic conductivity", "permeability" );
    model.Apply( avg );
    printRangeOfVariable( model, "hydraulic conductivity", true );
    // divide it by the node variable 'viscosity' to obtain 'hydraulic conductivity' on integration points
    cout <<"\nmain: dividing it by a constant viscosity:\n";    
    const double  viscosity(1.6e-3);
    ConstantFactor<DIM,divides>  Kdiv_mu( model.Database(), "hydraulic conductivity", "hydraulic conductivity", viscosity );
    model.Apply( Kdiv_mu );
    printRangeOfVariable( model, "hydraulic conductivity", true );
 
    // text files to examine constraint-point and regional variables
    // -------------------------------------------------------------
    TextInterface  txt_output;
    txt_output.OutputDataAsTextColumns( model, "hydraulic-conductivity", "hydraulic conductivity" );
 
    VTK_Interface<DIM>  vtk_output;
    //vtk_output.OutputNodeDataToVTK( model, "region1", "node_vars", 0 );
    vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
    vtk_output.OutputDataToVTK( model, "hydraulic-conductivity", "hydraulic conductivity", 0 );

    // ---------------------------------------------------------------------------------------------
    // 2. Working with regional variables
    // ---------------------------------------------------------------------------------------------
    cout <<"\nmain: creating a new property 'average permeability' and computing it with ArithmeticMean:\n";
    // region property
    PropertyHandle<DIM>  kavg( model, "average permeability", SCALAR, REGION );
    ArithmeticMean<DIM,ScalarVariable>  ravg( model.Database(), "average permeability", "permeability" );
    model.Apply( ravg );
    //printRangeOfVariable( model, "average permeability", true );
   
   cout <<"\nRegionProperties_Example: That's it..."<< endl;

   fs::current_path("../../example_inputs/");
  
} // end RegionProperties_Example::run

} // end namespace csmp


