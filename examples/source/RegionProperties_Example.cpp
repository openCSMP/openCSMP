#include "RegionProperties_Example.h" 
#include "vset_makers.h"
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
   AddRequirement( "'b25' .dat, .asc, -regions & -configuration .txt");
   AddRequirement( "example1_regions.txt" );
} 


void RegionProperties_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

    enum {DIM=2};
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
    fem_type.insert("QUAD_4");
    ModelTopology   mesh_topology( "variable-access-test model", isoparametric );
    vector<size_t>  elms; 
    elms.reserve(n_squares_on_side);
    for ( size_t i=0U; i<n_squares_on_side; i++ ) elms.push_back(i);
    mesh_topology.AddRegion( "region1", fem_type, elms );
    elms.erase( elms.begin(), elms.end() );
    elms.reserve(mesh_container.Elements()-n_squares_on_side);
    for ( size_t i=n_squares_on_side; i<mesh_container.Elements(); i++ ) elms.push_back(i);
    mesh_topology.AddRegion( "region2", fem_type, elms );
    elms.erase( elms.begin(), elms.end() );
    mesh_topology.Out();
 
    // constructing the model with the constructor for ANSYS meshes 
    Model<DIM>  model( mesh_topology, mesh_container, "example1.txt" );
    
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
    
    cout <<"\nmain: setting K to to permeability:\n";
    ArithmeticMean<DIM,ScalarVariable>  avg( model.Database(), "hydraulic conductivity", "permeability" );
    model.Apply( avg );
    printRangeOfVariable( model, "hydraulic conductivity", true );
    
    cout <<"\nmain: dividing it by a constant viscosity:\n";    
    const double64  viscosity(1.6e-3);
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
    printRangeOfVariable( model, "average permeability", true );
   
   cout <<"\nRegionProperties_Example: That's it..."<< endl;
  
} // end RegionProperties_Example::run

} // end namespace csmp


