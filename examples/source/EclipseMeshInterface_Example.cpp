#include "EclipseMeshInterface_Example.h"


#include "CSMP_definitions.h"
#include "VTU_Interface.h"
#include "InputDataManager.h"

#include "EclipseModel_UoM.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

  using namespace eclipse;
  
void EclipseMeshInterface_Example::Specifications()
{
   SetTitle( "EclipseMeshInterface_Example" );
   SetDifficulty( 1 );
   SetCategory( "Software Functionality" );
   AddAuthor( "Roman Manasipov" );
   AddDescription( "source in: EclipseMeshInterface_Example.cpp" );
   AddDescription( "Example show's how to use EclispeModel constructor." );
   AddRequirement( "One has to have corresponding mesh files to be present.");
   AddRequirement( "by default: NPD5.grdecl, Johansen Data Set ( http://www.sintef.no/projectweb/matmora/downloads/johansen/ )");
}

void EclipseMeshInterface_Example::Run()
{
   // Model name
    std::string model_name;
    cout <<"\nPlease input model (mesh file prefix) name: ";
    cin  >> model_name;

    // Load Model
    std::string variables_file( "CSMP_Eclipse_example-variables.txt" );

    // Setup mesh
    std::string regions_file( model_name );
    bool exclude_inactive_cells( true );
    bool tetra_mesh( false );
    const bool create_boundaries( false );
    EclipseModelSettings settings( model_name );
    settings.MeshSetup( regions_file,
                        exclude_inactive_cells,
                        tetra_mesh,
                        create_boundaries );
    // Setup properties
    // porosity
    std::string   porosity_name("porosity");
    VARIABLE_TYPE porosity_type( SCALAR );
    PLACEMENT     porosity_place( ELEMENT );
    settings.PoroPropertySetup( porosity_name,
                                porosity_type,
                                porosity_place );
    // permeability
    std::string   permeability_name("tensor permeability");
    VARIABLE_TYPE permeability_type( TENSOR );
    PLACEMENT     permeability_place( ELEMENT );
    std::string   permeability_unit("mD"); // other options m2,D
    settings.PermPropertySetup( permeability_name,
                                permeability_type,
                                permeability_place,
                                permeability_unit );
    // rock type
    std::string   rocktype_name("rock type");
    VARIABLE_TYPE rocktype_type( SCALAR );
    PLACEMENT     rocktype_place( ELEMENT );
    settings.RockNumPropertySetup( rocktype_name,
                                   rocktype_type,
                                   rocktype_place );
  
    EclipseModel modelOut( settings, model_name, variables_file );

    // Get Fault Regions
    std::vector<std::string> faults;
    modelOut.GetFaults( faults );

    // Get Well Regions
    std::vector<std::string> wells;
    modelOut.GetWells( wells );

  // We are in the process of rewriting the Eclipse interface, and the
  // following part is not yet fully ported.  - AJB
  
  VTK_Interface<3U>  vtk_output;
  vtk_output.OutputNodeDataToVTK( modelOut, "EclipseInterfaceExample", 0 );

  #if 0
    // 2. Assign Fault and Matrix properties
    std::string          vol_source_name ( "fluid volume source" );
    csmp::ScalarVariable vol_source_val  ( csmp::PLAIN, 0.0 );
    csmp::ScalarVariable porosity_val    ( csmp::PLAIN, 1.0 );
    csmp::TensorVariable<3U> permeability_val( csmp::PLAIN, 0.0 );
    permeability_val(0,0) = 8.5e-11;
    permeability_val(1,1) = 8.5e-11;
    permeability_val(2,2) = 0.5*8.5e-11;

    modelOut.InputPropertyValue( vol_source_name.c_str(), vol_source_val );
    if( !faults.empty() )
    {
        vol_source_val = -1.0;
        const size_t num_faults( faults.size() );
        for( size_t fid = 0; fid<num_faults; ++fid )
        {
            vol_source_val += 1.0e-3;
            modelOut.Region( faults[fid] ).InputPropertyValue( vol_source_name.c_str(),   vol_source_val,   csmp::COMPLETE );
            modelOut.Region( faults[fid] ).InputPropertyValue( porosity_name.c_str(),     porosity_val,     csmp::COMPLETE );
            modelOut.Region( faults[fid] ).InputPropertyValue( permeability_name.c_str(), permeability_val, csmp::COMPLETE );
        }
    }
    if( !wells.empty() )
    {
        vol_source_val = 1.0;
        const size_t num_wells( wells.size() );
        for( size_t wid = 0; wid<num_wells; ++wid )
        {
            vol_source_val -= 1.0e-3;
            modelOut.Region( wells[wid] ).InputPropertyValue( vol_source_name.c_str(),   vol_source_val,   csmp::COMPLETE );
            modelOut.Region( wells[wid] ).InputPropertyValue( porosity_name.c_str(),     porosity_val,     csmp::COMPLETE );
            modelOut.Region( wells[wid] ).InputPropertyValue( permeability_name.c_str(), permeability_val, csmp::COMPLETE );
        }
    }

    // 3. Output properties
    // TODO: only output properties found in the Eclipse file and not a pre-defined list of them!

    std::string problem_title       ( "EclipseMeshInterface_Example" );
    std::string output_file_name    ( model_name );
    std::string model_subdomain_part( "Model" );
    std::list<std::string> properties;
    properties.push_back( porosity_name.c_str() );
    properties.push_back( permeability_name.c_str() );
    properties.push_back( vol_source_name.c_str() );
    properties.push_back( rocktype_name.c_str() );
    size_t time_step( 0 );
    VTU_Interface<3U>  vtu_output1( modelOut, problem_title.c_str() );
    vtu_output1.OutputDataToVTU( output_file_name.c_str(), properties, model_subdomain_part.c_str(), time_step );

    // 4. Output to csmp native binary format
    modelOut.OutputToBinaryFile( model_name.c_str() );
  
    // 5. reading the model back into memory
    //Model<3U> modelIn( model_name.c_str() );
    //VTU_Interface<3U>  vtu_output2( modelIn, problem_title.c_str() );
    //vtu_output2.OutputDataToVTU( output_file_name.c_str(), properties, model_subdomain_part.c_str(), time_step );
#endif

    cout <<"\nEclipseMeshInterface_Example: That's it!\n";
}

} // csmp
