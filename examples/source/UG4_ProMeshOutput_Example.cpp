#include "UG4_ProMeshOutput_Example.h"

#include "UG4_UGX_FileExport.h"

#include "Model.h"
#include "Region.h"

// File I/O and Initialization
#include "TRIANGLE_Interface.h"
#include "VTU_Interface.h"
#include "InputDataManager.h"
#include "Standard_IO_Handler.h"
#include "Boundary.h"
#include "VSet.h"
#include "Box.h"

using namespace std;

namespace csmp {

void UG4_ProMeshOutput_Example::Specifications()
{
  SetTitle( "Output of CSMP Model to ProMesh Pre-Processor for UG4" );
  SetDifficulty( 1 );
  SetCategory( "Software Interfaces" );
  AddAuthor( "SKM" );
  AddDescription( "Demonstrates the use of interfaces/UG4_UGX_FileExport" );
  AddDescription( "source in: UG4_ProMeshOutput_Example.cpp" );
  AddRequirement( "input file set: 'blunt30deg.1', blunt30deg.1-configuration.txt, UG4_ProMeshOutput_Example-variables.txt or other files");
}



/** *****************************************************************************************

   Reads a model from the desired InterFace and outputs it as a '*.ugx' file for processinfg it with UG4's ProMesh graphical
   preprocessor..

   *****************************************************************************************

*/
void UG4_ProMeshOutput_Example::Run()
{
    // ---------------------------------------------------------------------------------------
    // 0. Import model from Shewchuk's Triangle mesher
    //    (for other input formats, see CSMPInterfaces_Example, or use CSMP
    // ---------------------------------------------------------------------------------------
    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;
    string              file_name;
    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
    cin >> file_name;
    const bool isoparametric{true};
    mesh_interface.ReadTriangle2DMesh( file_name.c_str(), mesh_container, isoparametric );

    // 'UG4_ProMeshOutput_Example-variables.txt' is the text file that defines the variables used in this example
    Model<2U>  model( mesh_container, "UG4_ProMeshOutput_Example-variables.txt" );
    string     config_file{ file_name };
    // give this model a name as none is supplied to the constructor
    model.Name( file_name.c_str() );
    mesh_container.Erase();
    printModelDimensions( model, true );
    printRangeOfVariable( model, "permeability", true );
    
    // 0.1 Fix potentially incorrect boundary flagging (if the model is not perfectly rectangular)
    // ------------------------------------------------------------------------------------------
    // create csmp::Boundary objects consisting of Face objects
    model.EstablishBoxBoundariesFromOrientation();
    // using the newly generated boundary faces, fix potentially broken corner and side flags
    recreateBoxBoundaryFlags( model );
    
    // 0.2 since Triangle interface assigns the values from the .ele file to the variable 'permeability' but we need 'rocktype'
    model.Region("Model").CopyReplace( "permeability", "rocktype" );
    
    
    // ---------------------------------------------------------------------------------------
    // 1. Define region names and assign properties to model subregions
    // ---------------------------------------------------------------------------------------
    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile( model, config_file.c_str() );
    // verifying that the regions were created correctly
    model.RegionsOut();
    // verifying the property distribution
    // defining input properties
    list<string> assigned_properties;
    assigned_properties.emplace_back( "rocktype" );
    assigned_properties.emplace_back( "porosity" );
    assigned_properties.emplace_back( "permeability" );
    //assigned_properties.emplace_back( "vertical permeability" );
    assigned_properties.emplace_back( "brooks corey parameter" );
    assigned_properties.emplace_back( "residual saturation aqueous phase" );
    assigned_properties.emplace_back( "residual saturation carbonic phase" );
    assigned_properties.emplace_back( "entry pressure" );
    assigned_properties.emplace_back( "fluid pressure" );
    const long time{0ul};

    VTU_Interface<2>  vtu(model);
    vtu.OutputDataToVTU( file_name, assigned_properties, "Model", time );


    // ---------------------------------------
    // 2. Output of model to UG4 (.ugx) format
    // ---------------------------------------
    UG4_UGX_FileExport<2> ugx_interface( model );

    ugx_interface.Write_UGX_FileASCII( model, model.Name() );
    
    
    // ---------------------------------------
    // 3. Output model to CSMP's native format
    // ---------------------------------------
    model.OutputToBinaryFile( file_name.c_str() );

    cout <<"\nmain: That's it..."<< endl;

} // end Run




/* TEST CODE FOR BOUNDARY FLAGS

    cout <<"\nmain: checking validity of boundary flags:"<< endl;
    for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); ++bit ) {
         cout <<"\n"<<"'"<< (*bit).first <<"', with boudary flags: "<< endl;
         for ( auto b : (*bit).second.NodeVector() )
           cout <<"  "<< parseBoundary(b->AtBoundary());
           
         if ( (*bit).first == "BOTTOM" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              if ( flags.find(CNR1) == flags.end() ) cout <<" CNR1 is missing. ";
              if ( flags.find(CNR2) == flags.end() ) cout <<" CNR2 is missing. ";
           }
         else if ( (*bit).first == "RIGHT" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              if ( flags.find(CNR2) == flags.end() ) cout <<" CNR2 is missing. ";
              if ( flags.find(CNR3) == flags.end() ) cout <<" CNR3 is missing. ";
           }
         else if ( (*bit).first == "TOP" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              if ( flags.find(CNR3) == flags.end() ) cout <<" CNR3 is missing. ";
              if ( flags.find(CNR4) == flags.end() ) cout <<" CNR4 is missing. ";
           }
         else if ( (*bit).first == "LEFT" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
             if ( flags.find(CNR1) == flags.end() ) cout <<" CNR1 is missing. ";
             if ( flags.find(CNR4) == flags.end() ) cout <<" CNR4 is missing. ";
           }
      }
    cout << endl;
    
*/


} // csmp
