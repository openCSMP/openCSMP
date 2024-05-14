#include "Averaging_Example.h"

#include "ANSYS_Model2D.h"
#include "Region.h"
#include "VTU_Interface.h"
#include "PL_Utilities.h"

using namespace std;

namespace csmp {

void Averaging_Example::Specifications()
  {
    SetTitle( "Extrapolation & Averaging" );
    SetDifficulty( 2 );
    SetCategory( "Software Functionality" );
    AddAuthor( "P. Lang" );
    AddDescription( "Comparison of element to node extrapolation approaches" );
    AddDescription( "Illustration of VTU Interface" );
    AddDescription( "source in: Averaging_Example.cpp" );
    AddRequirement( "file set: 'LeftRight'");
    AddRequirement( "CSMP-1phase-variables.txt");
  }



/**
    Visualisation of element-to-node extrapolation options within csmp.
*/
void Averaging_Example::Run()
{
  string model_name;

#ifdef USE_ANSYS_INPUT_FILE
  // choose an input model with large variations in element size to maximise the difference between the different
  // extrapolation approaches
  cout<< "\nPlease enter the name of ANSYS input model: "<<endl;
  cin >> model_name;
  // initializing model and properties
  ANSYS_Model2D model( model_name.c_str(), "CSMP-1phase-variables.txt", true );
#else
  cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'LeftRight':"<<endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "LeftRight";

  //find the name of current example source file
  string file_name = GetExampleFileName(__FILE__);
  string variable_file = "CSMP-1phase-variables.txt";
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file);
  //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model<2U>  model(model_name, variable_file);
#endif
  PropertyHandle<2U> sourceSink( model, "fluid volume source", SCALAR, ELEMENT );
  PropertyHandle<2U> nodalSourceSinkByDistance( model, "nodal fluid volume source distance", SCALAR, NODE );
  PropertyHandle<2U> nodalSourceSinkByVolume( model, "nodal fluid volume source volume", SCALAR, NODE );
  PropertyHandle<2U> nodalSourceSinkByCount( model, "nodal fluid volume source count", SCALAR, NODE );

  // resetting target variables
  nodalSourceSinkByDistance = 0.;
  nodalSourceSinkByVolume = 0.;
  nodalSourceSinkByCount = 0.;

  // setting source variables
  sourceSink = 1.;

  // setting up VTU interface for visualization
  VTU_Interface<2U> vtu( model, "Extrapolation Comparison" );

  // creating a reference to the model region
  Region<2U>& region( model.Region( "Model" ) );

  // using the ModelSubDomain methods to extrapolate single-valued field
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source distance" );
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source volume", false );

  // output
  vtu.OutputDataToVTU( "ExtrapolationByDistance", "nodal fluid volume source distance", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByVolume", "nodal fluid volume source volume", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByCount", "nodal fluid volume source count", "Model", static_cast<int>(0) );

  // creating a perturbed field and extrapolating this
  randomPerturb( model, "fluid volume source", 20. );

  // weight_element_influence_by_distance_of barycentre to node
  const bool by_distance{true};
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source distance", by_distance );
  // weighting by element area (opposite effect)
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source volume", false );

  vtu.OutputDataToVTU( "PerturbedFluidVolumeSource", "fluid volume source", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByDistance", "nodal fluid volume source distance", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByVolume", "nodal fluid volume source volume", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByCount", "nodal fluid volume source count", "Model", static_cast<int>(1) );

  filesystem::current_path("../../example_inputs/");

}

} // csmp
