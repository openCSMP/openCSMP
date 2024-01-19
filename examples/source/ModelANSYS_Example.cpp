#include "ModelANSYS_Example.h"

#include "ANSYS_Model3D.h"


using namespace std;

namespace csmp {

void ModelANSYS_Example::Specifications()
{
  SetTitle( "ANSYS_Model3D: loading mesh created by ANSYS(ICEM CFD Engineering Tetra)" );
  SetDifficulty( 1 );
  SetCategory( "Software Interfaces" );
  AddAuthor( "P. Lang" );
  AddDescription( "how to use CSMP's ANSYS interface" );
  AddDescription( "source in: ModelANSYS_Example.cpp" );
  AddRequirement( "a binary ansys mesh (.asc & .dat), e.g., fracs4");
  AddRequirement( "CSMP-1phase-variables.txt" );
} 


void ModelANSYS_Example::Run() {
  string model_name;
  cout << "\nPlease enter the name of input model, or press ENTER to use the default model 'fracs4':" << endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "fracs4";

  string variable_file("CSMP-1phase-variables.txt");

  //create a working directory with current example name, go into this directory, and copy input files into it.
  CopyInputFiles(model_name, variable_file);

  //build a model from the ANSYS binary files
  const bool irregular_mesh(true);
  ANSYS_Model3D model(model_name.c_str(), variable_file.c_str(), irregular_mesh);

  cout<<"\nANSYS model '"<<model_name<<"' is loaded successfully"<<endl;

  // now you need to onfigure the model from file and you can start your calculation

  fs::current_path("../../example_inputs/");

}


void ModelANSYS_Example::CopyInputFiles(std::string& model_name, std::string& variable_file) {
  //find the name of current example source file
  string example_name = GetExampleFileName(__FILE__);
  //create a working directory with the name of this example and go into it
  fs::create_directory("../example_outputs");
  fs::current_path("../example_outputs");
  if (fs::is_directory(example_name)) fs::remove_all(example_name); //if directory already exists, delete it
  fs::create_directory(example_name);
  fs::current_path(example_name);

  //copy ANSYS model files into working directory
  string input_directory = (fs::current_path().parent_path().parent_path()).string();
  input_directory += "/example_inputs/input_meshes/";

  string path = "../../example_inputs/input_meshes/";
  string name = model_name + ".asc";
  string file_name = path + name;
  if (fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    string error_message = "\n\nError: file '";
    error_message += (name + "' does not exist in directory " + input_directory);
    error_message += (", example cannot run, please copy input ANSYS model files (*.dat, *.asc) to this directory\n");
    throw std::runtime_error(error_message);
  }
  name = model_name + ".dat";
  file_name = path + name;
  if (fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    string error_message = "\n\nError: file '";
    error_message += (name + "' does not exist in directory " + input_directory);
    error_message += (", example cannot run, please copy input ANSYS model files (*.dat, *.asc) to this directory\n");
    throw std::runtime_error(error_message);
  }

  //copy variable file into working directory
  input_directory = (fs::current_path().parent_path().parent_path()).string();
  input_directory += "/example_inputs/variables_and_configuration_files/";

  path = "../../example_inputs/variables_and_configuration_files/";
  file_name = path + variable_file;
  if (fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    string error_message = "\n\nError: file '";
    error_message += (variable_file + "' does not exist in directory " + input_directory);
    error_message += (", example cannot run, please copy this file into this directory\n");
    throw std::runtime_error(error_message);
  }
}


} // csmp
