// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Visitor_Example.h"

#include "Model.h"
#include "Region.h"
#include "ANSYS_Interface.h"
#include "VSet.h"
#include "ModelTopology.h"
#include "VSetConverter.h"

#include "PressureSaturationInitializer.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp{

void Visitor_Example::Specifications()
{
  SetTitle( "Visitor" );
  SetDifficulty( 2 );
  SetCategory( "Software Functionality" );
  AddAuthor( "Shaho" );
  AddDescription( "source in: Visitor_Example.cpp" );
  AddDescription( "application and implementation of csmp::Visitor" );
  AddRequirement( "LeftRight (CSMP binary files)" );
  AddRequirement( "variable file (VisitorExample-var.txt)" );
}



void Visitor_Example::Run()
{
  /*
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  const bool          isoparametric(true);
  ANSYS_Interface     mesh_interface(isoparametric);
  VSet<2U>            mesh_container;
  ModelTopology       mesh_topology(isoparametric);
  VSetConverter<2U>   mesh_converter;

  string model_name("undefined");
  cout <<"\nEnter name of model: ";
  cin >> model_name;
  const bool binary_file( true );
  mesh_interface.Read_ANSYS_Mesh( model_name.c_str(), mesh_container, mesh_topology, binary_file, true );

  const bool get_domain_info_from_regions_file{true};
  Model<2U> reservoir_model( mesh_topology, mesh_container, "VisitorExample-var.txt", get_domain_info_from_regions_file );
  */

  string model_name;
  cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'LeftRight':"<<endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "LeftRight";

  //find the name of current example source file
  string file_name = GetExampleFileName(__FILE__);
  string variable_file = "VisitorExample-var.txt";
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file);
  //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model<2U>  reservoir_model(model_name, variable_file);

  PressureSaturationInitializer<2U> pressure_saturation_initializer ( reservoir_model,
                                                                      "saturation water", "saturation oil",
                                                                      "water pressure", "oil pressure",
                                                                      1000., 700., 3., 5., 1e5, 2, 1e4, 0., 0.);
  
  std::cout << "\nPressureSaturationInitializer::Visit(Node): Visiting reservoir model.\n";    
  reservoir_model.Accept( pressure_saturation_initializer );
  std::cout << "\nPressureSaturationInitializer::Visit(Node): Pressures and saturations initialized successfully!\n";
  VTK_Interface<2U>  vtk_output;

  vtk_output.OutputDataToVTK( reservoir_model, "oil-pressure", "oil pressure", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "water-pressure", "water pressure", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "saturation-oil", "saturation oil", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "saturation-water", "saturation water", 0 );

  filesystem::current_path("../../example_inputs/");
  
} // Run()

} // csmp
