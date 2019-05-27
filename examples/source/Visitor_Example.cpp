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
  AddRequirement( "LeftRight .asc,.dat,-regions.txt" );
  AddRequirement( "Visitor_Example-var.txt" );
}



void Visitor_Example::Run()
{
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
  const bool irregular_mesh( false );
  mesh_interface.Read_ANSYS_Mesh( model_name.c_str(), mesh_container, mesh_topology, binary_file, irregular_mesh );

  Model<2U> reservoir_model( mesh_topology, mesh_container, "VisitorExample-var.txt" );


  PressureSaturationInitializer<2U> pressure_saturation_initializer ( reservoir_model,
                                                                      "saturation water", "saturation oil",
                                                                      "water pressure", "oil pressure",
                                                                      1000., 700., 3., 5., 1e5, 2, 1e4, 0., 0.);


  pressure_saturation_initializer.Visit( &reservoir_model );
  VTK_Interface<2U>  vtk_output;

  vtk_output.OutputDataToVTK( reservoir_model, "oil-pressure", "oil pressure", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "water-pressure", "water pressure", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "saturation-oil", "saturation oil", 0 );
  vtk_output.OutputDataToVTK( reservoir_model, "saturation-water", "saturation water", 0 );
} // Run()

} // csmp
