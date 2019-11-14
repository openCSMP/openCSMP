#include "ModelANSYS_Example.h"

#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp {

void ModelANSYS_Example::Specifications()
{
  SetTitle( "ANSYS_Model3D: loading mesh created by ANSYS(ICEM CFD Engineering Tetra)" );
  SetDifficulty( 1 );
  SetCategory( "Software Functionality" );
  AddAuthor( "P. Lang" );
  AddDescription( "how to use CSMP's ANSYS interface" );
  AddDescription( "source in: ModelANSYS_Example.cpp" );
  AddRequirement( "a binary ansys mesh (.asc & .dat)");
  AddRequirement( "CSMP-1phase-variables.txt" );
} 


void ModelANSYS_Example::Run()
{
  string modelName;
  cout << "\nModel Name (ex: fracs4): ";
  cin >> modelName;
  cout << endl << endl;
  const bool irregular_mesh(true);
  ANSYS_Model3D model( modelName.c_str(), "CSMP-1phase-variables.txt", irregular_mesh );
  // now you need to onfigure the model from file and you can start your calculation
}

} // csmp
