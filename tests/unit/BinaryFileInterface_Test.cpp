#include "BinaryFileInterface_Test.h"
#include "BinaryFileInterface.h"
#include "ANSYS_Model3D.h"
#include "VSet.h"
#include "PropertyHandle.h"
#include "FEM_Data.h"

using namespace std;

namespace csmp{

/**
     Tests the copy constructor of the Model and the correctness of the file reading routines.
     
     @todo SKM: the logic of this test is not clear and it fails probably because it is muddled.
*/
void BinaryFileInterface_Test::run()
{
  // 1. creating model initialised with a scalar permeability value of 1D
  // --------------------------------------------------------------------
  ANSYS_Model3D model3D( "prism_test", "CSMP-2phase-variables.txt", true );
  // remove once model copy operator becomes available
  ANSYS_Model3D model3Dcopy( "prism_test", "CSMP-2phase-variables.txt", true );
  Region<3>& region = model3D.Region( "Model" );
  Index permKey = model3D.Database().StorageKey( "permeability" );
  const double64 permValue( 1.0E-12 );
  const vector<Element<3>*>::const_iterator elmtsEnd = region.ElementsEnd();
  PropertyHandle<3> permeability( model3D, "permeability", SCALAR, ELEMENT );
  permeability = permValue;

  // 2. writing VData from this model to binary file
  // -----------------------------------------------
  BinaryFileInterface<3> toBin;
  // binary ops depend on vset binary ops, complete test there
  _test( toBin.WriteConnectivityFile( model3D, "BinaryFileInterface_Test3" ) );
  VSet<3> vSet;
  double64 time( 0. );
  _test( toBin.ReadConnectivityFile( "BinaryFileInterface_Test3", vSet, time ) );
  
  // 3. Reading the model data back from binary file and constructing a new model
  // ----------------------------------------------------------------------------
  const bool isoparametric_elements(true);
  Model<3> model3DvsetCopy( vSet, "CSMP-2phase-variables.txt", isoparametric_elements, false );
  // passed stays true, but is never ever used afterwards?
  bool passed( true );
  for( vector<Element<3>*>::const_iterator it=region.ElementsBegin(); it!=elmtsEnd; ++it )
    passed = ( fabs((*it)->Read(permKey) - permValue) <= numeric_limits<double64>::epsilon() );
  _test( toBin.WriteDataTo( model3D, "BinaryFileInterface_Test3", "permeability", 0 ) );
  
 // FAILS HERE: instead of the variable name, the string returned is the header of the binary file
  _test( toBin.ReadVariableName( "BinaryFileInterface_Test3" ) == "permeability" );
  
  FEM_Data<ScalarVariable> femData;
  // FAILS again because of the faulty string is recovered a second time
  _test( toBin.ReadDataFrom<ScalarVariable>( "BinaryFileInterface_Test3", model3Dcopy.Database(), femData ) == "permeability" );
  ScalarVariable min,max;
  femData.MinMaxOf( min, max );
  _equal( min(), permValue, numeric_limits<double64>::epsilon() );
  _equal( max(), permValue, numeric_limits<double64>::epsilon() );
}

} // csmp
