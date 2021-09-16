#include "BinaryFileInterface_Test.h"
#include "BinaryFileInterface.h"
#include "ANSYS_Model3D.h"
#include "VSet.h"
#include "PropertyHandle.h"
#include "FEM_Data.h"

using namespace std;

namespace csmp {

/**
     Tests the copy constructor of the Model and the correctness of the file reading routines.
     
     SKM: refactored test 1/3/2017.
*/
void BinaryFileInterface_Test::run()
{
  // 1. creating model initialised with a scalar permeability value of 1D
  // --------------------------------------------------------------------
  ANSYS_Model3D model3D( "prism_test", "CSMP-1phase-variables.txt", true );

  // 2. writing mesh as VData from this model to binary file
  // -------------------------------------------------------
  BinaryFileInterface<3> toBin;
  // binary ops depend on vset binary ops, complete test there
  _test( toBin.WriteConnectivityFile( model3D, "BinaryFileInterface_Test3" ) );
  VSet<3> vSet;
  double64 time( 0. );
  // connectivity only, no data
  _test( toBin.ReadConnectivityFile( "BinaryFileInterface_Test3", vSet, time ) );
  
  // 3. Reading the mesh back from binary file and constructing a new model
  // ----------------------------------------------------------------------
  const bool isoparametric_elements(true);
  Model<3> model3DvsetCopy( vSet, "CSMP-1phase-variables.txt", isoparametric_elements );
  // assigning a permeability
  const ScalarVariable permValue(ANY,1.0e-12);
  model3DvsetCopy.InputPropertyValue( "permeability", permValue );

  _test( toBin.WriteDataTo( model3DvsetCopy, "BinaryFileInterface_Test4", "permeability", 0 ) );
  
 // FAILS HERE: instead of the variable name, the string returned is the header of the binary file
  _test( toBin.ReadVariableName( "BinaryFileInterface_Test4" ) == "permeability" );
  
  FEM_Data<ScalarVariable> femData;
  // FAILS again because of the faulty string is recovered a second time
  _test( toBin.ReadDataFrom<ScalarVariable>( "BinaryFileInterface_Test4", model3DvsetCopy.Database(), femData ) == "permeability" );
  ScalarVariable min,max;
  femData.MinMaxOf( min, max );
  _equal( min(), permValue(), numeric_limits<double64>::epsilon() );
  _equal( max(), permValue(), numeric_limits<double64>::epsilon() );
}

} // csmp
