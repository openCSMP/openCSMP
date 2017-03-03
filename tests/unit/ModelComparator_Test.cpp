#include "ModelComparator_Test.h"

#include "ANSYS_Model3D.h"
#include "PropertyHandle.h"
#include "ModelComparator.h"

namespace csmp{

ModelComparator_Test::ModelComparator_Test()
{
}

void ModelComparator_Test::run()
{
  ANSYS_Model3D model( "Cube", "CSMP-1phase-variables.txt", true );

  PropertyHandle<3U> permeability( model, "permeability", SCALAR, ELEMENT );
  PropertyHandle<3U> fluid_pressure( model, "fluid pressure", SCALAR, NODE );
  permeability = 1.0E-12;
  fluid_pressure = 1.0e+5;
  model.OutputToBinaryFile( "ComparatorTestVSet1" );
  permeability = 1.0E-13;
  fluid_pressure = 1.5e+5;
  model.OutputToBinaryFile( "ComparatorTestVSet2" );

  ModelComparator<3U> comparitor;
  double64 shouldBeZero( comparitor.CompareVSets( "ComparatorTestVSet1.vset", "ComparatorTestVSet1.vset",
                                                  "permeability", "permeability",
                                                  "CSMP-1phase-variables.txt","CSMP-1phase-variables.txt",
                                                  true ) );
  _equal( shouldBeZero, 0., 1.0E-3 );

  double64 shouldBeAlsoZero( comparitor.CompareVSetsAtPoints( "ComparatorTestVSet1.vset", "ComparatorTestVSet1.vset",
                                                  "fluid pressure", "fluid pressure",
                                                  "CSMP-1phase-variables.txt", "CSMP-1phase-variables.txt",
                                                   true ) );
  _equal( shouldBeAlsoZero, 0., 1.0E-3 );

}

} // csmp
