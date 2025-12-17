#include "ModelComparator_Test.h"

#include "ANSYS_Model3D.h"
#include "PropertyHandle.h"
#include "ModelComparator.h"

namespace csmp{

void ModelComparator_Test::run()
{
  ANSYS_Model3D model( "Cube", "CSMP-ModelComparator_Test-variables.txt", true );

  PropertyHandle<3U> permeability( model, "permeability", SCALAR, ELEMENT );
  PropertyHandle<3U> fluid_pressure( model, "fluid pressure", SCALAR, NODE );
  permeability = 1.0E-12;
  fluid_pressure = 1.0e+5;
  model.OutputToBinaryFile( "ComparatorTestVSet1" );
  permeability = 1.0E-13;
  fluid_pressure = 1.5e+5;
  model.OutputToBinaryFile( "ComparatorTestVSet2" );

  ModelComparator<3U> comparator;
  double shouldBeZero( comparator.CompareVSets( "ComparatorTestVSet1.vset", "ComparatorTestVSet1.vset",
                                                "permeability", "permeability",
                                                "CSMP-1phase-variables.txt", "CSMP-1phase-variables.txt" ) );
  _equal( shouldBeZero, 0., 1.0E-3 );

  double shouldBeAlsoZero( comparator.CompareVSetsAtPoints( "ComparatorTestVSet1.vset", "ComparatorTestVSet1.vset",
                                                            "fluid pressure", "fluid pressure",
                                                            "CSMP-1phase-variables.txt", "CSMP-1phase-variables.txt" ) );
  _equal( shouldBeAlsoZero, 0., 1.0E-3 );

}

} // csmp
