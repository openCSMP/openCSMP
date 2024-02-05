#ifndef TUTORIAL4_EXAMPLE_REVISITED_H
#define TUTORIAL4_EXAMPLE_REVISITED_H

#include "Example.h"

#include "Model.h"

namespace csmp {

/**
       As refactored by Hani Akbari (2018)
*/
class  Tutorial4_Example_revised : public Example
{

public:

  virtual void Run();
  virtual void Specifications();

private:

  void scaleRegion( Model<2U>& sg, double scale_factor );
  void constructVelocityVector( Model<2U>& mdl );
  void assignFluxToPointSource( Model<2U>& mdl, const char* flux );

};

} // csmp

#endif // TUTORIAL4_EXAMPLE_REVISITED_H
