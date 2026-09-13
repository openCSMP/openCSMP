#ifndef PERIODIC_BOUNDARY_CONDITION_EXAMPLE_H
#define PERIODIC_BOUNDARY_CONDITION_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  PeriodicBoundaryCondition_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // PERIODIC_BOUNDARY_CONDITION_EXAMPLE_H
