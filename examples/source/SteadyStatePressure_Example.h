#ifndef STEADY_STATE_PRESSURE_EXAMPLE_H
#define STEADY_STATE_PRESSURE_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  SteadyStatePressure_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // STEADY_STATE_PRESSURE_EXAMPLE_H
