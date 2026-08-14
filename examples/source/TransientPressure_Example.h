#ifndef TRANSIENTPRESSURE_EXAMPLE_H
#define TRANSIENTPRESSURE_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  TransientPressure_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // TRANSIENTPRESSURE_EXAMPLE_H
