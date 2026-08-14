#ifndef POLICYBASED_EXAMPLE_H
#define POLICYBASED_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  PolicyBased_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // POLICYBASED_EXAMPLE_H
