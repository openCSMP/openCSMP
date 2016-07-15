#ifndef MODELANSYS_EXAMPLE_H
#define MODELANSYS_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<size_t> class Model;

class  ModelANSYS_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // MODELANSYS_EXAMPLE_H
