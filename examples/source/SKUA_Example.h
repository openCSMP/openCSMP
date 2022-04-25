#ifndef SKUA_EXAMPLE_H
#define SKUA_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  SKUA_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  static void ListModels();
  static void RunSuperSimple();
};

} // csmp

#endif // SKUA_EXAMPLE_H

