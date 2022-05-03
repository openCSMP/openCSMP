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
  static void RunSKUA_model();
  static void RunSKUA_box_shaped_with_boundary();
  static void RunSKUA_split_boundary_layer();
};

} // csmp

#endif // SKUA_EXAMPLE_H

