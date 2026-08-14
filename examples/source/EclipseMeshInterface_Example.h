#ifndef ECLIPSE_MESH_INTERFACE_EXAMPLE_H
#define ECLIPSE_MESH_INTERFACE_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  EclipseMeshInterface_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();

private:
  void CopyInputFiles(std::string& model_name, std::string& variable_file);
};

} // csmp

#endif // ECLIPSE_MESH_INTERFACE_EXAMPLE_H
