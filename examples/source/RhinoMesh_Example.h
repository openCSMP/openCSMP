#ifndef RHINOMESH_EXAMPLE_H
#define RHINOMESH_EXAMPLE_H

#include "Example.h"
#include "CSMP_definitions.h"

namespace csmp {

  template<size_t> class Model;

class  RhinoMesh_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  void SideBoundaryConditions( Model<3U>& );
  void ConcentrationRectangle( Model<3U>&, double64 );
};

} // csmp

#endif // RHINOMESH_EXAMPLE_H
