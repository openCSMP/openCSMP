#ifndef DES_ADVECTIONOFTRACER3D_EXAMPLE_H
#define DES_ADVECTIONOFTRACER3D_EXAMPLE_H

#include "Example.h"

#include "CSMP_number_types.h"

namespace csmp {

  template<size_t dim> class Model;
  template<size_t dim> class VTK_Interface;

class  DESAdvectionOfTracer3D_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
};

} // csmp

#endif // DES_ADVECTIONOFTRACER3D_EXAMPLE_H
