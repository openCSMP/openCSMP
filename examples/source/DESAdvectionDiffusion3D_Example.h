#ifndef DES_ADVECTIONDIFFUSION3D_EXAMPLE_H
#define DES_ADVECTIONDIFFUSION3D_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<uint32_t dim> class Model;
  template<uint32_t dim> class VTK_Interface;

class  DESAdvectionDiffusion3D_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
};

} // csmp

#endif // DES_ADVECTIONDIFFUSION3D_EXAMPLE_H
