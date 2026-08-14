#ifndef CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H
#define CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// Stokes - Darcy flow through fracture network in a porous rock
class  ParallelPlateFracture_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H
