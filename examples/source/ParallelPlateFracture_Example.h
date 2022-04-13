#ifndef CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H
#define CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  ParallelPlateFracture_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // CSMP_PARALLEL_PLATE_FRACTURE_EXAMPLE_H
