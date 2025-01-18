#ifndef AVERAGING_EXAMPLE_H
#define AVERAGING_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// element to node interpolation using different approaches; see ErrorMetric_Example for a practical application
class Averaging_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
};

} // csmp


#endif // AVERAGING_EXAMPLE_H
