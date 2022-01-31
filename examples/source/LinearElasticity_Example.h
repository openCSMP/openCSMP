#ifndef LINEARELASTICITY_EXAMPLE_H
#define LINEARELASTICITY_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<size_t> class Model;

class  LinearElasticity_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
private:
  void SteadyStatePressure( Model<2U>& );
};

} // csmp

#endif // LINEARELASTICITY_EXAMPLE_H
