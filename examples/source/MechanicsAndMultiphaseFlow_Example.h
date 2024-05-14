#ifndef MECHANICS_AND_MULTIPHASE_FLOW_EXAMPLE_H
#define MECHANICS_AND_MULTIPHASE_FLOW_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<uint32_t> class Model;

class  MechanicsAndMultiphaseFlow_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  void SteadyStatePressure( Model<2U>& );
  void PrintModelProperties( const Model<2U>& );
};

} // csmp

#endif // MECHANICS_AND_MULTIPHASE_FLOW_EXAMPLE_H
