#ifndef TUTORIAL4_EXAMPLE_REVISITED_H
#define TUTORIAL4_EXAMPLE_REVISITED_H

#include "Example.h"

namespace csmp {

/// incomplete forward declaration
template<uint32_t> class Model;

/**
    Single-phase flow through the pore-space of a sandstone as modelled with the Stokes lubrication equation.
    @author Jan Zaretsky, model later refactored by Hani Akbari
    @date 2018
*/
class  Tutorial4_Example : public Example
{

public:

  virtual void Run();
  virtual void Specifications();

private:

  void scaleRegion( Model<2U>&, double scale_factor );
  void constructVelocityVector( Model<2U>& );
  void assignFluxToPointSource( Model<2U>&, const char* flux );

};

} // csmp

#endif // TUTORIAL4_EXAMPLE_REVISITED_H
