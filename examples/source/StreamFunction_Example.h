#ifndef STREAMFUNCTION_EXAMPLE_H
#define STREAMFUNCTION_EXAMPLE_H

#include "Example.h"

#include "Box.h"
#include "PDE_Integrator.h"
#include "Region.h"

namespace csmp {

  template<uint32_t> class Model;
  template<uint32_t> class Interrelation;
  class Standard_IO_Handler;


class  StreamFunction_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
private:
  void computeStreamFunction( Model<2U>& sg,
                              BOX_BOUNDARY boundary0, BOX_BOUNDARY boundary1,
                              double total_flux, const char* stream_func_var );

  void analyze_sensitivity( Model<2U>& sg, const char* group, Standard_IO_Handler& io,
                            Interrelation<2U>& itr, PDE_Integrator<2U,Region>& algo );

  template<uint32_t dim> double integrateDomainBoundaryFlux( Model<dim>& sg );


};

} // csmp

#endif // STREAMFUNCTION_EXAMPLE_H
