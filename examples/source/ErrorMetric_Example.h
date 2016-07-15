#ifndef ERRORMETRIC_EXAMPLE_H
#define ERRORMETRIC_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<size_t> class Model;

class  ErrorMetric_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
private:
  void assignLargestEigenValueOfTo( Model<3U>&, const char* of_var, const char* to_var );

  void discretizationError3D( Model<3U>&, const char* hessian_var,
                              const char* error_var, const char* sc_err_var );

};

} // csmp

#endif // ERRORMETRIC_EXAMPLE_H
