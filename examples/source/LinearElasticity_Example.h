// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_LINEAR_ELASTICITY_EXAMPLE_H
#define CSMP_LINEAR_ELASTICITY_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<uint32_t> class Model;

class  LinearElasticity_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  void SteadyStatePressure( Model<2U>& );
  void PrintModelProperties( const Model<2U>& );
};

} // csmp

#endif // CSMP_LINEAR_ELASTICITY_EXAMPLE_H
