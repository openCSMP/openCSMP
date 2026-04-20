#ifndef LINEAR_ELASTICITY_E_TEST_CASE_H
#define LINEAR_ELASTICITY_E_TEST_CASE_H

#include "Test.h"

namespace csmp {
  
  /** Linear Elasticity Test Case E
  =================================
  Mesh:       Linear Bars
  Model:      Beam 3 m 1D
  Test:       Tensile strain
  BC:         const displ / const displ
  Criterion:  analytical solution: strain
  =================================
  */
  class LinearElasticityE_VVCase : public Test
    {
    public:
      virtual void run();
    };

  } // csmp

#endif
