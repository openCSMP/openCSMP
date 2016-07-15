#ifndef LINEAR_ELASTICITY_C_TEST_CASE_H
#define LINEAR_ELASTICITY_C_TEST_CASE_H

#include "Test.h"

namespace csmp
  {
  
  /** Linear Elasticity Test Case C
  =================================
  Mesh:       Linear Bars
  Model:      Beam 3 m 1D
  Test:       Tensile strain
  BC:         const displ / const force
  Criterion:  analytical solution: strain
  =================================
  */
  class LinearElasticityC_VVCase : public Test
    {
    public:
      virtual void run();
    };

  } // csmp

#endif