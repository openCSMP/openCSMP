#ifndef LINEAR_ELASTICITY_B_TEST_CASE_H
#define LINEAR_ELASTICITY_B_TEST_CASE_H

#include "Test.h"

namespace csmp {
 
  /** LinearElasticity Test Case B
  ================================
  Mesh:       2D Linear Triangles
  Model:      Beam 40x2 m
  Test:       Deflection
  BC:         const displ / const force
  Criterion:  analytical solution: deflection
  ================================
  */
  class LinearElasticityB_VVCase : public Test
    {
    public:
      virtual void run();
    };

  } // csmp

#endif
