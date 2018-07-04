#ifndef LINEAR_ELASTICITY_D_TEST_CASE_H
#define LINEAR_ELASTICITY_D_TEST_CASE_H

#include "Test.h"

namespace csmp
  {
 
  /** Linear Elasticity Test Case D
  =================================
  Mesh:       Linear Triangles
  Model:      Beam 40x2 m 2D
  Test:       Tensile stress
  BC:         const displ / const disp
  Criterion:  analytical solution: strain, stress, poisson deformation
  =================================
  */
  class LinearElasticityD_VVCase : public Test
    {
    public:
      virtual void run();
    };

  } // csmp

#endif