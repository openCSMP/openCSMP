#ifndef LINEAR_ELASTICITY_A_VV_CASE_H
#define LINEAR_ELASTICITY_A_VV_CASE_H

#include "Test.h"

namespace csmp
  {
  
  /** Linear Elasticity Test Case A
  =================================
  Mesh:       Linear Triangles
  Model:      Beam 40x2 m 2D
  Test:       Tensile stress
  BC:         const displ / const force
  Criterion:  analytical solution: strain, stress, poisson deformation
  =================================
  */
  class LinearElasticityA_VVCase : public Test
    {
    public:
      LinearElasticityA_VVCase(const char* prefix);
      virtual void run();
    };

  } // csmp

#endif
