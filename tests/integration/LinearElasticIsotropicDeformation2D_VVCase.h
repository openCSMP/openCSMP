#ifndef LINEAR_ELASTIC_ISOTROPIC_DEFORMATION_2D_VV_CASE_H
#define LINEAR_ELASTIC_ISOTROPIC_DEFORMATION_2D_VV_CASE_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp
  {

  class LinearElasticIsotropicDeformation2D_VVCase: public Test
    {
    public:
      LinearElasticIsotropicDeformation2D_VVCase(const char* prefix);

      void SetupLinearTriangles();
      void SetupQuadraticTriangles();

      void SetupLinearQuadrilaterals();
      void SetupQuadraticQuadrilaterals();

      void Solve();
      virtual void run();

  private:
      bool isLinear_;
      bool isTriangle_;
      const char* prefix_;
      bool plane_strain_ = true;

    };

  } // csmp

#endif
