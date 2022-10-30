#ifndef LINEAR_ELASTIC_FRACTURE_APERTURE_2D_VVCASE_H
#define LINEAR_ELASTIC_FRACTURE_APERTURE_2D_VVCASE_H

#include "Test.h"
#include "stdlib.h"
#include "vector"
#include "SAMG_Settings.h"

namespace csmp
  {

  class LinearElasticFractureAperture2D_VVCase: public Test
    {
    public:
      LinearElasticFractureAperture2D_VVCase(const char* prefix);
      virtual void run();

      void set_min_avg_error(double epsilon);

      template<typename T>
      void WriteSolutionToFile( std::vector<std::vector<T>> data, std::string name);

      void SetSettings( SAMG_Settings& settings);

  private:
      double min_avg_error_ = 10.0;

    };

  } // csmp

#endif
