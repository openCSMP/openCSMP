#ifndef GEOTHERMAL_PSEUDO1D_VVCASE_H
#define GEOTHERMAL_PSEUDO1D_VVCASE_H

#include "Test.h"
#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;

  /** Geothermal pseudo 1D Test Case
  =================================
  Mesh:       Linear Triangles
  Model:      2000x1000 m 2D
  Test:       Temperature
  BC:         constant pressure, constant temperature
  Criterion:
  =================================
  */
  class Geothermal_pseudo1D_VVCase : public Test
    {
      public:
        Geothermal_pseudo1D_VVCase(const char* prefix);
        virtual void run();
      private:
        void outputToVTU( Model<2U>& model,std::string model_name, const std::list<std::string>& props, size_t timestep );
        void ComputeMassConductivity ( Model<2U>& model );
    bool Compare (Model<2U>& model, std::string file );
    };

  } // csmp

#endif // GEOTHERMAL_PSEUDO1D_VVCASE_H
