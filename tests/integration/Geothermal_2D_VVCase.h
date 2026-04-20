#ifndef GEOTHERMAL_2D_VVCASE_H
#define GEOTHERMAL_2D_VVCASE_H

#include "CSMP_definitions.h"
#include "Test.h"


using namespace std;
using namespace csmp;


namespace csmp
  {
  
template<uint32_t> class Model;

  /** Geothermal 2D Test Case
  =================================
  Mesh:       Linear Triangles
  Model:      2000x1000 m 2D
  Test:       Temperature
  BC:         constant pressure, constant temperature
  Criterion:
  =================================
  */
  class Geothermal_2D_VVCase : public Test
    {
      public:
        Geothermal_2D_VVCase(const char* prefix);
        virtual void run();
      private:
        void outputToVTU( Model<3U>& model,std::string model_name, const list<string>& props, size_t timestep );
        void ComputeMassConductivity (Model<3U>& model);
        void ComputeMassGravityTerm (Model<3U>& model);
    };

  } // csmp

#endif // GEOTHERMAL_2D_VVCASE_H
