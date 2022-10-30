#ifndef CSMP_GEOTHERMAL_1D_VVCASE_H
#define CSMP_GEOTHERMAL_1D_VVCASE_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

  /** Geothermal pseudo 1D Test Case
  =================================
  Mesh:       Linear Line Elements
  Model:      2000 m in horizontal direction
  Test:        Temperature
  BC:          Constant pressure, constant temperature
  Criterion: Comparison with TOUGH
  =================================
  */
  class Geothermal_1D_VVCase : public Test
    {
      public:
        Geothermal_1D_VVCase(); ///< sets name of Test
        virtual void run();
        
      private:
        void OutputToVTU( Model<1U>&, std::string model_name, const std::list<std::string>& props, size_t timestep );
        void ComputeMassConductivity( Model<1U>& );
		    bool Compare ( Model<1U>& model, std::string file );
      
        // VTU_Interface<1U> vtu(model);
    };

  } // csmp

#endif // GEOTHERMAL_PSEUDO1D_VVCASE_H
