#ifndef CSMP_FINITE_VOLUME_STENCIL_SPEED_TEST_H
#define CSMP_FINITE_VOLUME_STENCIL_SPEED_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp
  {
     template<uint32_t> class Model;
  /**
      @brief benchmarking the FiniteVolumeStencil container class underpinning CSMP's transport schemes.
      
       Compares the access speed of FiniteVolumeStencil with a refactored version based on DynamicArray2D and DynamicArray3D
       
       @attention a tetrahedral / triangular element mesh is used to to keep  variable access cost to a mininum
   */
  class FiniteVolumeStencilSpeed_Test : public Test
    {
       public:
         virtual void run();
         
         void AssignFlowProperties( Model<3U>& );
         void DivergenceFreeTotalVelocityField(  Model<3U>&, double delta_pf );
         void TestFlowThroughModel( Model<3U>& model, bool prescribed_velocity );
         
       private:
         const bool verbose_ = true;
    };

  } // csmp

#endif
