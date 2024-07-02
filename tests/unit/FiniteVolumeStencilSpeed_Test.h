#ifndef CSMP_FINITE_VOLUME_STENCIL_SPEED_TEST_H
#define CSMP_FINITE_VOLUME_STENCIL_SPEED_TEST_H

#include "Test.h"

namespace csmp
  {

  /**
      @brief benchmarking the FiniteVolumeStencil container class underpinning CSMP's transport schemes.
      
       Compares the access speed of FiniteVolumeStencil with a refactored version based on DynamicArray2D and DynamicArray3D
   */
  class FiniteVolumeStencilSpeed_Test : public Test
    {
    public:
      virtual void run();
    };

  } // csmp

#endif
