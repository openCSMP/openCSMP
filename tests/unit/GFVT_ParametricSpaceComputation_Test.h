//
//  GenericFiniteVolumeTransport_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/01/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_IN_PARAMETRIC_SPACE_TEST_H
#define CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_IN_PARAMETRIC_SPACE_TEST_H

#include "Test.h"

namespace csmp {

/**
      Compares computations done entirely in parametric space and then mapped to global space as
      volumes, areas, and vectors, with the physical space computations that are currently used in CSMP.
      
      @attention FiniteVolumePolicy and FiniteVolumeStencil have their own dedicated tests that pass without issues.
      This test was written with the goal to conduct the transport computations directly in parametric space, even using
      the distorted pressure gradients there (which are needed instead of velocities in the case of tensor properties).
      The advantage of this approach would be that quadrilateral facets are always planar, which leads to volume conservation issues
      with the current physical space computations.
      
      @todo currently, the parametric space computations are not correct and need fixing.
 */
class GFVT_ParametricSpaceComputation_Test : public Test {
  public:
    virtual void run();
  
    void TestBasics();
    void BenchmarkGlobalVersusParametricIntegration();
  
  private:
    constexpr static bool verbose_ = false;
};



} // end csmp


#endif /* CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_IN_PARAMETRIC_SPACE_TEST_H */
