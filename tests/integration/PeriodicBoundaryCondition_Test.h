#ifndef CSMP_PERIODIC_BOUNDARY_CONDITION_TEST_H
#define CSMP_PERIODIC_BOUNDARY_CONDITION_TEST_H

#include "Test.h"

namespace csmp {

/**
    Verifies PDE_Integrator's periodic ('master_node_id'-linked) DOF handling
    for steady-state Darcy flow through a heterogeneous fracture network,
    mirroring the two boundary condition setups offered by the
    'PeriodicBoundaryCondition_Example':

      1) single periodic direction (periodic top/bottom, Dirichlet left/right);
      2) double periodic (both pairs of edges periodic, with a pressure jump
         driving flow across one pair, and the four corners Dirichlet-pinned).
*/
class PeriodicBoundaryCondition_Test : public Test {
public:
  virtual void run();

private:
  void TestSinglePeriodicDirection();
  void TestDoublePeriodicDirections();
};

} // csmp

#endif // CSMP_PERIODIC_BOUNDARY_CONDITION_TEST_H
