#ifndef INCOMPRESSIBLE_SINGLE_PHASE_FLOW_TENSOR_PERM_3DFEM_VV_CASE
#define INCOMPRESSIBLE_SINGLE_PHASE_FLOW_TENSOR_PERM_3DFEM_VV_CASE

#include "Test.h"

namespace csmp
  {
  
/** Finite Element Incompressible Single Phase Porous Media Flow Test Case
  =================================
  Element type(s) : Linear Tets
  Target Mesh set : UnitCubeFine
  Test            : Steady State Pressure (Incompressible Single Phase Porous Media Flow)
                    1) Dirichlet (const Pressure) on bottom-left and top-right corners.
                    Criterion:  Comparison with analytical solution: pressure, flow velocity.

  =================================
  */
  class IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase : public Test
    {
    public:
      IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase(const char* prefix);
      virtual void run();
    };

  } // csmp

#endif
