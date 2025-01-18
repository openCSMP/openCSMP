#ifndef INCOMPRESSIBLE_SINGLE_PHASE_FLOW_FEM_VV_CASE
#define INCOMPRESSIBLE_SINGLE_PHASE_FLOW_FEM_VV_CASE

#include "Test.h"

namespace csmp
  {
  
/** Finite Element Incompressible Single Phase Porous Media Flow Test Case
=================================
Mesh:       Linear Triangles
Model:      40x2 m 2D
Test:       Steady State Pressure (Incompressible Single Phase Porous Media Flow)
BC:         1) Dirichlet (const Pressure) on RIGHT and LEFT boundaries. TOP and BOTTOM natural.
            2) Dirichlet (const Pressure) on bottom-left and top-right corners.
            3) Dirichlet (const Pressure) on bottom-left corner, NEUMANN (const rate)
Criterion:  Comparison with analytical solution: pressure, flow velocity.

=================================
*/
  class IncompressibleSinglePhaseFlowFEM_VVCase : public Test
    {
    bool verbose_;
    public:
      explicit IncompressibleSinglePhaseFlowFEM_VVCase(const char* prefix, bool verbose=true );
      virtual void run();
    };

  } // csmp

#endif  // INCOMPRESSIBLE_SINGLE_PHASE_FLOW_FEM_VV_CASE
