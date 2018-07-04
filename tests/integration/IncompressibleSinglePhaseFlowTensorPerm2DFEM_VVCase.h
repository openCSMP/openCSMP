#ifndef INCOMPRESSIBLE_SINGLE_PHASE_FLOW_TENSOR_PERM_FEM_VV_CASE
#define INCOMPRESSIBLE_SINGLE_PHASE_FLOW_TENSOR_PERM_FEM_VV_CASE

#include "Test.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

namespace csmp
{

/** Finite Element Incompressible Single Phase Porous Media Flow Test Case
  =================================
  Element type(s) : Linear Triangles
  Target Mesh set : simple_square_domain_unstructured_triangles
  Test            : Steady State Pressure (Incompressible Single Phase Porous Media Flow)
                    1) Dirichlet (const Pressure) on bottom-left and top-right corners.
                    Criterion:  Comparison with analytical solution: pressure, flow velocity.

  =================================
  */
class IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase : public Test
{
public:
    IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase(const char* prefix);
    virtual void run();
};

} // csmp

#endif  // INCOMPRESSIBLE_SINGLE_PHASE_FLOW_FEM_VV_CASE
