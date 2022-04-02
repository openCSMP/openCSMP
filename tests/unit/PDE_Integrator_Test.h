#ifndef CSMP_PDE_INTEGRATOR_TEST_H
#define CSMP_PDE_INTEGRATOR_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "PDE_Integrator.h"
#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t> class Region;
template<uint32_t> class Element;


/**

TODO

@todo test that the PDE operators are accumulated into the right places in the lefthand matrix
@todo test accumulation of boundary integrals
@todo test coupling of domains by SplitBoundary integrals
@todo test accumulation of a coupled system
@todo test that the contacts between domain and boundaries are identified correctly for the accumulation of surface integrals at the boundary
@todo move all the testing that relates to matrix inversion into Solver_Test and specific subclasses
@todo test repeated use of integrator in a time-dependent problem:  matrix retention vs. reconstruction, test change of DOF from step to step

*/
class PDE_Integrator_Test : public Test {
  public:
    explicit PDE_Integrator_Test( Model<2U>& model );
    explicit PDE_Integrator_Test( Model<3U>& model );
    ~PDE_Integrator_Test();
    void run();

  private:
    Model<2U>& model_;
    PDE_Integrator<2U, Region>* pde_reference_ = nullptr;         // for validating
    PDE_Integrator<2U, Region>* pde_test_      = nullptr;         // for testing
    
    const static bool verbose_ = true;

    /*
    =================> LHS:
    NumIntegral_BT_D_B_dV
    NumIntegral_BT_D_op_dV
    NumIntegral_dNT_dN_dV
    NumIntegral_dNT_mixed_op_dN_dV
    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV
    NumIntegral_dNT_op_dN_dV                                      1
    NumIntegral_dNT_op_dN_dV_NT_v_dN_dV
    NumIntegral_dNT_op_dN_NT_op_dop_dN_dV
    NumIntegral_DNT_op_DN_NT_v_DN_dV
    NumIntegral_NT_dNi_dV
    NumIntegral_NT_dNi_dV_sc
    NumIntegral_NT_lhsop_N_dV
    NumIntegral_PT_lhsop_P_dV        

    ==================> RHS:
    MathOperatorRHS
    NumIntegral_DNi_rhsop_dV
    NumIntegral_dNT_op_dV                                         1
    NumIntegral_DNT_rhsop_DN_dV
    NumIntegral_DNT_v_dV
    NumIntegral_NT_mixed_op_dNi_dV
    NumIntegral_NT_op_dNi_dV
    NumIntegral_NT_op_N_dS
    NumIntegral_NT_op_N_dV                                        1
    NumIntegral_NT_op1_op2_dNi_dV
    NumIntegral_op_NT_dN_orthogonal_dV
    NumIntegral_op_NT_N_dV
    NumIntegral_op_PT_P_dV
    NumIntegral_PT_op_dS
    NumIntegral_PT_op_dV
    NumIntegral_PT_op_P_dV
    PointSource_rhsop                                             1
    */

  private:
    void Reset();
    void TestSingleVariable() {}       // TODO: reinstate method
    void TestTwoScalarVariables() {}   // TODO: reinstate method
    void TestOutputSingleVariable() {} // TODO: reinstate method
  };

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TEST_H */
