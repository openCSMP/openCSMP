// NumIntegral_dNT_rhsop_dN_dV_Test.h
#ifndef CSMP_NUMINTEGRAL_DNT_RHSOP_DN_DV_TEST
#define CSMP_NUMINTEGRAL_DNT_RHSOP_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
 * Unit test for NumIntegral_DNT_rhsop_DN_dV — the RHS stiffness-vector operator.
 *
 * Computes: RHS[j] = integral( grad(N_j) . [D] . grad(N_k) * u_k ) dV
 *                  = sum_k K_jk * u_k
 *
 * where K_jk is the stiffness matrix and u_k are nodal values of the test variable.
 *
 * Tests:
 *   - scalarTest:      scalar diffusivity, non-uniform nodal values
 *   - multiplierTest:  scalar diffusivity with integral multiplier
 *   - uniformTest:     uniform nodal values -> RHS must be zero (row sum property)
 *   - ignoreTest:      IgnoreOperand(true) -> RHS must be zero
 *
 * Mesh: iso.1, 4 ISOPARAMETRIC_LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 * Stiffness matrix (scalar E={1,2,3,4}, s=0.5):
 *   K_00=5s,  K_11=3s,  K_22=5s,  K_33=7s,  K_44=20s
 *   K_04=K_40=-5s, K_14=K_41=-3s, K_24=K_42=-5s, K_34=K_43=-7s
 *   All other off-diagonal entries = 0
 */
class NumIntegral_dNT_rhsop_dN_dV_Test : public Test {
  public:
    explicit NumIntegral_dNT_rhsop_dN_dV_Test( bool verbose );
    ~NumIntegral_dNT_rhsop_dN_dV_Test();

    // non-copyable, non-movable
    NumIntegral_dNT_rhsop_dN_dV_Test( const NumIntegral_dNT_rhsop_dN_dV_Test& ) = delete;
    NumIntegral_dNT_rhsop_dN_dV_Test& operator=( const NumIntegral_dNT_rhsop_dN_dV_Test& ) = delete;
    NumIntegral_dNT_rhsop_dN_dV_Test( NumIntegral_dNT_rhsop_dN_dV_Test&& ) = delete;
    NumIntegral_dNT_rhsop_dN_dV_Test& operator=( NumIntegral_dNT_rhsop_dN_dV_Test&& ) = delete;

    void run() override;
    
    void scalarTest();    ///< scalar diffusivity, non-uniform nodal values
    void multiplierTest();///< scalar diffusivity with integral multiplier
    void uniformTest();   ///< uniform nodal values -> RHS must be zero

  private:
    void setNodeVariable(    const std::vector<double>& val, const char* var_name );
    void setElementScalar(   const std::vector<double>& val, const char* var_name );
    void calculateRHS( std::vector<double>& rhs, MathOperatorRHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_NUMINTEGRAL_DNT_RHSOP_DN_DV_TEST

