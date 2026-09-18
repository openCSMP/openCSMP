// NumIntegral_NT_rhsop_N_dV_Test.h
#ifndef CSMP_NUMINTEGRAL_NT_RHSOP_N_DV_TEST
#define CSMP_NUMINTEGRAL_NT_RHSOP_N_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
 * Unit test for NumIntegral_NT_rhsop_N_dV — the RHS mass vector operator.
 *
 * Tests the consistent and lumped formulations against analytically derived
 * reference values computed from the iso.1 mesh (4 linear triangles, 5 nodes).
 *
 * Mesh geometry (0-indexed nodes):
 *
 *   3 -------- 2
 *   |  \  e3  /|
 *   |   \    / |
 *   | e2 \  / e1
 *   |     \/   |
 *   |     4    |
 *   |     /\   |
 *   |    /  \  |
 *   |   / e0 \ |
 *   |  /      \|
 *   0 -------- 1
 *
 * Element connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
 * All triangles have area A = 62500.
 *
 * The operator computes:
 *
 *   {RHS}[j] = integral( N_j * op * N_k ) dV  (row-summed)
 *
 * which after row-sum assembly gives:
 *
 *   RHS[j] = sum_k M_jk   (consistent)
 *   RHS[j] = M_jj^lumped  (lumped)
 */
class NumIntegral_NT_rhsop_N_dV_Test : public Test {
  public:
    explicit NumIntegral_NT_rhsop_N_dV_Test( bool verbose );
    ~NumIntegral_NT_rhsop_N_dV_Test();

    // non-copyable, non-movable
    NumIntegral_NT_rhsop_N_dV_Test( const NumIntegral_NT_rhsop_N_dV_Test& ) = delete;
    NumIntegral_NT_rhsop_N_dV_Test& operator=( const NumIntegral_NT_rhsop_N_dV_Test& ) = delete;
    NumIntegral_NT_rhsop_N_dV_Test( NumIntegral_NT_rhsop_N_dV_Test&& ) = delete;
    NumIntegral_NT_rhsop_N_dV_Test& operator=( NumIntegral_NT_rhsop_N_dV_Test&& ) = delete;

    void run() override;
    void consistentTest();
    void lumpedTest();
    void lumpedEqualsConsistentRowSumTest();

  private:
    void setElementVariable( const std::vector<double>& var, const char* var_name );
    void calculateRHS( std::vector<double>& rhs, MathOperatorRHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp

#endif //

