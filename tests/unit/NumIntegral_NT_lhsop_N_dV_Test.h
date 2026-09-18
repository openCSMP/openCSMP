#ifndef CSMP_NUMINTEGRAL_NT_LHSOP_N_DV_TEST
#define CSMP_NUMINTEGRAL_NT_LHSOP_N_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "SparseMatrix.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "MathOperatorLHS.h"

namespace csmp {

/**
 * Unit test for NumIntegral_NT_lhsop_N_dV — the capacitance matrix operator.
 *
 * Tests the consistent and lumped formulations against analytically derived
 * reference values computed from the iso.1 mesh (4 linear triangles, 5 nodes).
 *
 * Mesh geometry (0-indexed nodes):
 *
 *   3 -------- 2
 *   |  \  e3  /|
 *   |   \      / |
 *   | e2 \  / e1
 *   |     \/   |
 *   |     4    |
 *   |     /\   |
 *   |    /  \  |
 *   |   / e0  \ |
 *   |  /         \|
 *   0 -------- 1
 *
 * Element connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
 * All triangles have area A = 62500.
 */
class NumIntegral_NT_lhsop_N_dV_Test : public Test {
  public:
    explicit NumIntegral_NT_lhsop_N_dV_Test( bool verbose );
    ~NumIntegral_NT_lhsop_N_dV_Test();

    // non-copyable, non-movable
    NumIntegral_NT_lhsop_N_dV_Test( const NumIntegral_NT_lhsop_N_dV_Test& ) = delete;
    NumIntegral_NT_lhsop_N_dV_Test& operator=( const NumIntegral_NT_lhsop_N_dV_Test& ) = delete;
    NumIntegral_NT_lhsop_N_dV_Test( NumIntegral_NT_lhsop_N_dV_Test&& ) = delete;
    NumIntegral_NT_lhsop_N_dV_Test& operator=( NumIntegral_NT_lhsop_N_dV_Test&& ) = delete;

    void run()              override;
    void consistentTest();
    void lumpedTest();
    void rowSumTest();
    void symmetryTest();

  private:
    void setNodeVariable(    const std::vector<double>& var, const char* var_name );
    void setElementVariable( const std::vector<double>& var, const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    static constexpr double tol_ = 1.0e-4;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_NUMINTEGRAL_NT_LHSOP_N_DV_TEST
