// NumIntegral_dNT_lhsop_dN_dV_Test.h
#ifndef CSMP_NUMINTEGRAL_DNT_OP_DN_DV_TEST
#define CSMP_NUMINTEGRAL_DNT_OP_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "SparseMatrix.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"

namespace csmp {

/**
 * Unit test for NumIntegral_dNT_lhsop_dN_dV — the stiffness matrix operator.
 *
 * Tests scalar, vector and tensor operand cases against analytically derived
 * reference values computed from the iso.1 mesh (4 linear triangles, 5 nodes).
 *
 * Mesh geometry (0-indexed nodes):
 *
 *   0(0,500) -------- 3(500,500)
 *   |        \  e3  / |
 *   |         \    /  |
 *   |    e2    \  / e1 |
 *   |           \/    |
 *   |           4     |
 *   |          (250,  |
 *   |           250)  |
 *   |           /\    |
 *   |          /  \   |
 *   |    e0   /    \  |
 *   |        /      \ |
 *   1(0,0) ---------- 2(500,0)
 *
 * Element connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
 * All triangles have area A = 62500.
 *
 * Shape function gradient vectors g_j = 2A * grad(N_j):
 *   e0: g0=(-250,250),  g1=(-250,-250), g4=(500,0)
 *   e1: g1=(-250,-250), g2=(250,-250),  g4=(0,500)
 *   e2: g2=(250,-250),  g3=(250,250),   g4=(-500,0)
 *   e3: g3=(250,250),   g0=(-250,250),  g4=(0,-500)
 *
 * Stiffness matrix entry: K_jk^(e) = E * (g_j . g_k) / (4A)
 * Base unit: s = 125000/(4A) = 0.5
 *
 * Global stiffness matrix (scalar E={1,2,3,4}):
 *   K_00=5s,  K_11=3s,  K_22=5s,  K_33=7s,  K_44=20s
 *   K_04=K_40=-5s, K_14=K_41=-3s, K_24=K_42=-5s, K_34=K_43=-7s
 *   All other off-diagonal entries = 0
 */
class NumIntegral_dNT_lhsop_dN_dV_Test : public Test {
  public:
    explicit NumIntegral_dNT_lhsop_dN_dV_Test( bool verbose );
    ~NumIntegral_dNT_lhsop_dN_dV_Test();

    // non-copyable, non-movable
    NumIntegral_dNT_lhsop_dN_dV_Test( const NumIntegral_dNT_lhsop_dN_dV_Test& ) = delete;
    NumIntegral_dNT_lhsop_dN_dV_Test& operator=( const NumIntegral_dNT_lhsop_dN_dV_Test& ) = delete;
    NumIntegral_dNT_lhsop_dN_dV_Test( NumIntegral_dNT_lhsop_dN_dV_Test&& ) = delete;
    NumIntegral_dNT_lhsop_dN_dV_Test& operator=( NumIntegral_dNT_lhsop_dN_dV_Test&& ) = delete;

    void run()                  override;
    void scalarTest();          ///< isotropic scalar diffusivity
    void vectorTest();          ///< diagonal anisotropic diffusivity
    void tensorTest();          ///< full tensor anisotropic diffusivity
    void symmetryTest();        ///< stiffness matrix must be symmetric
    void rowSumTest();          ///< row sums must be zero (no-flux property)

  private:
    void setElementScalar( const std::vector<double>& val,  const char* var_name );
    void setElementVector( const std::vector<std::array<double,2>>& val,
                           const char* var_name );
    void setElementTensor( const std::vector<std::array<double,4>>& val,
                           const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    // known connected pairs from mesh topology
    // used to avoid accessing non-existent sparse matrix entries
    static constexpr std::array<std::pair<size_t,size_t>,8> connected_pairs_{{
        {0U,4U}, {1U,4U}, {2U,4U}, {3U,4U},
        {4U,0U}, {4U,1U}, {4U,2U}, {4U,3U}
    }};

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_NUMINTEGRAL_DNT_OP_DN_DV_TEST
