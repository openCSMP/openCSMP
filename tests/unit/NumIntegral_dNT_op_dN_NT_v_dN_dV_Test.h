// NumIntegral_dNT_op_dN_NT_v_dN_dV_Test.h
#ifndef CSMP_NUMINTEGRAL_DNT_OP_DN_NT_V_DN_DV_TEST
#define CSMP_NUMINTEGRAL_DNT_OP_DN_NT_V_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "SparseMatrix.h"
#include "NumIntegral_DNT_op_DN_NT_v_DN_dV.h"
#include "MathOperatorLHS.h"

namespace csmp {

/**
 * Unit test for NumIntegral_DNT_op_DN_NT_v_DN_dV —
 * the advection-dispersion matrix operator.
 *
 * Assembles:
 *   A_jk = integral( (grad N_j)^T [D] grad N_k + N_j v . grad N_k ) dV
 *        = K_jk + V_jk
 *
 * where K_jk is the dispersion (stiffness) matrix and V_jk is the
 * advection matrix.
 *
 * Tests:
 *   - dispersionOnlyTest:  v=0, result must equal pure stiffness matrix
 *   - advectionOnlyTest:   D=0, result must equal pure advection matrix
 *   - combinedTest:        D=1, v=(1,0), full advection-dispersion matrix
 *   - nonSymmetryTest:     A_jk != A_kj for non-zero v
 *   - zeroVelocityTest:    v=0 -> A_jk = K_jk (symmetric)
 *
 * Mesh: iso.1, 4 ISOPARAMETRIC_LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 *
 * Key constants:
 *   s = 125000/(4A) = 0.5       (stiffness base unit)
 *   q = 250/6 = 41.667          (advection base unit)
 */
class NumIntegral_dNT_op_dN_NT_v_dN_dV_Test : public Test {
  public:
    explicit NumIntegral_dNT_op_dN_NT_v_dN_dV_Test( bool verbose );
    ~NumIntegral_dNT_op_dN_NT_v_dN_dV_Test();

    // non-copyable, non-movable
    NumIntegral_dNT_op_dN_NT_v_dN_dV_Test( const NumIntegral_dNT_op_dN_NT_v_dN_dV_Test& ) = delete;
    NumIntegral_dNT_op_dN_NT_v_dN_dV_Test& operator=( const NumIntegral_dNT_op_dN_NT_v_dN_dV_Test& ) = delete;
    NumIntegral_dNT_op_dN_NT_v_dN_dV_Test( NumIntegral_dNT_op_dN_NT_v_dN_dV_Test&& ) = delete;
    NumIntegral_dNT_op_dN_NT_v_dN_dV_Test& operator=( NumIntegral_dNT_op_dN_NT_v_dN_dV_Test&& ) = delete;

    void run()                override;
    void dispersionOnlyTest();  ///< v=0 -> pure stiffness matrix
    void advectionOnlyTest();   ///< D=0 -> pure advection matrix
    void combinedTest();        ///< D=1, v=(1,0) -> full matrix
    void nonSymmetryTest();     ///< A_jk != A_kj for non-zero v
    void zeroVelocityTest();    ///< v=0 -> symmetric result

  private:
    void setElementScalar( const std::vector<double>& val,  const char* var_name );
    void setElementVector( const std::array<double,2>& v,   const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_NUMINTEGRAL_DNT_OP_DN_NT_V_DN_DV_TEST

