// Integral_dNT_lhsop_dN_NT_v_dN_dV_Test.h
#ifndef CSMP_INTEGRAL_DNT_LHSOP_DN_NT_V_DN_DV_TEST
#define CSMP_INTEGRAL_DNT_LHSOP_DN_NT_V_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "SparseMatrix.h"
#include "Integral_dNT_op_dN_NT_v_dN_dV.h"
#include "MathOperatorLHS.h"

namespace csmp {

/**
 * Unit test for Integral_DNT_op_DN_NT_v_DN_dV —
 * the analytical advection-dispersion matrix operator.
 *
 * Uses exact analytical integration in global coordinates for straight-sided
 * LINEAR_TRIANGLE elements. Reference values are identical to the numerical
 * operator NumIntegral_DNT_op_DN_NT_v_DN_dV since both compute the same
 * integral exactly for linear triangles.
 *
 * Assembles:
 *   A_jk = integral( (grad N_j)^T [D] grad N_k + N_j v . grad N_k ) dV
 *        = K_jk + V_jk
 *
 * Mesh: iso.1, 4 LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 * s = 0.5 (stiffness base unit), q = 250/6 (advection base unit)
 */
class Integral_dNT_lhsop_dN_NT_v_dN_dV_Test : public Test {
  public:
    explicit Integral_dNT_lhsop_dN_NT_v_dN_dV_Test( bool verbose );
    ~Integral_dNT_lhsop_dN_NT_v_dN_dV_Test();

    Integral_dNT_lhsop_dN_NT_v_dN_dV_Test( const Integral_dNT_lhsop_dN_NT_v_dN_dV_Test& ) = delete;
    Integral_dNT_lhsop_dN_NT_v_dN_dV_Test& operator=( const Integral_dNT_lhsop_dN_NT_v_dN_dV_Test& ) = delete;
    Integral_dNT_lhsop_dN_NT_v_dN_dV_Test( Integral_dNT_lhsop_dN_NT_v_dN_dV_Test&& ) = delete;
    Integral_dNT_lhsop_dN_NT_v_dN_dV_Test& operator=( Integral_dNT_lhsop_dN_NT_v_dN_dV_Test&& ) = delete;

    void run() override;
    void dispersionOnlyTest();  ///< v=0 -> pure stiffness matrix
    void advectionOnlyTest();   ///< D=0 -> pure advection matrix
    void combinedTest();        ///< D=1, v=(1,0) -> full matrix
    void nonSymmetryTest();     ///< A_jk != A_kj for non-zero v
    void zeroVelocityTest();    ///< v=0 -> symmetric result
    void matchesNumericalTest();///< must match NumIntegral counterpart exactly

  private:
    void setElementScalar( const std::vector<double>& val,  const char* var_name );
    void setElementVector( const std::array<double,2>& v,   const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_TEST
