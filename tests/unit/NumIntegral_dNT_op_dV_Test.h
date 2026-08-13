// NumIntegral_dNT_op_dV_Test.h
#ifndef CSMP_NUMINTEGRAL_DNT_OP_DV_TEST
#define CSMP_NUMINTEGRAL_DNT_OP_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "NumIntegral_dNT_op_dV.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
 * Unit test for NumIntegral_dNT_op_dV — the gradient-vector RHS operator.
 *
 * Computes: RHS[j] = integral( grad(N_j) . v ) dV
 *                  = sum_e A * grad(N_j) . v
 *
 * where v is an element-placed vector operand (e.g. "gravity term").
 *
 * For linear triangles grad(N_j) is constant, so:
 *   RHS[j] = sum_{e containing j} (g_j . v) / 2
 *
 * where g_j = 2A * grad(N_j) are the gradient vectors:
 *   e0={0,1,4}: g0=(-250,250),  g1=(-250,-250), g4=(500,0)
 *   e1={1,2,4}: g1=(-250,-250), g2=(250,-250),  g4=(0,500)
 *   e2={2,3,4}: g2=(250,-250),  g3=(250,250),   g4=(-500,0)
 *   e3={3,0,4}: g3=(250,250),   g0=(-250,250),  g4=(0,-500)
 *
 * With v=(1,0):  RHS = {-250, -250,  250,  250, 0}
 * With v=(0,1):  RHS = { 250, -250, -250,  250, 0}
 * With v=(1,1):  RHS = {   0, -500,    0,  500, 0}
 *
 * Conservation: sum_j RHS[j] = 0 for any uniform v (closed domain)
 */
class NumIntegral_dNT_op_dV_Test : public Test {
  public:
    explicit NumIntegral_dNT_op_dV_Test( bool verbose );
    ~NumIntegral_dNT_op_dV_Test();

    // non-copyable, non-movable
    NumIntegral_dNT_op_dV_Test( const NumIntegral_dNT_op_dV_Test& ) = delete;
    NumIntegral_dNT_op_dV_Test& operator=( const NumIntegral_dNT_op_dV_Test& ) = delete;
    NumIntegral_dNT_op_dV_Test( NumIntegral_dNT_op_dV_Test&& ) = delete;
    NumIntegral_dNT_op_dV_Test& operator=( NumIntegral_dNT_op_dV_Test&& ) = delete;

    void run()                override;
    void xDirectionTest();    ///< v=(1,0) — isolates x-gradient components
    void yDirectionTest();    ///< v=(0,1) — isolates y-gradient components
    void diagonalTest();      ///< v=(1,1) — combined x+y
    void conservationTest();  ///< sum of RHS entries must be zero for any v

  private:
    void setElementVector( const std::array<double,2>& v, const char* var_name );
    void calculateRHS( std::vector<double>& rhs, MathOperatorRHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_NUMINTEGRAL_DNT_OP_DV_TEST

