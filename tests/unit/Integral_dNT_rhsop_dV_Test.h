// Integral_dNT_rhsop_dV_Test.h
#ifndef CSMP_INTEGRAL_DNT_RHSOP_DV_TEST
#define CSMP_INTEGRAL_DNT_RHSOP_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Model;
/**
 * Unit test for Integral_dNT_op_dV —
 * the analytical gradient-vector RHS operator.
 *
 * Computes: RHS[j] = integral( grad(N_j) . v ) dV
 *                  = sum_{e containing j} (g_j . v) / 2
 *
 * Uses exact analytical integration in global coordinates for straight-sided
 * LINEAR_TRIANGLE elements. Reference values are identical to the numerical
 * operator NumIntegral_dNT_op_dV since both compute the same integral
 * exactly for linear triangles.
 *
 * Mesh: iso.1, 4 LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 */
class Integral_dNT_rhsop_dV_Test final : public Test {
  public:
    explicit Integral_dNT_rhsop_dV_Test( bool verbose );
    ~Integral_dNT_rhsop_dV_Test();

    Integral_dNT_rhsop_dV_Test( const Integral_dNT_rhsop_dV_Test& ) = delete;
    Integral_dNT_rhsop_dV_Test& operator=( const Integral_dNT_rhsop_dV_Test& ) = delete;
    Integral_dNT_rhsop_dV_Test( Integral_dNT_rhsop_dV_Test&& ) = delete;
    Integral_dNT_rhsop_dV_Test& operator=( Integral_dNT_rhsop_dV_Test&& ) = delete;

    void run() override;
    
    void xDirectionTest();      ///< v=(1,0)
    void yDirectionTest();      ///< v=(0,1)
    void diagonalTest();        ///< v=(1,1)
    void conservationTest();    ///< sum of RHS must be zero
    void matchesNumericalTest();///< must match NumIntegral counterpart

  private:
    void setElementVector( const std::array<double,2>& v, const char* var_name );
    void calculateRHS( std::vector<double>& rhs, MathOperatorRHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_INTEGRAL_DNT_OP_DV_TEST
