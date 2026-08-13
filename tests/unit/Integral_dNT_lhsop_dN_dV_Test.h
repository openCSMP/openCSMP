// Integral_dNT_lhsop_dN_dV_Test.h
#ifndef CSMP_INTEGRAL_DNT_LHSOP_DN_DV_TEST
#define CSMP_INTEGRAL_DNT_LHSOP_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "MathOperatorLHS.h"

namespace csmp {

class SparseMatrix;
template<uint32_t> class Model;

/**
 * Unit test for Integral_dNT_op_dN_dV — the analytical stiffness matrix operator.
 *
 * Uses exact analytical integration in global coordinates for straight-sided
 * LINEAR_TRIANGLE elements. Reference values are identical to the numerical
 * operator NumIntegral_dNT_lhsop_dN_dV since both compute the same integral
 * exactly for linear triangles.
 *
 * Mesh: iso.1, 4 LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 * s = 125000/(4A) = 0.5
 */
class Integral_dNT_lhsop_dN_dV_Test : public Test {
  public:
    explicit Integral_dNT_lhsop_dN_dV_Test( bool verbose );
    ~Integral_dNT_lhsop_dN_dV_Test();

    Integral_dNT_lhsop_dN_dV_Test( const Integral_dNT_lhsop_dN_dV_Test& ) = delete;
    Integral_dNT_lhsop_dN_dV_Test& operator=( const Integral_dNT_lhsop_dN_dV_Test& ) = delete;
    Integral_dNT_lhsop_dN_dV_Test( Integral_dNT_lhsop_dN_dV_Test&& ) = delete;
    Integral_dNT_lhsop_dN_dV_Test& operator=( Integral_dNT_lhsop_dN_dV_Test&& ) = delete;

    void run() override;
    void scalarTest();
    void symmetryTest();
    void rowSumTest();

  private:
    void setElementScalar( const std::vector<double>& val, const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_INTEGRAL_DNT_OP_DN_DV_TEST

