// Integral_NT_lhsop_N_dV_Test.h
#ifndef CSMP_INTEGRAL_NT_LHSOP_N_DV_TEST
#define CSMP_INTEGRAL_NT_LHSOP_N_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "SparseMatrix.h"
#include "Integral_NT_lhsop_N_dV.h"

namespace csmp {

/**
 * Unit test for Integral_NT_lhsop_N_dV — the analytical mass matrix operator.
 *
 * Uses exact analytical integration in global coordinates for straight-sided
 * LINEAR_TRIANGLE elements.
 *
 * Mesh: iso.1, 4 LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 * A/6=10416.67, A/12=5208.33
 */
class Integral_NT_lhsop_N_dV_Test : public Test {
  public:
    explicit Integral_NT_lhsop_N_dV_Test( bool verbose );
    ~Integral_NT_lhsop_N_dV_Test();

    Integral_NT_lhsop_N_dV_Test( const Integral_NT_lhsop_N_dV_Test& ) = delete;
    Integral_NT_lhsop_N_dV_Test& operator=( const Integral_NT_lhsop_N_dV_Test& ) = delete;
    Integral_NT_lhsop_N_dV_Test( Integral_NT_lhsop_N_dV_Test&& ) = delete;
    Integral_NT_lhsop_N_dV_Test& operator=( Integral_NT_lhsop_N_dV_Test&& ) = delete;

    void run() override;
    
    void consistentTest();
    void lumpedTest();
    void rowSumTest();
    void symmetryTest();

  private:
    void setElementScalar( const std::vector<double>& val, const char* var_name );
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_INTEGRAL_NT_LHSOP_N_DV_TEST
