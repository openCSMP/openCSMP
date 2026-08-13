// Integral_dNT_rhsop_dN_dV_Test.h
#ifndef CSMP_INTEGRAL_DNT_RHSOP_DN_DV_TEST
#define CSMP_INTEGRAL_DNT_RHSOP_DN_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "Integral_DNT_rhsop_DN_dV.h"
#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
 * Unit test for Integral_DNT_rhsop_DN_dV —
 * the analytical stiffness-vector RHS operator.
 *
 * Computes: RHS[j] = integral( grad(N_j)^T [D] grad(N_k) * u_k ) dV
 *                  = sum_k K_jk * u_k
 *
 * Uses exact analytical integration in global coordinates for straight-sided
 * LINEAR_TRIANGLE elements. Reference values are identical to the numerical
 * operator NumIntegral_DNT_rhsop_DN_dV since both compute the same integral
 * exactly for linear triangles.
 *
 * Mesh: iso.1, 4 LINEAR_TRIANGLE elements, 5 nodes
 * Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}, A=62500
 * s = 0.5 (stiffness base unit)
 */
class Integral_dNT_rhsop_dN_dV_Test : public Test {
  public:
    explicit Integral_dNT_rhsop_dN_dV_Test( bool verbose );
    ~Integral_dNT_rhsop_dN_dV_Test();

    Integral_dNT_rhsop_dN_dV_Test( const Integral_dNT_rhsop_dN_dV_Test& ) = delete;
    Integral_dNT_rhsop_dN_dV_Test& operator=( const Integral_dNT_rhsop_dN_dV_Test& ) = delete;
    Integral_dNT_rhsop_dN_dV_Test( Integral_dNT_rhsop_dN_dV_Test&& ) = delete;
    Integral_dNT_rhsop_dN_dV_Test& operator=( Integral_dNT_rhsop_dN_dV_Test&& ) = delete;

    void run() override;
    
    void scalarTest();          ///< RHS = K * u, non-uniform nodal values
    void uniformTest();         ///< uniform u -> RHS must be zero
    void matchesNumericalTest();///< must match NumIntegral counterpart

  private:
    void setNodeVariable(    const std::vector<double>& val, const char* var_name );
    void setElementScalar(   const std::vector<double>& val, const char* var_name );
    void calculateRHS( std::vector<double>& rhs, MathOperatorRHS<2U>& oper );

    double     tol_;
    Model<2U>* sg_{ nullptr };
    const bool verbose_;
};

} // namespace csmp
#endif // CSMP_INTEGRAL_DNT_RHSOP_DN_DV_TEST

