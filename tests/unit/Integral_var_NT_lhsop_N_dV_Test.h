#ifndef CSMP_INTEGRAL_VAR_NT_LHSOP_N_DV_TEST
#define CSMP_INTEGRAL_VAR_NT_LHSOP_N_DV_TEST

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"

#include "Integral_var_NT_lhsop_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"

namespace csmp {

template<uint32_t> class Model;

class Integral_var_NT_lhsop_N_dV_Test : public Test {
  public:
    explicit Integral_var_NT_lhsop_N_dV_Test( bool verbose );
    ~Integral_var_NT_lhsop_N_dV_Test();
    // rule of five — non-copyable, non-movable (owns raw Model pointer)
    Integral_var_NT_lhsop_N_dV_Test( const Integral_var_NT_lhsop_N_dV_Test& ) = delete;
    Integral_var_NT_lhsop_N_dV_Test& operator=( const Integral_var_NT_lhsop_N_dV_Test& ) = delete;
    Integral_var_NT_lhsop_N_dV_Test( Integral_var_NT_lhsop_N_dV_Test&& ) = delete;
    Integral_var_NT_lhsop_N_dV_Test& operator=( Integral_var_NT_lhsop_N_dV_Test&& ) = delete;

    void run() override;
    void valueTest();
    void compareConsistentTest();
    void compareLumpedTest();
    void lumpedTest();
    void rowSumTest();

  private:
    void compareTest(bool lumped);

    void setNodeVariable( const std::vector<double>& var, const char* var_name);
    void showNodeVariable(const char* var_name) const;
    void setElementVariable( const std::vector<double>& var, const char* var_name);
    void calculateGlobalMatrix(SparseMatrix& sm, MathOperatorLHS<2U>& oper);

    static constexpr double tol_ = 1.0e-5;

    Model<2U>* sg_ = nullptr;
    const bool verbose_;
};

} // end namespace csmp
#endif
