#ifndef Integral_var_NT_rhsop_N_dV_Test_h
#define Integral_var_NT_rhsop_N_dV_Test_h

#include "Test.h"
#include "Model.h"
#include "VSet.h"
#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "PDE_Integrator.h"
#include "Integral_var_NT_lhsop_N_dV.h"
#include "Integral_var_NT_rhsop_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"
#include "Concatenate.h"

namespace csmp {

  class Integral_var_NT_rhsop_N_dV_Test : public Test {
  public:
    Integral_var_NT_rhsop_N_dV_Test( bool verbose );
    ~Integral_var_NT_rhsop_N_dV_Test();
    void run();
    void valueTest();
    void compareConsistentTest();
    void compareLumpedTest();

  private:
    void compareTest(bool lumped);

    void setNodeVariable( std::vector<double64>& var, const char* var_name);
    void showNodeVariable( const char* var_name);
    void setElementVariable( std::vector<double64>& var, const char* var_name);
    void calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper);
    void calculateGlobalRHS( std::vector<double64>& rhs, MathOperatorRHS<2U>& oper);

    Model<2U>*    sg_;
    const double  tol_;
    const bool    verbose_;
  };

} // end namespace csmp

#endif
