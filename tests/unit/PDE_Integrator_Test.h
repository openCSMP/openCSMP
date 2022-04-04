#ifndef CSMP_PDE_INTEGRATOR_TEST_H
#define CSMP_PDE_INTEGRATOR_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "MathOperatorLHS.h"
#include "MathOperatorRHS.h"
#include "PropertyDatabase.h"

namespace csmp {

// LHS operator that puts fixed values in the matrix
template<uint32_t dim,class CELL=Element<dim> >
class LHS_FixedValueMatrix : public MathOperatorLHS<dim> {
  public:
    LHS_FixedValueMatrix( const PropertyDatabase<dim>& p, const char* oper, const char* basic, const char* test, double value )
      : MathOperatorLHS<dim>(p,oper,basic,test), value_for_matrix_{value} {}
    virtual ~LHS_FixedValueMatrix() {}
    virtual void GetOperands( const CELL& ) {}
    virtual void ComputeContribution( const CELL& e ) {
         MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
         MathOperatorLHS<dim>::LHS = -value_for_matrix_; // off-diagonal is negative
         MathOperatorLHS<dim>::LHS.AssignToDiagonal( value_for_matrix_ ); // diagonal is negative
      }
    virtual LHS_FixedValueMatrix<dim,CELL >* clone() const { return new LHS_FixedValueMatrix<dim,CELL >(*this); }
  private:
    double value_for_matrix_;
};

// RHS operator that puts fixed values in the matrix
template<uint32_t dim,class CELL=Element<dim> >
class RHS_FixedValueMatrix : public MathOperatorRHS<dim> {
  public:
    RHS_FixedValueMatrix( const PropertyDatabase<dim>& p, const char* basic, const char* test, double value )
     : MathOperatorRHS<dim>(p,basic,test), value_for_vector_{value} {}
    virtual ~RHS_FixedValueMatrix() {}
    virtual void GetOperands( const CELL& ) {}
    virtual void ComputeContribution( const CELL& e )
      { this->RHS.resize(e.Nodes()); for ( auto i{0U}; i<dim; i++ ) this->RHS[i] = value_for_vector_; } // off-diagonal is negative
    virtual RHS_FixedValueMatrix<dim,CELL >* clone() const { return new RHS_FixedValueMatrix<dim,CELL >(*this); }
  private:
    double value_for_vector_;
};


template<uint32_t> class Model;

/**

TODO

@todo test that the PDE operators are accumulated into the right places in the lefthand matrix
@todo test that late accumulate methods do the right things
@todo test accumulation of boundary integrals
@todo test coupling of domains by SplitBoundary integrals
@todo test accumulation of a coupled system
@todo test that the contacts between domain and boundaries are identified correctly for the accumulation of surface integrals at the boundary
@todo move all the testing that relates to matrix inversion into Solver_Test and specific subclasses
@todo test repeated use of integrator in a time-dependent problem:  matrix retention vs. reconstruction, test change of DOF from step to step

*/
class PDE_Integrator_Test : public Test {
  public:
    explicit PDE_Integrator_Test( Model<2U>& model );
    explicit PDE_Integrator_Test( Model<3U>& model );
    ~PDE_Integrator_Test();
    
    void run();

  private:
    const static bool verbose_ = true;

    /*
    =================> LHS:
    NumIntegral_BT_D_B_dV
    NumIntegral_BT_D_op_dV
    NumIntegral_dNT_dN_dV
    NumIntegral_dNT_mixed_op_dN_dV
    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV
    NumIntegral_dNT_op_dN_dV                                      1
    NumIntegral_dNT_op_dN_dV_NT_v_dN_dV
    NumIntegral_dNT_op_dN_NT_op_dop_dN_dV
    NumIntegral_DNT_op_DN_NT_v_DN_dV
    NumIntegral_NT_dNi_dV
    NumIntegral_NT_dNi_dV_sc
    NumIntegral_NT_lhsop_N_dV
    NumIntegral_PT_lhsop_P_dV        

    ==================> RHS:
    MathOperatorRHS
    NumIntegral_DNi_rhsop_dV
    NumIntegral_dNT_op_dV                                         1
    NumIntegral_DNT_rhsop_DN_dV
    NumIntegral_DNT_v_dV
    NumIntegral_NT_mixed_op_dNi_dV
    NumIntegral_NT_op_dNi_dV
    NumIntegral_NT_op_N_dS
    NumIntegral_NT_op_N_dV                                        1
    NumIntegral_NT_op1_op2_dNi_dV
    NumIntegral_op_NT_dN_orthogonal_dV
    NumIntegral_op_NT_N_dV
    NumIntegral_op_PT_P_dV
    NumIntegral_PT_op_dS
    NumIntegral_PT_op_dV
    NumIntegral_PT_op_P_dV
    PointSource_rhsop                                             1
    */

  private:
    Model<2U>& model_; // 2D is fully sufficient
    
    /// checks that element matrices and vectors end up in the right place
    void TestAssembly();
    
    // methods from Luat
    void Reset();
    void TestSingleVariable() {}       // TODO: reinstate method
    void TestTwoScalarVariables() {}   // TODO: reinstate method
    void TestOutputSingleVariable() {} // TODO: reinstate method
  };

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TEST_H */
