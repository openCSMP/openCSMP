#ifndef CSMP_PDE_INTEGRATOR_TEST_H
#define CSMP_PDE_INTEGRATOR_TEST_H

#include "Test.h"
#include "Attorney.h"
#include "PDE_Integrator.h"
#include "Box.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;


/** Attorney design pattern gives access to protected / private member variables and methods of PDE_Integrator

template<uint32_t dim>
class PDE_Integrator_Attorney : public Attorney<class PDE_Integrator<dim>> {
public:
    using Attorney<PDE_Integrator<dim>>::Attorney; // Inherit Attorney constructor
    using PDE_Integrator<dim>::Accumulate;
    using PDE_Integrator<dim>::EstablishMatrixSetup;
    using PDE_Integrator<dim>::lhs_operators_;
    using PDE_Integrator<dim>::rhs_operators_;
    using PDE_Integrator<dim>::G_;
    using PDE_Integrator<dim>::rh_;
};
*/

template<uint32_t dim,
         template<uint32_t> class CELLTYPE = Element, 
         class MATRIXTYPE = CompressedRowMatrix>
class PDE_Integrator_Attorney : public Attorney<PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>> {
    
    // Convenience alias to avoid repeating long template signatures
    using TargetIntegrator = PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>;

public:
    using Attorney<TargetIntegrator>::Attorney; // Inherit Attorney constructor

    // Expose member functions
    using TargetIntegrator::Accumulate;
    using TargetIntegrator::EstablishMatrixSetup;
    using TargetIntegrator::EliminateEssentialConditions;
    using TargetIntegrator::AssignEssentialConditions;
    using TargetIntegrator::AssignInitialConditions;

    // Expose member variables (operators and matrices)
    using TargetIntegrator::lhs_operators_;
    using TargetIntegrator::rhs_operators_;
    using TargetIntegrator::G_;
    using TargetIntegrator::rh_;
};


/** LHS operator that puts fixed values in the matrix */
template<uint32_t dim,class CELL=Element<dim> >
class LHS_FixedValueMatrix : public MathOperatorLHS<dim> {
  public:
    LHS_FixedValueMatrix( const PropertyDatabase<dim>& p, const char* oper, const char* basic, const char* test, double value )
      : MathOperatorLHS<dim>(p,oper,basic,test), value_label_for_matrix_entry_{value} {
          MathOperatorLHS<dim>::Name( "LHS_FixedValueMatrix", oper, basic, test );
          MathOperatorLHS<dim>::MultiplyBy(value); // sets factor_ for retrieval by MultiplyBy
    }
    virtual ~LHS_FixedValueMatrix() {}
    // use method of base class: virtual void GetOperands( const CELL& );
    virtual void ComputeContribution( const CELL& );
    virtual LHS_FixedValueMatrix<dim,CELL >* clone() const { return new LHS_FixedValueMatrix<dim,CELL >(*this); }
  private:
    double value_label_for_matrix_entry_;
};



/** RHS operator that puts fixed values in the matrix */
template<uint32_t dim,class CELL=Element<dim> >
class RHS_FixedValueMatrix : public MathOperatorRHS<dim> {
  public:
    RHS_FixedValueMatrix( const PropertyDatabase<dim>& p, const char* oper, const char* test, double value )
     : MathOperatorRHS<dim>(p,oper,test), value_label_for_vector_entry_{value} {
          MathOperatorRHS<dim>::Name( "RHS_FixedValueMatrix", test );
          MathOperatorRHS<dim>::MultiplyBy(value); // for retrieval
    }
    virtual ~RHS_FixedValueMatrix() {}
    virtual void ComputeContribution( const CELL& );
    virtual RHS_FixedValueMatrix<dim,CELL >* clone() const { return new RHS_FixedValueMatrix<dim,CELL >(*this); }
  private:
    double value_label_for_vector_entry_;
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
    PDE_Integrator_Test();
    explicit PDE_Integrator_Test( Model<2U>& model );
    ~PDE_Integrator_Test();
    
    void run();

  private:
    const static bool verbose_ = false;

  private:
    Model<2U>* model_{nullptr}; // any dimension is fully sufficient
    bool delete_model_{false}; // in case an existing model was passed in
    
    /// checks that element matrices and vectors end up in the right place
    void TestAssembly();
    void TestAssemblySingleScalarNoDirichlet( bool debug );
    void TestAssemblySingleScalarDirichlet( bool debug );
    void TestAssemblySingleVectorNoDirichlet( bool debug );
    void TestAssemblySingleVectorDirichlet( bool debug );
    
    // SKM test development
    void TestAssemblyTwoScalarVariablesNoDirichlet( bool debug );
    void TestAssemblyTwoScalarVariablesDirichlet( bool debug );
    
    void TestAssemblyScalarAndVectorVariableNoDirichlet( bool debug );
    void TestAssemblyScalarAndVectorVariableDirichlet( bool debug );
    
    /// Biot-type poromechanics (u,p) - same TINY 2D model with 2 triangles and 1 quadrilateral
    void TestBlockStructuredCoupledFEM_AssemblyWithDirichletElimination( bool debug );
    
    /// Test dataset for model TINY
    void GenerateTinyPoromechanicsTestSystem( CompressedRowMatrix& A_reduced, std::vector<double>& b_reduced );
    void PrintTinyTestSystem( const CompressedRowMatrix& A_reduced, const std::vector<double>& b_reduced );
    
    // remaining methods from Luat
    void ResetModelProperties();
    void TestSingleVariable() {}       // TODO: reinstate method
    void TestTwoScalarVariables() {}   // TODO: reinstate method
    void TestOutputSingleVariable() {} // TODO: reinstate method
    
    /// testing accumulation inside of the PDE_Integrator
    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestMatrix(const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                    const MeshManager<dim>& mesh,
                    const std::set<BOX_BOUNDARY>& dirich,
                    const VARIABLE_TYPE variable_type);
    
    /// for fixed 'val' material operators
    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestRHSVector(const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                       const MeshManager<dim>& mesh,
                       const double val,
                       const std::set<BOX_BOUNDARY>& dirich );

    /// for fixed material operators read from model using 'material operand'
    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestRHSVector(const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                       const MeshManager<dim>& mesh );

    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestMultiVariableMatrix( const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                                  const MeshManager<dim>& mesh,
                                  const std::vector<size_t>& DOF_indexes );

    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestMultiVariableRHSVector( const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                                     const std::map<Parameter, size_t>& test_operands,
                                     const std::vector<size_t>& DOF_indexes,
                                     const ModelSubDomain<dim, CELLTYPE>& );
    template<uint32_t dim,
             template<uint32_t> class CELLTYPE = Element,
             class MATRIXTYPE = CompressedRowMatrix>
    void TestPivotVector( const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                          const MeshManager<dim>& mesh,
                          const double mat_val,
                          const double dirich_val,
                          const std::set<BOX_BOUNDARY>& dirich,
                          const VARIABLE_TYPE variable_typev );
                       
    /// Test whether the matrix is symmetric and likely to be positive definite (SPD)
    void TestMatrixIsSPD( const CompressedRowMatrix& );
    void TestMatrixIsSPD( const SparseMatrix& );

 };
     
 
} // end namespace csmp


    /* INTEGRALS TO TEST SEPARATELLY
    
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
    PointSource_rhsop
                                               1
    */

#endif /* CSMP_PDE_INTEGRATOR_TEST_H */
