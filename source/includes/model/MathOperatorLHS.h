#ifndef CSMP_MATH_OPERATOR_LHS_H
#define CSMP_MATH_OPERATOR_LHS_H

#include "Parameter.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"


namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class PropertyDatabase;

/**
@brief Base class for FE or FVM integrals for the lefthandside (matrix).

@attention In a linear system A x = b, every single row of that matrix A represents one equation.
Each equation corresponds to one specific Test Function basis.

Therefore, the rows are determined by the Test Operand.

The columns represent the degrees of freedom (DOF) of the unknown variable 'u''.
These correspond to the Basic (Trial) Functions.
Therefore, the columns are determined by the Basic Operand.
 
@author S.K. Matthai
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

To be able to write PDE equations in finite element form, using the 
PDE_Integrator class. 

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class MathOperatorLHS {
  public:

    MathOperatorLHS( const PropertyDatabase<dim>&,
                     const char* basic, 
                     const char* test );

    MathOperatorLHS( const PropertyDatabase<dim>&,
                     const char* oper, 
                     const char* basic, 
                     const char* test );
    
    MathOperatorLHS( const MathOperatorLHS& mo );
    
    MathOperatorLHS& operator=( const MathOperatorLHS& );
    
    virtual ~MathOperatorLHS() = default;
    

    /// Operand Functions

    std::string   Name() const noexcept;
    void          Name( const char* mathoperator_name, const char* bopname, const char* topname ) noexcept;
    void          Name( const char* mathoperator_name, const char* opname,  const char* bopname, const char* topname ) noexcept;

    virtual void  Out() const;

    const Parameter& MaterialOperand()          const noexcept;
    const Index&  MaterialOperandKey()          const noexcept;
    std::string   MaterialOperandName()         const noexcept;
    VARIABLE_TYPE MaterialOperandType()         const noexcept;
    PLACEMENT     MaterialOperandPlacement()    const noexcept;
    uint32_t      MaterialOperandDataDepth()    const noexcept;

    const Parameter& BasicOperand()             const noexcept;
    const Index&  BasicOperandKey()             const noexcept;
    std::string   BasicOperandName()            const noexcept;
    VARIABLE_TYPE BasicOperandType()            const noexcept;
    PLACEMENT     BasicOperandPlacement()       const noexcept;
    uint32_t      BasicOperandDataDepth()       const noexcept;
    size_t        BasicOperandOffset()          const noexcept;
    void          BasicOperandOffset( size_t ) noexcept;

    const Parameter& TestOperand()              const noexcept;
    const Index&  TestOperandKey()              const noexcept;
    std::string   TestOperandName()             const noexcept;
    VARIABLE_TYPE TestOperandType()             const noexcept;
    PLACEMENT     TestOperandPlacement()        const noexcept;
    uint32_t      TestOperandDataDepth()        const noexcept;
    size_t        TestOperandOffset()           const noexcept;
    void          TestOperandOffset( size_t ) noexcept;

    /// Accumulation Process Settings

    /// Get flags and properties
    bool          Add()                         const noexcept;
    bool          Subtract()                    const noexcept;
    bool          AddLater()                    const noexcept;
    bool          SubtractLater()               const noexcept;
    bool          Multiply()                    const noexcept;
    bool          LumpedFormulation()           const noexcept;
    uint32_t      ApplicationCycle()            const noexcept;
    uint32_t      ApplicationCycles()           const noexcept;
    /// returns 'factor_' used for time multiplication
    [[nodiscard]] double MultiplyBy()           const noexcept;
    bool          MultiplyWithTimeIncrement()   const noexcept;
    bool          DivideByTimeIncrement()       const noexcept;

    /// Set flags and properties
    void          AddAccumulate() noexcept;
    void          SubtractAccumulate() noexcept;
    void          AddAccumulateLater() noexcept;
    void          SubtractAccumulateLater() noexcept;
    void          MultiplyAccumulate() noexcept;
    void          LumpedFormulation ( bool ) noexcept;
    void          ApplicationCycle  ( uint32_t ) noexcept;
    void          ApplicationCycles ( uint32_t ) noexcept;
    void          MultiplyBy( double factor ) noexcept;
    void          MultiplyWithTimeIncrement( bool multiply ) noexcept;
    void          DivideByTimeIncrement( bool divide ) noexcept;

    /// interpolation of property if isoparametric elements are used
    void          PropertyAtIntegrationPoint( const CELL<dim>&,
                                              const csmp::Index&,
                                              uint32_t ip,
                                              DenseMatrix<dim>& );
                                              
    /// getting data from the Element, Face, InterFace
    virtual void  GetOperands( const CELL<dim>& );

    /// integration performed on Element, Face, InterFace
    virtual void  ComputeContribution( const CELL<dim>& );

    /// writing data to the Element, Face, InterFace
    virtual void  WriteOperands( CELL<dim>& );

    /// multiply with time increment if this is desired
    void MultiplyWithTimeFactor( double dt ) noexcept { LHS *= dt; }

    /// assigment to the left hand side global matrix (in case where no full elimination of Dirichlet constraints is carried out )
    void AssignToGlobal( const CELL<dim>&, SparseMatrix& );

    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector
    void AssignToGlobal( const CELL<dim>&, SparseMatrix&, std::vector<double>&, const std::vector<size_t>& );

    /// used by PDE_Integrator_CRM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions)
    void AssignToGlobal( const CELL<dim>&, CompressedRowMatrix&, std::vector<double>&, const std::vector<size_t>& );

    virtual MathOperatorLHS<dim,CELL>* clone() const = 0;

  protected:

    /// prevent default construction, yet set defaults when used
    MathOperatorLHS();

    std::string                      name_;   ///< name of operator

    Parameter                        op;      ///< material property operand
    std::pair<Parameter, size_t>     bop;     ///< basic function operand and calculation offset
    std::pair<Parameter, size_t>     top;     ///< test function operand and calculation offset

    DenseMatrix<DM_MIN>              LHS;     ///< solution matrix to be accumulated
    DenseMatrix<DM_MIN>              DERIV;   ///< shape function derivative matrix TODO: not thread safe
    std::vector<double>              IPOL;    ///< shape function vector

    std::vector<DenseMatrix<dim> >   MTRL;    ///< material property matrix(es) needed for PDE operand
    double                           factor_; ///< constant scaling factor

    // constraints on the accumulation procedure
    bool                             add_accumulate_;
    bool                             subtract_accumulate_;
    bool                             add_accumulate_later_;
    bool                             subtract_accumulate_later_;
    bool                             multiply_accumulate_;
    bool                             lump_matrices_;
    uint32_t                         application_cycles_;
    uint32_t                         application_cycle_;
    // time-dependent multipliers
    bool                             time_multiply_;
    bool                             time_divide_;

  private:

    /**
     * @brief Reads prescribed Dirichlet values of the basic (column) operand
     *        into vals, sized to IDB.size() = LHS.Cols().
     *
     * Used by both AssignToGlobal overloads to populate the pivot vector
     * when basic DOFs are Dirichlet-constrained.
     *
     * @param e        The element being processed
     * @param n_nodes  Number of nodes on the element
     * @param vals     Output vector, must be pre-sized to IDB.size()
     */
    void ReadBasicOperandValues( const CELL<dim>& e,
                                 uint32_t         n_nodes,
                                 vector<double>&  vals ) const
    {
        switch ( this->BasicOperandType() )
        {
            case SCALAR:
                for ( uint32_t nIdx{0U}; nIdx < n_nodes; ++nIdx )
                    vals[nIdx] = e.N(nIdx)->Read( this->BasicOperandKey() );
                break;

            case VECTOR:
            {
                VectorVariable<dim> var;
                for ( uint32_t nIdx{0U}; nIdx < n_nodes; ++nIdx ) {
                    e.N(nIdx)->Read( this->BasicOperandKey(), var );
                    for ( uint32_t i{0U}; i < dim; ++i )
                        vals[nIdx*dim + i] = var.Component(i);
                }
                break;
            }

            case TENSOR:
            {
                constexpr uint32_t dim2 = dim * dim;
                TensorVariable<dim> var;
                for ( uint32_t nIdx{0U}; nIdx < n_nodes; ++nIdx ) {
                    e.N(nIdx)->Read( this->BasicOperandKey(), var );
                    for ( uint32_t i{0U}; i < dim; ++i )
                        for ( uint32_t j{0U}; j < dim; ++j )
                            vals[nIdx*dim2 + i*dim + j] = var.Component(i*dim + j);
                }
                break;
            }

            case ARRAY:
            {
                const auto depth( this->BasicOperandDataDepth() );
                ArrayVariable var;
                for ( uint32_t nIdx{0U}; nIdx < n_nodes; ++nIdx ) {
                    e.N(nIdx)->Read( this->BasicOperandKey(), var );
                    for ( uint32_t i{0U}; i < depth; ++i )
                        vals[nIdx*depth + i] = var.Component(i);
                }
                break;
            }

            case FLAGGEDARRAY:
            {
                const auto depth( this->BasicOperandDataDepth() );
                FlaggedArrayVariable var;
                for ( uint32_t nIdx{0U}; nIdx < n_nodes; ++nIdx ) {
                    e.N(nIdx)->Read( this->BasicOperandKey(), var );
                    for ( uint32_t i{0U}; i < depth; ++i )
                        vals[nIdx*depth + i] = var.Component(i);
                }
                break;
            }

            default:
              throw logic_error( "ERROR, MathOperatorLHS::ReadBasicOperandValues: Undefined basic operand variable type" );
        }
    }
};


} // csmp

#endif
