#ifndef CSMP_MATH_OPERATOR_LHS_H
#define CSMP_MATH_OPERATOR_LHS_H

#include "Parameter.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"
//#include "CompressedSparseRowMatrix.h"

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
    
    virtual ~MathOperatorLHS();
    

    /// Operand Functions

    std::string   Name() const;
    void          Name( const char* mathoperator_name, const char* bopname, const char* topname );
    void          Name( const char* mathoperator_name, const char* opname,  const char* bopname, const char* topname );

    virtual void  Out() const;

    const Parameter& MaterialOperand()          const;
    const Index&  MaterialOperandKey()          const;
    std::string   MaterialOperandName()         const;
    VARIABLE_TYPE MaterialOperandType()         const;
    PLACEMENT     MaterialOperandPlacement()    const;
    uint32_t      MaterialOperandDataDepth()    const;

    const Parameter& BasicOperand()             const;
    const Index&  BasicOperandKey()             const;
    std::string   BasicOperandName()            const;
    VARIABLE_TYPE BasicOperandType()            const;
    PLACEMENT     BasicOperandPlacement()       const;
    uint32_t      BasicOperandDataDepth()       const;
    size_t        BasicOperandOffset()          const;
    void          BasicOperandOffset( size_t );

    const Parameter& TestOperand()              const;
    const Index&  TestOperandKey()              const;
    std::string   TestOperandName()             const;
    VARIABLE_TYPE TestOperandType()             const;
    PLACEMENT     TestOperandPlacement()        const;
    uint32_t      TestOperandDataDepth()        const;
    size_t        TestOperandOffset()           const;
    void          TestOperandOffset( size_t );

    /// Accumulation Process Settings

    /// Get flags and properties
    bool          Add()                         const;
    bool          Subtract()                    const;
    bool          AddLater()                    const;
    bool          SubtractLater()               const;
    bool          Multiply()                    const;
    bool          LumpedFormulation()           const;
    uint32_t      ApplicationCycle()            const;
    uint32_t      ApplicationCycles()           const;
    double        MultiplyBy()                  const;
    bool          MultiplyWithTimeIncrement()   const;
    bool          DivideByTimeIncrement()       const;

    /// Set flags and properties
    void          AddAccumulate();
    void          SubtractAccumulate();
    void          AddAccumulateLater();
    void          SubtractAccumulateLater();
    void          MultiplyAccumulate();
    void          LumpedFormulation ( bool );
    void          ApplicationCycle  ( uint32_t );
    void          ApplicationCycles ( uint32_t );
    void          MultiplyBy( double factor );
    void          MultiplyWithTimeIncrement( bool multiply );
    void          DivideByTimeIncrement( bool divide );

    /// interpolation of property if isoparametric elements are used
    void          PropertyAtIntegrationPoint( const CELL<dim>&,
                                              const csmp::Index&,
                                              uint32_t ip,
                                              DenseMatrix<DM_MIN>& );
                                              
    /// getting data from the Element, Face, InterFace
    virtual void  GetOperands( const CELL<dim>&  );

    /// integration performed on Element, Face, InterFace
    virtual void  ComputeContribution( const CELL<dim>&  );

    /// writing data to the Element, Face, InterFace
    virtual void  WriteOperands( CELL<dim>& );

    /// multiply with time increment if this is desired
    virtual void  MultiplyWithTimeFactor( double dt );

    /// assigment to the left hand side global matrix (after everything was calculated )
    virtual void  AssignToGlobal( const CELL<dim>&, SparseMatrix& );

    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector
    virtual void AssignToGlobal( const CELL<dim>&, SparseMatrix&, std::vector<double>&, const std::vector<size_t>& );

    /// used by PDE_Integrator_CRM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions)
    virtual void AssignToGlobal( const CELL<dim>&, CompressedRowMatrix&, std::vector<double>&, const std::vector<size_t>& );

    virtual MathOperatorLHS<dim,CELL>* clone() const = 0;

  protected:

    /// prevent default construction
    MathOperatorLHS();

    std::string                         name_;   ///< name of operator

    Parameter                           op;      ///< material property operand
    std::pair<Parameter, size_t>        bop;     ///< basic function operand and calculation offset
    std::pair<Parameter, size_t>        top;     ///< test function operand and calculation offset

    DenseMatrix<DM_MIN>                 LHS;     ///< solution matrix to be accumulated

    std::vector<DenseMatrix<DM_MIN> >   MTRL;    ///< material property matrix(es) needed for PDE operand
    DenseMatrix<DM_MIN>                 DERIV;   ///< shape function derivative matrix
    std::vector<double>                 IPOL;    ///< shape function vector
    double                              factor_; ///< constant scaling factor

    // constraints on the accumulation procedure
    bool                                add_accumulate_;
    bool                                subtract_accumulate_;
    bool                                add_accumulate_later_;
    bool                                subtract_accumulate_later_;
    bool                                multiply_accumulate_;
    bool                                lump_matrices_;
    uint32_t                            application_cycles_;
    uint32_t                            application_cycle_;
    // time-dependent multipliers
    bool                                time_multiply_;
    bool                                time_divide_;
};



} // csmp

#endif
