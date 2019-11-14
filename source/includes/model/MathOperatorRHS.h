#ifndef CSMP_MATH_OPERATOR_RHS_H
#define CSMP_MATH_OPERATOR_RHS_H

#include "Parameter.h"

#include "DenseMatrix.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class PropertyDatabase;

/**
@brief Base class for FE or FV integrals accumulated into righthand vector.

@author S.K. Matthai
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

To be able to write PDE equations in finite element form, using the
CSMP Algorithm class.

RULES for using righthand Mathoperators

1. supplied data:

    - Property vectors contain data from: nodes, integration points, or elements
 
    - node and integration point data are supplied together in a numbering as
      specified for the particular element type.
 
    - if a node or a integration point has multiple degrees of freedom, then
      the node property vector<double64> will be of dimension (nds+ips) * dof
      in this case, the entries will be 'dof' per node or constraint point in the
      order as given above.

*/
template<size_t dim>
class MathOperatorRHS {

  protected:

    // prevent default construction
    MathOperatorRHS();

  public:

    MathOperatorRHS( const MathOperatorRHS& );

    MathOperatorRHS( const PropertyDatabase<dim>&,
                     const char* test );

    MathOperatorRHS( const PropertyDatabase<dim>&,
                     const char* oper,
                     const char* test );

    MathOperatorRHS& operator=( const MathOperatorRHS& );

    virtual ~MathOperatorRHS();

    /// Operand Functions
    std::string   Name() const ;
    void          Name( const char*, const char* topname );
    void          Name( const char*, const char* opname, const char* topname );

    virtual void  Out() const;

    /// op
    const Parameter& MaterialOperand()          const;
    const Index&  MaterialOperandKey()          const;
    std::string   MaterialOperandName()         const;
    VARIABLE_TYPE MaterialOperandType()         const;
    PLACEMENT     MaterialOperandPlacement()    const;
    size_t        MaterialOperandDataDepth()    const;

    /// bop
    const Parameter& BasicOperand()             const;
    const Index&  BasicOperandKey()             const;
    std::string   BasicOperandName()            const;
    VARIABLE_TYPE BasicOperandType()            const;
    PLACEMENT     BasicOperandPlacement()       const;
    size_t        BasicOperandDataDepth()       const;
    size_t        BasicOperandOffset()          const;
    void          BasicOperandOffset( size_t );

    /// top
    const Parameter& TestOperand()              const;
    const Index&  TestOperandKey()              const;
    std::string   TestOperandName()             const;
    VARIABLE_TYPE TestOperandType()             const;
    PLACEMENT     TestOperandPlacement()        const;
    size_t        TestOperandDataDepth()        const;
    size_t        TestOperandOffset()           const;
    void          TestOperandOffset( size_t );

    /// Accumulation Process Settings

    /// Get Flags and Properties
    bool          Add()                         const;
    bool          Subtract()                    const;
    bool          AddLater()                    const;
    bool          SubtractLater()               const;
    bool          Multiply()                    const;
    bool          LumpedFormulation()           const;
    size_t        ApplicationCycle()            const;
    size_t        ApplicationCycles()           const;
    double64      MultiplyBy()                  const;
    bool          MultiplyWithTimeIncrement()   const;
    bool          DivideByTimeIncrement()       const;

    /// Set Flags and Properties
    void          AddAccumulate();
    void          SubtractAccumulate();
    void          AddAccumulateLater();
    void          SubtractAccumulateLater();
    void          MultiplyAccumulate();
    void          LumpedFormulation ( bool );
    void          ApplicationCycle  ( size_t        );
    void          ApplicationCycles ( size_t        );
    void          MultiplyBy( double64 integral_mult_factor );
    void          MultiplyWithTimeIncrement( bool multiply );
    void          DivideByTimeIncrement( bool divide );

    /// interpolation of property if isoparametric elements are used
    void          PropertyAtIntegrationPoint( const Element<dim>&,
                                              const csmp::Index&,
                                              size_t ip, DenseMatrix<DM_MIN>& );
    void          PropertyAtIntegrationPoint( const Face<dim>&,
                                              const csmp::Index&,
                                              size_t ip, DenseMatrix<DM_MIN>& );

    /// getting data from the Element, Face, InterFace
    virtual void  GetOperands        ( Element<dim>& );
    virtual void  GetOperands        ( Face<dim>& );
    virtual void  GetOperands        ( InterFace<dim>& );

    /// writing data to the Element, Face, InterFace
    virtual void  WriteOperands      ( Element<dim>& );
    virtual void  WriteOperands      ( Face<dim>&  );
    virtual void  WriteOperands      ( InterFace<dim>& );

    /// integration performed on Element, Face, InterFace
    virtual void  ComputeContribution( Element<dim>& );
    virtual void  ComputeContribution( Face<dim>& );
    virtual void  ComputeContribution( InterFace<dim>& );
  
    /// if so specified multiply with time increment
    virtual void  MultiplyWithTimeFactor( double64 dt );

    /// assigment to the right hand side global vector (after everything was calculated )
    virtual void  AssignToGlobal( const Element<dim>&,   std::vector<double64>& rhs );
    virtual void  AssignToGlobal( const Face<dim>&,      std::vector<double64>& rhs );
    virtual void  AssignToGlobal( const InterFace<dim>&, std::vector<double64>& rhs );

    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions, Luat Khoa Tran)
    virtual void  AssignToGlobal(const Element<dim>&, std::vector<double64>& rhs, const std::vector<size_t>&  );
    virtual void  AssignToGlobal(const Face<dim>&, std::vector<double64>& rhs, const std::vector<size_t>&  );
    virtual void  AssignToGlobal(const InterFace<dim>&, std::vector<double64>& rhs, const std::vector<size_t>& );

    virtual MathOperatorRHS<dim>* clone() const = 0;

  protected:
    std::string                         name_;               ///< name of operator

    Parameter                           op;                  ///< material property operand
    std::pair<Parameter, size_t >       top;                 ///< test function operand

    std::vector<double64>               RHS;                 ///< solution vector<double64> to be accumulated
    std::vector<size_t>                 IDT;                 ///< node-ID & global constraint points vector ( test operand )

    std::vector<DenseMatrix<DM_MIN> >   MTRL;                ///< material property matrix(es) needed for PDE operand
    DenseMatrix<DM_MIN>                 DERIV;               ///< shape function derivative matrix
    std::vector<double64>               IPOL;                ///< shape function vector
    std::vector<ScalarVariable >        SC;                  ///< node property vector<double64> of scalars
    std::vector<VectorVariable<dim> >   VC;                  ///< vectors
    std::vector<TensorVariable<dim> >   TS;                  ///< tensors
    std::vector<ArrayVariable >         AR;                  ///< arrays
    std::vector<FlaggedArrayVariable >  FR;                  ///< flagged arrays

    double64                            factor_;             ///< constant factor

    /// specifies accumulation procedure
    bool        add_accumulate_;
    bool        subtract_accumulate_;
    bool        add_accumulate_later_;
    bool        subtract_accumulate_later_;
    bool        multiply_accumulate_;
    bool        lump_matrices_;

    size_t      application_cycles_;
    size_t      application_cycle_;

    /// time-dependent multipliers
    bool        time_multiply_;
    bool        time_divide_;

};

// for multi-dimensional solution variables
void transformNodeIndexVector( size_t dim, const csmp::Index&, std::vector<size_t>& );



} // csmp



#endif



















