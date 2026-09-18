// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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

template<uint32_t> class Element;
template<uint32_t> class PropertyDatabase;

/**
@brief Base class for FE or FV integrals accumulated into righthand vector.

@attention In a linear systemA x = b, every single row of that matrix A represents one equation.
Each equation corresponds to one specific Test Function basis.

Therefore, the rows are determined by the Test Operand.

The columns represent the degrees of freedom (DOF) of your unknown variable 'u''.
These correspond to the Basic (Trial) Functions.
Therefore, the columns are determined by the Basic Operand and do not appear in the righthand side.

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
      the node property vector<double> will be of dimension (nds+ips) * dof
      in this case, the entries will be 'dof' per node or constraint point in the
      order as given above.

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class MathOperatorRHS {
    MathOperatorRHS();

  public:

    MathOperatorRHS( const MathOperatorRHS& );

    MathOperatorRHS( const PropertyDatabase<dim>&,
                     const char* test );

    MathOperatorRHS( const PropertyDatabase<dim>&,
                     const char* oper,
                     const char* test );

    MathOperatorRHS& operator=( const MathOperatorRHS& );

    virtual ~MathOperatorRHS() = default;

    /// Operand Functions
    std::string   Name() const noexcept;
    void          Name( const char*, const char* topname ) noexcept;
    void          Name( const char*, const char* opname, const char* topname ) noexcept;

    virtual void  Out() const;

    /// op
    const Parameter& MaterialOperand()          const noexcept;
    const Index&  MaterialOperandKey()          const noexcept;
    std::string   MaterialOperandName()         const noexcept;
    VARIABLE_TYPE MaterialOperandType()         const noexcept;
    PLACEMENT     MaterialOperandPlacement()    const noexcept;
    uint32_t      MaterialOperandDataDepth()    const noexcept;

    /// bop
    const Parameter& BasicOperand()             const noexcept;
    const Index&  BasicOperandKey()             const noexcept;
    std::string   BasicOperandName()            const noexcept;
    VARIABLE_TYPE BasicOperandType()            const noexcept;
    PLACEMENT     BasicOperandPlacement()       const noexcept;
    uint32_t      BasicOperandDataDepth()       const noexcept;
    size_t        BasicOperandOffset()          const noexcept;
    void          BasicOperandOffset( size_t ) noexcept;

    /// top
    const Parameter& TestOperand()              const noexcept;
    const Index&  TestOperandKey()              const noexcept;
    std::string   TestOperandName()             const noexcept;
    VARIABLE_TYPE TestOperandType()             const noexcept;
    PLACEMENT     TestOperandPlacement()        const noexcept;
    uint32_t      TestOperandDataDepth()        const noexcept;
    size_t        TestOperandOffset()           const noexcept;
    void          TestOperandOffset( size_t ) noexcept;

    /// Accumulation Process Settings

    /// Get Flags and Properties
    bool          Add()                         const noexcept;
    bool          Subtract()                    const noexcept;
    bool          AddLater()                    const noexcept;
    bool          SubtractLater()               const noexcept;
    bool          Multiply()                    const noexcept;
    bool          LumpedFormulation()           const noexcept;
    uint32_t      ApplicationCycle()            const noexcept;
    uint32_t      ApplicationCycles()           const noexcept;
    [[nodiscard]] double MultiplyBy()           const noexcept;
    bool          MultiplyWithTimeIncrement()   const noexcept;
    bool          DivideByTimeIncrement()       const noexcept;

    /// Set Flags and Properties
    void          AddAccumulate() noexcept;
    void          SubtractAccumulate() noexcept;
    void          AddAccumulateLater() noexcept;
    void          SubtractAccumulateLater() noexcept;
    void          MultiplyAccumulate() noexcept;
    void          LumpedFormulation ( bool ) noexcept;
    void          ApplicationCycle  ( uint32_t ) noexcept;
    void          ApplicationCycles ( uint32_t ) noexcept;
    void          MultiplyBy( double integral_mult_factor ) noexcept;
    void          MultiplyWithTimeIncrement( bool multiply ) noexcept;
    void          DivideByTimeIncrement( bool divide ) noexcept;

    /// interpolation of property if isoparametric elements are used
    void          PropertyAtIntegrationPoint( const CELL<dim>&,
                                              const csmp::Index&,
                                              uint32_t ip, DenseMatrix<dim>& mtrl );

    /// getting data from the Element, Face, InterFace
    virtual void  GetOperands( const CELL<dim>& );

    /// integration performed on Element, Face, InterFace
    virtual void  ComputeContribution( const CELL<dim>& );
  
    /// writing data to the Element, Face, InterFace
    virtual void  WriteOperands( CELL<dim>& );

    /// if so specified multiply with time increment
    virtual void  MultiplyWithTimeFactor( double dt ) noexcept;

    /// assigment to the right hand side global vector (after everything was calculated )
    void  AssignToGlobal( const CELL<dim>&, std::vector<double>& rhs );

    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions, Luat Khoa Tran)
    void  AssignToGlobal( const CELL<dim>&, std::vector<double>& rhs, const std::vector<size_t>& );

    virtual MathOperatorRHS<dim,CELL>* clone() const = 0;

  protected:
    std::string                   name_;               ///< name of operator

    Parameter                      op;                  ///< material property operand
    std::pair<Parameter, size_t>   top;                 ///< test function operand

  // TODO: remove these non members to avoid race conditions during parallel accumulation
    std::vector<DenseMatrix<dim> > MTRL;                ///< material property matrix(es) needed for PDE operand
    std::vector<double>            RHS;                 ///< righthand vector<double> to be accumulated

    double                         factor_;             ///< constant factor

    /// specifies accumulation procedure
    bool        add_accumulate_;
    bool        subtract_accumulate_;
    bool        add_accumulate_later_;
    bool        subtract_accumulate_later_;
    bool        multiply_accumulate_;
    bool        lump_matrices_;

    uint32_t    application_cycles_;
    uint32_t    application_cycle_;

    /// time-dependent multipliers
    bool        time_multiply_;
    bool        time_divide_;

};

/// expands element dof vec in var1_comp0, var1_comp2... form
void transformNodeIndexVector( const csmp::Index& test_function_operand_idx, std::vector<size_t>& ) noexcept;

} // csmp



#endif



















