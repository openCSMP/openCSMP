#ifndef CSMP_MATRIX_OPERATOR_H
#define CSMP_MATRIX_OPERATOR_H

#include "GenericTransportScheme.h"

namespace csmp {

/**

\brief     Matrix operator
\details   Part of the Colleoli transport scheme.
\version   0a
\date      6/12/2017
\pre       base class for matrix operators
\warning
\copyright The University of Melbourne

@section motivation Motivation


@section design Design Intent


@section applicability Applicability


@section collaborations Collaborations


@section implementation Implementation


@section examples Application Examples

@code

@endcode

*/
template<size_t dim>
class MatrixOperator
{
public:
    size_t Stage() const { return stage_; }
  
    void Stage( size_t stage ) { stage_ = stage; }

    bool MultiplyWithTimeIncrement() const { return multiply_with_dt_; }
    void MultiplyWithTimeIncrement( bool multiply_with_dt ) { multiply_with_dt_ = multiply_with_dt; }

    void TimeIncrement( double64 dt ) { dt_ = dt; }
    virtual void AccumulateStencil( Element<dim>& fe, SparseMatrix& lhs ) const = 0;
    virtual void AccumulateFiniteVolume( Node<dim>& fv, SparseMatrix& lhs ) const = 0;

    virtual ~MatrixOperator() { }

protected:
    explicit MatrixOperator( size_t stage )
      : stage_(stage), multiply_with_dt_(false), dt_(std::numeric_limits<double64>::quiet_NaN())
    {
    }
  
    size_t stage_;
    bool multiply_with_dt_;
    double64 dt_;
};


} // end csmp


#endif /* CSMP_MATRIX_OPERATOR_H */
