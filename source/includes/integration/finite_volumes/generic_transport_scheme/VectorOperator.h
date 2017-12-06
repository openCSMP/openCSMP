#ifndef CSMP_VECTOR_OPERATOR_H
#define CSMP_VECTOR_OPERATOR_H

#include "GenericTransportScheme.h"

namespace csmp {

/**
\brief     Vector operator
\details   Part of the Colleoli transport scheme.
\author    
\version   0a
\date      6/12/2017
\pre       base class for matrix operators
\bug
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
class VectorOperator
{
public:
  ACCUMULATION_MODE AccumulationMode() const { return mode_; }
  
  void AccumulationMode( ACCUMULATION_MODE mode ) { mode_ = mode; }
  
  bool MultiplyWithTimeIncrement() const { return multiply_with_dt_; }
  void MultiplyWithTimeIncrement( bool multiply_with_dt ) { multiply_with_dt_ = multiply_with_dt; }
  
  void TimeIncrement( double64 dt ) { dt_ = dt; }
  virtual void AccumulateStencil( const Element<dim>* fe, std::vector<double64>& rhs ) const = 0;
  virtual void AccumulateFiniteVolume( const Node<dim>* fv, std::vector<double64>& rhs ) const = 0;
  
  virtual ~VectorOperator() { }
  
protected:
  VectorOperator( ACCUMULATION_MODE mode )
  : mode_(mode), multiply_with_dt_(false), dt_(std::numeric_limits<double64>::quiet_NaN())
  {
  }
  
  ACCUMULATION_MODE mode_;
  bool multiply_with_dt_;
  double64 dt_;
};



} // end csmp


#endif /* CSMP_VECTOR_OPERATOR_H */
