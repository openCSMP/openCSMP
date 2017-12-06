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

    virtual void AccumulateStencil( const Element<dim>* fe, std::vector<double64>& mat ) const = 0;
    virtual void AccumulateFiniteVolume( const Node<dim>* fv, std::vector<double64>& mat ) const = 0;

    virtual ~VectorOperator() { }

protected:
    VectorOperator( ACCUMULATION_MODE mode )
       : mode_(mode)
    {
    }

    ACCUMULATION_MODE mode_;    
};



} // end csmp


#endif /* CSMP_VECTOR_OPERATOR_H */
