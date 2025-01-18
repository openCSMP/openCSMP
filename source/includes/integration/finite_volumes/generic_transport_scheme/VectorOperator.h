#ifndef CSMP_VECTOR_OPERATOR_H
#define CSMP_VECTOR_OPERATOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Node;

/**
\brief     Vector operator
\details   Part of the Colleoli transport scheme.
\author    Stephan Matthai
\version   0a
\date      6/12/2017
\pre       base class for matrix operators
\copyright The University of Melbourne

*/
template<uint32_t dim>
class VectorOperator {
  public:
    VectorOperator() : factor_(1.) {}
    virtual ~VectorOperator() {}
    
    virtual void AccumulateFiniteVolume( const Node<dim>&, std::vector<double>& rhs ) const = 0;

    virtual void AccumulateStencil( const Element<dim>&, std::vector<double>& rhs ) const = 0;
//    virtual void AccumulateStencil( Face<dim>&, std::vector<double>& rhs ) const = 0;
//    virtual void AccumulateStencil( InterFace<dim>&, std::vector<double>& rhs ) const = 0;

    /// includes the time increment in the multiplication factor for this operator
    void MultiplyWithTimeIncrement( double dt ) { factor_ *= dt; }

    /// set factor to achieve multiplication with time increment (fac=dt), subtraction (fac=-1), multiplication or division (fac=1/value)
    void     Factor( double value ) { factor_ = value; }
    double Factor() const { return factor_; }
    
    // virtual void Out() const = 0; TODO: use verbose function to print term to be accumulated

  private:
    double factor_ = 1.; ///<  1=add, -1=subtract, factor=val = multiply, factor=1/val divide,  factor=time_increment if so needed
};



} // end csmp


#endif /* CSMP_VECTOR_OPERATOR_H */
