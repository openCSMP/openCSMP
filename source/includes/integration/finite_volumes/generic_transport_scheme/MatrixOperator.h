#ifndef CSMP_MATRIX_OPERATOR_H
#define CSMP_MATRIX_OPERATOR_H

#include "CSMP_definitions.h"
#include "SparseMatrix.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Node;

/**

\brief     Matrix operator
\details   Part of the Colleoli transport scheme.
\version   0a
\date      6/12/2017
\pre       base class for matrix operators

\copyright The University of Melbourne

TODO: move MultiplyWithDt etc to operation
TODO: find way of treating terms in equation as groups in a specific way
TODO: perhaps chain operations like mathematic expressions
TODO: perhaps perform multiple sequential operations in an element by element fashion

*/
template<size_t dim>
class MatrixOperator {
  public:
    MatrixOperator() : factor_(1.) {}
    virtual ~MatrixOperator() {}

    virtual void AccumulateFiniteVolume( const Node<dim>& fv, SparseMatrix& lhs ) const = 0;

    virtual void AccumulateStencil( const Element<dim>& fe, SparseMatrix& lhs ) const = 0;
//    virtual void AccumulateStencil( Face<dim>& fe, SparseMatrix& lhs ) const = 0;

    /// includes the time increment in the multiplication factor for this operator
    void MultiplyWithTimeIncrement( double64 dt ) { factor_ *= dt; }

    /// set factor to achieve multiplication with time increment (fac=dt), subtraction (fac=-1), multiplication or division (fac=1/value)
    void     Factor( double64 value ) { factor_ = value; }
    double64 Factor() const { return factor_; }

    // virtual void Out() const = 0; TODO: rather use verbose function to print term to be accumulated

  private:
    double64 factor_ = 1.; ///<  1=add, -1=subtract, factor=val = multiply, factor=1/val divide,  factor=time_increment if so needed
};


} // end csmp


#endif /* CSMP_MATRIX_OPERATOR_H */
