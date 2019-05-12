#ifndef POINT_SOURCE_LHSOP_H
#define POINT_SOURCE_LHSOP_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/** @brief Assignment of scalar nodal (Neumann) source or sink terms.

    This operator applies (absolute) nodal source terms to the LHS matrix diagonal as they are.
    The nodal values must be precomputed elsewhere.

    @attention the nodal term gets applied as many times as the node is a member of an element, i.e.
    n.Parents() times. Therefore the value should be divided by the number of parent elements of the node.
    This is already accounted for by the method GetOperands().

@date 12/5/2019
@author Stephan K. Matthai

*/      
template<size_t dim,class SIMPLEX=Element<dim> >
class PointSource_lhsop : public MathOperatorLHS<dim> {
  public:
    PointSource_lhsop( const PropertyDatabase<dim>&, const char* nodal_src, const char* basic, const char* test );
    virtual void GetOperands( SIMPLEX& );
    virtual void ComputeContribution( SIMPLEX& );
    virtual PointSource_lhsop<dim,SIMPLEX>* clone() const { return new PointSource_lhsop<dim,SIMPLEX> (*this); }
  private:
    std::vector<ScalarVariable>  SRC_;
};

} // csmp

#endif
















