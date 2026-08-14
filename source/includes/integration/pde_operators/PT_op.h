#ifndef PT_OP_H
#define PT_OP_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
 
@brief Point sources, loads etc. for vector solution variable.
In mechanics, this operator is known as nodal force vector.

@author S.K. Matthai
@author S. Roberts
@date 2000

Use this operand to accumulate nodal vector-type properties flagged 
Dirichlet into the righthand side of the global matrix equation.

The forces are distributed onto the Dirichlet nodes such that nodes which
are shared by two elements are accumulated correctly:

F(midpoint node) = (F / boundary-nodes)
F(corner node)   = (F / boundary-nodes) / 2

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class PT_op : public MathOperatorRHS<dim,CELL> {
  public:
    PT_op( const PropertyDatabase<dim>&, const char* oper, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    PT_op<dim,CELL>* clone() const override { return new PT_op<dim,CELL> (*this); }
    
  private:
    std::vector<VectorVariable<dim> >  NODAL_FORCE;
};


} // csmp

#endif
















