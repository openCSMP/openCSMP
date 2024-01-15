#ifndef NUMINTEGRAL_SET_RHS_TO_ONE_H
#define NUMINTEGRAL_SET_RHS_TO_ONE_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// for a scalar solution variable this operator creates a unit vector for accumulation into righthand side vector
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_SetRHS_to_One : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_SetRHS_to_One( const PropertyDatabase<dim>&,
                               const char* test );
    
    virtual ~NumIntegral_SetRHS_to_One();
    
    /// since there is no material Operand nothing needs to be done
    virtual void GetOperands( const CELL<dim>& ) {}

    virtual void ComputeContribution( const CELL<dim>& );

    virtual NumIntegral_SetRHS_to_One<dim,CELL>* clone() const
      { return new NumIntegral_SetRHS_to_One<dim,CELL> (*this); }
};

} // csmp

#endif
















