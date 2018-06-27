#ifndef NUMINTEGRAL_SET_RHS_TO_ONE_H
#define NUMINTEGRAL_SET_RHS_TO_ONE_H

#include "MathOperatorRHS.h"

namespace csmp {

/// for a scalar solution variable this operator creates a unit vector for accumulation into righthand side vector
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_SetRHS_to_One : public MathOperatorRHS<dim> {
  public:
    NumIntegral_SetRHS_to_One( const PropertyDatabase<dim>& p,
                               const char* test );
    
    virtual ~NumIntegral_SetRHS_to_One();
    
    /// since there is no material Operand nothing needs to be done
    virtual void GetOperands( CELL& ) {}

    virtual void ComputeContribution( CELL& e );

    virtual NumIntegral_SetRHS_to_One<dim,CELL>* clone() const
      { return new NumIntegral_SetRHS_to_One<dim,CELL> (*this); }
};

} // csmp

#endif
















