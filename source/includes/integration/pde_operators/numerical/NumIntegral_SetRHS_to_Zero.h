#ifndef NUMINTEGRAL_SET_RHS_TO_ZERO_H
#define NUMINTEGRAL_SET_RHS_TO_ZERO_H

#include "MathOperatorRHS.h"

namespace csmp {

/// invokes accumulation of zero right-hand side = "homogeneous" boundary conditions
template<uint32_t dim,class CELL=Element<dim> >
class NumIntegral_SetRHS_to_Zero : public MathOperatorRHS<dim> {
  public:
    NumIntegral_SetRHS_to_Zero( const PropertyDatabase<dim>&,
                                const char* test );
    
    virtual ~NumIntegral_SetRHS_to_Zero();
    
    /// since there is no material Operand nothing needs to be done
    virtual void GetOperands( const CELL& ) {}

    virtual void ComputeContribution( const CELL& );

    virtual NumIntegral_SetRHS_to_Zero<dim,CELL>* clone() const
      { return new NumIntegral_SetRHS_to_Zero<dim,CELL> (*this); }
};

} // csmp

#endif
















