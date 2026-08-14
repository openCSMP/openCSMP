#ifndef NUMINTEGRAL_SET_RHS_TO_ZERO_H
#define NUMINTEGRAL_SET_RHS_TO_ZERO_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// invokes accumulation of zero right-hand side = "homogeneous" boundary conditions
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_SetRHS_to_Zero : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_SetRHS_to_Zero( const PropertyDatabase<dim>&,
                                const char* test );
    
    /// since there is no material Operand nothing needs to be done
    void GetOperands( const CELL<dim>& ) override final { /* no operands need to be collected */ }

    void ComputeContribution( const CELL<dim>& ) override final;

    NumIntegral_SetRHS_to_Zero<dim,CELL>* clone() const override final { return new NumIntegral_SetRHS_to_Zero<dim,CELL> (*this); }
};

} // csmp

#endif
















