#ifndef NUMINTEGRAL_SET_RHS_TO_ZERO_H
#define NUMINTEGRAL_SET_RHS_TO_ZERO_H

#include "MathOperatorRHS.h"

namespace csmp {

/// invokes accumulation of zero right-hand side = "homogeneous" boundary conditions
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_SetRHS_to_Zero : public MathOperatorRHS<dim> {
  public:
    NumIntegral_SetRHS_to_Zero( const PropertyDatabase<dim>& p,
                                const char* test );
    
    virtual ~NumIntegral_SetRHS_to_Zero();
    
    virtual void ComputeContribution( SIMPLEX& );

    virtual void GetOperands( SIMPLEX& );
    virtual NumIntegral_SetRHS_to_Zero<dim,SIMPLEX>* clone() const
      { return new NumIntegral_SetRHS_to_Zero<dim,SIMPLEX> (*this); }
};

// since there is no material Operand nothing needs to be done
template<size_t dim,class SIMPLEX>
inline void NumIntegral_SetRHS_to_Zero<dim,SIMPLEX>::GetOperands( SIMPLEX& ) {}

} // csp

#endif
















