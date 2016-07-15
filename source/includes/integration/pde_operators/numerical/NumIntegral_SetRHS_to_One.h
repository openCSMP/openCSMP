#ifndef NUMINTEGRAL_SET_RHS_TO_ONE_H
#define NUMINTEGRAL_SET_RHS_TO_ONE_H

#include "MathOperatorRHS.h"

namespace csmp {

/// scalar solution variable: invokes accumulation of a right-hand side equal to unity
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_SetRHS_to_One : public MathOperatorRHS<dim> {
  public:
    NumIntegral_SetRHS_to_One( const PropertyDatabase<dim>& p,
                               const char* test );
    
    virtual ~NumIntegral_SetRHS_to_One();
    
    virtual void ComputeContribution( SIMPLEX& e );

    virtual void GetOperands( SIMPLEX& );

    virtual NumIntegral_SetRHS_to_One<dim,SIMPLEX>* clone() const
      { return new NumIntegral_SetRHS_to_One<dim,SIMPLEX> (*this); }
};

// since there is no material Operand nothing needs to be done
template<size_t dim,class SIMPLEX>
inline void NumIntegral_SetRHS_to_One<dim,SIMPLEX>::GetOperands( SIMPLEX& ) {}

} // csmp

#endif
















