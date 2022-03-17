#ifndef INTEGRAL_NT_N_DV_H
#define INTEGRAL_NT_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// integral test function products (no variable is actually used)
template<uint32_t dim,class CELL=Element<dim> >
class Integral_NT_N_dV : public MathOperatorLHS<dim> {
  public:
    // variable must be the same as in the RHS (it is not used to form integral)
    Integral_NT_N_dV( const PropertyDatabase<dim>& p, const char* test_variable );
    
    virtual void GetOperands( const CELL& ) {}
    virtual void ComputeContribution( const CELL& );
};



/**
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */



} // csmp

#endif
















