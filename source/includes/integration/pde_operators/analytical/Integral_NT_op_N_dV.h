#ifndef INTEGRAL_NT_OP_N_DV_H
#define INTEGRAL_NT_OP_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {
  /**
  @author S.K. Matthaei
  @author S. Roberts
  @date 1999 */

template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_NT_op_N_dV : public MathOperatorRHS<dim> {
  public:
    Integral_NT_op_N_dV( const PropertyDatabase<dim>& p, 
                         const char* oper, const char* test );
    
    virtual void GetOperands( const SIMPLEX& );
    // integration etc.
    virtual void ComputeContribution( const SIMPLEX& );
    
    virtual Integral_NT_op_N_dV<dim,SIMPLEX>* clone() const { return new Integral_NT_op_N_dV<dim,SIMPLEX> (*this); }
  private:
    DenseMatrix<DM_MIN>  INN;
    ScalarVariable      sc;
};



/**
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















