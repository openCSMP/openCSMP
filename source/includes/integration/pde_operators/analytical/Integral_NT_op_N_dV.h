#ifndef INTEGRAL_NT_OP_N_DV_H
#define INTEGRAL_NT_OP_N_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

  /**
  @author S.K. Matthaei
  @author S. Roberts
  @date 1999 */

template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_NT_op_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_NT_op_N_dV( const PropertyDatabase<dim>&,
                         const char* oper, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_NT_op_N_dV<dim,CELL>* clone() const override final
      { return new Integral_NT_op_N_dV<dim,CELL> (*this); }
    
  private:
    DenseMatrix<DM_MIN>  INN;
    ScalarVariable       sc;
};



/**
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















