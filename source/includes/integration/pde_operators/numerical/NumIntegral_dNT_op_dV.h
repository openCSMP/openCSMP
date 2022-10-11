#ifndef NUM_INTEGRAL_DNT_OP_DV_H
#define NUM_INTEGRAL_DNT_OP_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/// to integrate over a gradient represented by a vector property
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNT_op_dV( const PropertyDatabase<dim>&,
                           const char* oper,    // (vector) gradient property, e.g., rho g grad z
                           const char* test );  // scalar, for instance fluid pressure
    
    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual NumIntegral_dNT_op_dV<dim,CELL>* clone() const { return new NumIntegral_dNT_op_dV<dim,CELL> (*this); }
    
  private:
    DenseMatrix<DM_MIN>  B_, BT_; 
};

} // csmp

#endif
