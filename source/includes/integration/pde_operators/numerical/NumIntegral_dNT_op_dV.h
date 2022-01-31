#ifndef NUM_INTEGRAL_DNT_OP_DV_H
#define NUM_INTEGRAL_DNT_OP_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/// to integrate over a gradient represented by a vector property
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_dNT_op_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_dNT_op_dV( const PropertyDatabase<dim>& pref, 
                           const char*             oper,    // (vector) gradient property, e.g., rho g grad z
                           const char*             test );  // scalar, for instance fluid pressure
    
    virtual void ComputeContribution( const CELL& );
    
    virtual NumIntegral_dNT_op_dV<dim,CELL>* clone() const { return new NumIntegral_dNT_op_dV<dim,CELL> (*this); }
    
  private:
    DenseMatrix<DM_MIN>  B, BT; 
};

} // csmp

#endif
