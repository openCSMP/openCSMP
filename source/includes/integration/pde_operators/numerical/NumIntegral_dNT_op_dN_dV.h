#ifndef NUM_INTEGRAL_DNT_OP_DN_DN_H
#define NUM_INTEGRAL_DNT_OP_DN_DN_H

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/// PDE operator:  oper div^2 N - interpolation functions squared.
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_op_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>&, 
                              const char* oper,
                              const char* basic,
                              const char* test );
                              
    virtual ~NumIntegral_dNT_op_dN_dV();
    
    virtual void ComputeContribution( const CELL<dim>& );
  
    virtual NumIntegral_dNT_op_dN_dV<dim,CELL>* clone() const { return new NumIntegral_dNT_op_dN_dV<dim,CELL>(*this); }
    
  protected:
    DenseMatrix<DM_MIN>  B_, BT_;
};

} // csmp

#endif
