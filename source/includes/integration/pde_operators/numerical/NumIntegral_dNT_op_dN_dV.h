#ifndef NUM_INTEGRAL_DNT_OP_DN_DN_H
#define NUM_INTEGRAL_DNT_OP_DN_DN_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  oper div^2 N - interpolation functions squared.
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_dNT_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                        const char* oper, 
                        const char* basic, 
                        const char* test );
    
    virtual void ComputeContribution( SIMPLEX& e );
  
    virtual NumIntegral_dNT_op_dN_dV<dim,SIMPLEX >* clone() const { return new NumIntegral_dNT_op_dN_dV<dim,SIMPLEX >(*this); }
    
  protected:
    DenseMatrix<DM_MIN>  B, BT; 
};

} // csmp

#endif
