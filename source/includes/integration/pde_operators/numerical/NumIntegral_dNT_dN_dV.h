#ifndef NUM_INTEGRAL_DNT_DN_DV_H
#define NUM_INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  div^2 N = interpolation function derivate matrix squared.
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_dNT_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_dN_dV( const PropertyDatabase<dim>& pref, 
                           const char* basic, 
                           const char* test );
    
    virtual void GetOperands( SIMPLEX& );
    virtual void ComputeContribution( SIMPLEX& );
    virtual NumIntegral_dNT_dN_dV<dim,SIMPLEX >* clone() const { return new NumIntegral_dNT_dN_dV<dim,SIMPLEX >(*this); }
  private:
    DenseMatrix<DM_MIN>  B, BT; 
};



/// since there is no material Operand nothing needs to be done
template<size_t dim,class SIMPLEX>
inline void NumIntegral_dNT_dN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& )
 {
 }

} // csmp

#endif
