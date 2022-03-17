#ifndef NUM_INTEGRAL_DNT_DN_DV_H
#define NUM_INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  div^2 N = interpolation function derivate matrix squared.
template<uint32_t dim,class CELL=Element<dim> >
class NumIntegral_dNT_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_dN_dV( const PropertyDatabase<dim>& pref, 
                           const char* basic, 
                           const char* test );
  
    /// no operands need to be fetched from computational domain
    virtual void GetOperands( const CELL& ) {}
  
    virtual void ComputeContribution( const CELL& );
  
    virtual NumIntegral_dNT_dN_dV<dim,CELL >* clone() const { return new NumIntegral_dNT_dN_dV<dim,CELL >(*this); }
    
  private:
    DenseMatrix<DM_MIN>  B, BT; 
};

} // csmp

#endif
