#ifndef LHS_INTEGRAL_DNT_DN_DV_H
#define LHS_INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

/**
    Interpolation function derivative matrix squared.
    
    @author SKM
    @date 30/01/2018
*/
template<size_t dim,class SIMPLEX=Element<dim> >
class LHS_Integral_dNT_dN_dV : public MathOperatorLHS<dim> {
  public:
    /// @param basic weighting function operand, @param test interpolation function operand
    LHS_Integral_dNT_dN_dV( const PropertyDatabase<dim>&,
                            const char* basic,
                            const char* test );
  
    virtual void GetOperands( const SIMPLEX& e ) { /* noting to do here since there is no material operand */ }
    
    virtual void ComputeContribution( const SIMPLEX& );
  
    virtual LHS_Integral_dNT_dN_dV<dim,SIMPLEX>* clone() const { return new LHS_Integral_dNT_dN_dV<dim,SIMPLEX> (*this); }
  
  private:
    DenseMatrix<DM_MIN>          DN, DNT;
    const std::vector<double>  unity;
};

} // csmp

#endif
