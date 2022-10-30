#ifndef LHS_INTEGRAL_DNT_DN_DV_H
#define LHS_INTEGRAL_DNT_DN_DV_H

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    Interpolation function derivative matrix squared.
    
    @author SKM
    @date 30/01/2018
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class LHS_Integral_dNT_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    /// @param basic weighting function operand, @param test interpolation function operand
    LHS_Integral_dNT_dN_dV( const PropertyDatabase<dim>&,
                            const char* basic,
                            const char* test );
  
    /// @note  to do here since there is no material operand
    virtual void GetOperands( const CELL<dim>& ) override {}
    virtual void ComputeContribution( const CELL<dim>& ) override;
  
    virtual LHS_Integral_dNT_dN_dV<dim,CELL>* clone() const { return new LHS_Integral_dNT_dN_dV<dim,CELL> (*this); }
  
  private:
    DenseMatrix<DM_MIN>          DN, DNT;
    const std::vector<double>  unity;
};

} // csmp

#endif
