#ifndef INTEGRAL_DNT_DN_DV_H
#define INTEGRAL_DNT_DN_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
    Interpolation function derivatives squared.
    No material operand required.
    
    @author S.K. Matthai
    @author S. Roberts
    @date 1999

    RHS PDE operator representing the divergence squared of the dependent
    variable.

*/
template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_dNT_dN_dV : public MathOperatorRHS<dim> {
  public:
    Integral_dNT_dN_dV( const PropertyDatabase<dim>&, const char* test );
                        
    virtual ~Integral_dNT_dN_dV() {}

    virtual void ComputeContribution( const SIMPLEX& );
  
    virtual Integral_dNT_dN_dV<dim,SIMPLEX>* clone() const { return new Integral_dNT_dN_dV<dim,SIMPLEX> (*this); }

  private:
    DenseMatrix<DM_MIN>  DN, DNT, UNITY; 
};


/**
 
@class Integral_dNT_dN_dV Integral_dNT_dN_dV "pde_operators/Integral_dNT_dN_dV.h"
*/

} // csmp

#endif
