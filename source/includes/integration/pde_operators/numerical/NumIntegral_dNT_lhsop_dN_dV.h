#ifndef NUM_INTEGRAL_DNT_LHSOP_DN_DN_H
#define NUM_INTEGRAL_DNT_LHSOP_DN_DN_H

#include "MathOperatorLHS.h"
#include "FiniteElement.h"
#include <Eigen/Dense>

namespace csmp {

template<uint32_t> class Element;

/// PDE operator:  oper div^2 N - interpolation functions squared = conductance matrix
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_lhsop_dN_dV final : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_lhsop_dN_dV( const PropertyDatabase<dim>&, 
                              const char* oper,
                              const char* basic,
                              const char* test );
                              
    void ComputeContribution( const CELL<dim>& ) override final;
  
    NumIntegral_dNT_lhsop_dN_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNT_lhsop_dN_dV<dim,CELL>(*this); }
};

} // csmp

#endif
