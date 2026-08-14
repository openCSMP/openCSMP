#ifndef Upwind_Integral_dNT_op_dN_dV_h
#define Upwind_Integral_dNT_op_dN_dV_h

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/// Known as: element conductance matrix or K div^2 P
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Upwind_Integral_dNT_op_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>&,
                           		  const char* oper,
                           		  const char* basic,
                           		  const char* test,
                           		  const char* upwind,
                           		  const char* trigger,
                           		  const double prefactor = 1. );
    
    void ComputeContribution( const CELL<dim>& ) override final;
    void GetOperands( const CELL<dim>& ) override final;

     Upwind_Integral_dNT_op_dN_dV<dim,CELL>* clone() const override final
      { return new Upwind_Integral_dNT_op_dN_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  DN, DNT;
    
    Index uvar_;
    Index tvar_;
    
    std::vector<ScalarVariable > el_uvar;
    std::vector<ScalarVariable > el_tvar;
    
    const double prefactor_;
};


/**
 
@class Upwind_Integral_dNT_op_dN_dV Upwind_Integral_dNT_op_dN_dV "pde_operators/Upwind_Integral_dNT_op_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999

@section motivation Motivation
 
PDE operator representing the divergence squared of the dependent 
variable.
*/

} // csmp

#endif
