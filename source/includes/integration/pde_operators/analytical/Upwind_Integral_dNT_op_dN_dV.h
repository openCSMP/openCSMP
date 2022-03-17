#ifndef Upwind_Integral_dNT_op_dN_dV_h
#define Upwind_Integral_dNT_op_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "DenseMatrix.h"

namespace csmp {

/// Known as: element conductance matrix or K div^2 P
template<uint32_t dim,class CELL=Element<dim> >
class Upwind_Integral_dNT_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                           		  const char* oper,
                           		  const char* basic,
                           		  const char* test,
                           		  const char* upwind,
                           		  const char* trigger,
                           		  const double prefactor = 1. );
    
    void ComputeContribution( const CELL& );
    void GetOperands( const CELL& );

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
