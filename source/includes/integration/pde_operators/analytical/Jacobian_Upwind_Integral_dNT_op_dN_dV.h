#ifndef Jacobian_Upwind_Integral_dNT_op_dN_dV_h
#define Jacobian_Upwind_Integral_dNT_op_dN_dV_h

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/// Known as: jacobian of element conductance matrix or K div^2 P
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Jacobian_Upwind_Integral_dNT_op_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    Jacobian_Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>&,
                                            const char* oper,
                                            const char* basic,
                                            const char* test,
                                            const char* test_orig,
                                            const char* upwind,
                                            const char* d_upwind,
                                            const char* trigger,
                                            const double delta,
                                            const double prefactor = 1.0);
    
    void ComputeContribution( const CELL<dim>& );
    void GetOperands( const CELL<dim>& );

  private:
    DenseMatrix<DM_MIN> DN, DNT;
    
    Index upwind_;
    Index trigger_;
    Index test_orig_; // fluid pressure
    Index d_upwind_; // lambda (S+dS)
    
    std::vector<ScalarVariable > el_upwind;
    std::vector<ScalarVariable > el_trigger;
    std::vector<ScalarVariable > el_test_orig;
    std::vector<ScalarVariable > el_d_upwind;
    
    const double delta_;
    const double prefactor_;
};


/**
@class Jacobian_Upwind_Integral_dNT_op_dN_dV Jacobian_Upwind_Integral_dNT_op_dN_dV "pde_operators/Jacobian_Upwind_Integral_dNT_op_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999
 

PDE operator for the tangent stiffness matrix, using an upwind formulation.
It differs from the usual Upwind_Integral_dNT_op_dN_dV in the test-function,
since the placement in the global matrix is different.
*/

} // csmp

#endif
