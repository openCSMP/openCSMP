#ifndef Jacobian_Upwind_Integral_dNT_op_dN_dV_h
#define Jacobian_Upwind_Integral_dNT_op_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "DenseMatrix.h"

namespace csmp {

/// Known as: jacobian of element conductance matrix or K div^2 P
template<size_t dim,class SIMPLEX=Element<dim> >
class Jacobian_Upwind_Integral_dNT_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    Jacobian_Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                           		  const char* oper,
                           		  const char* basic,
                           		  const char* test,
                           		  const char* test_orig,
                           		  const char* upwind,
                           		  const char* d_upwind,
                           		  const char* trigger,
                           		  const double64 delta,
                           		  const double64 prefactor = 1.0);
    
    void ComputeContribution( SIMPLEX& e );
    void GetOperands( SIMPLEX& e );

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
    
    const double64 delta_;
    const double64 prefactor_;
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
