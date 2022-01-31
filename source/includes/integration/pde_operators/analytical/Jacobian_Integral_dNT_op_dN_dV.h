#ifndef Jacobian_Integral_dNT_op_dN_dV_h
#define Jacobian_Integral_dNT_op_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "DenseMatrix.h"

namespace csmp {

/// Known as: jacobian of element conductance matrix or K div^2 P
template<size_t dim,class SIMPLEX=Element<dim> >
class Jacobian_Integral_dNT_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    Jacobian_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                           		  const char* oper,
                           		  const char* basic,
                           		  const char* test,
                           		  const char* test_orig,
                           		  const char* lambda,
                           		  const char* d_lambda,
                           		  const double delta,
                           		  const double prefactor = 1.0);
    
    void ComputeContribution( const SIMPLEX& );
    void GetOperands( const SIMPLEX& );

  private:
    DenseMatrix<DM_MIN> DN, DNT;
    DenseMatrix<DM_MIN> res_;
    
    Index lambda_;
    Index test_orig_; // fluid pressure
    Index d_lambda_; // lambda (S+dS)
    
    std::vector<ScalarVariable > el_lambda;
    std::vector<ScalarVariable > el_test_orig;
    std::vector<ScalarVariable > el_d_lambda;
    
    const double delta_;
    const double prefactor_;
};


/**
@class Jacobian_Integral_dNT_op_dN_dV Jacobian_Integral_dNT_op_dN_dV "pde_operators/Jacobian_Integral_dNT_op_dN_dV.h"
@author S.K. Matthaei
@author S. Roberts
@date 1999
 
@section motivation Motivation
 
PDE operator for the tangent stiffness matrix, using an upwind formulation.
It differs from the usual Upwind_Integral_dNT_op_dN_dV in the test-function,
since the placement in the global matrix is different.
 */

} // csmp

#endif
