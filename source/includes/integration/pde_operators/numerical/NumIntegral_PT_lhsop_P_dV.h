#ifndef NUM_INTEGRAL_PT_LHSOP_P_DV_H
#define NUM_INTEGRAL_PT_LHSOP_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"


namespace csmp {
/**
@author S.K. Matthaei
@author S. Geiger
@date 2000 */

/// vector solution variable: "mass matrix" for vector dependent variables
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_PT_lhsop_P_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_PT_lhsop_P_dV( const PropertyDatabase<dim>& p, 
                               const char* oper, const char* oper2, 
                               const char* basic, const char* test );
    
    virtual void GetOperands( SIMPLEX& e );
    virtual void ComputeContribution( SIMPLEX& e );

  private:
    int                 nodal_degrees_of_freedom;
    csmp::Index         phi_key;
    ScalarVariable      phi;
    DenseMatrix<DM_MIN> PT, P;
    
    void N_to_P( const std::vector<double64>& N, DenseMatrix<DM_MIN>& P );
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp


#endif
















