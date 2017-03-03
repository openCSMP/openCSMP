#ifndef CVFE_NUM_INTEGRAL_DNT_OP_DN_DN_H
#define CVFE_NUM_INTEGRAL_DNT_OP_DN_DN_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

template<size_t> class TwoPhaseModel;

/// PDE operator:  oper div^2 N
template<size_t dim,class SIMPLEX>
class CVFE_NumIntegral_dNT_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    CVFE_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                   TwoPhaseModel<dim>&, 
                                   const char* permeability, 
                                   const char*, 
                                   const char* );

    virtual void GetOperands( SIMPLEX& e );
    virtual void ComputeContribution( SIMPLEX& e );

  private:
    DenseMatrix<DM_MIN>  DN_, LK_, FN_; 
    TwoPhaseModel<dim>&  kri_;
};

} // csmp

#endif
