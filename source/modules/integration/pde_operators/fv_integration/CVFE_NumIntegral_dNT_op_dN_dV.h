#ifndef CVFE_NUM_INTEGRAL_DNT_OP_DN_DN_H
#define CVFE_NUM_INTEGRAL_DNT_OP_DN_DN_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class TwoPhaseModel;

/// PDE operator:  oper div^2 N
template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFE_NumIntegral_dNT_op_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    CVFE_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                   TwoPhaseModel<dim>&, 
                                   const char* permeability, 
                                   const char*, 
                                   const char* );

    virtual void GetOperands( CELL<dim>& );
    virtual void ComputeContribution( CELL<dim>& );

  private:
    DenseMatrix<DM_MIN>  DN_, LK_, FN_; 
    TwoPhaseModel<dim>&  kri_;
};

} // csmp

#endif
