#ifndef NUM_INTEGRAL_OP_PT_P_DV_H
#define NUM_INTEGRAL_OP_PT_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Geiger
@date 2000 */

/// vector solution variable: "mass matrix" for computation of vector properties (2-3DOF) and constant coefficients
template<uint32_t dim,class CELL=Element<dim> >
class NumIntegral_op_PT_P_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_op_PT_P_dV( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* test );
    
    virtual void ComputeContribution( const CELL& );

  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















