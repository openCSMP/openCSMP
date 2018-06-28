#ifndef NUM_INTEGRAL_PT_OP_P_DV_H
#define NUM_INTEGRAL_PT_OP_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/**
@author S.K. Matthaei
@author S. Geiger
@date 2000 */

/// vector solution variable: equivalent of mass matrix
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_PT_op_P_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_PT_op_P_dV( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* test );

    virtual void ComputeContribution( CELL& e );

  private:
    size_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















