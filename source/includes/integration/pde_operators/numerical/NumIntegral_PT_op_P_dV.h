#ifndef NUM_INTEGRAL_PT_OP_P_DV_H
#define NUM_INTEGRAL_PT_OP_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**

Vector solution variable: equivalent of mass matrix

@author S.K. Matthaei
@author S. Geiger
@date 2000

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_PT_op_P_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_PT_op_P_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* test );

    void ComputeContribution( const CELL<dim>& ) override final;

  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















