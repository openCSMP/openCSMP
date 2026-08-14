#ifndef NUM_INTEGRAL_OP_PT_P_DV_H
#define NUM_INTEGRAL_OP_PT_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**
vector solution variable: "mass matrix" for computation of vector properties (2-3DOF) and constant coefficients

@author S.K. Matthai
@author S. Geiger
@date 2000

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_op_PT_P_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_op_PT_P_dV( const PropertyDatabase<dim>&, 
                            const char* oper, const char* test );

    void ComputeContribution( const CELL<dim>& ) override final;

  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















