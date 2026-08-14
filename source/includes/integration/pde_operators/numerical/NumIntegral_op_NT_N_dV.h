#ifndef NUM_INTEGRAL_OP_NT_N_DV_H
#define NUM_INTEGRAL_OP_NT_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**  "mass matrix", "fluid sources or sinks", or "capacitance matrix" for constant coefficients

@author S.K. Matthai
@author S. Geiger
@date 2000

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_op_NT_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_op_NT_N_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* test );

    void GetOperands( const CELL<dim>& ) override final { /* do not load any data */ }
    void ComputeContribution( const CELL<dim>& ) override final;
    
  private:
    uint32_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger
} // csmp

#endif
















