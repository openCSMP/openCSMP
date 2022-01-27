#ifndef NUM_INTEGRAL_OP_NT_N_DV_H
#define NUM_INTEGRAL_OP_NT_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/**
@author S.K. Matthaei
@author S. Geiger
@date 2000 */

/// "mass matrix", "fluid sources or sinks", or "capacitance matrix" for constant coefficients
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_op_NT_N_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_op_NT_N_dV( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* test );

    virtual void ComputeContribution( const CELL& );
    
  private:
    size_t  nodal_degrees_of_freedom;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger
} // csmp

#endif
















