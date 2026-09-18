// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_OP_PT_P_DV_H
#define NUM_INTEGRAL_OP_PT_P_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**

Vector solution variable: "mass matrix" for computation of vector properties (2-3DOF) and constant coefficients

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

} // csmp

#endif
















