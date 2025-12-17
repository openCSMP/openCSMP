#ifndef NumIntegral_NT_op_N_dV_h
#define NumIntegral_NT_op_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief Known as: "mass matrix", "fluid sources or sinks", or "capacitance matrix"

@author S.K. Matthaei
@author S. Roberts
@date 1999

@note use only for scalar-type dependent variables (1 DOF per node)

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_op_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_op_N_dV( const PropertyDatabase<dim>&,
                            const char* oper, const char* test );

    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_NT_op_N_dV<dim,CELL>* clone() const override final { return new NumIntegral_NT_op_N_dV<dim,CELL> (*this); }
    
  private:
    uint32_t  nodal_degrees_of_freedom;

    DenseMatrix<DM_MIN>  N_, NT_,
                         RHS_TEMP_;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















