#ifndef NumIntegral_NT_op_N_dV_h
#define NumIntegral_NT_op_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/**

@brief Known as: "mass matrix", "fluid sources or sinks", or "capacitance matrix"

@author S.K. Matthaei
@author S. Roberts
@date 1999

@note use only for scalar-type dependent variables (1 DOF per node)

*/
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_NT_op_N_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_NT_op_N_dV( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* test );

    virtual void ComputeContribution( SIMPLEX& e );
    virtual NumIntegral_NT_op_N_dV<dim,SIMPLEX>* clone() const { return new NumIntegral_NT_op_N_dV<dim,SIMPLEX> (*this); }
  private:
    size_t  nodal_degrees_of_freedom;

    DenseMatrix<DM_MIN>  N, NT, 
                         RHS_TEMP;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















