#ifndef NUM_INTEGRAL_NT_OP_N_DS_H
#define NUM_INTEGRAL_NT_OP_N_DS_H

#include "MathOperatorRHS.h"

namespace csmp {

/**

@brief Known as: "mass matrix", "fluid sources or sinks", or "capacitance matrix"

@author S.K. Matthaei
@author S. Roberts
@date 1999

@note use only for scalar-type dependent variables (1 DOF per node)

@attention to be used to assign boundary integrals to higher-dimensional domains

*/
template<size_t dim, template<size_t> class CELL>
class NumIntegral_NT_op_N_dS : public MathOperatorRHS<dim> {
  public:
    NumIntegral_NT_op_N_dS( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* test );

    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual NumIntegral_NT_op_N_dS<dim,CELL>* clone() const { return new NumIntegral_NT_op_N_dS<dim,CELL> (*this); }

  private:
    size_t  nodal_degrees_of_freedom;

    DenseMatrix<DM_MIN>  N, NT, 
                         RHS_TEMP;
};

// copyright (c) 2000 by Stephan K. Matthai & Sebastian Geiger

} // csmp

#endif
















