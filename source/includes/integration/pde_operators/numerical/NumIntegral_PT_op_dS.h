#ifndef NUM_INTEGRAL_PT_OP_DS_H
#define NUM_INTEGRAL_PT_OP_DS_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {
  
/** Scalars, or vectors acting on model boundary surfaces as tractions/stresses

@date 2000-2013
@author Dr. Stephan K. Matthai
@author Stephen G. Roberts
@author M. Nejati
@author P.S. Lang

Operator to be applied to boundary, with NEUMANN variables considered only. Scalars are considered positiv if directed
towards inside of boundary.

*/
template<size_t dim>
class NumIntegral_PT_op_dS : public MathOperatorRHS<dim> {
  public:
    NumIntegral_PT_op_dS( const PropertyDatabase<dim>& p, 
                          const char* oper,    // VECTOR/SCALAR variable on FACE   
                          const char* test );  // VECTOR variable on NODE
    
    virtual void GetOperands( Face<dim>& f );
    virtual void ComputeContribution( Face<dim>& f );
  
    virtual ::csmp::NumIntegral_PT_op_dS<dim>* clone() const { return new NumIntegral_PT_op_dS<dim> (*this); }

  private:
    VectorVariable<dim>  oper_;            ///< Face variable value
    ScalarVariable       scalar_;          ///< Face scalar value
}; 


} // csmp

#endif
















