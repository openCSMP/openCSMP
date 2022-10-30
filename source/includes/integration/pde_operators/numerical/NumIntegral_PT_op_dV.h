#ifndef NUM_INTEGRAL_PT_OP_DV_H
#define NUM_INTEGRAL_PT_OP_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/** vector solution variable: integration of 'body forces', e.g., action of gravity
 */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_PT_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_PT_op_dV( const PropertyDatabase<dim>&, const char* oper, const char* test );
    virtual ~NumIntegral_PT_op_dV() {}
    
    virtual void GetOperands( const CELL<dim>& );
 
    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual NumIntegral_PT_op_dV<dim,CELL>* clone() const { return new NumIntegral_PT_op_dV<dim,CELL> (*this); }
  private:
    std::vector<double>  BFORCE;
};

/**
 
@class NumIntegral_PT_op_dV NumIntegral_PT_op_dV "pde_operators/NumIntegral_PT_op_dV.h"

@author S.K. Matthaei
@date 2000 

Add contributions to the righthandside of a matrix equation which arise
due to forces acting on the mass represented by each element. These 
forces are specified as nodal vector<double> variables and NumIntegral_PT_op_dV
distributes them evenly over the element. 
 
 */

} // csmp

#endif
