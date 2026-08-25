#ifndef NUM_INTEGRAL_PT_OP_DV_H
#define NUM_INTEGRAL_PT_OP_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@class NumIntegral_PT_op_dV NumIntegral_PT_op_dV "pde_operators/NumIntegral_PT_op_dV.h"

Add contributions to the righthandside of a matrix equation which arise
due to forces acting on the mass represented by each element. These 
forces are specified as nodal vector<double> variables and NumIntegral_PT_op_dV
distributes them evenly over the element.

  Vector solution variable (test): integration of 'body forces', e.g., action of gravity
  
@author S.K. Matthaei
@date 2000 

 */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_PT_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_PT_op_dV( const PropertyDatabase<dim>&, const char* oper, const char* test );
    ~NumIntegral_PT_op_dV() = default;
    
    void GetOperands( const CELL<dim>& ) override final;
 
    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_PT_op_dV<dim,CELL>* clone() const override final { return new NumIntegral_PT_op_dV<dim,CELL> (*this); }
  private:
    std::vector<double>  BFORCE;
};


} // csmp

#endif
