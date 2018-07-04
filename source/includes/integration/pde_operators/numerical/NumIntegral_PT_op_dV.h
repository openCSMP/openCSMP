#ifndef NUM_INTEGRAL_PT_OP_DV_H
#define NUM_INTEGRAL_PT_OP_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

/// vector solution variable: integration of 'body forces', e.g., action of gravity
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_PT_op_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_PT_op_dV( const PropertyDatabase<dim>& pref, const char* oper, const char* test );
    
    virtual void GetOperands( CELL& e );
 
    virtual void ComputeContribution( CELL& e );
    virtual NumIntegral_PT_op_dV<dim,CELL>* clone() const { return new NumIntegral_PT_op_dV<dim,CELL> (*this); }
  private:
    std::vector<double64>  BFORCE;
};

/**
 
@class NumIntegral_PT_op_dV NumIntegral_PT_op_dV "pde_operators/NumIntegral_PT_op_dV.h"

@author S.K. Matthaei
@date 2000 

Add contributions to the righthandside of a matrix equation which arise
due to forces acting on the mass represented by each element. These 
forces are specified as nodal vector<double64> variables and NumIntegral_PT_op_dV
distributes them evenly over the element. 
 
 */

} // csmp

#endif
