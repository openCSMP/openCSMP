#ifndef NUM_INTEGRAL_PT_OP_DV_H
#define NUM_INTEGRAL_PT_OP_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    @brief Body force RHS (vector test function, volume integral).

    f_j = ∫_Ω Pⱼᵀ [f] dV

    where P is the vector interpolation matrix and [f] is a body force vector (e.g. ρg for gravitational loading in mechanics).

    Operand: Vector — element-placed.

    Test variable: Vector (displacement), node-placed.

    Application: Gravitational body force in linear elasticity and geomechanics.

    Vector solution variable (test): integration of 'body forces', e.g., action of gravity
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
