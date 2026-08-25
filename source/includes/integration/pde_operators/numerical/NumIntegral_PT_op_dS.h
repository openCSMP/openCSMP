#ifndef NUM_INTEGRAL_PT_OP_DS_H
#define NUM_INTEGRAL_PT_OP_DS_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Face;
  
/**
   @brief Neumann traction RHS (surface integral)
   
    f_j = ∫_Γ Pⱼᵀ t dS  ≈  Σᵢ wᵢ |Jᵢ| N_j(ξᵢ) · t
   
   where t is the applied traction vector (or scalar pressure) on the boundary face Γ, wᵢ are quadrature weights, and |Jᵢ| is the surface Jacobian determinant.

   Operand: Vector traction or scalar pressure — face-placed, flagged NEUMANN.
   Only faces where the operand carries the NEUMANN status flag are accumulated; all others are skipped.

   Test variable: Vector (displacement), node-placed.

  Application:
  - Far-field in-situ stress conditions on outer model boundaries
  - Mud pressure on the borehole wall
  - Any distributed surface traction in mechanics
*/
template<uint32_t dim>
class NumIntegral_PT_op_dS : public MathOperatorRHS<dim,Face> {
  public:
    NumIntegral_PT_op_dS( const PropertyDatabase<dim>& p, 
                          const char* oper,    // VECTOR/SCALAR variable on FACE   
                          const char* test );  // VECTOR variable on NODE
    
    void GetOperands( const Face<dim>& ) override final;
    void ComputeContribution( const Face<dim>& ) override final;
  
    ::csmp::NumIntegral_PT_op_dS<dim>* clone() const override { return new NumIntegral_PT_op_dS<dim> (*this); }

  private:
    VectorVariable<dim>  oper_;            ///< Face variable value
    ScalarVariable       scalar_;          ///< Face scalar value
}; 


} // csmp

#endif
















