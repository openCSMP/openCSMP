#ifndef NumIntegral_BT_op_dV_h
#define NumIntegral_BT_op_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    @brief Initial stress / pore pressure body force RHS

    f_j = ∫_Ω Bⱼᵀ [σ₀] dV

    where [σ₀] is a known initial stress or pore-pressure-induced stress tensor. B is the strain-displacement matrix.
    This integral adds the equivalent nodal forces produced by a pre-existing stress state or fluid pressure field.

    Operand: Scalar (isotropic, e.g. fluid pressure) or tensor (anisotropic initial stress) — element or integration point-placed.

    Test variable: Vector (displacement), node-placed.

    Application:

    Pore pressure contribution to effective stress in poromechanics (Biot coupling term on the RHS)
    Initial stress loading in geomechanics.
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_BT_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_BT_op_dV( const PropertyDatabase<dim>&,
                          const char* oper, const char* test );

    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_BT_op_dV<dim,CELL>* clone() const override { return new NumIntegral_BT_op_dV<dim,CELL> (*this); }
    
  private:  
    DenseMatrix<DM_MIN>  B, BT, STR; 
};

} // csmp

#endif
