#ifndef NumIntegral_op_NT_dN_orthogonal_dV_h
#define NumIntegral_op_NT_dN_orthogonal_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**
@brief Stream function orthogonal gradient RHS

f_k = -∫_Ω Nₖ φ (∇⊥N)ᵀ dV

where ∇⊥ denotes the 90°-rotated gradient operator:

∇⊥N = ( -∂N/∂y,  ∂N/∂x )

applied component-wise to the shape functions, and φ is the scalar stream function interpolated from nodal values φⱼ.
The negative sign reflects the convention that flow is directed down the pressure gradient.

The rotation transforms the stream function gradient into the orthogonal velocity field: if ψ is the stream function then

u = ∂ψ/∂y,   v = -∂ψ/∂x

so the velocity vector is perpendicular to the isolines of ψ.

Operand: Scalar stream function — node-placed. Nodal values are read directly via Node::Read into NPROP.

Test variable: Scalar, node-placed.

Spatial dimension: 2D only. The 90° rotation is only defined in 2D — the operator should not be used for dim == 3.

Key distinction from other gradient operators:

NumIntegral_dNT_op_dV contracts (∇N)ᵀ with a pre-computed element body force vector
NumIntegral_dNi_rhsop_dV computes the gradient of a nodal scalar in one coordinate direction
This operator rotates the gradient by 90° before contracting,
converting stream function gradients into orthogonal velocity contributions

Application:
- Stream function formulation of 2D incompressible flow
- Post-processing of stream function fields to recover velocity
- Linear elasticity computations involving rotated gradient fields

@note see streamfunction example

 */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_op_NT_dN_orthogonal_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_op_NT_dN_orthogonal_dV( const PropertyDatabase<dim>&,
                                        const char* basic,            // e.g., fluid pressure
                                        const char* test );           // streamfunction
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>* clone() const override final { return new NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  M, DNORTHO, NT; 
    std::vector<double>  NPROP, IPOL, UNITY, RES;
};

} // csmp

#endif

























