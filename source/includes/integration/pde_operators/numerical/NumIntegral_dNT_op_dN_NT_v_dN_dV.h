#ifndef CSMP_NUM_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_H
#define CSMP_NUM_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_H

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
    Advection-dispersion matrix (FEM integral), see Istok (1989), suitable for steady-state solutions for small Peclet number flows
    
    K_jk = ∫_Ω (∇N_j)ᵀ [σ] ∇N_k dV  +  ∫_Ω N_j (v · ∇N_k) dV
    
    The first term is the diffusive contribution (identical to NumIntegral_dNT_lhsop_dN_dV).
    The second term is the advective contribution weighted by the velocity field v.

    Operands: Diffusivity tensor [σ] and advection velocity v — both element or integration point-placed.
    Test variable: Scalar, node-placed.
    Application: Advection-diffusion transport of scalar quantities (concentration, temperature, tracer).
 */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_op_dN_NT_v_dN_dV final : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>&, 
                                      const char* diffusion_oper,   ///< element prop, for instance thermal conductivity
                                      const char* advection_oper,   ///< element prop, for instance heat transport velocity
                                      const char* basic,            ///< e.g., fluid pressure
                                      const char* test );           ///< e.g., fluid pressure
                            
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

    NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>(*this); }

  private:
    DenseMatrix<DM_MIN>  DN, DNT, NT3;
    DenseMatrix<DM12>    VIP;
    std::vector<double>  IPOL;    ///< basis function values (at integration point)
    VectorVariable<dim>  velo_;   ///< Darcy flow velocity
    csmp::Index          adv_key; ///< index of the variable that shall be advected
};

} // csmp

#endif /* CSMP_NUM_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_H */

























