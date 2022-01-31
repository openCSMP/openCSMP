#ifndef NUM_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_H
#define NUM_INTEGRAL_DNT_OP_DN_NT_V_DN_DV_H

#include "MathOperatorLHS.h"
#include "VectorVariable.h"

namespace csmp {

template<size_t> class PropertyDatabase;

/// advection-dispersion integral, see Istok (1989), use to obtain steady-state solutions for small Peclet number flows
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_DNT_op_DN_NT_v_DN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_DNT_op_DN_NT_v_DN_dV( const PropertyDatabase<dim>& pref, 
                                      const char* diffusion_oper,   ///< element prop, for instance thermal conductivity
                                      const char* advection_oper,   ///< element prop, for instance heat transport velocity
                                      const char* basic,            ///< e.g., fluid pressure
                                      const char* test );           ///< e.g., fluid pressure
                            
    virtual void GetOperands( const CELL& );
    virtual void ComputeContribution( const CELL& );

  private:
    DenseMatrix<DM_MIN>    DN, DNT, VIP, NT3; 
    std::vector<double>  IPOL;    ///< basis function values (at integration point)
    VectorVariable<dim>    velo_;   ///< Darcy flow velocity
    csmp::Index            adv_key; ///< index of the variable that shall be advected
};

} // csmp

#endif

























