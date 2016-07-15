#ifndef Integral_dNT_op_dN_NT_v_dN_dV_h
#define Integral_dNT_op_dN_NT_v_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// Diffusion - advection lefthandside PDE operator; works for Peclet numbers <1
template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_dNT_op_dN_NT_v_dN_dV : public MathOperatorLHS<dim> {
  public:
    Integral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref, 
                                    const char* oper,  const char* velo, 
                                    const char* basic, const char* test );
    
    virtual void GetOperands( SIMPLEX& e );
    
    virtual void ComputeContribution( SIMPLEX& e );
  
  private:
    DenseMatrix<DM_MIN>  B, BT;
    csmp::Index          velo_key;
    VectorVariable<dim>  VXYZ;
};

/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
