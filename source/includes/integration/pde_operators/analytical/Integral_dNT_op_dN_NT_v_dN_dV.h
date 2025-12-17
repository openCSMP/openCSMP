#ifndef Integral_dNT_op_dN_NT_v_dN_dV_h
#define Integral_dNT_op_dN_NT_v_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// Diffusion - advection lefthandside PDE operator; works for Peclet numbers <1
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_op_dN_NT_v_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    Integral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref, 
                                    const char* oper,  const char* velo, 
                                    const char* basic, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
  
  private:
    DenseMatrix<DM_MIN>  B, BT;
    csmp::Index          velo_key;
    VectorVariable<dim>  VXYZ;
};

/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
