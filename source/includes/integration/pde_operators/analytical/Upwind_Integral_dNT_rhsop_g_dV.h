#ifndef Upwind_Integral_dNT_rhsop_g_dV_h
#define Upwind_Integral_dNT_rhsop_g_dV_h

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// buoyancy for instance
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Upwind_Integral_dNT_rhsop_g_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Upwind_Integral_dNT_rhsop_g_dV(const PropertyDatabase<dim>& p, 
                                   const char* oper,
                                   const char* test,
                                   const char* upwind,
                                   const char* trigger,
                                   const double prefactor = 1. );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    void SpatialDerivative( uint32_t xyz=2 );
  
  private:
    DenseMatrix<DM_MIN>  DN, DNT, coords;
    const double         gravity;   // acceleration of gravity
    uint32_t             xyz;       // 1=x, 2=y, 3=z
    
    Parameter upwind_;
    Parameter trigger_;
    std::vector<ScalarVariable > upwind_var_;
    std::vector<ScalarVariable > trigger_var_;
    const double prefactor_;
};



/**
 
@class Upwind_Integral_dNT_rhsop_g_dV
 
@author S.K. Matthaei
@author S. Roberts
@date 1999
 */


} // csmp

#endif
















