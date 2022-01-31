#ifndef Upwind_Integral_dNT_rhsop_g_dV_h
#define Upwind_Integral_dNT_rhsop_g_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/// buoyancy for instance
template<size_t dim,class SIMPLEX=Element<dim> >
class Upwind_Integral_dNT_rhsop_g_dV : public MathOperatorRHS<dim> {
  public:
    Upwind_Integral_dNT_rhsop_g_dV(const PropertyDatabase<dim>& p, 
                                   const char* oper,
                                   const char* test,
                                   const char* upwind,
                                   const char* trigger,
                                   const double prefactor = 1. );
    
    virtual void GetOperands( const SIMPLEX& );
    /// integration etc.
    virtual void ComputeContribution( const SIMPLEX& );
    
    void SpatialDerivative( size_t xyz=2 );
  
  private:
    DenseMatrix<DM_MIN>  DN, DNT, coords;
    double                          gravity;   // acceleration of gravity
    size_t                   xyz;       // 1=x, 2=y, 3=z
    
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
















