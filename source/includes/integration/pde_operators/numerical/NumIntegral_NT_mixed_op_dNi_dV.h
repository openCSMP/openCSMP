#ifndef NumIntegral_NT_mixed_op_dNi_dV_h
#define NumIntegral_NT_mixed_op_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// for instance for calculation of hydrostatic gradient
template<size_t dim,class SIMPLEX>
class NumIntegral_NT_mixed_op_dNi_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_NT_mixed_op_dNi_dV( const PropertyDatabase<dim>& p, 
                              const char* nodal_mtrl_multiplier, // e.g., viscosity
                              const char* oper,                  // e.g., fluid density
                              const char* mtrl,                  // e.g., permeability
                              const char* test );                // e.g., fluid pressure
    
    virtual void GetOperands( SIMPLEX& e );
    virtual void ComputeContribution( SIMPLEX& e );
    
    void SpatialDerivative( size_t xyz=2 );
    
    /// do not use MultiplyWithTimeIncrement() here since this would multiply
    /// the whole contribution
    void MaterialPropertyTimeMultiplier( double64 time_increment );
  
  private:
    std::vector<double64>         IPOL, DNI;
    DenseMatrix<DM_MIN>           DN;
    Index                         mtrl_key, nmult_key;
    ScalarVariable                oper_eprop, eprop;
    const double64                gravity;   // acceleration of gravity
    size_t                        xyz; // 1=x, 2=y, 3=z
    double64                      mtrl_time_multiplier;
    std::vector<ScalarVariable >  oper_nprop, mtrl_nprop;
};



/**
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

}


#endif
















