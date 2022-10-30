#ifndef NumIntegral_NT_mixed_op_dNi_dV_h
#define NumIntegral_NT_mixed_op_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// for instance for calculation of hydrostatic gradient
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_mixed_op_dNi_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_mixed_op_dNi_dV( const PropertyDatabase<dim>&, 
                              const char* nodal_mtrl_multiplier, // e.g., viscosity
                              const char* oper,                  // e.g., fluid density
                              const char* mtrl,                  // e.g., permeability
                              const char* test );                // e.g., fluid pressure
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );
    
    void SpatialDerivative( uint32_t xyz=2 );
    
    /// do not use MultiplyWithTimeIncrement() here since this would multiply
    /// the whole contribution
    void MaterialPropertyTimeMultiplier( double time_increment );
  
  private:
    std::vector<double>           IPOL, DNI;
    DenseMatrix<DM_MIN>           DN;
    Index                         mtrl_key, nmult_key;
    ScalarVariable                oper_eprop, eprop;
    const double                  gravity;   // acceleration of gravity
    uint32_t                      xyz; // 1=x, 2=y, 3=z
    double                        mtrl_time_multiplier;
    std::vector<ScalarVariable >  oper_nprop, mtrl_nprop;
};



/**
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

}


#endif
















