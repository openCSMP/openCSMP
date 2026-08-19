#ifndef NumIntegral_NT_op_dNi_dV_h
#define NumIntegral_NT_op_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief For the calculation of directional integrals in direction xyz;
example hydrostatic pf-gradient;
basic operands can be node or element variables.

@attention operator used in TopographyDrivenFlow_Example

@author S.K. Matthai
@date 2005

*/
enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_op_dNi_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>&,
                              const char* oper,    // e.g., fluid density
                              const char* test,    // e.g., fluid pressure
                              double acc_gravity=9.8601 ); // m/s2

    NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>&,
                              const char* oper,    // e.g., fluid density
                              const char* mtrl,    // e.g., conductivity
                              const char* test,    // e.g., fluid pressure
                              double acc_gravity=9.8601 );
    
    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
    
    void SpatialDerivative( uint32_t xyz );
    
    NumIntegral_NT_op_dNi_dV<dim,CELL>* clone() const override final { return new NumIntegral_NT_op_dNi_dV<dim,CELL> (*this); }
  
  private:
    std::vector<double>           IPOL;
    DenseMatrix<DM_MIN>           DN;
    Index                         mtrl_key;
    double                        oper_eprop, eprop;
    const double                  gravity;   ///< acceleration of gravity (m/s2)
    uint32_t                      xyz;       ///< 0=x, 1=y, 2=z
    double                        mtrl_time_multiplier;
    std::vector<ScalarVariable >  oper_nprop;
};

} // csmp

#endif
















