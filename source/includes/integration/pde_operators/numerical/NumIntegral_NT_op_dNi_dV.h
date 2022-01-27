#ifndef NumIntegral_NT_op_dNi_dV_h
#define NumIntegral_NT_op_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**

@brief For the calculation of directional integrals in direction xyz;
example hydrostatic pf-gradient;
basic operands can be node or element variables.

@author S.K. Matthai
@date 2005

*/
enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

template<size_t dim,class CELL=Element<dim> >
class NumIntegral_NT_op_dNi_dV : public MathOperatorRHS<dim> {
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
    
    virtual void GetOperands( const CELL& );

    virtual void ComputeContribution( const CELL& );
    
    void SpatialDerivative( size_t xyz );
    
    virtual NumIntegral_NT_op_dNi_dV<dim,CELL>* clone() const { return new NumIntegral_NT_op_dNi_dV<dim,CELL> (*this); }
  
  private:
    std::vector<double>         IPOL;
    DenseMatrix<DM_MIN>           DN;
    Index                         mtrl_key;
    double                      oper_eprop, eprop;
    const double                gravity;   ///< acceleration of gravity (m/s2)
    size_t                        xyz;       ///< 0=x, 1=y, 2=z
    double                      mtrl_time_multiplier;
    std::vector<ScalarVariable >  oper_nprop;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















