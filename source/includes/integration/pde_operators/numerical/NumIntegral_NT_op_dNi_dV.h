#ifndef NumIntegral_NT_op_dNi_dV_h
#define NumIntegral_NT_op_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**

@brief For the calculation of directional integrals in direction xyz;
example hydrostatic pf-gradient;
conductivity can be node or element property but is not computed
from basic variables. 

@author S.K. Matthai
@date 2005

*/
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_NT_op_dNi_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>& p, 
                              const char* oper,    // e.g., fluid density
                              const char* test,    // e.g., fluid pressure
                              double64 acc_gravity=9.8601 ); // m/s2

    NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>& p, 
                              const char* oper,    // e.g., fluid density
                              const char* mtrl,    // e.g., conductivity
                              const char* test,    // e.g., fluid pressure
                              double64 acc_gravity=9.8601 );
    
    virtual void GetOperands( SIMPLEX& e );

    virtual void ComputeContribution( SIMPLEX& e );
    
    void SpatialDerivative( size_t xyz );
    virtual NumIntegral_NT_op_dNi_dV<dim,SIMPLEX>* clone() const { return new NumIntegral_NT_op_dNi_dV<dim,SIMPLEX> (*this); }
  
  private:
    std::vector<double64>         IPOL;
    DenseMatrix<DM_MIN>           DN;
    Index                         mtrl_key;
    double64                      oper_eprop, eprop;
    const double64                gravity;   ///< acceleration of gravity (m/s2)
    size_t                        xyz;       ///< 0=x, 1=y, 2=z
    double64                      mtrl_time_multiplier;
    std::vector<ScalarVariable >  oper_nprop;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















