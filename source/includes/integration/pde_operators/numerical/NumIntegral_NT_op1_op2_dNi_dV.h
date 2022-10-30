#ifndef NumIntegral_NT_op1_op2_dNi_dV_h
#define NumIntegral_NT_op1_op2_dNi_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@date 2005

gravity term in transient flow: -S / dt  +  K     g   delta_rho grad Z
                              mtrl1 * dt   mtrl2  oper      dN (pf)

@note only use in BE scheme where storage term is divided by time-increment
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_op1_op2_dNi_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_op1_op2_dNi_dV( const PropertyDatabase<dim>&,
                                   const char* oper,          // e.g., fluid density
                                   const char* mtrl1,         // e.g., storativity
                                   const char* mtrl2,         // e.g., conductivity
                                   const char* test );        // e.g., fluid pressure
    
    virtual void GetOperands( const CELL<dim>& );

    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual void MultiplyWithTimeFactor( double dt );
  
  private:
    std::vector<double>           IPOL;
    DenseMatrix<DM_MIN>           DN;
    csmp::Index                   mtrl1_key, mtrl2_key;
    ScalarVariable                oper_eprop, mtrl1_prop, mtrl2_prop;
    const double                  gravity;   // acceleration of gravity
    const uint32_t                xyz; // 1=x, 2=y
    std::vector<ScalarVariable >  oper_nprop;
};



} // csmp

#endif
















