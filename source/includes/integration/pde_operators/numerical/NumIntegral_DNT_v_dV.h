#ifndef NumIntegral_DNT_v_dV_h
#define NumIntegral_DNT_v_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// for instance for calculation of hydrostatic gradient
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_DNT_v_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_DNT_v_dV( const PropertyDatabase<dim>& p, 
                          const char* oper,                  // e.g., Darcy velocity
                          const char* r_factor,
                          const char* dens,                  // e.g., fluid density
                          const char* test );                // e.g., streaming potential
    
    virtual void GetOperands( CELL& e );

    virtual void ComputeContribution( CELL& e );
  
  private:
    std::vector<double64>         IPOL;
    DenseMatrix<DM_MIN>           DN, DNT, OPMAT;
    Index                         rho_key,   rfac_key;
    ScalarVariable                rfac;
    VectorVariable<dim>           oper_eprop;
    std::vector<ScalarVariable >  oper_nprop, dens_nprop;
};



/**

copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















