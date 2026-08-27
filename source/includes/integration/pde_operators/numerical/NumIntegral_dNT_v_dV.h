#ifndef CSMP_NUM_INTEGRAL_DNT_V_DV_H
#define CSMP_NUM_INTEGRAL_DNT_V_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// Ffor calculation of hydrostatic and other gradient
template<uint32_t dim, template<uint32_t> class CELL>
class NumIntegral_dNT_v_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNT_v_dV( const PropertyDatabase<dim>&, 
                          const char* oper,                  // e.g., Darcy velocity
                          const char* r_factor,
                          const char* dens,                  // e.g., fluid density
                          const char* test );                // e.g., streaming potential
    
    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
  
  private:
    std::vector<double>           IPOL;
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
















