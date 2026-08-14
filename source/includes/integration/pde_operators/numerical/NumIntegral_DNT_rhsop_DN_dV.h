#ifndef NumIntegral_dNT_rhsop_dN_dV_h
#define NumIntegral_dNT_rhsop_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief RHS = div^2 . operand - interpolation function derivative matrix collapsed 
into righthand side vector.

@author S.K. Matthai
@author S. Roberts
@date 1999 

@copyright 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_rhsop_dN_dV final : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& p, 
                                 const char* oper,          
                                 const char* test );       

    NumIntegral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& p, 
                                 const char* integral_multiplier,
                                 const char* oper,          
                                 const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
  
    NumIntegral_dNT_rhsop_dN_dV<dim,CELL>* clone() const override final
      { return new NumIntegral_dNT_rhsop_dN_dV<dim,CELL> (*this); }
      
  private:
    DenseMatrix<DM_MIN>           DN, DNT, OPMAT, TEMP;
    std::vector<ScalarVariable >  noperand;
    ScalarVariable                eoperand, multiplier;
    bool                          has_multiplier_;
    csmp::Index                   mult_key;
};


} // csmp

#endif
















