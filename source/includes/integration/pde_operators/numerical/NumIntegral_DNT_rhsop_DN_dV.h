#ifndef NumIntegral_DNT_rhsop_DN_dV_h
#define NumIntegral_DNT_rhsop_DN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**

@brief RHS = div^2 . operand - interpolation function derivative matrix collapsed 
into righthand side vector.

@author S.K. Matthai
@author S. Roberts
@date 1999 

@copyright 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts

*/
template<uint32_t dim,class CELL=Element<dim> >
class NumIntegral_DNT_rhsop_DN_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& p, 
                                 const char* oper,          
                                 const char* test );       

    NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& p, 
                                 const char* integral_multiplier,
                                 const char* oper,          
                                 const char* test );       
    
    virtual void GetOperands( const CELL& );

    virtual void ComputeContribution( const CELL& );
    
    void IgnoreOperand( bool ignore );
    virtual NumIntegral_DNT_rhsop_DN_dV<dim,CELL>* clone() const
      { return new NumIntegral_DNT_rhsop_DN_dV<dim,CELL> (*this); }
      
  private:
    DenseMatrix<DM_MIN>        DN, DNT, OPMAT;
    std::vector<ScalarVariable >  noperand;
    ScalarVariable                eoperand, multiplier;
    bool                              ignore_operand;
    csmp::Index                         mult_key;
};


} // csmp

#endif
















