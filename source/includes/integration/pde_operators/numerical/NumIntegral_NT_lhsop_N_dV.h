#ifndef NumIntegral_NT_lhsop_N_dV_h
#define NumIntegral_NT_lhsop_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1997 */

/// Known as: capacitance matrix
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_NT_lhsop_N_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>& p, 
                               const char* oper, const char* basic, const char* test );
    
    virtual void ComputeContribution( CELL& e );
    virtual NumIntegral_NT_lhsop_N_dV<dim,CELL >* clone() const { return new NumIntegral_NT_lhsop_N_dV<dim,CELL >(*this); }
  private:
    DenseMatrix<DM_MIN>  TEMP;
};

/* copyright (c) 1997 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















