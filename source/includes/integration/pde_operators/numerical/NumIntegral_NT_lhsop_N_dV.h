#ifndef NumIntegral_NT_lhsop_N_dV_h
#define NumIntegral_NT_lhsop_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;
/**
@author S.K. Matthaei
@author S. Roberts
@date 1997 */

/// Known as: capacitance matrix
template<uint32_t dim, template<uint32_t> class CELL=csmp::Element>
class NumIntegral_NT_lhsop_N_dV : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>&,
                               const char* oper, const char* basic, const char* test );
    
    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_NT_lhsop_N_dV<dim,CELL>* clone() const override { return new NumIntegral_NT_lhsop_N_dV<dim,CELL >(*this); }
    
  private:
    DenseMatrix<DM_MIN>  TEMP;
};

/* copyright (c) 1997 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















