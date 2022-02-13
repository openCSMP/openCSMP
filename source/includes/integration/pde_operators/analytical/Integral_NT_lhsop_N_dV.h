#ifndef INTEGRAL_NT_LHSOP_N_DV_H
#define INTEGRAL_NT_LHSOP_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {
  /**
  @author S.K. Matthaei
  @author S. Roberts
  @date 1999 */

/// known as mass or capacitance matrix
template<uint32_t dim,class CELL=Element<dim> >
class Integral_NT_lhsop_N_dV : public MathOperatorLHS<dim> {
  public:
    Integral_NT_lhsop_N_dV( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* basic, const char* test );

    virtual void GetOperands( const CELL& e );

    virtual void ComputeContribution( const CELL& e );
    
    virtual Integral_NT_lhsop_N_dV<dim,CELL>* clone() const { return new Integral_NT_lhsop_N_dV<dim,CELL> (*this); }

  private:
    double scalar_value_;
};



/**

copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















