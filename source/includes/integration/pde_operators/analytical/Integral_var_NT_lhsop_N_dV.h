#ifndef Integral_var_NT_lhsop_N_dV_h
#define Integral_var_NT_lhsop_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// known as mass or capacitance matrix
template<uint32_t dim,class CELL=Element<dim> >
class Integral_var_NT_lhsop_N_dV : public MathOperatorLHS<dim> {
  public:
    Integral_var_NT_lhsop_N_dV(const PropertyDatabase<dim>& p, 
                               const char* oper,
                               const char* basic,
                               const char* test,
                               const char* var,
                               const double prefactor = 1. );
    
    virtual void GetOperands( const CELL& );

    virtual void ComputeContribution( const CELL& );
    
    virtual Integral_var_NT_lhsop_N_dV<dim,CELL>* clone() const { return new Integral_var_NT_lhsop_N_dV<dim,CELL> (*this); }
  private:
    void ComputeIntegral( const CELL& );
    
    ScalarVariable                op_;
    Parameter                     var_;
    std::vector<ScalarVariable >  vvar_;
    const double                prefactor_;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















