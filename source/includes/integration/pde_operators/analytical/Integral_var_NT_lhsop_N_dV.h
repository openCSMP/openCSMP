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
template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_var_NT_lhsop_N_dV : public MathOperatorLHS<dim> {
  public:
    Integral_var_NT_lhsop_N_dV(const PropertyDatabase<dim>& p, 
                               const char* oper,
                               const char* basic,
                               const char* test,
                               const char* var,
                               const double prefactor = 1.0);
    
    virtual void GetOperands( SIMPLEX& e );

    virtual void ComputeContribution( SIMPLEX& e );
    virtual Integral_var_NT_lhsop_N_dV<dim,SIMPLEX>* clone() const { return new Integral_var_NT_lhsop_N_dV<dim,SIMPLEX> (*this); }
  private:
    void ComputeIntegral( SIMPLEX& e );
    
    ScalarVariable                op_;
    Parameter                     var_;
    std::vector<ScalarVariable >  vvar_;
    const double                prefactor_;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















