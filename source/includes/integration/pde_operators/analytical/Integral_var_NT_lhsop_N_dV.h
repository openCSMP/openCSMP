#ifndef Integral_var_NT_lhsop_N_dV_h
#define Integral_var_NT_lhsop_N_dV_h

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// known as mass or capacitance matrix
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_var_NT_lhsop_N_dV : public MathOperatorLHS<dim,CELL> {
  public:
    Integral_var_NT_lhsop_N_dV(const PropertyDatabase<dim>&,
                               const char* oper,
                               const char* basic,
                               const char* test,
                               const char* var,
                               const double prefactor = 1. );
    
    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_var_NT_lhsop_N_dV<dim,CELL>* clone() const override final { return new Integral_var_NT_lhsop_N_dV<dim,CELL> (*this); }
    
  private:
    void ComputeIntegral( const CELL<dim>& );
    
    ScalarVariable                op_;
    Parameter                     var_;
    std::vector<ScalarVariable >  vvar_;
    const double                  prefactor_;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















