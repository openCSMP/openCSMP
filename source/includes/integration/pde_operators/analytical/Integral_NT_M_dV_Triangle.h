#ifndef Integral_NT_M_dV_Triangle_h
#define Integral_NT_M_dV_Triangle_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// to map FV solution to FEM framework (Bijective Mapping)
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_NT_M_dV_Triangle : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_NT_M_dV_Triangle( const PropertyDatabase<dim>&, const char* mapped_property );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

     virtual Integral_NT_M_dV_Triangle<dim,CELL>* clone() const override final
        { return new Integral_NT_M_dV_Triangle<dim,CELL> (*this); }

  private:
    std::vector<ScalarVariable > NPROP;
    double                       vol_div3;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















