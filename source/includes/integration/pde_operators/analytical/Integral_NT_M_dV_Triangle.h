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
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );
  
  private:
    std::vector<ScalarVariable > NPROP;
    double                       vol_div3;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















