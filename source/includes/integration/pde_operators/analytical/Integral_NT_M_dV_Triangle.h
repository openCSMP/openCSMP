#ifndef Integral_NT_M_dV_Triangle_h
#define Integral_NT_M_dV_Triangle_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// to map FV solution to FEM framework (Bijective Mapping)
template<size_t dim,class SIMPLEX=Element<dim> >
class Integral_NT_M_dV_Triangle : public MathOperatorRHS<dim> {
  public:
    Integral_NT_M_dV_Triangle( const PropertyDatabase<dim>& p, const char* mapped_property );
    
    virtual void GetOperands( SIMPLEX& e );
    virtual void ComputeContribution( SIMPLEX& e );
  
  private:
    std::vector<ScalarVariable > NPROP;
    std::vector<double64>        xyz, ctr, IPOL;
    double64                     vol_div3;
};



/**
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















