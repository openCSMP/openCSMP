#ifndef IAPWS_H2O_CONDUCTIVITY_H
#define IAPWS_H2O_CONDUCTIVITY_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class IAPWS_H2O_Conductivity : public Interrelation<dim> {
    Operand<dim>&  k;     
    Operand<dim>&  K;  
    Operand<dim>&  v;  
    ScalarVariable  perm, viscos;

  public: // 					                   
    IAPWS_H2O_Conductivity( const PropertyDatabase<dim>& p );
    ~IAPWS_H2O_Conductivity() {};
    void Calculate();
};

}

#endif

