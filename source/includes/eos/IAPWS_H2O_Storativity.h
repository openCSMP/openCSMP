#ifndef IAPWS_H2O_STORATIVITY_H
#define IAPWS_H2O_STORATIVITY_H

#include "CSMP_definitions.h"
#include "Interrelation.h"
#include "thirdparty/iaps_h2o_eos/steam4.h"

namespace csmp {

template<size_t dim>
class IAPWS_H2O_Storativity : public Interrelation<dim> {
  public:
    IAPWS_H2O_Storativity( const PropertyDatabase<dim>& p, double64 z_rock=1.0e-12 );
    ~IAPWS_H2O_Storativity() {};
    void Calculate();
  
  private:
    Operand<dim>&  B;  // fluid compressibility
    Operand<dim>&  X;  // porosity
    Operand<dim>&  S;  // storativity
    ScalarVariable beta, phi;
    double64       Z, Z_rock;
};

} // end namespace csmp

#endif

