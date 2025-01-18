#ifndef IAPWS_H2O_HEAT_TRANSFER_VELOCITY
#define IAPWS_H2O_HEAT_TRANSFER_VELOCITY

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class IAPWS_H2O_HeatTransferVelocity : public Interrelation<dim> {
    Operand<dim>&  V;     
    Operand<dim>&  VH;  
    Operand<dim>&  P;  
    Operand<dim>&  CPF;  
    Operand<dim>&  RF;  
    Operand<dim>&  CPR;  
    Operand<dim>&  RR;  
    ScalarVariable      cpf, rf, phi,  rr, cpr;
    VectorVariable<dim> vel;
    double factor;

  public: // 					                               
    IAPWS_H2O_HeatTransferVelocity( const PropertyDatabase<dim>& p );
    ~IAPWS_H2O_HeatTransferVelocity() {};
    void Calculate();
};

}

#endif

