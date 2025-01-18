#ifndef GROUND_WATER_DARCY_VELOCITY_H
#define GROUND_WATER_DARCY_VELOCITY_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class GroundwaterDarcyVelocity : public Interrelation<dim> {
    Operand<dim>&  DH;   // hydraulic head gradient
    Operand<dim>&  K;    // hydraulic conductivity
    Operand<dim>&  V;    // Darcy velocity
    VectorVariable<dim>  head_grad;
    ScalarVariable       hcond;
    
  public:
    GroundwaterDarcyVelocity( const PropertyDatabase<dim>& p );
    ~GroundwaterDarcyVelocity() {};
    void Calculate();
};

} 

#endif

