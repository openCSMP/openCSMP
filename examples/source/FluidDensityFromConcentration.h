#ifndef FLUID_DENSITY_FROM_CONCENTRATION_H
#define FLUID_DENSITY_FROM_CONCENTRATION_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class FluidDensityFromConcentration : public Interrelation<dim> {
    Operand<dim>&   DENS;    
    Operand<dim>&   CONC;    
    ScalarVariable  conc;
    const double64  rho_zero;
    double64        increment;
    
  public:
    FluidDensityFromConcentration( const PropertyDatabase<dim>& p, double64 rho_max );
    ~FluidDensityFromConcentration() {};
    void Calculate();
};

}

#endif

