#ifndef FLUID_DENSITY_FROM_CONCENTRATION_H
#define FLUID_DENSITY_FROM_CONCENTRATION_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class FluidDensityFromConcentration : public Interrelation<dim> {
    Operand<dim>&   DENS;    
    Operand<dim>&   CONC;    
    ScalarVariable  conc;
    const double  rho_zero;
    double        increment;
    
  public:
    FluidDensityFromConcentration( const PropertyDatabase<dim>& p, double rho_max );
    ~FluidDensityFromConcentration() = default;
    
    void Calculate() override final;
};

}

#endif

