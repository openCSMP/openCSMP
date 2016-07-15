#ifndef HYDROSTATIC_PRESSURE_H
#define HYDROSTATIC_PRESSURE_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class HydrostaticPressure : public Interrelation<dim> {
    Operand<dim>&    P;    // absolute fluid pressure
    Operand<dim>&    E;    // vertical elevation (zero at model base)
    ScalarVariable  height, pres;
    const double64  rho0, g, zmax; // g = acceleration of gravity
    
  public:
    HydrostaticPressure( const PropertyDatabase<dim>& p,
                         double64 highest_elevation, 
                         double64 ref_density=1000. );
    ~HydrostaticPressure() {};
    void Calculate();
};

}

#endif

