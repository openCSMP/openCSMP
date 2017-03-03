#ifndef DYNAMIC_VISCOSITY_H
#define DYNAMIC_VISCOSITY_H

#include "CSMP_definitions.h"

namespace csmp {

class DynamicViscosity {
  
  public:
    DynamicViscosity();
    ~DynamicViscosity();
    double64  ViscosityFromTemperatureAndDensity( double64 t, double64 rho );

  private:
    double64 ak[4], bij[6][5];
    const double64 kelvin, rhostar, tstar;

};

}

#endif
