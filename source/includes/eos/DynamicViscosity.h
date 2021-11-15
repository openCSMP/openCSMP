#ifndef DYNAMIC_VISCOSITY_H
#define DYNAMIC_VISCOSITY_H

#include "CSMP_definitions.h"

namespace csmp {

class DynamicViscosity {
  
  public:
    DynamicViscosity();
    ~DynamicViscosity();
    double  ViscosityFromTemperatureAndDensity( double t, double rho );

  private:
    double ak[4], bij[6][5];
    const double kelvin, rhostar, tstar;

};

}

#endif
