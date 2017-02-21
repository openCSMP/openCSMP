#ifndef  IAPWS_H20_DENSITY_H
#define  IAPWS_H20_DENSITY_H

#include "CSMP_definitions.h"
#include "iaps_h2o_eos/steam4.h"

namespace csmp {

class IAPWS_H2O_Density {
    Prop           *props, *lprops, *sprops;   // PROST output structures
    double64       t, p, d, h;
    const double64 dp;
  public:
    IAPWS_H2O_Density();
    ~IAPWS_H2O_Density();
    double64 Density( double64 T, double64 P );
    double64 Viscosity( double64 T, double64 P );
    double64 Density_P_MPa( double64 T, double64 P );
    double64 VapourDensity( double64 T, double64 P );
    double64 LiquidDensity( double64 T, double64 P );
    double64 Pressure( double64 T, double64 rho );
    double64 HeatCapacity( double64 T, double64 P );
    double64 HeatCapacityDensTemp( double64 rho, double64 T );
    double64 EnthalpyDensTemp( double64 rho, double64 T );
    double64 Enthalpy( double64 P, double64 T );
    double64 LiquidSaturation( double64 T, double64 P );
    double64 SaturationTemperature( double64 P );
    double64 SaturationPressure( double64 T );
};

} // end namespace csmp

#endif
