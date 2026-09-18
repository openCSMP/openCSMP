// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef  IAPWS_H20_DENSITY_H
#define  IAPWS_H20_DENSITY_H

#include "CSMP_definitions.h"
#include "steam4.h"

namespace csmp {

class IAPWS_H2O_Density {
    Prop           *props, *lprops, *sprops;   // PROST output structures
    double       t, p, d, h;
    const double dp;
  public:
    IAPWS_H2O_Density();
    ~IAPWS_H2O_Density();
    double Density( double T, double P );
    double Viscosity( double T, double P );
    double Density_P_MPa( double T, double P );
    double VapourDensity( double T, double P );
    double LiquidDensity( double T, double P );
    double Pressure( double T, double rho );
    double HeatCapacity( double T, double P );
    double HeatCapacityDensTemp( double rho, double T );
    double EnthalpyDensTemp( double rho, double T );
    double Enthalpy( double P, double T );
    double LiquidSaturation( double T, double P );
    double SaturationTemperature( double P );
    double SaturationPressure( double T );
};

} // end namespace csmp

#endif
