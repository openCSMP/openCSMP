// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IAPWS_H2O_Density.h"

using namespace std;

namespace csmp {

IAPWS_H2O_Density::IAPWS_H2O_Density()
 : dp(1.0e-08) // convergence criterion 
 {
    props  = newProp('p', 't', 2);
    lprops = newProp('p', 't', 2);
    sprops = newProp('p', 't', 2);
 }
 
IAPWS_H2O_Density::~IAPWS_H2O_Density() {}





/// outputs density (kg m-3) for temperature (oC) and pressure (Pa)
double IAPWS_H2O_Density::Density( double T, double P )
{
    if ( P < 101325.0 )
      {
         cout <<"\n IAPWS_H2O_Density::Density: detected sub-atmospheric pressure: ";
         cout << P <<" at ToC: "<< T << endl;
         P = 101325.0;
      }
    T += 273.15;
    water_tp(T, P, d, dp, props);
    return props->d;
}





/// outputs density (kg m-3) for temperature (oC) and pressure (MPa)
double IAPWS_H2O_Density::Density_P_MPa( double T, double P )
 {
    T += 273.15;
    P *= 1.0e+06;
    water_tp(T, P, d, dp, props);
    return props->d;
 }




/// outputs pressure (Pa) for temperature (oC) and density (kg m-3)
double IAPWS_H2O_Density::Pressure( double T, double rho )
 {
    T += 273.15;
    water_td(T, rho, props);
    return props->p;
 }



/// outputs heat capacity (J K-1 kg-1) for temperature (oC) and pressure (Pa)
double IAPWS_H2O_Density::HeatCapacity( double T, double P )
 {
    T += 273.15;
    water_tp(T, P, d, dp, props);
    return props->cp;
 }




/// outputs heat capacity (J K-1 kg-1) from density (oC) and temperature (oC)
double IAPWS_H2O_Density::HeatCapacityDensTemp( double rho, double T )
 {
    T += 273.15;
    water_td(T, rho, props);
    return props->cp;
 }




/// outputs enthalpy=heat content (J K-1 kg-1) from density (oC) and temperature (oC)
double IAPWS_H2O_Density::EnthalpyDensTemp( double rho, double T )
 {
    T += 273.15;
    water_td(T, rho, props);
    return props->h;
 }

double IAPWS_H2O_Density::Enthalpy( double P, double T )
 {
    T += 273.15;
    water_tp(T, P, d, dp, props);
    return props->h;
 }

double IAPWS_H2O_Density::Viscosity( double P, double T )
 {
    T += 273.15;
    water_tp(T, P, d, dp, props);
    return viscos(props);
 }





/// outputs  vapour density (kg m-3) for temperature (oC) and pressure (Pa)
double IAPWS_H2O_Density::LiquidSaturation( double T, double P )
 {
    sat_p( P, lprops, sprops ); 
    T +=273.15;
    water_tp(T, P, d, dp, props);
    // all liquid
    if      ( T < lprops->T ) return 1.0;
    // all steam
    else if ( T > lprops->T ) return 0.0;
    // liquid + steam
    else                      return (1.0 - props->x);
    
 }


}








