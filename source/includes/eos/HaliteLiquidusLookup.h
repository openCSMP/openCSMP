// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef   HALITELIQUIDUSLOOKUP_H
#define   HALITELIQUIDUSLOOKUP_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

#include "VLH_LiquidLookup.h"
#include "NaClMeltingCurveLiquidLookup.h"
#include "States.h"

namespace csmp
{

  class HaliteLiquidusLookup 
  {

  public:
    HaliteLiquidusLookup(  const double& externaltemperature,
                           const double& externalpressure);
    ~HaliteLiquidusLookup();

    double                      Temperature();
    double                      Pressure();
    double                      MassFractionNaCl();
    // these may be abandoned upon code clean-up
    double                      Density();
    double                      Enthalpy();
    double                      HeatCapacity();
    double                      Compressibility();
    double                      Viscosity();
    double                      DCompositionDT();
    double                      DEnthalpyDT();
    double                      DSaltMassFractionDT();
    double                      ReportComposition();
    double                      ReportEnthalpy();
    double                      ValueOf(const int& property_index);
 
  private:

    const double&               temperature;
    const double&               pressure;

    double                      tcurrent;
    double                      pcurrent;
    double                      xcurrent;
    double                      tdummy;
    double                      pdummy;
    double                      t_res;
    double                      p_res;
    double                      tnorm;
    double                      pnorm;
    double                      pvlh;
    double                      t_iA;
    double                      t_iB;
    double                      t_iC;
    double                      t_iD;
    double                      p_iA;
    double                      p_iB;
    double                      p_iC;
    double                      p_iD;
    double                      value_iA;
    double                      value_iB;
    double                      value_iC;
    double                      value_iD;
    double                      value_before;
    double                      value_behind;
    double                      value_bottom;
    double                      value_top;
    double                      value_interpolated;
    double                      value_vlh;
    double                      value_vlh_behind;
    double                      value_vlh_before;
    double                      tcrit;
    double                      tvlh;
    double                      x_high;

    VLH_LiquidLookup              vlh_liquid;
    NaClMeltingCurveLiquidLookup  naclmelt_liquid;

    const double                vlh_tmax;
    const double                vlh_pmax;

    long                          it;
    long                          ip;
    long                          t_dim;
    long                          p_dim;
    long                          iA;
    long                          iB;
    long                          iC;
    long                          iD;
    long                          i_dummy;
    long                          it_p_max;
    long                          ip_p_max;

    States                        state;
    States                        state_iA;
    States                        state_iB;
    States                        state_iC;
    States                        state_iD;

    std::vector<double>         storage_vector;    // stores data in sequence t-p-x-rho-h at each Lookup point
    std::vector<States>           state_vector;         // stores fluid state at each Lookup point

    double                      NormalInterpolation( const int& property_index );
    double                      NearVLHInterpolationLowT( const int& property_index );
    double                      NearVLHInterpolationHighT( const int& property_index );
    double                      NearNaClMeltInterpolation( const int& property_index );
    double                      NearVLHMaxInterpolation(const int& property_index);
  
    void                          SetTemperatureAndPressure();
    void                          GetIndex_iA(const int& property_index);
    void                          GetTemperatureIndex(const double& t);
    void                          GetPressureIndex(const double& p);

    ErrorHandler&                 csmp_error;
   };

}// namespace csmp
#endif
