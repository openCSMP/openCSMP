#ifndef   TWOPHASELIQUIDLOOKUP_H
#define   TWOPHASELIQUIDLOOKUP_H

#include "CSMP_definitions.h"

#include "CriticalCurveLookup.h"
#include "CriticalPointH2O.h"
#include "VLH_LiquidLookup.h"
#include "H2OLookup.h"
#include "States.h"
#include "ErrorHandler.h"

namespace csmp{

  class TwophaseLiquidLookup {

  public:
    TwophaseLiquidLookup(const double& externaltemperature, 
                         const double& externalpressure);
    ~TwophaseLiquidLookup();

    double              Temperature();
    double              Pressure();
    double              MassFractionNaCl();
    double              Density();
    double              Enthalpy();
    double              HeatCapacity();
    double              Compressibility();
    double              Viscosity();
    double              DCompositionDT();
    double              DSaltMassFractionDT();
    double              DEnthalpyDT();
    double              ReportMassFractionNaCl();
    double              ReportEnthalpy();
    double              MinXResolution();
    long                  EqType();

  private:

    const double&       temperature;
    const double&       pressure;

    double              tcurrent;
    double              pcurrent;
    double              xcurrent;
    double              tdummy;
    double              pdummy;
    double              t_res;
    double              p_res;
    double              tnorm;
    double              pnorm;
    double              pvlh;
    double              ph2o;
    double              t_iA;
    double              t_iB;
    double              t_iC;
    double              t_iD;
    double              p_iA;
    double              p_iB;
    double              p_iC;
    double              p_iD;
    double              value_iA;
    double              value_iB;
    double              value_iC;
    double              value_iD;
    double              value_bottom;
    double              value_top;
    double              value_before;
    double              value_behind;
    double              value_interpolated;
    double              value_crit;
    double              value_crit_behind;
    double              value_vlh;
    double              value_vlh_behind;
    double              value_vlh_before;
    double              value_boil;
    double              value_boil_behind;
    double              pboil_iA;
    double              pboil_iB;
    double              tboil;
    double              tcrit;
    double              tvlh;
    double              x_low;
    double              x_high;
    double              pcrit_iB;
    double              pcrit_iA;

    long                   it;
    long                   ip;
    long                   t_dim;
    long                   p_dim;
    long                   iA;
    long                   iB;
    long                   iC;
    long                   iD;
    long                   i_dummy;
    long                   eqtype;
    long                   it_p_max;
    long                   ip_p_max;

    States                 state;
    States                 state_iA;
    States                 state_iB;
    States                 state_iC;
    States                 state_iD;

    std::vector<double> storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point
    std::vector<States>    state_vector;         // stores fluid state at each Lookup point

    double              Interpolate(const int& property_index);
    double              NormalInterpolation( const int& property_index );
    double              NearCritpointInterpolation( const int& property_index );
    double              NearCritcurveInterpolation( const int& property_index );
    double              NearBoilingCurveInterpolation( const int& property_index );
    double              NearVLHInterpolationLowT( const int& property_index );
    double              NearVLHInterpolationHighT( const int& property_index );
    double              InterpolateBetweenBoilingCurveAndVLH( const int& property_index);
    double              NearVLHMaxInterpolation(const int& property_index);
   
    void                   SetTemperatureAndPressure();
    void                   GetIndex_iA(const int& property_index);
    void                   GetTemperatureIndex(const double& t);
    void                   GetPressureIndex(const double& p);

    CriticalPointH2O       cp_h2o;
    CriticalCurveLookup    critcurve;
    VLH_LiquidLookup       vlh_liquid;
    H2OLookup              water;
 
    const double        vlh_tmax;
    const double        vlh_pmax;

    ErrorHandler&         csmp_error;
  };

  inline long TwophaseLiquidLookup::EqType(){ return eqtype;}

}// namespace csmp
#endif
