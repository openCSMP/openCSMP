#ifndef __LH_HALITELOOKUP_H_
#define __LH_HALITELOOKUP_H_

#include "CSMP_definitions.h"

#include "VLH_HaliteLookup.h"
#include "NaClMeltingCurveHaliteLookup.h"
#include "States.h"
#include "ErrorHandler.h"

namespace csmp
{

  class LH_HaliteLookup 
  {

  public:
    LH_HaliteLookup(  const double& externaltemperature, 
                      const double& externalpressure);
    ~LH_HaliteLookup();

    double               Temperature();
    double               Pressure();
    double               MassFractionNaCl();
    // these may be abandoned upon code-cleanup
    double               Density();
    double               Enthalpy();
    double               HeatCapacity();
    double               Compressibility();
    double               Viscosity();
    double               DCompositionDT();
    double               DEnthalpyDT();
    double               DSaltMassFractionDT();
    double               ReportComposition();
    double               ReportEnthalpy();
    double               ValueOf(const int& property_index);
 
  private:

    const double&        temperature;
    const double&        pressure;

    double               tcurrent;
    double               pcurrent;
    double               xcurrent;
    double               tdummy;
    double               pdummy;
    double               t_res;
    double               p_res;
    double               tnorm;
    double               pnorm;
    double               pvlh;
    double               t_iA;
    double               t_iB;
    double               t_iC;
    double               t_iD;
    double               p_iA;
    double               p_iB;
    double               p_iC;
    double               p_iD;
    double               value_iA;
    double               value_iB;
    double               value_iC;
    double               value_iD;
    double               value_before;
    double               value_behind;
    double               value_bottom;
    double               value_top;
    double               value_interpolated;
    double               value_vlh;
    double               value_vlh_behind;
    double               value_vlh_before;
    double               tcrit;
    double               tvlh;
    double               x_high;

    long                   it;
    long                   ip;
    long                   t_dim;
    long                   p_dim;
    long                   iA;
    long                   iB;
    long                   iC;
    long                   iD;
    long                   i_dummy;
    long                   it_p_max;
    long                   ip_p_max;

    States                 state;
    States                 state_iA;
    States                 state_iB;
    States                 state_iC;
    States                 state_iD;

    std::vector<double>  storage_vector;    // stores data in sequence t-p-x-rho-h at each Lookup point
    std::vector<States>    state_vector;         // stores fluid state at each Lookup point

    double               NormalInterpolation( const int& property_index );
    double               NearVLHInterpolationLowT( const int& property_index );
    double               NearVLHInterpolationHighT( const int& property_index );
    double               NearNaClMeltInterpolation( const int& property_index );
    double               NearVLHMaxInterpolation(const int& property_index);
  
    void                   SetTemperatureAndPressure();
    void                   GetIndex_iA(const int& property_index);
    void                   GetTemperatureIndex(const double& t);
    void                   GetPressureIndex(const double& p);

    VLH_HaliteLookup                vlh_halite;
    NaClMeltingCurveHaliteLookup    naclmelt_halite;

    const double        vlh_tmax;
    const double        vlh_pmax;

    ErrorHandler&         csmp_error;
   };

}// namespace csmp
#endif
