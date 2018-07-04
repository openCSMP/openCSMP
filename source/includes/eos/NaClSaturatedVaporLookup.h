#ifndef   NACLSATURATEDVAPORLOOKUP_H
#define   NACLSATURATEDVAPORLOOKUP_H

#include "CSMP_definitions.h"
#include "H2OLookup.h"
#include "CriticalPointH2O.h"
#include "VLH_VaporLookup.h"
#include "States.h"
#include "ErrorHandler.h"

#include <cmath>

namespace csmp
{

  class NaClSaturatedVaporLookup
  {

  public:
    NaClSaturatedVaporLookup(const double64& externaltemperature, const double64& externalpressure);
    ~NaClSaturatedVaporLookup();
    NaClSaturatedVaporLookup(const NaClSaturatedVaporLookup&);
    NaClSaturatedVaporLookup& operator=(const NaClSaturatedVaporLookup&);

    double64                Temperature();
    double64                Pressure();
    double64                MassFractionNaCl();
    double64                Density();
    double64                Enthalpy();
    double64                HeatCapacity();
    double64                Compressibility();
    double64                Viscosity();
    double64                DCompositionDT();
    double64                DMassFractionNaClDT();
    double64                ReportMassFractionNaCl();
    double64                ReportEnthalpy();
    double64                DEnthalpyDT();
    double64                ValueOf(const int& property_index);

  private:

    const double64&         temperature;
    const double64&         pressure;

    double64                tcurrent;
    double64                pcurrent;
    double64                xcurrent;
    double64                tdummy;
    double64                t_res;
    double64                p_res;
    double64                tnorm;
    double64                pnorm;
    double64                t_iA;
    double64                t_iB;
    double64                t_iC;
    double64                t_iD;
    double64                p_iA;
    double64                p_iB;
    double64                p_iC;
    double64                p_iD;
    double64                pvlh_iB;
    double64                pvlh_iA;
    double64                value_bottom;
    double64                value_top;
    double64                value_interpolated;
    double64                value_iA;
    double64                value_iB;
    double64                value_iC;
    double64                value_iD;
    double64                value_before;
    double64                value_behind;
    double64                value_vlh;
    double64                value_vlh_behind;
    double64                value_vlh_before;
    double64                tvlh;
    double64                x_low;
    double64                x_high;
    double64                pmax_low;
    double64                pmax_high;

    long                    it;
    long                    ip;
    long                    t_dim;
    long                    p_dim;
    long                    iA;
    long                    iB;
    long                    iC;
    long                    iD;
    long                    i_dummy;
    long                    it_p_max;
    long                    ip_p_max_low;
    long                    ip_p_max_high;
    States                  state;
    States                  state_iA;
    States                  state_iB;
    States                  state_iC;
    States                  state_iD;
    
    std::vector<double64>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point
    std::vector<States>     state_vector;         // stores fluid state at each Lookup point

    double64 NormalInterpolation( const int& property_index );
    double64 NearVLHInterpolationLowT( const int& property_index );
    double64 NearVLHInterpolationHighT( const int& property_index );
    double64 NearVLHMaxInterpolation( const int& property_index );
   
    void SetTemperatureAndPressure();
    void GetIndex_iA(const int& property_index);
    void GetTemperatureIndex(const double64& t);
    void GetPressureIndex(const double64& p);

    CriticalPointH2O               cp_h2o;
    VLH_VaporLookup                vlh_vapor;
    H2OLookup                      water;

    ErrorHandler&                  csmp_error;
  };

}// namespace csmp
#endif
