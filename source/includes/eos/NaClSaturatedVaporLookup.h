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
    NaClSaturatedVaporLookup(const double& externaltemperature, const double& externalpressure);
    ~NaClSaturatedVaporLookup();
    NaClSaturatedVaporLookup(const NaClSaturatedVaporLookup&);
    NaClSaturatedVaporLookup& operator=(const NaClSaturatedVaporLookup&);

    double                Temperature();
    double                Pressure();
    double                MassFractionNaCl();
    double                Density();
    double                Enthalpy();
    double                HeatCapacity();
    double                Compressibility();
    double                Viscosity();
    double                DCompositionDT();
    double                DMassFractionNaClDT();
    double                ReportMassFractionNaCl();
    double                ReportEnthalpy();
    double                DEnthalpyDT();
    double                ValueOf(const int& property_index);

  private:

    const double&         temperature;
    const double&         pressure;

    double                tcurrent;
    double                pcurrent;
    double                xcurrent;
    double                tdummy;
    double                t_res;
    double                p_res;
    double                tnorm;
    double                pnorm;
    double                t_iA;
    double                t_iB;
    double                t_iC;
    double                t_iD;
    double                p_iA;
    double                p_iB;
    double                p_iC;
    double                p_iD;
    double                pvlh_iB;
    double                pvlh_iA;
    double                value_bottom;
    double                value_top;
    double                value_interpolated;
    double                value_iA;
    double                value_iB;
    double                value_iC;
    double                value_iD;
    double                value_before;
    double                value_behind;
    double                value_vlh;
    double                value_vlh_behind;
    double                value_vlh_before;
    double                tvlh;
    double                x_low;
    double                x_high;
    double                pmax_low;
    double                pmax_high;

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
    
    std::vector<double>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point
    std::vector<States>     state_vector;         // stores fluid state at each Lookup point

    double NormalInterpolation( const int& property_index );
    double NearVLHInterpolationLowT( const int& property_index );
    double NearVLHInterpolationHighT( const int& property_index );
    double NearVLHMaxInterpolation( const int& property_index );
   
    void SetTemperatureAndPressure();
    void GetIndex_iA(const int& property_index);
    void GetTemperatureIndex(const double& t);
    void GetPressureIndex(const double& p);

    CriticalPointH2O               cp_h2o;
    VLH_VaporLookup                vlh_vapor;
    H2OLookup                      water;

    ErrorHandler&                  csmp_error;
  };

}// namespace csmp
#endif
