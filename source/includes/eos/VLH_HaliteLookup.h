#ifndef   VLH_HALITELOOKUP_H_
#define   VLH_HALITELOOKUP_H_

#include "ErrorHandler.h"
#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"
#include "States.h"

namespace csmp
{
  
  class VLH_HaliteLookup 
  {
    
  public:
    VLH_HaliteLookup(const double64& externaltemperature);
    ~VLH_HaliteLookup();
    
    double64                Temperature();
    double64                Pressure();
    double64                MassFractionNaCl();
    double64                Density();
    double64                Enthalpy();
    double64                HeatCapacity();
    double64                Compressibility();
    double64                Viscosity();
    double64                Pmax();
    double64                Tmax();
    double64                DEnthalpyDT();
    double64                TfromP(const double64& press, const double64& t_estimate); 
    double64                ValueOf(const int& property_index);

    std::vector<double64>   properties_at_tmax;

  private:

    const double64&         temperature;

    double64                tcurrent;
    double64                pcurrent;
    double64                xcurrent;
    double64                t_res;
    double64                tnorm;
    double64                pnorm;
    double64                dp;
    double64                tmax;
    double64                pmax;

    long                    it;
    long                    t_dim;
    long                    i_guess;
    long                    i_max;
    long                    i_min;
    long                    it_p_max;
    States                  state;
    
    std::vector<double64>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point

    TriplePointNaCl         tp_nacl;

    void GetTemperatureIndex(const double64& myt);

    ErrorHandler&           csmp_error;
  };

  inline double64 VLH_HaliteLookup::Tmax(){ return tmax; }
  inline double64 VLH_HaliteLookup::Pmax(){ return pmax; }
}// namespace csmp
#endif






















