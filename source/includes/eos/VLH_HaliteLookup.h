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
    VLH_HaliteLookup(const double& externaltemperature);
    ~VLH_HaliteLookup();
    
    double                Temperature();
    double                Pressure();
    double                MassFractionNaCl();
    double                Density();
    double                Enthalpy();
    double                HeatCapacity();
    double                Compressibility();
    double                Viscosity();
    double                Pmax();
    double                Tmax();
    double                DEnthalpyDT();
    double                TfromP(const double& press, const double& t_estimate); 
    double                ValueOf(const int& property_index);

    std::vector<double>   properties_at_tmax;

  private:

    const double&         temperature;

    double                tcurrent;
    double                pcurrent;
    double                xcurrent;
    double                t_res;
    double                tnorm;
    double                pnorm;
    double                dp;
    double                tmax;
    double                pmax;

    long                    it;
    long                    t_dim;
    long                    i_guess;
    long                    i_max;
    long                    i_min;
    long                    it_p_max;
    States                  state;
    
    std::vector<double>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point

    TriplePointNaCl         tp_nacl;

    void GetTemperatureIndex(const double& myt);

    ErrorHandler&           csmp_error;
  };

  inline double VLH_HaliteLookup::Tmax(){ return tmax; }
  inline double VLH_HaliteLookup::Pmax(){ return pmax; }
}// namespace csmp
#endif






















