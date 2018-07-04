#ifndef   VLH_VAPORLOOKUP_H
#define   VLH_VAPORLOOKUP_H

#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"
#include "ErrorHandler.h"
#include "States.h"

namespace csmp
{

  class VLH_VaporLookup 
  {

  public:
    VLH_VaporLookup(const double64& externaltemperature);
    ~VLH_VaporLookup();
 
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
    double64                TfromP(const double64& press, const double64& t_estimate); 
    double64                DPressureDT();
    double64                DEnthalpyDT();
    double64                ValueOf(const int& property_index);

    std::vector<double64>   properties_at_tmax;

  private:
    const double64&         temperature;

    double64                tcurrent;
    double64                pcurrent;
    double64                xcurrent;
    double64                tmax;
    double64                pmax;
    double64                t_res;
    double64                tnorm;
    double64                pnorm;
    double64                dp;

    long                    it;
    long                    t_dim;
    long                    i_guess;
    long                    i_max;
    long                    i_min;
    long                    it_p_max;
    long                    tdim_times_pindex;

    States                  state;

    std::vector<double64>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point

    TriplePointNaCl                      tp_nacl;

    void GetTemperatureIndex(const double64& t);

    ErrorHandler&           csmp_error;
  };

  inline double64 VLH_VaporLookup::Tmax(){ return tmax; }
  inline double64 VLH_VaporLookup::Pmax(){ return pmax; }
}// namespace csmp
#endif
