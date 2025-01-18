#ifndef NACLMELTINGCURVEHALITELOOKUP_H
#define NACLMELTINGCURVEHALITELOOKUP_H

#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"


namespace csmp
{

  class NaClMeltingCurveHaliteLookup
  {

  public:
    NaClMeltingCurveHaliteLookup(const double& externaltemperature,
                                 const double& externalpressure);
    ~NaClMeltingCurveHaliteLookup();

    double               TmeltFromP();
    double               PmeltFromT();
    double               MassFractionNaCl();
    double               Density();
    double               Enthalpy();
    double               HeatCapacity();
    double               Compressibility();
    double               Viscosity();
    double               TfromP(const double& press); 
    double               ValueOf(const int& property_index);

  private:

    const double&        temperature;
    const double&        pressure;

    double               tcurrent;
    double               pcurrent;
    double               xcurrent;
    double               t_res;
    double               tnorm;
    double               pnorm;
    double               tmin;

    long                   it;
    long                   t_dim;
    long                   i_guess;
    long                   i_max;
    long                   i_min;
    long                   it_max;
    long                   it_min;
    int                    state;

    std::vector<double>  storage_vector; // stores data in sequence t-p-x-rho-h at each lookup point

    TriplePointNaCl        tp_nacl;

    void GetTemperatureIndex(const double& myt);

  };

}// namespace csmp
#endif






















