#ifndef   NACLMELTINGCURVELIQUIDLOOKUP_H
#define   NACLMELTINGCURVELIQUIDLOOKUP_H

#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"
#include "States.h"
#include <cmath>

namespace csmp
{

  class NaClMeltingCurveLiquidLookup
  {

  public:

    NaClMeltingCurveLiquidLookup(const double64& externaltemperature,
                                 const double64& externalpressure);
    ~NaClMeltingCurveLiquidLookup();

    double64                TmeltFromP();
    double64                PmeltFromT();
    double64                MassFractionNaCl();
    double64                Density();
    double64                Enthalpy();
    double64                HeatCapacity();
    double64                Compressibility();
    double64                Viscosity();
    double64                TfromP(const double64& press); 
    double64                ValueOf(const int& property_index);

  private:

    const double64&         temperature;
    const double64&         pressure;

    double64                tcurrent;
    double64                pcurrent;
    double64                xcurrent;
    double64                t_res;
    double64                tnorm;
    double64                pnorm;
    double64                tmin;

    long                    it;
    long                    t_dim;
    long                    i_guess;
    long                    i_max;
    long                    i_min;
    long                    it_min;
    long                    it_max;
    States                  state;


    std::vector<double64>   storage_vector; // stores data in sequence t-p-x-rho-h at each lookup point

    TriplePointNaCl         tp_nacl;

    void                    GetTemperatureIndex(const double64& myt);

  };

}// namespace csmp
#endif
