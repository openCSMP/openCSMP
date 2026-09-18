// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
    VLH_VaporLookup(const double& externaltemperature);
    ~VLH_VaporLookup();
 
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
    double                TfromP(const double& press, const double& t_estimate); 
    double                DPressureDT();
    double                DEnthalpyDT();
    double                ValueOf(const int& property_index);

    std::vector<double>   properties_at_tmax;

  private:
    const double&         temperature;

    double                tcurrent;
    double                pcurrent;
    double                xcurrent;
    double                tmax;
    double                pmax;
    double                t_res;
    double                tnorm;
    double                pnorm;
    double                dp;

    long                    it;
    long                    t_dim;
    long                    i_guess;
    long                    i_max;
    long                    i_min;
    long                    it_p_max;
    long                    tdim_times_pindex;

    States                  state;

    std::vector<double>   storage_vector; // stores data in sequence t-p-x-rho-h at each Lookup point

    TriplePointNaCl                      tp_nacl;

    void GetTemperatureIndex(const double& t);

    ErrorHandler&           csmp_error;
  };

  inline double VLH_VaporLookup::Tmax(){ return tmax; }
  inline double VLH_VaporLookup::Pmax(){ return pmax; }
}// namespace csmp
#endif
