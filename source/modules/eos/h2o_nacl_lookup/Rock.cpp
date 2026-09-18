// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <iostream>
#include <cmath>
#include "Rock.h"
#include "CSMP_mathUtilities.h"


using namespace std;

namespace csmp
{

  Rock::Rock()
    :
    csmp_error( ErrorHandler::Instance() ),
    cp           (880),
    t_dependent  (true),

    b            (1.),
    nu           (1.78),
    sigma1       (1.),
    h_fusion     (300000.)
  {
    cout << "\nRock constructor:\n";
    cout << " Creating a Rock with";
    cout << " temperature dependent heat capacity\n";
    cout << " with a baseline (low temperature) value of: " << cp << "J/kg/K\n\n";
  }

  Rock::Rock( double heatcapacity )
    :
    csmp_error( ErrorHandler::Instance() ),
    cp (heatcapacity),
    t_dependent (false),

    b            (1.),
    nu           (1.78),
    sigma1       (1.),
    h_fusion     (300000.)
  {
    cout << "Rock::Rock( const double& temperature, double heatcapacity )\n";
    cout << "constructor:\n";
    cout << "      creating Rock with temperature-independent\n";
    cout << "      heat capacity of " << cp << "J/kg/K\n\n";
  }

  Rock::~Rock()
  {
  }

  double Rock::HeatCapacity( double t )
  {
    if (!t_dependent)
      {
        return cp;
      }

    else
      {
        if ( t < 750.)
          return cp;
        else if ( t < 800.)
          return cp + (t - 750.) / 50.*cp; // to avoid a jump function
        else
          return 2.0 * cp; // doubled
      }
  }

  double Rock::HeatCapacity( double t, double tl, double ts)
  {
    if (tl < 0. or ts < 0.)
      {
        cerr << "\nRock::HeatCapacity( double t, double tl, double ts) called with tl = " << tl << " and ts = " << ts;
        csmp_error.Note( FATAL_ERROR, "Rock::HeatCapacity( double t, double tl, double ts)",
                         "The variable tl (liquidus temperature) and/or ts (solidus temperature) has not been assigned a proper value" );
      }

    if (tl < ts)
      {
        cerr << "\nRock::HeatCapacity( double t, double tl, double ts) called with tl = " << tl << " and ts = " << ts;
        csmp_error.Note( FATAL_ERROR, "Rock::HeatCapacity( double t, double tl, double ts)",
                         "The rock liquidus temperature is inferior to the solidus temperature!" );
      }

    if (crystallization_curve == "power law")
      {
        if ( t <= ts )
          return cp;
        else if ( t <= tl )
          return cp * (nu + (1 - nu) * (1 - pow(((t - ts) / (tl - ts)), b)));
        else
          return nu * cp;
      }

    else if (crystallization_curve == "marxer ulmer")
      {
        // heat capacity as per crystallization curve parameterization Marxer Ulmer 2019;
        double a(2.74619571e-12), b(-1.08683749e-08), c(1.69816241e-05), d(-1.30731969e-02), e(4.94565516e+00), f(-7.32026831e+02);

        if (t <= ts)
          return cp;
        else if (t <= tl)
          return cp * (nu + (1 - nu) * (a * pow(t, 5) + b * pow(t, 4) + c * pow(t, 3) + d * pow(t, 2) + e * (t) + f));
        else
          return nu * cp;
      }

    else if (crystallization_curve == "error function")
      {
        if ( t <= ts )
          return cp;
        else if ( t <= tl )
          return cp * (nu + (1 - nu) * (0.5 * (erfc(sqrt(30) * sigma1 * ((t - ts) / (tl - ts) - 0.5)))));
        else
          return nu * cp;
      }
    else
      {
        throw csmp::Exception( FATAL_ERROR, "Rock::HeatCapacity",
                               "Invalid crystallization curve: " + crystallization_curve +
                               "\nAvailable options: 'error function', 'power law', 'marxer ulmer'." );
      }
  }

  double Rock::MinimumHeatCapacity()
  {
    return cp;
  }

  double Rock::Enthalpy( double t )
  {
    if (!t_dependent)
      {
        return cp * t;
      }

    else
      {
        if (t < 750.)
          return cp * t;
        else if (t < 800.)
          return
            t * cp // base value
            + 0.5 * cp / 50. * (t - 750.) * (t - 750.); // 0.5*slope*delta_t^2-0.5*slope*0^2;
        else
          return
            t * cp
            + 0.5 * cp * 50.
            + (t - 800.) * cp;
      }
  }

  double Rock::Enthalpy(double t, double tl, double ts)
  {
    if (tl < 0. or ts < 0.)
      {
        cerr << "\nRock::Enthalpy(double t, double tl, double ts) called with tl = " << tl << " and ts = " << ts;
        csmp_error.Note( FATAL_ERROR, "Rock::HeatCapacity( double t, double tl, double ts)",
                         "The variable tl (liquidus temperature) and/or ts (solidus temperature) has not been assigned a proper value" );
      }

    if (tl < ts)
      {
        cerr << "\nRock::HeatCapacity( double t, double tl, double ts) called with tl = " << tl << " and ts = " << ts;
        csmp_error.Note( FATAL_ERROR, "Rock::HeatCapacity( double t, double tl, double ts)",
                         "The rock liquidus temperature is inferior to the solidus temperature!" );
      }

    if (crystallization_curve == "power law")
      {
        if (t <= ts)
            return cp * t;
        else if(t <= tl)
            return cp * (t-(1-nu)*(tl-ts)/(b+1.)*pow(((t-ts)/(tl-ts)),(b+1.)))+h_fusion*pow(((t-ts)/(tl-ts)),b);
        else
            return cp*(tl-(1-nu)*(tl-ts)/(b+1)+nu*(t-tl))+h_fusion;
      }

      else if (crystallization_curve == "marxer ulmer")
      {
          // enthalpy as per crystallization curve parameterization Marxer Ulmer 2019;
          double a(2.74619571e-12), b(-1.08683749e-08), c(1.69816241e-05), d(-1.30731969e-02), e(4.94565516e+00), f(-7.32026831e+02);

          if (t <= ts)
              return cp * t;
          else if (t <= tl)
              return cp * nu * (t - ts) + cp * ts  + cp * (1 - nu) *
                                                        (a * (tl - ts) / 6 * pow((t - ts) / (tl - ts), 6) +
                                                         b * (tl - ts) / 5 * pow((t - ts) / (tl - ts), 5) +
                                                         c * (tl - ts) / 4 * pow((t - ts) / (tl - ts), 4) +
                                                         d * (tl - ts) / 3 * pow((t - ts) / (tl - ts), 3) +
                                                         e * (tl - ts) / 2 * pow((t - ts) / (tl - ts), 2) +
                                                         f * (t - ts) / (tl - ts)) +
                     h_fusion * (1 - (a * pow(t, 5) + b * pow(t, 4) + c * pow(t, 3) + d * pow(t, 2) + e * t + f));
          else
              return cp * (nu * (tl - ts) + (1 - nu) *
                                                (a / 6 * (tl - ts) +
                                                 b / 5 * (tl - ts) +
                                                 c / 4 * (tl - ts) +
                                                 d / 3 * (tl - ts) +
                                                 e / 2 * (tl - ts) +
                                                 f) ) +
                     cp * ts + h_fusion * (1 - (a * pow(tl, 5) + b * pow(tl, 4) + c * pow(tl, 3) + d * pow(tl, 2) + e * tl + f)) + nu * cp * (t - tl) ;
      }

    else if (crystallization_curve == "error function")
      {
        alpha = sqrt(30.) * sigma1 / (tl - ts);
        beta = -sqrt(30.) * sigma1 * ts / (tl - ts) - sqrt(30.) * sigma1 * 0.5;

        if (t <= ts)
          return cp * t;
        else if (t <= tl)
            return cp * ts + cp / 2.*((1 + nu) * (t - ts) - (1 - nu) * ((beta / alpha + t) * erf(alpha * t + beta) - (beta / alpha + ts) * erf(alpha * ts + beta) +
                                                                        1. / (sqrt(PI) * alpha) * (exp(-pow((alpha * t + beta), 2.)) - exp(-pow((alpha * ts + beta), 2.))))) +
                                                                        h_fusion * (1 - 0.5 * (erfc(alpha * t + beta)));

        else
            return cp * ts + cp / 2.*((1 + nu) * (tl - ts) - (1 - nu) * ((beta / alpha + tl) * erf(alpha * tl + beta) - (beta / alpha + ts) * erf(alpha * ts + beta) + 1. /
                                                                        (sqrt(PI) * alpha) * (exp(-pow((alpha * tl + beta), 2.)) - exp(-pow((alpha * ts + beta), 2.))))) +
                                                                        nu * cp * (t - tl) + h_fusion;
      }

      else
      {
          throw csmp::Exception( FATAL_ERROR, "Rock::Enthalpy",
                          "Invalid crystallization curve: " + crystallization_curve +
                              "\nAvailable options: 'error function', 'power law', 'marxer ulmer'." );
      }
  }


  void Rock::SetRockHeatCapacity(double mini_cp)
  {
    cp          = mini_cp;
    cout << "\nSetting the baseline low temperature rock heat capacity to: " << cp << endl;
  }

  void Rock::SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                         double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_)
  {
      b           = b_coefficient;
      nu          = nu_coefficient;
      sigma1      = sigma1_coefficient;
      h_fusion    = latent_heat_of_fusion;
      crystallization_curve = crystallization_curve_;
      cout << "\nSetting the rock melting/crystallization curve parameters";
  }

}
