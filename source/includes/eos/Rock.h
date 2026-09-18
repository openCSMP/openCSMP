// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef   ROCK_H
#define   ROCK_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

/**
     @class Rock Rock.h

     @author Philipp Weis, Benoit Lamy-Chappuis, ETH Zuerich

     @section motivation Motivation

     The rock class is used to compute the rock heat capacity and enthalpy at a given temperature.
     The rock heat capacity can either be a constant or be temperature dependent.
     As of now, the sole motivation to make the heat capacity temperature dependent is to account for the latent heat of fusion/crystallization of the rock.

     In HeatCapacity( double t); the heat capacity is doubled in a linear fashion in the 750-800C temperature window.
     This was used in several papers from Philipp Weis.

     With HeatCapacity(double t, double tl, double ts) one can set specific liquidus and solidus temperatures as well as various parameters that
     define the crystallization curve between the solidus and liquidus: it can either be a power law function or a "logistic function".
     nu; ratio of melt to rock heat capacity
     b;  exponent used for power law option
     sigma1; a parameter of the logistic function
     h_fusion; latent heat of fusion (J/Kg)
     crystallization_curve; string defines option to use "power law", "error function", or "marxer ulmer" parametrisation to calculate HeatCapacity and Enthalpy.

     You can find a more detailed description of those parameters in Lamy-Chappuis et al. 2020 +supp materials and references therein.

     @section usage Usage
     Create an instance of the rock class: Rock rock, it will be initialized with default static values;
     Rock rock(cp) will create a rock with constant heat capacity of value "cp".
     Optional: modify the parameters of the rock with Set_Rock or Set_Rock_crystallization_curve

     @section issues Known issues
     The way we use this class in the CVFEM_PHX_Scheme is fairly limited. A single rock "type" is ever used.
     If one want to use various rock types in a simulation, it is possible to create several rock instances, however one would then need to query heat capacity
     and enthalpy values for these various rocks on a region by region basis in the CVFEM_PHX_Scheme (+dependencies), it would need to be implemented.

  */

namespace csmp
{

class Rock
{

public:

    Rock();                        // to be used for T-dependent heat capacity, various parameters can be set with Set_Rock_heat_capacity and Set_Rock_crystallization_curve
    Rock( double heatcapacity ); // to be used with constant heat capacity

    ~Rock();

    double HeatCapacity( double t);
    double HeatCapacity(double t, double tl, double ts);
    double MinimumHeatCapacity();
    double Enthalpy( double t);
    double Enthalpy(double t, double tl, double ts);


    void SetRockHeatCapacity( double mini_cp );
    void SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                       double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve);

private:

    double alpha;
    double beta;
    ErrorHandler &csmp_error;
    double cp;
    bool   t_dependent;

    double b;
    double nu;
    double sigma1;
    double h_fusion;

    std::string crystallization_curve;
};


}

#endif
