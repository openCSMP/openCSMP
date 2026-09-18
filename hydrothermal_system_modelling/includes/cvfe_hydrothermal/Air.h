// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef   Air_H
#define   Air_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

namespace csmp
{

class Air
{

public:

    Air();
    ~Air();

    double HeatCapacity( double t );
    double Enthalpy( double t );
    double Density( double t, double p );
    double Viscosity( double t );
    double Compressibility( double p );

private:
    ErrorHandler& csmp_error;
};

/**
     @class Air Air.h

     @author Benoit Lamy-Chappuis, ETH Zuerich

     @section motivation Motivation

     The Air class is used to compute the Air heat capacity and enthalpy at a given temperature.
     The Air heat capacity can either be a constant or be temperature dependent.


     @section usage Usage
     Create an instance of the Air class: Air Air, it will be initialized with default static values;
     Air Air(cp) will create a Air with constant heat capacity of value "cp".

     @section issues Known issues

  */

}

#endif
