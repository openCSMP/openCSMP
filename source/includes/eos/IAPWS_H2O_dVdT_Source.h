// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef IAPWS_H2O_DVDT_SOURCE_H
#define IAPWS_H2O_DVDT_SOURCE_H

#include "CSMP_definitions.h"
#include "Interrelation.h"
#include "steam4.h"

namespace csmp {

template<size_t dim>
class IAPWS_H2O_dVdT_Source : public Interrelation<dim> {
    Operand<dim>&  T;     // temperature (oC)
    Operand<dim>&  preT;  // previous temperature (oC)
    Operand<dim>&  X;     // porosity
    Operand<dim>&  Q;     // thermal expansion
    Operand<dim>&  A;     // fluid expansivity
    ScalarVariable alpha, temperature, T_prev, phi;
    double       eT, eT_rock, // thermal expansivity (m3 K-1)
                   deltaT,      // temperature change
                   delta_t;     // change in time

  public: // 					                     length-expansion coeff. for rock (concrete)
    IAPWS_H2O_dVdT_Source( const PropertyDatabase<dim>& p, double dt, double eL_rock=1.2e-5 );
    ~IAPWS_H2O_dVdT_Source() {};
    void Calculate();
    void SetTimeIncrement( double dt );
};

}

#endif

