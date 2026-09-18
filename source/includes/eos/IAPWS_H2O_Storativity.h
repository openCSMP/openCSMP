// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef IAPWS_H2O_STORATIVITY_H
#define IAPWS_H2O_STORATIVITY_H

#include "CSMP_definitions.h"
#include "Interrelation.h"
#include "steam4.h"

namespace csmp {

template<size_t dim>
class IAPWS_H2O_Storativity : public Interrelation<dim> {
  public:
    IAPWS_H2O_Storativity( const PropertyDatabase<dim>& p, double z_rock=1.0e-12 );
    ~IAPWS_H2O_Storativity() {};
    void Calculate();
  
  private:
    Operand<dim>&  B;  // fluid compressibility
    Operand<dim>&  X;  // porosity
    Operand<dim>&  S;  // storativity
    ScalarVariable beta, phi;
    double       Z, Z_rock;
};

} // end namespace csmp

#endif

