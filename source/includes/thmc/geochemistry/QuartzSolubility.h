// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef QUARTZ_SOLUBILITY_H
#define QUARTZ_SOLUBILITY_H

#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class QuartzSolubility : public Interrelation<dim> {
    Operand<dim>&  T; /// < temperature (oC)
    Operand<dim>&  R; /// < fluid density (kg m-3)
    Operand<dim>&  S; /// < quartz solubility (kg silica / kg fluid)
    ScalarVariable  sc_tC, rho;
    double tC, v, d1, d2, d3, d4, sol;

  public:
    QuartzSolubility( const PropertyDatabase<dim>& );

    void Calculate() override final;
};

} // csmp

#endif
















