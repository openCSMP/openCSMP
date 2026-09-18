// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HYDROSTATIC_PRESSURE_H
#define HYDROSTATIC_PRESSURE_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

/// Interrelation subclass that computes gravity-induced vertical variation of pore-pressure in a rock sequence
template<uint32_t dim>
class HydrostaticPressure : public Interrelation<dim> {
    Operand<dim>&    P;    // absolute fluid pressure
    Operand<dim>&    E;    // vertical elevation (zero at model base)
    ScalarVariable  height, pres;
    const double  rho0, g, zmax; // g = acceleration of gravity
    
  public:
    HydrostaticPressure( const PropertyDatabase<dim>&,
                         double highest_elevation,
                         double ref_density=1000. );
                         
    ~HydrostaticPressure() = default;
    
    /// uses p_h(z) = rho g z + patm for the calculation
    void Calculate() override final;
};

}

#endif

