// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_INTERRELATION_STREAMING_POTENTIAL_H
#define CSMP_INTERRELATION_STREAMING_POTENTIAL_H

#include "Interrelation.h"
#include "PropertyHandle.h"

namespace csmp {

/// algorithm and interrelation required to compute streaming potential
template<uint32_t dim>
class StreamingPotential {
  public:
    StreamingPotential( Model<dim>& );
    ~StreamingPotential();
    
    void EvaluatePotential( Model<dim>&,
                            const char* potential="streaming potential" );

  private:
    PropertyHandle<dim>  ones,    /// <  potential source
                         ccoeff;  /// < coupling coefficient
 };
 
/// to compute the input coefficients
template<uint32_t dim>
class PotentialSource : public Interrelation<dim> {
    Operand<dim>&  mu;     // dynamic viscosity [Pa s-1]
    Operand<dim>&  k;      // permeability    [m2]
    Operand<dim>&  rhof;   // fluid density   [kg m-3]
    Operand<dim>&  R;      // R factor
    
    ScalarVariable  temp1, temp2;
    
  public:
    PotentialSource( const PropertyDatabase<dim>& );

    void Calculate() override final;
};
 
 
} // csmp 

#endif // CSMP_INTERRELATION_STREAMING_POTENTIAL_H


