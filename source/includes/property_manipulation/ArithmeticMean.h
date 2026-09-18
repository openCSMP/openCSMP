// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ARITHMETIC_AVERAGE_H
#define ARITHMETIC_AVERAGE_H

#include "Interrelation.h"

namespace csmp {

/// hualp! this one does only the arithmetic average, else you need a visitor!
//enum AVERAGE { WEIGHTED_AVERAGE, ARITHMETIC_MEAN, HARMONIC_MEAN, GEOMETRIC_MAIN };

/// the usual kinds of averages
template<uint32_t dim, typename var>
class ArithmeticMean : public Interrelation<dim> {
    Operand<dim>&  I_;   // property to average
    Operand<dim>&  O_;   // property into which to write the result
    var  prop_value_;  
    ArithmeticMean();     
    
  public:
    ArithmeticMean( const PropertyDatabase<dim>& p, const char* res_prop, const char* prop_to_avg );
    
    void Calculate() override final;
};


} // csmp

#endif // ARITHMETIC_AVERAGE_H

