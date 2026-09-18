// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Limiter.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim>
Limiter<dim>::Limiter( const PropertyDatabase<dim>& p,
					    const char* prop, double min_value, double max_value )
      : Interrelation<dim>(p),
        op( Interrelation<dim>::GlobalProperty(prop) ),
        min(min_value),
        max(max_value)
 {
    Interrelation<dim>::Name("Limiter");
    Interrelation<dim>::OutputCondition( op, INIT_COND );
    Interrelation<dim>::ResultProperty(prop);
 }

/** Limit values to specified range between min_value and max_value */
template<uint32_t dim>
void Limiter<dim>::Calculate()
 {

    if ( op < min )
     {
      op = min;
     }
    if ( op > max )
     {
      op = max;
     }
 }

template class Limiter<1U>;
template class Limiter<2U>;
template class Limiter<3U>;

} // csmp
