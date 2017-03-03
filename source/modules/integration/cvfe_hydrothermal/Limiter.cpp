#include "Limiter.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<size_t dim>
Limiter<dim>::Limiter( const PropertyDatabase<dim>& p,
					    const char* prop, double64 min_value, double64 max_value )
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
template<size_t dim>
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
