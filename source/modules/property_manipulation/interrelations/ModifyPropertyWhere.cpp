#include "ModifyPropertyWhere.h"

using namespace std;


namespace csmp {
/** calculates 'hydraulic conductivity'

the hydraulic conductivity of element i is calculated from the permeability of
the element and the average of the nodal fluid viscosities.
This conductivity is multiplied with the nodal fluid density average since
the goal is to calculate the mass flux.

@note If the property has the value of 'where', it is changed to the value of 'to'
*/

template<size_t dim>
ModifyPropertyWhere<dim>::ModifyPropertyWhere( const PropertyDatabase<dim>& p, 
                                               const char* prop, 
                                               double64 where, double64 to ) 
      : Interrelation<dim>(p),
        k( Interrelation<dim>::GlobalProperty(prop) ),
        if_value(where),
        to_value(to)
 {
    Interrelation<dim>::Name("ModifyPropertyWhere");
    Interrelation<dim>::OutputCondition( k, INIT_COND );
    Interrelation<dim>::ResultProperty(prop);
 }

template class ModifyPropertyWhere<1U>;
template class ModifyPropertyWhere<2U>;
template class ModifyPropertyWhere<3U>;

} // csmp
