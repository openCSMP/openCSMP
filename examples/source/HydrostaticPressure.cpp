#include "HydrostaticPressure.h"

using namespace std;

namespace csmp {

template<size_t dim>
HydrostaticPressure<dim>::HydrostaticPressure( const PropertyDatabase<dim>& p,
                                               double highest_elevation, 
                                               double ref_density ) 
      : Interrelation<dim>(p),
        E( Interrelation<dim>::GlobalProperty("elevation") ),
        P( Interrelation<dim>::GlobalProperty("absolute fluid pressure") ),
        g(9.80665),              // m N
        rho0(ref_density),       // kg m-3
        zmax(highest_elevation)  // m above base of model
 {
    Interrelation<dim>::Name("HydrostaticPressure");
    Interrelation<dim>::OutputCondition( P, PLAIN );
    Interrelation<dim>::ResultProperty("fluid pressure");
 }


template<size_t dim>
void HydrostaticPressure<dim>::Calculate()
 {
    E.AssignTo( height );
    
    P = (zmax - height()) * rho0 * g;
 } 

template class HydrostaticPressure<1U>;
template class HydrostaticPressure<2U>;
template class HydrostaticPressure<3U>;

}
