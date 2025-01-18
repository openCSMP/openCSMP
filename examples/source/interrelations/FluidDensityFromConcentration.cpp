#include "FluidDensityFromConcentration.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
FluidDensityFromConcentration<dim>::FluidDensityFromConcentration( const PropertyDatabase<dim>& p, double rho_max )
      : Interrelation<dim>(p),
        CONC( Interrelation<dim>::GlobalProperty("concentration") ),
        DENS( Interrelation<dim>::GlobalProperty("fluid density") ),
        rho_zero( 1000.0 ),
        increment( rho_max - rho_zero )
 {
    Interrelation<dim>::Name("FluidDensityFromConcentration");
    Interrelation<dim>::OutputCondition( DENS, PLAIN );
    Interrelation<dim>::ResultProperty("fluid density");
 }


template<uint32_t dim>
void FluidDensityFromConcentration<dim>::Calculate()
 {
    CONC.AssignTo( conc );
    DENS = rho_zero + conc * increment;
 } 

template class FluidDensityFromConcentration<2U>;
template class FluidDensityFromConcentration<3U>;

}
