#include "HydraulicHead.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
HydraulicHead<dim>::HydraulicHead( const PropertyDatabase<dim>& p, double ref_density ) 
      : Interrelation<dim>(p),
        P( Interrelation<dim>::GlobalProperty("absolute fluid pressure") ),
        E( Interrelation<dim>::GlobalProperty("elevation") ),
        H( Interrelation<dim>::GlobalProperty("hydraulic head") ),
        g(9.80665),        // m N
        rho0(ref_density)  // kg m-3
 {
    Interrelation<dim>::Name("HydraulicHead");
    Interrelation<dim>::OutputCondition( H, PLAIN );
    Interrelation<dim>::ResultProperty("hydraulic head");
 }
 
 
 
template<uint32_t dim>
inline void HydraulicHead<dim>::Calculate()
 {
    P.AssignTo( pres );
    E.AssignTo( height );
    
    // hydraulic head, e.g., Bear (1972), 
    // Freeze & Cherry (1979, eq. 2.16, p. 21) is:
    //  elevation head + piezometric head 
    H = height() + pres() / (rho0 * g);
 } 



template class HydraulicHead<1U>;
template class HydraulicHead<2U>;
template class HydraulicHead<3U>;

} // csmp

