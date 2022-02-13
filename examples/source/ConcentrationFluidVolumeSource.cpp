#include "ConcentrationFluidVolumeSource.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ConcentrationFluidVolumeSource<dim>::ConcentrationFluidVolumeSource( const PropertyDatabase<dim>& p,
                                                                     double max_rho ) 
      : Interrelation<dim>(p),
        CONCN( Interrelation<dim>::GlobalProperty("concentration") ),
        CONCP( Interrelation<dim>::GlobalProperty("previous concentration") ),
        Q( Interrelation<dim>::GlobalProperty("fluid volume source") ),
        PHI( Interrelation<dim>::GlobalProperty("porosity") ),
        dt(1.0),
        rho_zero(1000.0),
        rho_increment( max_rho - 1000.0 )
 {
    Interrelation<dim>::Name("ConcentrationFluidVolumeSource");
    Interrelation<dim>::OutputCondition( Q, PLAIN );
    Interrelation<dim>::ResultProperty("fluid volume source");
 }


template<uint32_t dim>
void ConcentrationFluidVolumeSource<dim>::Calculate()
 {
    CONCN.AssignTo( concn );
    CONCP.AssignTo( concp );
    PHI.AssignTo( phi );
    Q  = ( rho_increment* ( concn() - concp() ) );
    Q /= ( rho_zero + rho_increment * concp() );
    Q *= -1.0;
    Q /= dt;
    Q *= phi;
 } 

template class ConcentrationFluidVolumeSource<2U>;
template class ConcentrationFluidVolumeSource<3U>;

}
