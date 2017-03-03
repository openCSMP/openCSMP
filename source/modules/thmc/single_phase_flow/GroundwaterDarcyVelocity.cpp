#include "GroundwaterDarcyVelocity.h"

using namespace std;

namespace csmp {

template<size_t dim>
GroundwaterDarcyVelocity<dim>::GroundwaterDarcyVelocity( const PropertyDatabase<dim> &p )
      : Interrelation<dim>(p),
        DH( Interrelation<dim>::GlobalProperty("hydraulic head gradient") ),
        K( Interrelation<dim>::GlobalProperty("conductivity") ),
        V( Interrelation<dim>::GlobalProperty("velocity") )
 {
    Interrelation<dim>::Name("GroundwaterDarcyVelocity");
    Interrelation<dim>::OutputCondition( V, PLAIN );
    Interrelation<dim>::ResultProperty("velocity");
 }


template<size_t dim>
void GroundwaterDarcyVelocity<dim>::Calculate()
 {
    DH.AssignTo( head_grad );
    K.AssignTo( hcond );
   
    head_grad *= hcond();
    head_grad *= -1.0;    // since flow is down the gradient
    
    V = head_grad;
 } 

template class GroundwaterDarcyVelocity<1U>;
template class GroundwaterDarcyVelocity<2U>;
template class GroundwaterDarcyVelocity<3U>;

}
