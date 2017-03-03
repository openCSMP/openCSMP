#include "IAPWS_H2O_Conductivity.h"

using namespace std;

namespace csmp {

template<size_t dim>
IAPWS_H2O_Conductivity<dim>::IAPWS_H2O_Conductivity( const PropertyDatabase<dim>& pref ) 
      : Interrelation<dim>(pref),
        k( Interrelation<dim>::GlobalProperty("permeability") ),
        K( Interrelation<dim>::GlobalProperty("conductivity") ),
        v( Interrelation<dim>::GlobalProperty("fluid viscosity") )
 {
    Interrelation<dim>::Name("IAPWS_H2O_Conductivity");
    Interrelation<dim>::OutputCondition( K, PLAIN );
    Interrelation<dim>::ResultProperty("conductivity");
 }


template<size_t dim>
void IAPWS_H2O_Conductivity<dim>::Calculate()
 {
    k.AssignTo( perm );
    v.AssignTo( viscos );
    K = perm / viscos; // permeability / fluid viscosity
 } // end Calculate
 
 
template class IAPWS_H2O_Conductivity<1U>;
template class IAPWS_H2O_Conductivity<2U>;
template class IAPWS_H2O_Conductivity<3U>;

} // end namespace csmp
