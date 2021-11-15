#include "IAPWS_H2O_Storativity.h"

using namespace std;

namespace csmp {

template<size_t dim>
IAPWS_H2O_Storativity<dim>::IAPWS_H2O_Storativity( const PropertyDatabase<dim>& pref, double z_rock ) 
      : Interrelation<dim>(pref),
        B( Interrelation<dim>::GlobalProperty("fluid compressibility") ),
        X( Interrelation<dim>::GlobalProperty("porosity") ),
        S( Interrelation<dim>::GlobalProperty("storativity") ),
        Z_rock(z_rock)
 {
    Interrelation<dim>::Name("IAPWS_H2O_Storativity");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty("storativity");

 }



template<size_t dim>
void IAPWS_H2O_Storativity<dim>::Calculate()
 {
    B.AssignTo( beta );
    X.AssignTo( phi );
    S = phi() * beta() + (1.0 - phi()) * Z_rock;

 } // end Calculate
 
 
template class IAPWS_H2O_Storativity<1U>;
template class IAPWS_H2O_Storativity<2U>;
template class IAPWS_H2O_Storativity<3U>;

} // end namespace csmp

 
