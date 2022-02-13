#include "IAPWS_H2O_dVdT_Source.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
IAPWS_H2O_dVdT_Source<dim>::IAPWS_H2O_dVdT_Source( const PropertyDatabase<dim>& pref,
                                                   double dt, double eL_rock ) 
      : Interrelation<dim>(pref),
        T( Interrelation<dim>::GlobalProperty("temperature") ),
        preT( Interrelation<dim>::GlobalProperty("previous temperature") ),
        X( Interrelation<dim>::GlobalProperty("porosity") ),
        Q( Interrelation<dim>::GlobalProperty("thermal expansion") ),
        A( Interrelation<dim>::GlobalProperty("fluid expansivity") ),
        // the change of volume, dV, of rock is calculated from the isotropic change in length, a
        // as dV = V1 3 a dT
        eT_rock(3 * eL_rock),
        delta_t(dt)
 {
    Interrelation<dim>::Name("IAPWS_H2O_dVdT_Source");
    Interrelation<dim>::OutputCondition( Q, PLAIN );
    Interrelation<dim>::ResultProperty("thermal expansion");
    
 }


template<uint32_t dim>
void IAPWS_H2O_dVdT_Source<dim>::SetTimeIncrement( double dt ) { delta_t = dt; }


template<uint32_t dim>
void IAPWS_H2O_dVdT_Source<dim>::Calculate()
 {
    T.AssignTo( temperature );
    preT.AssignTo( T_prev );
    A.AssignTo( alpha );
    X.AssignTo( phi );
    
    // assert that temperatures are non-negative
    assert( temperature() >= 0.0 );
    assert( T_prev()      >= 0.0 );
    
    // temperature change
    deltaT  = temperature() - T_prev();

    Q = ( ((phi() * alpha()) + (1.0 - phi()) * eT_rock) * deltaT ) / delta_t;
      
 } // end Calculate

 
template class IAPWS_H2O_dVdT_Source<1U>; 
template class IAPWS_H2O_dVdT_Source<2U>;
template class IAPWS_H2O_dVdT_Source<3U>;

} // end namespace csp
