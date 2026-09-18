// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IAPWS_H2O_HeatTransferVelocity.h"

using namespace std;

namespace csmp {

template<size_t dim>
IAPWS_H2O_HeatTransferVelocity<dim>::IAPWS_H2O_HeatTransferVelocity( const PropertyDatabase<dim>& pref ) 
      : Interrelation<dim>(pref),
        V( Interrelation<dim>::GlobalProperty("velocity") ),
        P( Interrelation<dim>::GlobalProperty("porosity") ),
        VH( Interrelation<dim>::GlobalProperty("heat transfer velocity") ),
        CPF( Interrelation<dim>::GlobalProperty("fluid heat capacity") ),
        RF( Interrelation<dim>::GlobalProperty("fluid density") ),
        CPR( Interrelation<dim>::GlobalProperty("rock heat capacity") ),
        RR( Interrelation<dim>::GlobalProperty("rock density") )
 {
    Interrelation<dim>::Name("IAPWS_H2O_HeatTransferVelocity");
    Interrelation<dim>::OutputCondition( VH, PLAIN );
    Interrelation<dim>::ResultProperty("heat transfer velocity");
 }


template<size_t dim>
void IAPWS_H2O_HeatTransferVelocity<dim>::Calculate()
 {
    V.AssignTo( vel );
    CPF.AssignTo( cpf );
    RF.AssignTo( rf );
    P.AssignTo( phi );
	  CPR.AssignTo( cpr );
	  RR.AssignTo( rr );
    
    // solving for 
    // [ ( 1.0 - phi ) cp_r * rho_r + phi * cp_f * rho_f ] * dT/dt = Grad [ vel * cp_f * rho_f * T ]
    // and computing the thermal retardation factor
    // note, the darcy and not the pore velocity must be taken!
    factor  = ( rf() * cpf() );
    factor /= ( phi() * cpf() * rf() + ( 1.0 - phi() ) * rr() * cpr() );
    VH = vel * factor; 
 
 } // end Calculate
 
 
template class IAPWS_H2O_HeatTransferVelocity<1U>;
template class IAPWS_H2O_HeatTransferVelocity<2U>;
template class IAPWS_H2O_HeatTransferVelocity<3U>;

} // end namespace csp
