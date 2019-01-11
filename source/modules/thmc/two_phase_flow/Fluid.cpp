//
//  Fluid.cpp
//  
//
//  Created by Stephan Matthai on 19/01/2018.
//
//

#include "Fluid.h"
#include "FlowFunctionsModule.h"

namespace csmp {

template<size_t dim, template<size_t> class USER>
Fluid<dim,USER>::Fluid()
 // testing the key that is used to indicate the phase state
 : state_key_( User()->key_k )
 {
     assert( state_key_.place == NODE );
     assert( state_key_.type == SCALAR );
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( const Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U ); 
    /*direct read/interpolation*/
    if ( phase == 0U ) return n->Read( User()->key_muH2O );
    return n->Read( User()->key_muCO2 );     
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( const Element<dim>* const e, size_t node, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    assert( node < e->Nodes() );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->N(node)->Read( User()->key_muH2O );
    return e->N(node)->Read( User()->key_muCO2 );
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->PropertyValueAtBaryCenter( User()->key_muH2O );
    return e->PropertyValueAtBaryCenter( User()->key_muCO2 );      
 }

 

  

 
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( const Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*direct read/interpolation*/
    if ( phase == 0U ) return n->Read( User()->key_rhoH2O );
    return n->Read( User()->key_rhoCO2 );
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( const Element<dim>* const e, size_t node, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    assert( node < e->Nodes() );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->N(node)->Read( User()->key_rhoH2O );
    return e->N(node)->Read( User()->key_rhoCO2 );
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->PropertyValueAtBaryCenter( User()->key_rhoH2O );
    return e->PropertyValueAtBaryCenter( User()->key_rhoCO2 );     
 }

  



/**
     aqueous phase viscosity / carbonic phase viscosity
*/
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::ViscosityRatio( const Node<dim>* const n ) const
 {
    assert( n != nullptr );
    return n->Read(User()->key_muH2O) / n->Read(User()->key_muCO2);
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::ViscosityRatio( const Element<dim>* const e, size_t node ) const
 {
    assert( e != nullptr );
    const double64 muw = e->N(node)->Read( User()->key_muH2O );
    const double64 mun = e->N(node)->Read( User()->key_muCO2 );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    return muw / mun;
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::ViscosityRatio( const Element<dim>* const e ) const
 {
    assert( e != nullptr );
    const double64 muw = e->PropertyValueAtBaryCenter( User()->key_muH2O );
    const double64 mun = e->PropertyValueAtBaryCenter( User()->key_muCO2 );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    return muw / mun;
 }



template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::MixtureDensity( const Node<dim>* const n ) const
 {
    assert( n != nullptr );
    const double64 sw = n->Read( User()->key_sH2O );
    return sw * n->Read(User()->key_rhoH2O) + (1.-sw) * n->Read(User()->key_rhoCO2);
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::MixtureDensity( const Element<dim>* const e, size_t node ) const
 {
    const double64 sw = e->N(node)->Read( User()->key_sH2O );
    assert( e != nullptr );
    const double64 rhow = e->N(node)->Read( User()->key_rhoH2O );
    const double64 rhon = e->N(node)->Read( User()->key_rhoCO2 );
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    return sw * rhow + (1. - sw) * rhon;
 }


template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::MixtureDensity( const Element<dim>* const e ) const
 {
    const double64 sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert( e != nullptr );
    const double64 rhow = e->PropertyValueAtBaryCenter( User()->key_rhoH2O );
    const double64 rhon = e->PropertyValueAtBaryCenter( User()->key_rhoCO2 );
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    return sw * rhow + (1. - sw) * rhon;
 }




template class Fluid<1U,FlowFunctionsModule1>;
template class Fluid<2U,FlowFunctionsModule1>;
template class Fluid<3U,FlowFunctionsModule1>;

template class Fluid<1U,FlowFunctionsModule2>;
template class Fluid<2U,FlowFunctionsModule2>;
template class Fluid<3U,FlowFunctionsModule2>;

template class Fluid<1U,FlowFunctionsModule3>;
template class Fluid<2U,FlowFunctionsModule3>;
template class Fluid<3U,FlowFunctionsModule3>; 

template class Fluid<1U,FlowFunctionsModule4>;
template class Fluid<2U,FlowFunctionsModule4>;
template class Fluid<3U,FlowFunctionsModule4>; 

template class Fluid<1U,FlowFunctionsModule5>;
template class Fluid<2U,FlowFunctionsModule5>;
template class Fluid<3U,FlowFunctionsModule5>; 

template class Fluid<1U,FlowFunctionsModule6>;
template class Fluid<2U,FlowFunctionsModule6>;
template class Fluid<3U,FlowFunctionsModule6>; 




// unit conversions


const double64 molarMassH2o ( 18.01528e-3);  // Kilograms per mole
const double64 molarMassCo2 ( 44.010e-3 );    // Kilograms per mole
const double64 molarMassNacl ( 58.443e-3 );    // Kilograms per mole


/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double64 massFracNaClToMolalNaClInAqueousPhase( double64 massFracSalt )
{
    // mass fraction in % weight
    return 1. * massFracSalt / ( molarMassNacl * ( 100 - massFracSalt ) );
}



/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double64 molalNaClToMassFracNaClInAqueousPhase( double64 mSalt)
{
    double64 dummy = mSalt * molarMassNacl;
    return dummy / (1 + dummy) * 100.; // in % weight
}



double64 massFracNaClToMolarFracNaClInAqueousPhase( double64 massFracSalt)
{
    double64 molarFrac =  massFracSalt * molarMassH2o;

    double64 dummy = massFracSalt * (molarMassH2o - molarMassNacl) + 100 * molarMassNacl;

    return molarFrac / dummy;
}




double64 molalNaClToMolarFracNaClInAqueousPhase( double64 mSalt)
{
    double64 molarFrac = mSalt / (mSalt + 55.508);

    return molarFrac;
}



double64 molalNaClToPpmInAqueousPhase( double64 mSalt )
{
    double64 conversion( molalNaClToMassFracNaClInAqueousPhase( mSalt) );

    conversion  *=  1.e-2*1.e6; //  ppm \in [0, 1.e6] ,1.e-2 because mass fraction is in %

    return      conversion;
}



double64 psiToPa( double64 pressureInPsi )
{
    return pressureInPsi * 6894.75729;
}


double64 paToPsi( double64 pressureInPa )
{
    return pressureInPa * 0.000145037738;
}



double64 paTobar( double64 pressureInPa )
{
    return  pressureInPa * 1.0e-5;
}



double64 barTopa( double64 pressureInbar )
{
    return  pressureInbar * 1.0e+5;
}




double64 ppmNaClToMolalNaClInAqueousPhase( double64 ppmSalt )
{
    return massFracNaClToMolalNaClInAqueousPhase(ppmSalt*1.e-6*1.e2); //conversion from massfraction to molality
}



double64 degreeCToKelvin( double64 temperatureInC )
{
   return  temperatureInC + 273.15;
}



double64 KelvinTodegreeC( double64 temperatureInK )
{
   return  temperatureInK - 273.15;
}




} // end csmp
