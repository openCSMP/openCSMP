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
 {
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Viscosity( Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U ); 
    /*direct read/interpolation*/
    if ( phase == 0U ) return n->Read( User()->key_muH2O );
    return n->Read( User()->key_muCO2 );     
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Viscosity( Element<dim>* const e, size_t node, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    assert( node < e->Nodes() );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->N(node)->Read( User()->key_muH2O );
    return e->N(node)->Read( User()->key_muCO2 );
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Viscosity( Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->PropertyValueAtBaryCenter( User()->key_muH2O );
    return e->PropertyValueAtBaryCenter( User()->key_muCO2 );      
 }

 

  

 
template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Density( Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*direct read/interpolation*/
    if ( phase == 0U ) return n->Read( User()->key_rhoH2O );
    return n->Read( User()->key_rhoCO2 );
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Density( Element<dim>* const e, size_t node, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    assert( node < e->Nodes() );
    /*direct read/interpolation*/
    if ( phase == 0U ) return e->N(node)->Read( User()->key_rhoH2O );
    return e->N(node)->Read( User()->key_rhoCO2 );
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::Density( Element<dim>* const e, size_t phase ) const
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
double Fluid<dim,USER>::ViscosityRatio( Node<dim>* const n ) const
 {
    assert( n != nullptr );
    return n->Read(User()->key_muH2O) / n->Read(User()->key_muCO2);
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::ViscosityRatio( Element<dim>* const e, size_t node ) const
 {
    assert( e != nullptr );
    const double muw = e->N(node)->Read( User()->key_muH2O );
    const double mun = e->N(node)->Read( User()->key_muCO2 );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    return muw / mun;
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::ViscosityRatio( Element<dim>* const e ) const
 {
    assert( e != nullptr );
    const double muw = e->PropertyValueAtBaryCenter( User()->key_muH2O );
    const double mun = e->PropertyValueAtBaryCenter( User()->key_muCO2 );
    assert( !isnan(muw) );
    assert( !isnan(mun) );
    return muw / mun;
 }



template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::MixtureDensity( Node<dim>* const n ) const
 {
    assert( n != nullptr );
    const double sw = n->Read( User()->key_sH2O );
    return sw * n->Read(User()->key_rhoH2O) + (1.-sw) * n->Read(User()->key_rhoCO2);
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::MixtureDensity( Element<dim>* const e, size_t node ) const
 {
    const double sw = e->N(node)->Read( User()->key_sH2O );
    assert( e != nullptr );
    const double rhow = e->N(node)->Read( User()->key_rhoH2O );
    const double rhon = e->N(node)->Read( User()->key_rhoCO2 );
    assert( !isnan(rhow) );
    assert( !isnan(rhon) );
    return sw * rhow + (1. - sw) * rhon;
 }


template<size_t dim, template<size_t> class USER>
double Fluid<dim,USER>::MixtureDensity( Element<dim>* const e ) const
 {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert( e != nullptr );
    const double rhow = e->PropertyValueAtBaryCenter( User()->key_rhoH2O );
    const double rhon = e->PropertyValueAtBaryCenter( User()->key_rhoCO2 );
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

template class Fluid<1U,FlowFunctionsModule7>;
template class Fluid<2U,FlowFunctionsModule7>;
template class Fluid<3U,FlowFunctionsModule7>;



// unit conversions


const double molarMassH2o ( 18.01528e-3);  // Kilograms per mole
const double molarMassCo2 ( 44.010e-3 );    // Kilograms per mole
const double molarMassNacl ( 58.443e-3 );    // Kilograms per mole


/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double massFracNaClToMolalNaClInAqueousPhase( double massFracSalt )
{
    // mass fraction in % weight
    return 1. * massFracSalt / ( molarMassNacl * ( 100 - massFracSalt ) );
}



/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double molalNaClToMassFracNaClInAqueousPhase( double mSalt)
{
    double dummy = mSalt * molarMassNacl;
    return dummy / (1 + dummy) * 100.; // in % weight
}



double massFracNaClToMolarFracNaClInAqueousPhase( double massFracSalt)
{
    double molarFrac =  massFracSalt * molarMassH2o;

    double dummy = massFracSalt * (molarMassH2o - molarMassNacl) + 100 * molarMassNacl;

    return molarFrac / dummy;
}




double molalNaClToMolarFracNaClInAqueousPhase( double mSalt)
{
    double molarFrac = mSalt / (mSalt + 55.508);

    return molarFrac;
}



double molalNaClToPpmInAqueousPhase( double mSalt )
{
    double conversion( molalNaClToMassFracNaClInAqueousPhase( mSalt) );

    conversion  *=  1.e-2*1.e6; //  ppm \in [0, 1.e6] ,1.e-2 because mass fraction is in %

    return      conversion;
}



double psiToPa( double pressureInPsi )
{
    return pressureInPsi * 6894.75729;
}


double paToPsi( double pressureInPa )
{
    return pressureInPa * 0.000145037738;
}



double paTobar( double pressureInPa )
{
    return  pressureInPa * 1.0e-5;
}



double barTopa( double pressureInbar )
{
    return  pressureInbar * 1.0e+5;
}




double ppmNaClToMolalNaClInAqueousPhase( double ppmSalt )
{
    return massFracNaClToMolalNaClInAqueousPhase(ppmSalt*1.e-6*1.e2); //conversion from massfraction to molality
}



double degreeCToKelvin( double temperatureInC )
{
   return  temperatureInC + 273.15;
}



double KelvinTodegreeC( double temperatureInK )
{
   return  temperatureInK - 273.15;
}




} // end csmp
