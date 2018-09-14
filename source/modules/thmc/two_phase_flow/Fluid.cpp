//
//  Fluid.cpp
//  
//
//  Created by Stephan Matthai on 19/01/2018.
//
//

#include "Fluid.h"
#include "CO2H2O_FunctionsModule1.h"
// this specific incarnation
#include "EOS_CO2H2ONaCl_Spycher2004.h"

namespace csmp {

EOS_CO2H2ONaCl_Spycher04  eos;

template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( const Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U ); 
    /*uncomment to compute from PTX properties*/
    if ( phase == 0U )
      return eos.mu_brine( n->Read(User()->key_pf), n->Read(User()->key_T), n->Read(User()->key_NaClaq) );
   
    return eos.mu_CarbonicPhase( n->Read(User()->key_pf), n->Read(User()->key_T) );
 }



template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    // TODO: replace with single interpolation loop
    const double64 pf = e->PropertyValueAtBaryCenter(User()->key_pf);
    const double64 T  = e->PropertyValueAtBaryCenter(User()->key_T);
   
    /*uncomment to compute from PTX properties*/
    if ( phase == 0U ) {
         const double64 NaCl_aq = e->PropertyValueAtBaryCenter(User()->key_NaClaq);
         return eos.mu_brine( pf, T, NaCl_aq );
      }
   
    return eos.mu_CarbonicPhase( pf, T );
 }

 

  

 
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( const Node<dim>* const n, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );
    /*uncomment to compute from PTX properties*/
    if ( phase == 0U ) return eos.Rho_brine( n->Read(User()->key_pf), n->Read(User()->key_T), n->Read(User()->key_NaClaq) );
      return eos.Rho_CarbonicPhase( n->Read(User()->key_pf), n->Read(User()->key_T) );
    
    /*direct interpolation*/
    if ( phase == 0U ) return n->Read( User()->key_rhoH2O );
    return n->Read( User()->key_rhoCO2 );
 }



template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( const Element<dim>* const e, size_t phase ) const
 {
    assert( phase == 0U or phase == 1U );

    // TODO: replace with single interpolation loop
    const double64 pf = e->PropertyValueAtBaryCenter(User()->key_pf);
    const double64 T  = e->PropertyValueAtBaryCenter(User()->key_T);

    /*uncomment to compute from PTX properties*/
    if ( phase == 0U ) {
         const double64 NaCl_aq = e->PropertyValueAtBaryCenter(User()->key_NaClaq);
         return eos.Rho_brine( pf, T, NaCl_aq );
      }

    return eos.Rho_CarbonicPhase( pf, T );
 }

  



 
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::DensityMixture( const Node<dim>* const n, double64 msalt ) const
 {
    return n->Read( User()->key_sH2O ) * eos.Rho_brine( n->Read(User()->key_pf), n->Read(User()->key_T), n->Read(User()->key_NaClaq) ) +
          (1. - n->Read(User()->key_sH2O)) * eos.Rho_CarbonicPhase( n->Read(User()->key_pf), n->Read(User()->key_T) );
 }




/**
     aqueous phase viscosity / carbonic phase viscosity
*/
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::ViscosityRatio( const Node<dim>* const n, double64 salinity ) const
 {
    return eos.mu_AqueousPhase( n->Read(User()->key_pf), n->Read(User()->key_T), n->Read(User()->key_NaClaq) ) / eos.mu_CarbonicPhase( n->Read(User()->key_pf), n->Read(User()->key_T) );
 }




/*
    // versions that account for dissolved CO2; TODO: check where the XCO2 and YH2) can be taken into account
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Viscosity( double64 pf, double64 T, double64 msalt, double64 XCO2, double64 YH2O, size_t phase ) const
 {
    assert( phase <= 1U );
    if ( phase == 0U ) return eos.mu_AqueousPhase( Pressure(p), n->Read(User()->key_T), msalt );
    return eos.mu_CarbonicPhase( Pressure(p), n->Read(User()->key_T) );
 }
 
 
 
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::Density( double64 pf, double64 T, double64 msalt, double64 XCO2, double64 YH2O, size_t phase ) const
 {
    assert( phase <= 1U );  // TODO: check how to get XCO2, YH2O accounted for?
    if ( phase == 0U ) return eos.Rho_AqueousPhase( pf, T, msalt );
    return eos.Rho_CarbonicPhase( pf, T );
 }
 
 
 
template<size_t dim, template<size_t> class USER>
double64 Fluid<dim,USER>::DensityMixture( double64 pf, double64 T, double64 sw, double64 msalt, double64 XCO2, double64 YH2O ) const
 {
    // TODO: check how to get XCO2, YH2O accounted for?
    return sw * eos.Rho_AqueousPhase( pf, T, msalt ) + (1. - sw) * eos.Rho_CarbonicPhase( pf, T );
 }
*/
 

template class Fluid<1U,CO2H2O_FunctionsModule1>;
template class Fluid<2U,CO2H2O_FunctionsModule1>;
template class Fluid<3U,CO2H2O_FunctionsModule1>;

template class Fluid<1U,CO2H2O_FunctionsModule2>;
template class Fluid<2U,CO2H2O_FunctionsModule2>;
template class Fluid<3U,CO2H2O_FunctionsModule2>;






// conversions


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
