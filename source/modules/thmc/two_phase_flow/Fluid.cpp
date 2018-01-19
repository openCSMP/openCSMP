//
//  Fluid.cpp
//  
//
//  Created by Stephan Matthai on 19/01/2018.
//
//

#include "Fluid.h"
#include "CSMP_definitions.h"

namespace csmp {

// conversions


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
    double64
            conversion( molalNaClToMassFracNaClInAqueousPhase( mSalt) );

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
