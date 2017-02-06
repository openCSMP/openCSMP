//
//  MechanicalProperties.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 6/23/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "MechanicalProperties.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

/// granite
MechanicalProperties::MechanicalProperties()
 : E(5.5e10),        // Young's modulus
   B(5.0e10),        // drained bulk modulus (K_dry)
   nu(0.2),          // Poisson's ratio
   TS(5.0e7),        // tensile strength
   mu(0.6),          // friction coefficient from 0 to >1 (some rubbers)
   G(0.5*E/(1+0.2)), // from relations between elastic moduli
   alpha(0.03),      // Biot coefficient phi <= alpha <= 1
   UCS(1.5e8),       // Hoek-Brown criterion inferred value
   pstar(0.6*E),     // about 5-7 times of the UCS of the rock
   C(0.2*E)          // cohesion = inherent shear strength
 {
 }


MechanicalProperties::MechanicalProperties( double64 E, double64 beta, double64 nu,
                                            double64 tensile_strength, double64 fric_coeff )
 : E(E),    // Young's modulus
   B(beta), // drained bulk modulus (K_dry)
   nu(nu),  // Poisson's ratio
   TS(tensile_strength),
   mu(fric_coeff),  // friction coefficient from 0 to >1 (some rubbers)
   // plausible values for sedimentary rocks
   G(0.5*E/(1+nu)), // from relations between elastic moduli
   alpha(0.6),      // Biot coefficient phi <= alpha <= 1
   UCS(0.1*E),      // Hoek-Brown criterion inferred value
   pstar(0.6*E),    // about 5-7 times of the UCS of the rock
   C(1.0e4)         // cohesion = inherent shear strength
 {
 }


/**
    Computes the friction angle (in degrees) from the friction coefficient
    
    The friction coefficient is:  mu = tan(theta),
 
    where theta is the friction angle in radians (degrees * 180./PI).
    Thus, theta as derived from the friction coefficient is:
    
    theta = arctan(mu)
*/
double64  frictionAngle( double64 fric_coeff )
 {
    const double64 radiansToDeg(180./CSMP_PI);
    return atan(fric_coeff) * radiansToDeg;
 }



/** 
    Tensile strength TS is calculated from the friction coefficient and
    the UCS = C0 using a modified Griffith criterian from Brace, 1960.
    Fjaer et al. 08', eqn. 2.27, p. 66
*/
double64  tensileStrengthFromUCS_Griffith( double64 UCS, double64 fric_coeff )
 {
    const double64 mu(fric_coeff);
    return 0.25 * UCS * sqrt(mu*mu - mu + 1.);
 }



void MechanicalProperties::Out(std::ostream& os) const
 {
    os <<"\nMechanicalProperties::Out:";
    os <<"\nE:     "<< E <<" Young's modulus (Pa)";
    os <<"\nB:     "<< B <<" bulk modulus = compressibility of rock skeleton = K_dry (Pa)";
    os <<"\nG:     "<< G <<" shear modulus = modulus of rigidity (Pa)";
    os <<"\nalpha: "<< alpha <<" Biot coefficient alpha = 1 - K_dry/K_grain (dimensionless ratio)";
    os <<"\nTS:    "<< TS <<" tensile strength (Pa)";
    os <<"\nUCS:   "<< UCS <<" unconfined compressive strength (Pa)";
    os <<"\npstar: "<< pstar <<" crushing pressure of the rock, see Fjaer et al. 08', p. 68 (Pa)";
    os <<"\nmu:    "<< mu <<" friction coefficient (dimensionless)";
    os <<"\nnu:    "<< nu <<" Poisson's ratio (dimensionless ratio)";
    os <<"\nC:     "<< C <<" cohesive strength = inherent shear strength (Pa).\n";
    os.flush();
 }




}  // end namespace csmp
