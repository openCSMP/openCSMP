// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  MechanicalProperties.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 6/23/13.
//

#include "MechanicalProperties.h"
#include "CSMP_physical_constants.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/// granite
MechanicalProperties::MechanicalProperties()
 : E_(5.5e10),        // Young's modulus
   B_(5.0e10),        // drained bulk modulus (K_dry)
   nu_(0.2),          // Poisson's ratio
   TS_(5.0e7),        // tensile strength
   mu_(0.6),          // friction coefficient from 0 to >1 (some rubbers)
   G_(0.5*E_/(1+0.2)), // from relations between elastic moduli
   alpha_(0.03),      // Biot coefficient phi <= alpha <= 1
   UCS_(1.5e8),       // Hoek-Brown criterion inferred value
   pstar_(0.6*E_),    // about 5-7 times of the UCS of the rock
   C_(0.2*E_)         // cohesion = inherent shear strength
 {
 }




MechanicalProperties::MechanicalProperties( const std::vector<double>& p )
 {
   if ( p.size() < 10 )
     throw csmp::Exception( ERROR, "MechanicalProperties::ctor(vector):", "input vector too small to initialise values." );
   E_     = p[0];  // Young's modulus
   B_     = p[1];  // drained bulk modulus (K_dry)
   nu_    = p[2];  // Poisson's ratio
   TS_    = p[3];  // tensile strength
   mu_    = p[4];  // friction coefficient from 0 to >1 (some rubbers)
   G_     = p[5];  // from relations between elastic moduli
   alpha_ = p[6];  // Biot coefficient phi <= alpha <= 1
   UCS_   = p[7];  // Hoek-Brown criterion inferred value
   pstar_ = p[8];  // about 5-7 times of the UCS of the rock
   C_     = p[9];  // cohesion = inherent shear strength
 }




MechanicalProperties::MechanicalProperties( double E, double beta, double nu,
                                            double tensile_strength, double fric_coeff )
 : E_(E),    // Young's modulus
   B_(beta), // drained bulk modulus (K_dry)
   nu_(nu),  // Poisson's ratio
   TS_(tensile_strength),
   mu_(fric_coeff),  // friction coefficient from 0 to >1 (some rubbers)
   // plausible values for sedimentary rocks
   G_(0.5*E/(1+nu)), // from relations between elastic moduli
   alpha_(0.6),      // Biot coefficient phi <= alpha <= 1
   UCS_(0.1*E),      // Hoek-Brown criterion inferred value
   pstar_(0.6*E),    // about 5-7 times of the UCS of the rock
   C_(1.0e4)         // cohesion = inherent shear strength
 {
 }


/**
    Computes the friction angle (in degrees) from the friction coefficient
    
    The friction coefficient is:  mu = tan(theta),
 
    where theta is the friction angle in radians (degrees * 180./PI).
    Thus, theta as derived from the friction coefficient is:
    
    theta = arctan(mu)
*/
double  frictionAngle( double fric_coeff )
 {
    const double radiansToDeg(180./CSMP_PI);
    return atan(fric_coeff) * radiansToDeg;
 }



/** 
    Tensile strength TS is calculated from the friction coefficient and
    the UCS = C0 using a modified Griffith criterian from Brace, 1960.
    Fjaer et al. 08', eqn. 2.27, p. 66
*/
double  tensileStrengthFromUCS_Griffith( double UCS, double fric_coeff )
 {
    const double mu(fric_coeff);
    return 0.25 * UCS * sqrt(mu*mu - mu + 1.);
 }



void MechanicalProperties::Out() const
 {
    cout <<"\nMechanicalProperties::Out:";
    cout <<"\n\tE:     "<< E_ <<" Young's modulus (Pa)";
    cout <<"\n\tB:     "<< B_ <<" bulk modulus = compressibility of rock skeleton = K_dry (Pa)";
    cout <<"\n\tG:     "<< G_ <<" shear modulus = modulus of rigidity (Pa)";
    cout <<"\n\talpha: "<< alpha_ <<" Biot coefficient alpha = 1 - K_dry/K_grain (dimensionless ratio)";
    cout <<"\n\tTS:    "<< TS_ <<" tensile strength (Pa)";
    cout <<"\n\tUCS:   "<< UCS_ <<" unconfined compressive strength (Pa)";
    cout <<"\n\tpstar: "<< pstar_ <<" crushing pressure of the rock, see Fjaer et al. 08', p. 68 (Pa)";
    cout <<"\n\tmu:    "<< mu_ <<" friction coefficient (dimensionless)";
    cout <<"\n\tnu:    "<< nu_ <<" Poisson's ratio (dimensionless ratio)";
    cout <<"\n\tC:     "<< C_ <<" cohesive strength = inherent shear strength (Pa).\n";
 }




}  // end namespace csmp
