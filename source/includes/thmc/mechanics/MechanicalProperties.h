//
//  MechanicalProperties.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 6/23/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#ifndef MECHANICAL_PROPERTIES_H
#define MECHANICAL_PROPERTIES_H

#include "CSMP_definitions.h"

namespace csmp {
  
/// returns arctan(mu) in degrees, @test OK
double64  frictionAngle( double64 friction_coeff );
  
/// modified Griffith criterion from Brace, 1960
double64  tensileStrengthFromUCS_Griffith( double64 UCS, double64 fric_coeff );


struct MechanicalProperties {
   MechanicalProperties();
   MechanicalProperties( double64 E, double64 K, double64 nu,
                         double64 tensile_strength, double64 fric_coeff );
  
   double64  E;          ///< Young's modulus
   double64  B;          ///< bulk modulus = compressibility of rock skeleton = K_dry
   double64  G;          ///< shear modulus = modulus of rigidity
   double64  alpha;      ///< Biot coefficient alpha = 1 - K_dry/K_grain
   double64  TS;         ///< tensile strength
   double64  UCS;        ///< unconfined compressive strength
   double64  pstar;      ///< crushing pressure of the rock, see Fjaer et al. 08', p. 68
   double64  mu;         ///< friction coefficient
   double64  nu;         ///< Poisson's ratio
   double64  C;          ///< cohesive strength = inherent shear strength
  
   void Out(std::ostream& os) const;
};

}  // end namespace csmp

#endif
