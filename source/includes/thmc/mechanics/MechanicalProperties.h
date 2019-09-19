//
//  MechanicalProperties.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 6/23/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MECHANICAL_PROPERTIES_H
#define CSMP_MECHANICAL_PROPERTIES_H

#include "CSMP_definitions.h"

namespace csmp {
  
/// returns arctan(mu) in degrees, @test OK
double64  frictionAngle( double64 friction_coeff );
  
/// modified Griffith criterion from Brace, 1960
double64  tensileStrengthFromUCS_Griffith( double64 UCS, double64 fric_coeff );


/**
    Regional properties shared by multiple elements, but have variations
    that can be calculated knowing other local properties discretized
    on the elements.
 
    Example: Young's modulus can be calculated from the value stored
    here with corrections based on the deviation of the local porosity
    from that for which E was determined.
 
    Such calculations must be informed by petrophysical correlations
    and the results might have to be perturbed by some randon variable
    to reflect natural variations.
 
*/
struct MechanicalProperties {
   /// granite properties as default
   MechanicalProperties();
  
   /// from property values supplied as a vector
   explicit MechanicalProperties( const std::vector<double64>& properties );
  
   /// partial initialisation
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
  
   // TODO: here we need the correlations that relate the properties to variations in the parent rock type
   // derived properties for fault rocks: dilatation, compaction etc.
  
   void Out() const;
};

}  // end namespace csmp

#endif /* CSMP_MECHANICAL_PROPERTIES_H */
