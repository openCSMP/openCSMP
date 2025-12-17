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
double  frictionAngle( double friction_coeff );
  
/// modified Griffith criterion from Brace, 1960
double  tensileStrengthFromUCS_Griffith( double UCS, double fric_coeff );


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
   explicit MechanicalProperties( const std::vector<double>& properties );
  
   /// partial initialisation
   MechanicalProperties( double E, double K, double nu,
                         double tensile_strength, double fric_coeff );
  
   double  E_;          ///< Young's modulus
   double  B_;          ///< bulk modulus = compressibility of rock skeleton = K_dry
   double  G_;          ///< shear modulus = modulus of rigidity
   double  alpha_;      ///< Biot coefficient alpha = 1 - K_dry/K_grain
   double  TS_;         ///< tensile strength
   double  UCS_;        ///< unconfined compressive strength
   double  pstar_;      ///< crushing pressure of the rock, see Fjaer et al. 08', p. 68
   double  mu_;         ///< friction coefficient
   double  nu_;         ///< Poisson's ratio
   double  C_;          ///< cohesive strength = inherent shear strength
  
   // TODO: here we need correlations that relate properties to variations in the parent rock type
   // derived properties for fault rocks: dilatation, compaction etc.
  
   void Out() const;
};

}  // end namespace csmp

#endif /* CSMP_MECHANICAL_PROPERTIES_H */
