//
//  StressInvariants.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 10/25/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef STRESS_INVARIANTS_H
#define STRESS_INVARIANTS_H

#include "TensorVariable.h"

namespace csmp {

/// standard decompositions of stress tensor, see Smith & Griffith
class StressInvariants {
  public:
    /// factors the upper triagonal of the symmetric tensor in terms of s, t, and theta
    explicit StressInvariants( const TensorVariable<2U>& );
    explicit StressInvariants( const TensorVariable<3U>& );
  
    /// theta angle in Mohr-Coulomb diagram that at which shear stresses are maximized
    double LodeAngle( const TensorVariable<2U>& ) const;
    double LodeAngle( const TensorVariable<3U>& ) const;
  
    /// mean stress = average of principal stresses
    double MeanStress() const;
  
    /// maximum of the deviatoric stresses = nonisostatic stresses
    double DeviatoricStress() const;
  
    /// sigma1 = maximum principal stress (Pa)
    double /* sigma1 */ MaximumPrincipalStress1() const;
    /// sigma2
    double /* sigma2 */ IntermediatePrincipalStress2() const;
    /// sigma3
    double /* sigma3 */ LeastPrincipalStress3() const;

    /// K(theta), see Zienkiewitz, volume 2, chapter 4.5.1, p. 233, alpha is in degrees
    double MohrCoulombYieldEnvelope( double friction_angle_alpha ) const;

    /// K(theta), see Zienkiewitz, volume 2, chapter 4.5.11, alpha is in degrees
    double SmoothMohrCoulombYieldEnvelope( double friction_angle_alpha ) const;

  private:
    const double sqrt3_, sqrt32_,
                   degrees_to_radians_,
                   s_, t_, theta_;
};

/**
  @class StressInvariants StressInvariants.h "main_library/constitutive_relationships/mechanics/"

  @author S.K. Matthaei

  @section motivation Motivation

  To extract the standard engineering invariants s, t, and theta from the upper 
  diagonal of a (symmetric diagonally dominant) Cartesian stress tensor,
  calculate the principal stresses, and evaluate a series of failure criteria.
  
  @section applicability Applicability

  CSMP geomechanics models.
  
  @section collaborations Collaborations


  @section implementation Implementation

  Follows Zienkiewitz volume 2 and Smith and Griffith.

  @section examples Application Examples

  To get an Index for a variable of interest in order to perform
  a computation, do the following:

  @test

*/

} // end csmp

#endif /* defined(STRESS_INVARIANTS_H) */
