//
//  StressInvariants.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 10/25/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include <limits>
#include "StressInvariants.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

/**
    factors the upper triagonal of the symmetric tensor in terms of s, t, and theta

    upper diagonal matrix:
 
    sigmaxx = Cartesian_stress_(0,0);
    sigmayy = Cartesian_stress_(1,1);
    sigmaxy = Cartesian_stress_(0,1);
*/
StressInvariants::StressInvariants( const TensorVariable<2U>& ts )
 : sqrt3_(sqrt(3.)),
   sqrt32_(sqrt(3./2.)),
   degrees_to_radians_((CSMP_PI)/180.),
   s_((ts(0,0) + ts(1,1)) / sqrt3_),
   t_(sqrt((ts(0,0)-ts(1,1))*(ts(0,0)-ts(1,1)) + 3.*(ts(0,1) * ts(0,1)))/sqrt3_),
   theta_(LodeAngle(ts))
 {
 } // end (2D constructor)


/**
    factors the upper triagonal of the symmetric tensor in terms of s, t, and theta

    upper diagonal matrix:
 
    sigmaxx = Cartesian_stress_(0,0);
    sigmayy = Cartesian_stress_(1,1);
    sigmazz = Cartesian_stress_(2,2);
    sigmaxy = Cartesian_stress_(0,1);
    sigmaxz = Cartesian_stress_(0,2);
    sigmayz = Cartesian_stress_(1,2);
*/
StressInvariants::StressInvariants( const TensorVariable<3U>& ts )
 : sqrt3_(sqrt(3.)),
   sqrt32_(sqrt(3./2.)),
   degrees_to_radians_((CSMP_PI)/180.),
   s_((ts(0,0) + ts(1,1) + ts(2,2)) / sqrt3_),
   t_(sqrt((ts(0,0)-ts(1,1))*(ts(0,0)-ts(1,1)) +
      (ts(1,1)-ts(2,2))*(ts(1,1)-ts(2,2)) +
      (ts(2,2)-ts(0,0))*(ts(2,2)-ts(0,0)) +
      6.*(ts(0,1) * ts(0,1)) + 6.*(ts(1,2) * ts(1,2)) + 6.*(ts(0,2) * ts(0,2)))/sqrt3_),
   theta_(LodeAngle(ts))
 {
 } // end (3D constructor)


/**
    theta angle in Mohr-Coulomb diagram that at which shear stresses are maximized

    @attention original angle(theta), is in radians (Smith & Griffiths, eq. 6.3. p. 227
 
*/
double StressInvariants::LodeAngle( const TensorVariable<3U>& ts ) const
 {
    const double sx  = (2. * ts(0,0) - ts(1,1) - ts(2,2)) / 3.;
    const double sy  = (2. * ts(1,1) - ts(2,2) - ts(0,0)) / 3.;
    const double sz  = (2. * ts(2,2) - ts(0,0) - ts(1,1)) / 3.;
    double J3(sx * sy * sz);
    J3 -= sx * (ts(1,2) * ts(1,2));
    J3 -= sy * (ts(0,2) * ts(0,2));
    J3 -= sz * (ts(0,1) * ts(0,1));
    J3 += 2. * ts(0,1) * ts(0,2) * ts(1,2);
   
    // theta
    return 1./3. * asin( (-3.* sqrt(6.) * J3) / (t_ * t_ * t_) );
   
 } // end LodeAngle

/// 2D version
double StressInvariants::LodeAngle( const TensorVariable<2U>& ts ) const
 {
    const double sx  = (2. * ts(0,0) - ts(1,1)) / 2.;
    const double sy  = (2. * ts(1,1) - ts(2,2)) / 2.;
    const double J3(sx * sy * ts(0,1));
   
    // theta
    return 1./3. * asin( (-3.* sqrt(6.) * J3) / (t_ * t_ * t_) );
   
 } // end LodeAngle


/// mean stress = average of principal stresses
double StressInvariants::MeanStress() const
 {
    return s_ / sqrt3_;
 }
  
/// maximum of the deviatoric stresses = nonisostatic stresses
double StressInvariants::DeviatoricStress() const
 {
    return t_ * sqrt32_;
 }


/**
    sigma1 = maximum principal stress (Pa)
    as obtained from the modified deviatoric and mean stresses, see Smith & Griffith, p. 234
*/
double StressInvariants::MaximumPrincipalStress1() const
 {
    return MeanStress() + 2./3. * DeviatoricStress() * sin( theta_ - (2.*CSMP_PI)/3. );
 }

/// sigma2
double StressInvariants::IntermediatePrincipalStress2() const
 {
    return MeanStress() + 2./3. * DeviatoricStress() * sin( theta_ );
 }

/// sigma3
double StressInvariants::LeastPrincipalStress3() const
 {
    return MeanStress() + 2./3. * DeviatoricStress() * sin( theta_ + (2.*CSMP_PI)/3. );
 }



/**
    hexagonal K(theta), see Zienkiewitz, volume 2, chapter 4.5.1, p. 233
    
    @attention no material parameters (other than the friction angle)
    are contained in this assessment. However, these need to enter failure
    calculations that are handled in other classes, such as BrittleFailure.
*/
double StressInvariants::MohrCoulombYieldEnvelope( double alpha ) const
 {
    double K_OfTheta = -sin(alpha * degrees_to_radians_);
    K_OfTheta *= sin(theta_);
    K_OfTheta /= sqrt3_;
    K_OfTheta += cos(theta_);
    K_OfTheta /= sqrt3_;
   
    return K_OfTheta;
 }


/**
   Computing K_OfTheta = 1/G_OfTheta ( modified Mohr-Coulomb envelope with smooth boundaries ):
   Zienkiewitz, Finite element method for solid and structural mechanics,
   Chapter 4. Inelastic and non-linear materials
   4.11. Non-uniqueness and localization in elasto-plastic deformations

    @attention no material parameters (other than the friction angle)
    are contained in this assessment. However, these need to enter failure
    calculations that are handled in other classes, such as BrittleFailure.
*/
double StressInvariants::SmoothMohrCoulombYieldEnvelope( double alpha ) const
 {
    double K(sin(alpha * degrees_to_radians_));
    K = (3. - K) / (3. + K);
    const double G_OfTheta = (2.*K) / ((1+K) - sin(3. * theta_) * (1-K));
    const double K_OfTheta = 1.0/ G_OfTheta;

    return K_OfTheta;
 }

} // end csmp
