// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DeformationInducedPermeability.cpp
//  CSMP_ReservoirSimulator
//
//

#include <cassert>
#include "CSMP_mathUtilities.h"
#include "DeformationInducedPermeability.h"

namespace csmp {

DeformationInducedPermeability::DeformationInducedPermeability( double max_aperture )
 : max_aperture_(max_aperture)
 {
 }


DeformationInducedPermeability::~DeformationInducedPermeability()
 {
 }







// inline functions

/**
    calculates aperture from fracture intensity (n / m), unless m=0. In this case, method returns 0.
    Also, if the average aperture exceeds the maximum aperture specified when the class was constructed,
    the former is returned.
*/
double  DeformationInducedPermeability::AverageAperture( double phi, double fracIntensity ) const
{
  assert( phi > 0. );
  assert( fracIntensity >= 0. );
  if ( fracIntensity <= std::numeric_limits<double>::epsilon() ) return 0.;
  
  // if calculated aperture would be greater than maximum aperture, the latter is returned
  if ( phi / fracIntensity > max_aperture_ ) return max_aperture_;
  
  return phi / fracIntensity;
}



/// from multiple fractures which all have the same aperture
double DeformationInducedPermeability::keff_UniformAperture( double phi, double matrix_k, double frac_intensity ) const
{
  assert( matrix_k > 0. );
  assert( phi > 0. );
  assert( frac_intensity >= 0. );
  // total porosity is used to create uniform aperture fractures
  const double a_avg = AverageAperture(phi, frac_intensity);
  const double a_avg_squared(a_avg * a_avg);
  //           keff                                    
  const double keff = matrix_k + frac_intensity * a_avg * (a_avg_squared / 12.);
  assert( keff > 0. );
  return keff;
}



/// Kozeny-Carman permeability at high porosity from mean breccia fragment size and tortuosity
double DeformationInducedPermeability::kKC_BrecciaFromFragSize( double phi, double matrix_k, double meanFragSize, double tau ) const 
{
  assert( phi > 0. );
  assert( matrix_k > 0. );
  assert( meanFragSize > 0. );
  assert( meanFragSize < 10. ); // fragments must be smaller than 10 meters
  assert( tau >= 1. );
  const double phi2(phi * phi);
  const double frag_sz2(meanFragSize * meanFragSize);
  const double t6(1. - phi);
  const double t7 = t6 * t6;
  // Kozeny-Carman version for spherical fragments with mean size
  const double k_breccia = matrix_k + (1./0.72e2) * ((phi2 * phi * frag_sz2) / (tau * t7));
  assert( k_breccia > 0. );
  return k_breccia;
}



/**  k_FaultCoreBreccia

     function generates k values as function of max-breccia k and distance from fault core.
     In core and plateau region, the thickness of which can be controlled by 
     @param plateau_factor
     k = kmax, but at a greater distance k descreases expontentially up to zero at 
     periphery
*/
double DeformationInducedPermeability::k_FaultCoreBreccia( double halfThickness, double x,
                                                           double plateau_factor, double kmax ) const
{
  assert( halfThickness > 0.01 );
  assert( x <= halfThickness );
  assert( plateau_factor > 0. );
  assert( kmax > 0. );
  const double t4 = std::exp(-plateau_factor * halfThickness / std::max(x,0.001) );
  const double result_k( kmax * (1. - t4) );
  assert( result_k > 0. );
  return result_k;
}




/**  k_DamageZoneFracture

     Generates from kmax, k values that follow a Gaussian distribution with 
     increasing distance from the fault zone. 
     This distribution is skewed towards the center of the fault 
     by a factor that can be controlled with
     @param skewness_factor which becomes unnoticeable at values < 0.01,
     recommended values ~5-50, the larger the values, the smaller the 
     zone at which k values near kmax are realized. 

*/
double DeformationInducedPermeability::k_DamageZoneFracture( double halfThickness, double x,
                                                             double skewness_factor, double kmax ) const
{
  assert( halfThickness > 0.01 );
  assert( x <= halfThickness );
  assert( skewness_factor > 0.001 );
  assert( kmax > 0. );
  const double t3 = 1. / halfThickness;
  const double t5 = std::exp(-skewness_factor * x * t3);
  const double t8 = std::sin(x * PI * t3);
  double result_k = kmax * skewness_factor * t5 * t8;
// SKM fix
result_k = std::max( result_k, 1.0e-18 );
  assert( result_k > 0. );
  return result_k;
}




} // end csmp
