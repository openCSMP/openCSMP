/*
 *  ExponentialFractureMatrixTransfer.cpp
 *  csmp_examples
 *
 *  Created by Stephan Matthai on 8/10/10.
 *  Copyright 2010 SKM. All rights reserved.
 *
 */

#include "ExponentialFractureMatrixTransfer.h"

#include <limits>

namespace csmp {

ExponentialFractureMatrixTransfer::ExponentialFractureMatrixTransfer( double swi,
                                                                      double swc,
                                                                      double snr,
                                                                      double lambda,
                                                                      double pd,
                                                                      double km,
                                                                      double phi,
                                                                      double muw,
                                                                      double mun )
 : swi_(swi), swc_(swc), snr_(snr), 
   lambda_(lambda), pd_(pd), km_(km), phi_(phi),
   muw_(muw), mun_(mun)
 {
 }
 
/// T(t,sw_avg)                                                                             
double ExponentialFractureMatrixTransfer::CurrentTransferRate( double t, 
                                                               double sw_matrix_avg,
                                                               double sw_frac ) const
 {
   if( t == 0 )
    return std::numeric_limits<double>::max();

   const double t1 = CapillaryDiffusivity( sw_matrix_avg );
   const double t2 = CurrentSaturationGradientAtFractureMatrixInterface( t, t1 );
   const double t4 = CapillaryPressureDerivative( sw_matrix_avg );
   
   return (t1 * t2 * t4);

 } // end CurrentTransferRate
  
  
    
double ExponentialFractureMatrixTransfer::CapillaryDiffusivity( double sw_matrix_avg ) const
 {
    const double t5 = (0. >= ((sw_matrix_avg - swc_) / (1. - swc_ - snr_)) ? 0. : ((sw_matrix_avg - swc_) / (1. - swc_ - snr_)));
    const double t6 = (1. <= t5 ? 1. : t5);
    const double t8 = std::pow(1. - t6, 2.);
    const double t10 = 1. / lambda_;
    const double t12 = std::pow(t6, ((2. + lambda_) * t10));
    const double t14 = t8 * (1. - t12);
    const double t15 = 1. / mun_;
    const double t19 = std::pow(t6, ((2. + 3. * lambda_) * t10));
    const double t22 = 1. / muw_;
    
    return t14 * t15 * t19 * t22 / (t14 * t15 + t19 * t22) * km_ / phi_;

 } // end CapillaryDiffusivity

void  ExponentialFractureMatrixTransfer::UpdateMatrixParameters( double swi, double swc, double snr, double lambda,
                                                                 double pd, double km, double phi,  double muw, double mun )
{
  swi_ = swi; swc_ = swc; snr_ = snr; lambda_ = lambda;
  pd_ = pd; km_ = km; phi_ = phi; muw_ = muw; mun_ = mun;
}


} // end csmp

