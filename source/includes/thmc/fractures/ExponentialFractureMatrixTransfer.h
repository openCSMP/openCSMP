// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EXPONENTIAL_FRACTURE_MATRIX_TRANSFER_H
#define EXPONENTIAL_FRACTURE_MATRIX_TRANSFER_H

#include <cmath>

namespace csmp {

/// documentation in 'SKM_linear_diffusion_based_transfer_function.mw'
class ExponentialFractureMatrixTransfer{
  public:
    ExponentialFractureMatrixTransfer( double swi    = 0.05,
                                       double swc    = 0.05,
                                       double snr    = 0.14,
                                       double lambda = 2.5,
                                       double pd     = 3000.,
                                       double km     = 1.E-15,
                                       double phi    = 0.25,
                                       double muw    = 1.E-3,
                                       double mun    = 5.E-3 );
                                       
    double CurrentTransferRate( double t, double sw_matrix_avg, double sw_frac ) const;

    void UpdateMatrixParameters( double swi    = 0.05,
                                 double swc    = 0.05,
                                 double snr    = 0.14,
                                 double lambda = 2.5,
                                 double pd     = 3000.,
                                 double km     = 1.E-15,
                                 double phi    = 0.25,
                                 double muw    = 1.E-3,
                                 double mun    = 5.E-3 );
    
  private:
    ExponentialFractureMatrixTransfer();

    double CapillaryDiffusivity( double sw_matrix_avg ) const;
    
    double CurrentSaturationGradientAtFractureMatrixInterface( double t, 
                                                               double kappa_BC ) const;
                                                               
    double CapillaryPressureDerivative( double sw_matrix_avg ) const;
  
  private:
    double swi_,    ///< initial saturation in matrix
           swc_,    ///< connate water saturation
           snr_,    ///< residual non-wetting phase saturation
           lambda_, ///< Brooks-Corey parameter (water-wet system)
           pd_,     ///< capillary entry pressure
           km_,     ///< matrix permeability
           phi_,    ///< matrix porosity
           muw_,    ///< viscosity water
           mun_;    ///< viscosity non-wetting phase
};


// definition of inline methods

/// ds/dt|_fracture-matrix interface
inline double ExponentialFractureMatrixTransfer::CurrentSaturationGradientAtFractureMatrixInterface( 
                                                                       double t, 
                                                                       double kappa_BC ) const
 {
    return -(1 - swi_) / std::sqrt(0.3141592654e1) / std::sqrt(kappa_BC * t);

 } // end CurrentSaturationGradientAtFractureMatrixInterface

                                                               
/// dpc/ds (sw_matrix_average)                                                               
inline double ExponentialFractureMatrixTransfer::CapillaryPressureDerivative( double sw_matrix_avg ) const
 {
    const double t1(1. / lambda_);
    return -pd_ * std::pow(sw_matrix_avg, -t1) * t1 / sw_matrix_avg;

 } // end CapillaryPressureDerivative


} // end csmp

#endif

