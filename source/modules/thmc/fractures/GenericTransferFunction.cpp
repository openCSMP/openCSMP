// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "GenericTransferFunction.h"

namespace csmp {

GenericTransferFunction::GenericTransferFunction( double phim,
                                                  double pd, double BC_param,
                                                  double muw, double mun )
 : phim_(phim),
   pd_(pd), BC_param_(BC_param), 
   muw_(muw), mun_(mun) 
 {
    assert( phim_ > 0. and phim_ <= 1. );
    assert( BC_param_ > 0.2 and BC_param_ <= 6. );
    assert( muw_ > 0. );
    assert( mun_ > 0. );
 }
 
 
GenericTransferFunction::~GenericTransferFunction()
 {
 }
 

/**

Computes transfer rate assuming that fracture saturation of the invading phase
is fixed=1. This approximation has to be seen in the context of the Af(sw)
scaling which is subsequently applied to the transfer rate in the block.
    
@return flux (volume * saturation per second) of the defending phase
*/
double GenericTransferFunction::Transfer( double sw_block,  // average saturation in gridblock
                                            double sw_initial,
                                            double sw_farfield,
                                            double km,
                                            double block_radius ) const
 {
    // compute beta parameter but without gravity term ((km delta rho g) / (phi_m block_radius))
    // SKM interpretation: lambda_overbar * km * dpcds * ((pc(sw_block) - pc_frac(sw=1)) / block_radius)
//    double beta = km * LambdaOverbar( sw_block ) * (Pc(sw_farfield) / block_radius);
    double beta = km * (Pc(sw_farfield) / block_radius);
    
    // Vermuelen correction for early time behaviour (see Zimmermann)
    double B = (std::fabs(sw_farfield-sw_initial) + std::fabs(sw_block - sw_initial)) / 
                  std::max( std::numeric_limits<double>::epsilon(), 2. * (sw_block - sw_initial) );
    
    // composing the function: beta * B * phi_m * (sw_farfield - sw_block)
    return beta * B * phim_ * (sw_farfield - sw_block);
 }
 



double  GenericTransferFunction::LambdaOverbar( double sw_eff ) const
{
    if ( sw_eff <= 0. ) sw_eff = 0.;
    if ( sw_eff >= 1. ) sw_eff = 1.;

    double krw = std::pow( sw_eff, (2. + 3. * BC_param_) / BC_param_ ); 
   	double krn = ((1. - sw_eff) * (1. - sw_eff)) *  std::pow( 1. - sw_eff, (2. + BC_param_) / BC_param_ ); 
   	double lw  = krw / muw_;
   	double ln  = krn / mun_;
   	
   	return (lw * ln) / (lw + ln);
}

} // end csmp
