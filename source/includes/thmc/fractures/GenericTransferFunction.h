// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GENERIC_TRANSFER_FUNCTION_H
#define GENERIC_TRANSFER_FUNCTION_H

#include "CSMP_definitions.h"

namespace csmp {

/// GTF, see doi 0.2118/121244-MS
class GenericTransferFunction {
  public:
    GenericTransferFunction( double phim,
                             double pd, double BC_param,
                             double muw, double mun );
                             
    ~GenericTransferFunction();
    
    double   Transfer( double sw_block,  
                         double sw_initial,
                         double sw_farfield,
                         double km,
                         double block_radius ) const;
    
    // Brooks-Corey parameters
    double  krw_BC( double sw_eff ) const;
    double  krn_BC( double sw_eff ) const;
    
    double  Lambda_w_BC( double sw_eff ) const;
    double  Lambda_n_BC( double sw_eff ) const;
    
    double  LambdaTotal( double sw_eff ) const;
    double  LambdaOverbar( double sw_eff ) const;
                                                    
    double  Pc( double sw ) const;
    double  dPcds( double sw ) const;
  
  private:
    double phim_, pd_, BC_param_, muw_, mun_;
};









// inline functions
inline double  GenericTransferFunction::krw_BC( double sw_eff ) const
 {
    if ( sw_eff <= 0. ) return static_cast<double>(0.);
    if ( sw_eff >= 1. ) return static_cast<double>(1.);

    return std::pow( sw_eff, (2. + 3. * BC_param_) / BC_param_ ); 
 }
 
 
inline double  GenericTransferFunction::krn_BC( double sw_eff ) const
 {
    if ( sw_eff <= 0. ) return static_cast<double>(1.);
    if ( sw_eff >= 1. ) return static_cast<double>(0.);
    
   	return ((1. - sw_eff) * (1. - sw_eff)) * (1. - std::pow( sw_eff, (2. + BC_param_) / BC_param_ )); 
 }


inline double  GenericTransferFunction::Lambda_w_BC( double sw_eff ) const
 {
    return krw_BC( sw_eff ) / muw_;
 }


inline double  GenericTransferFunction::Lambda_n_BC( double sw_eff ) const
 {
    return krn_BC( sw_eff ) / mun_;
 }


inline double  GenericTransferFunction::LambdaTotal( double sw_eff ) const
 {
    return krw_BC( sw_eff ) / muw_ + krn_BC( sw_eff ) / mun_;
 }



inline double GenericTransferFunction::Pc( double sw ) const
 {
    assert( sw >= 0. and sw <= 1. );
    return pd_ * std::pow( sw, -1. / BC_param_ );
 }


inline double GenericTransferFunction::dPcds( double sw ) const
 {
    assert( sw >= 0. and sw <= 1. );

    // computing the unlimited capillary pressure derivative
    return pd_ * std::pow( sw, -1. / BC_param_ ) / (sw * BC_param_);
 }

} // end csmp

#endif
