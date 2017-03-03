#ifndef GENERIC_TRANSFER_FUNCTION_H
#define GENERIC_TRANSFER_FUNCTION_H

#include "CSMP_definitions.h"

namespace csmp {

/// GTF, see doi 0.2118/121244-MS
class GenericTransferFunction {
  public:
    GenericTransferFunction( double64 phim,
                             double64 pd, double64 BC_param,
                             double64 muw, double64 mun );
                             
    ~GenericTransferFunction();
    
    double64   Transfer( double64 sw_block,  
                         double64 sw_initial,
                         double64 sw_farfield,
                         double64 km,
                         double64 block_radius ) const;
    
    // Brooks-Corey parameters
    double64  krw_BC( double64 sw_eff ) const;
    double64  krn_BC( double64 sw_eff ) const;
    
    double64  Lambda_w_BC( double64 sw_eff ) const;
    double64  Lambda_n_BC( double64 sw_eff ) const;
    
    double64  LambdaTotal( double64 sw_eff ) const;
    double64  LambdaOverbar( double64 sw_eff ) const;
                                                    
    double64  Pc( double64 sw ) const;
    double64  dPcds( double64 sw ) const;
  
  private:
    double64 phim_, pd_, BC_param_, muw_, mun_;
};









// inline functions
inline double64  GenericTransferFunction::krw_BC( double64 sw_eff ) const
 {
    if ( sw_eff <= 0. ) return static_cast<double64>(0.);
    if ( sw_eff >= 1. ) return static_cast<double64>(1.);

    return std::pow( sw_eff, (2. + 3. * BC_param_) / BC_param_ ); 
 }
 
 
inline double64  GenericTransferFunction::krn_BC( double64 sw_eff ) const
 {
    if ( sw_eff <= 0. ) return static_cast<double64>(1.);
    if ( sw_eff >= 1. ) return static_cast<double64>(0.);
    
   	return ((1. - sw_eff) * (1. - sw_eff)) * (1. - std::pow( sw_eff, (2. + BC_param_) / BC_param_ )); 
 }


inline double64  GenericTransferFunction::Lambda_w_BC( double64 sw_eff ) const
 {
    return krw_BC( sw_eff ) / muw_;
 }


inline double64  GenericTransferFunction::Lambda_n_BC( double64 sw_eff ) const
 {
    return krn_BC( sw_eff ) / mun_;
 }


inline double64  GenericTransferFunction::LambdaTotal( double64 sw_eff ) const
 {
    return krw_BC( sw_eff ) / muw_ + krn_BC( sw_eff ) / mun_;
 }



inline double64 GenericTransferFunction::Pc( double64 sw ) const
 {
    assert( sw >= 0. and sw <= 1. );
    return pd_ * std::pow( sw, -1. / BC_param_ );
 }


inline double64 GenericTransferFunction::dPcds( double64 sw ) const
 {
    assert( sw >= 0. and sw <= 1. );

    // computing the unlimited capillary pressure derivative
    return pd_ * std::pow( sw, -1. / BC_param_ ) / (sw * BC_param_);
 }

} // end csmp

#endif
