#include "SaturationFunction.h"
#include "FlowFunctions.h"

using namespace std;

namespace csmp {

/**
 
The effective saturation is always computed w.r.t. the wetting phase. 
(compare with Helmig 1997, p. 71 and note that phase 1 is the wetting 
phase; note also that Initialize() must be called first to get input 
parameters like residual saturations).  
*/
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::EffectiveSaturation( double64 sw ) const
 {
    return std::min( std::max( (sw - User()->swr) / (1. - User()->swr - User()->snr), 0. ), 1. );
    
 } // end EffectiveSaturation


/// return wetting-phase saturation based on effective saturation
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::SeffToSw( double64 seff ) const
 {
    return  seff * (1. - User()->swr - User()->snr) + User()->swr;

 }

 

/// relative permeabilities
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::krw( double64 sw ) const
 {
    cout <<"\nSaturationFunction<dim,USER>::krw (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::krw_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }



template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::krn( double64 sw ) const
 {
    cout <<"\nSaturationFunction<"<< dim <<">::krn (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::krn_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }

/// derivatives of relative permeabilities
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dkrwds( double64 sw, bool evaluate_numerically ) const
 {
    if ( evaluate_numerically ) return dkrwds_Numerical(sw);
    // linear model
    return EffectiveSaturation(sw);
 }



template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dkrnds( double64 sw, bool evaluate_numerically ) const
 {
    if ( evaluate_numerically ) return dkrnds_Numerical(sw);
    // linear model
    return 1. - EffectiveSaturation(sw);
 }


/// capillary pressure
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::pc( double64 sw ) const
 {
    cout <<"\nSaturationFunction<"<< dim <<">::pc_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::pc_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }



/// inverse capillary pressure function
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::Sw( double64 ) const
 {
    cout <<"\nSaturationFunction<"<< dim <<">::SwFromPc(base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::SwFromPc: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }



/// derivatives of inverse capillary pressure function
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dpcds( double64 sw, bool evaluate_numerically ) const
 {
    cout <<"\nSaturationFunction<"<<  dim <<">::dpcds (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::dpcdsw_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }



/// derivatives of inverse capillary pressure function
template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dsdpc( double64 pc, bool evaluate_numerically ) const
 {
    cout <<"\nSaturationFunction<"<<  dim <<">::dpcdsw_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("SaturationFunction<dim,USER>::dpcdsw_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }




/// Numerical derivatives

template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dkrwds_Numerical( double64 sw, double64 h ) const
{
  //const double64 dSedSw( 1.0/ (1.0 - Read( User()->swr_key ) - Read( User()->snr_key ) ) );
  const double64 dSedSw(EffectiveSaturation(1.+User()->swr));
  const double64 seff(EffectiveSaturation(sw));

  if( seff < 0.+h )
      return (krw( seff + h ) - krw( seff  ) ) / h * dSedSw;
  if( seff > 1.-h )
      return ( krw( seff ) - krw( seff - h  ) ) / h * dSedSw;

  return ( krw( seff+h  ) - krw( seff - h ) ) / (2.0 * h) * dSedSw;
}



template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dkrnds_Numerical( double64 sw, double64 h ) const
{
  const double64 dSedSw(EffectiveSaturation(1.+User()->swr));
  const double64 seff(EffectiveSaturation(sw));

  if( seff < 0.+h )
      return ( krn( seff + h ) - krn( seff ) ) / h * dSedSw;
  if( seff > 1.-h )
      return ( krn( seff ) - krn( seff - h ) ) / h * dSedSw;

  return ( krn( seff+h ) - krn( seff-h ) )/ (2.0 * h) * dSedSw;

}




template<size_t dim, template<size_t> class USER>
double64 SaturationFunction<dim,USER>::dpcds_Numerical( double64 sw, double64 h ) const
{
  double64 dpcds;

  if( sw < 0.+h )
      dpcds = ( pc( sw + h ) - pc( sw ) ) / h;
  else if( sw > 1.-h )
      dpcds = ( pc( sw ) - pc( sw - h ) ) / h;
  else
      dpcds = ( pc( sw + h ) - pc( sw - h ) )/ (2.0 * h);

  //if( dpcds < -MAX_CAPILLARY_PRESSURE_SLOPE_ )
  //    return -MAX_CAPILLARY_PRESSURE_SLOPE_;

  return dpcds;

}




/// General accessory functions



template<size_t dim, template<size_t> class USER>
void SaturationFunction<dim,USER>::Out() const
 { 
    cout <<"\nSaturationFunction<" << dim <<">::Out: ";
    cout <<"\nelement properties:";

    cout <<"\n                    maximum capillary pressure: "<< pc_max_;
    cout <<"\n         maximum capillary pressure derivative: "<< dpcds_max_;
    // former data members
    /*
    cout <<"\n       irreducible saturation of wetting phase: "<< swr_;
    cout <<"\n   irreducible saturation of non-wetting phase: "<< snr_;
    cout <<"\n                                  permeability: "<< k_;
    cout <<"\nnode properties:";
    cout <<"\ninterfacial tension (surface tension of fluid): "<< ift_;
    cout <<"\n                      saturation wetting phase: "<< sat_;
    cout <<"\n            effective saturation wetting phase: "<< seff_;
    cout <<"\n                                         pc(2): "<< pc_Phase();
    cout <<"\n                                      dpcdS(2): "<< dpcds_Phase();
    cout << endl << endl;
    */
}


template class SaturationFunction<1U,FlowFunctions>;
template class SaturationFunction<2U,FlowFunctions>;
template class SaturationFunction<3U,FlowFunctions>;


} // end namespace csp





