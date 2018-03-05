#include "BrooksCoreySaturationFunctions.h"
#include "FlowFunctions.h"

using namespace std;

namespace csmp {


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation( TARGET_PLACEMENT& p ) const
 {
    // seff = (sw - swr) / (1 - swr - snr)
    return (p.Interpolate(User()->key_sw) - p.Interpolate(User()->key_swr)) /
           (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr));
 }




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krw( TARGET_PLACEMENT& p ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Interpolate(User()->key_bcp) == static_cast<double64>(0.) )
      return EffectiveSaturation(p);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( EffectiveSaturation(p), 2. / p.Interpolate(User()->key_bcp) + 3. );
 }



/// as in Helmig 97, eqn.2.57, p. 75
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krn( TARGET_PLACEMENT& p ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Interpolate(User()->key_bcp) == static_cast<double64>(0.) )
      return 1. - EffectiveSaturation(p);

    const double64  seffn(1. - EffectiveSaturation(p));
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation(p), 2./ p.Interpolate(User()->key_bcp) + 1.) );
 }



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds( TARGET_PLACEMENT& p ) const
 {
    const double64 seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );
    //  switch to linear relperm model if lambda = 0
    if ( p.Interpolate(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ p.Interpolate(User()->key_bcp) + 3. ) * std::pow( TwoPhaseModel<dim>::seff_, 2./ p.Interpolate(User()->key_bcp) + 2. ) * seff_mult;
 }



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds( TARGET_PLACEMENT& p ) const
 {
    const double64 seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );

    //  switch to linear relperm model if lambda = 0
    if ( p.Interpolate(User()->key_bcp) == 0. ) return -seff_mult;

    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Interpolate(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;

 }


/**
    bounded pressure of the non-wetting phase Pc = Pnw - Pw

    @attention note that pc must be defined over the full saturation range, pc is capped based on maximum dpcds
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::pc( TARGET_PLACEMENT& p ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Interpolate(User()->key_bcp));
    const double64  entry_pressure(p.Interpolate(User()->key_pd));

    // linear relperm model
    if ( bcp == 0. ) {
         if ( seff <= 0. ) return MaxCapillaryPressure();
         if ( seff >= 1. ) return entry_pressure;
         if ( entry_pressure == MaxCapillaryPressure() ) return entry_pressure;
         return entry_pressure + ( 1. - seff ) * ( MaxCapillaryPressure() - entry_pressure );
      }

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure == 0. ) return 0.;

   const double64 seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative()/seff_mult ) ),
                                           ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( MaxCapillaryPressure() / entry_pressure, -bcp ) )
     {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double64 Se_min =  pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative()/seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - TwoPhaseModel<dim>::seff_ ) * MaxCapillaryPressureDerivative()/seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
}





/// dpcdS covers the full saturation range, dpcds is capped based on maximum dpcds
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds( TARGET_PLACEMENT& p ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Interpolate(User()->key_bcp));
    const double64  entry_pressure(p.Interpolate(User()->key_pd));
    const double64 seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );

    // linear relperm model
     if ( bcp == 0. ) {
          if ( ( seff < 0. ) || ( seff > 1. ) ) return 0.;
		      if ( entry_pressure == MaxCapillaryPressure() ) return 0.;

          return -( MaxCapillaryPressure() - entry_pressure ) * seff_mult;
       }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double64 Se_min = pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative()/seff_mult ) ), ( bcp / (1. + bcp) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( TwoPhaseModel<dim>::seff_ <= Se_min ) return -MaxCapillaryPressureDerivative();

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;
}




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::Sw( TARGET_PLACEMENT& p, double64 pc ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Interpolate(User()->key_bcp));
    const double64  entry_pressure(p.Interpolate(User()->key_pd));
    const double64  seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );

    // for linear relperm model
    if ( bcp == 0. )
      {
         // not a unique solution
         if ( MaxCapillaryPressure() == entry_pressure )
           return std::min( std::max( p.Interpolate(User()->key_sw), p.Interpolate(User()->key_swr) ), 1. - p.Interpolate(User()->key_snr) );

         if ( pc >= MaxCapillaryPressure() ) return TwoPhaseModel<dim>::SeffToSw(0.);
         if ( pc <= entry_pressure ) return SeffToSw(1.);

         return SeffToSw(seff);
      }

   // for zero entry pressure capillary pressure always is zero, so the saturation is not unique
   if ( entry_pressure == 0. ) return p.Interpolate(User()->key_sw);

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative()/seff_mult ) ),
                                          ( -1. / ( 1. + bcp ) ) );

   if ( pc >= MaxCapillaryPressure() ) {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
        const double64 Se_min =  pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative()/seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
        // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
        const double64 seff_dash = Se_min - ( pc - pcmax ) / MaxCapillaryPressureDerivative() * seff_mult;

        return TwoPhaseModel<dim>::SeffToSw(seff_dash);
     }

   const double64 seff_dash = std::pow( pc /entry_pressure, -bcp);
   return TwoPhaseModel<dim>::SeffToSw(seff_dash);
  
} // Sw(pc) inverse function




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dsdpc( TARGET_PLACEMENT& p, double64 pc ) const
{
    const double64  bcp(p.Interpolate(User()->key_bcp));
    const double64  entry_pressure(p.Interpolate(User()->key_pd));
    const double64  seff_mult( 1./ (1. - p.Interpolate(User()->key_swr) - p.Interpolate(User()->key_snr)) );

    // linear relperm model
    if ( bcp == 0. )
      {
         if ( MaxCapillaryPressure() == entry_pressure ) return 0.;
         return -1. / (MaxCapillaryPressure() - entry_pressure) / seff_mult;
      }
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double64 Se_min =  pow( ( entry_pressure / ( bcp * MaxCapillaryPressureDerivative() / seff_mult ) ),
                                ( bcp / ( 1. + bcp ) ) );

    const double64 Seff = (Sw( p, pc ) - p.Interpolate(User()->key_swr)) * seff_mult;

    if( Seff<= Se_min ) return -1. / MaxCapillaryPressureDerivative();

    return -bcp / entry_pressure * pow( pc / entry_pressure, -bcp - 1. ) / seff_mult;

} // end dsdpc









template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
void BrooksCoreySaturationFunctions<dim,USER>::Out( TARGET_PLACEMENT& p ) const
 {
    cout <<"\nBrooksCoreySaturationFunctions<"<< dim << ">::Out: Additional properties: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n       capillary entry pressure, pd: "<< p.Interpolate(User()->key_pd);
    cout <<"\n      Brooks-Corey lambda parameter: "<< p.Interpolate(User()->key_bcp);
    cout <<"\n             MAX_CAPILLARY_PRESSURE: "<< MaxCapillaryPressure() << endl << endl;
    cout <<"\n  MAX_CAPILLARY_PRESSURE_DERIVATIVE: "<< MaxCapillaryPressureDerivative() << endl << endl;

 } // end Out
 
 


//template class BrooksCoreySaturationFunctions<1U>;
//template class BrooksCoreySaturationFunctions<2U>;
template class BrooksCoreySaturationFunctions<3U,FlowFunctions>;

} // end namespace csmp


























