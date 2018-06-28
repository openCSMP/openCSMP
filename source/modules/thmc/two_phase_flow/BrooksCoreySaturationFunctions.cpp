#include "BrooksCoreySaturationFunctions.h"
#include "FlowFunctions.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class USER>
BrooksCoreySaturationFunctions<dim,USER>::BrooksCoreySaturationFunctions()
{
}


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation( TARGET_PLACEMENT& p ) const
 {
    // seff = (sw - swr) / (1 - swr - snr)
    double64 seff =  (p.Obtain(User()->key_sH2O) - p.Obtain(User()->key_srH2O)) /
           (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<1U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<2U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<3U,NODE>& ) const;

template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<3U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation( TARGET_PLACEMENT& p, double64 sw ) const
 {
    // seff = (sw - swr) / (1 - swr - snr)
    double64 seff = (sw - p.Obtain(User()->key_srH2O)) /
           (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::EffectiveSaturation( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 sw ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::SeffToSw( TARGET_PLACEMENT& p, double64 seff ) const
 {
    return seff * (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) + p.Obtain(User()->key_srCO2);
 }

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::SeffToSw( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::SeffToSw( FiniteElementPlacement<3U,NODE>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::SeffToSw( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::SeffToSw( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::SeffToSw( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krw( TARGET_PLACEMENT& p ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == static_cast<double64>(0.) )
      return EffectiveSaturation(p);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( EffectiveSaturation(p), 2. / p.Obtain(User()->key_bcp) + 3. );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krw( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krw( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krw( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krw( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krw( FiniteVolumePlacement<1U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krw( FiniteVolumePlacement<2U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteVolumePlacement<3U,NODE>& ) const;

template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krw( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


//Using prescribed sw value, instead of intepolated value.
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krw( TARGET_PLACEMENT& p, double64 sw ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == static_cast<double64>(0.) )
      return EffectiveSaturation(p,sw);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( EffectiveSaturation(p,sw), 2. / p.Obtain(User()->key_bcp) + 3. );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krw( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krw( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krw( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;



/// as in Helmig 97, eqn.2.57, p. 75
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krn( TARGET_PLACEMENT& p ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == static_cast<double64>(0.) )
      return 1. - EffectiveSaturation(p);

    const double64  seffn(1. - EffectiveSaturation(p));
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation(p), 2./ p.Obtain(User()->key_bcp) + 1.) );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krn( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krn( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krn( FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krn( FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krn( FiniteVolumePlacement<1U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krn( FiniteVolumePlacement<2U,NODE>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteVolumePlacement<3U,NODE>& ) const;

template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krn( FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


//Using prescribed sw value, instead of intepolated value.
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::krn( TARGET_PLACEMENT& p, double64 sw ) const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == static_cast<double64>(0.) )
      return 1. - EffectiveSaturation(p,sw);

    const double64  seffn(1. - EffectiveSaturation(p,sw));
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation(p,sw), 2./ p.Obtain(User()->key_bcp) + 1.) );
 }

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::krn( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::krn( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::krn( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 sw ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds( TARGET_PLACEMENT& p ) const
 {
    const double64 seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );
    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ p.Obtain(User()->key_bcp) + 3. ) *
           std::pow( EffectiveSaturation(p), 2./ p.Obtain(User()->key_bcp) + 2. ) * seff_mult;
 }

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds( TARGET_PLACEMENT& p ) const
 {
    const double64 seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );

    //  switch to linear relperm model if lambda = 0
    if ( p.Obtain(User()->key_bcp) == 0. ) return -seff_mult;

    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Obtain(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;

 }

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



/**
    bounded pressure of the non-wetting phase Pc = Pnw - Pw

    @attention note that pc must be defined over the full saturation range, pc is capped based on maximum dpcds
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::pc( TARGET_PLACEMENT& p ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Obtain(User()->key_bcp));
    const double64  entry_pressure(p.Obtain(User()->key_pd));

    // linear relperm model
    if ( bcp == 0. ) {
         if ( seff <= 0. ) return MaxCapillaryPressure();
         if ( seff >= 1. ) return entry_pressure;
         if ( entry_pressure == MaxCapillaryPressure() ) return entry_pressure;
         return entry_pressure + ( 1. - seff ) * ( MaxCapillaryPressure() - entry_pressure );
      }

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure == 0. ) return 0.;

   const double64 seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );

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
         return pcmax + ( Se_min - seff ) * MaxCapillaryPressureDerivative() / seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
}

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::pc( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::pc( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::pc( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::pc( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


// helper for numerical differentiation
double64 pcBC( double64 sw, double64 swr, double64 snr, double64 bcp, double64 pd )
  {
     const double64 pc_max(1.0e7), dpcds_max(1.0e6);
     const double64 seff( sw - swr / (1. - swr -snr) );
    
     if ( bcp == 0. ) {
          if ( seff <= 0. ) return pc_max;
          if ( seff >= 1. ) return pd;
          if ( pd == pc_max ) return pd;
          return pd + ( 1. - seff ) * ( pc_max - pd );
       }
   // for zero entry pressure capillary pressure always is zero
   if ( pd == 0. ) return 0.;

   const double64 seff_mult( 1./ (1. - swr - snr) );
   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = pd * pow( ( pd / ( bcp * dpcds_max/seff_mult ) ),  ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( pc_max / pd, -bcp ) ) {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double64 Se_min =  pow( ( pd / ( bcp * dpcds_max/seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * dpcds_max / seff_mult;
     }
   return pd * std::pow( seff, -1. / bcp );
}


/// dpcdS covers the full saturation range, dpcds is capped based on maximum dpcds
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds( TARGET_PLACEMENT& p ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Obtain(User()->key_bcp));
    const double64  entry_pressure(p.Obtain(User()->key_pd));
    const double64 seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );

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
    if ( seff <= Se_min ) return -MaxCapillaryPressureDerivative();

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;
}

template double64 BrooksCoreySaturationFunctions<1U,FlowFunctions>::dpcds( FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<2U,FlowFunctions>::dpcds( FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds( FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds( FiniteElementPlacement<3U,ELEMENT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::Sw( TARGET_PLACEMENT& p, double64 pc ) const
{
    const double64  seff(EffectiveSaturation(p));
    const double64  bcp(p.Obtain(User()->key_bcp));
    const double64  entry_pressure(p.Obtain(User()->key_pd));
    const double64  seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );

    // for linear relperm model
    if ( bcp == 0. )
      {
         // not a unique solution
         if ( MaxCapillaryPressure() == entry_pressure )
           return std::min( std::max( p.Obtain(User()->key_sH2O), p.Obtain(User()->key_srH2O) ), 1. - p.Obtain(User()->key_srCO2) );

         if ( pc >= MaxCapillaryPressure() ) return SeffToSw(p,0.);
         if ( pc <= entry_pressure ) return SeffToSw(p,1.);

         return SeffToSw(p,seff);
      }

   // for zero entry pressure capillary pressure always is zero, so the saturation is not unique
   if ( entry_pressure == 0. ) return p.Obtain(User()->key_sH2O);

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

        return SeffToSw(p,seff_dash);
     }

   const double64 seff_dash = std::pow( pc /entry_pressure, -bcp);
   return SeffToSw(p,seff_dash);
  
} // Sw(pc) inverse function

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::Sw( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::Sw( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::Sw( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::Sw( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dsdpc( TARGET_PLACEMENT& p, double64 pc ) const
{
    const double64  bcp(p.Obtain(User()->key_bcp));
    const double64  entry_pressure(p.Obtain(User()->key_pd));
    const double64  seff_mult( 1./ (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)) );

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

    const double64 Seff = (Sw( p, pc ) - p.Obtain(User()->key_srH2O)) * seff_mult;

    if( Seff<= Se_min ) return -1. / MaxCapillaryPressureDerivative();

    return -bcp / entry_pressure * pow( pc / entry_pressure, -bcp - 1. ) / seff_mult;

} // end dsdpc

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dsdpc( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dsdpc( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dsdpc( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dsdpc( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



// barebone functions
double64 kr_w( double64 seff, double64 bcp )  { return std::pow( seff, 2. / bcp + 3. ); }

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
{
   const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
   const double64 seff(EffectiveSaturation(p));
   const double64 bcp(p.Obtain(User()->key_bcp));

  if ( seff < 0.+h )
    return ( kr_w(seff + h, bcp) - kr_w(seff, bcp) ) / h * dSedSw;
  if ( seff > 1.-h )
    return ( kr_w(seff,bcp) - kr_w(seff - h, bcp) ) / h * dSedSw;

  return ( kr_w(seff+h,bcp) - kr_w( seff - h, bcp ) ) / (2. * h) * dSedSw;
}

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrwds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;



double64 kr_n( double64 seff, double64 bcp ) {
     const double64  seffn(1. - seff);
     return (seffn * seffn) * (1. - pow(  seff, 2./ bcp + 1.) );
 }

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
{
   const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
   const double64 seff(EffectiveSaturation(p));
   const double64 bcp(p.Obtain(User()->key_bcp));

  if ( seff < 0.+h )
    return ( kr_n( seff + h, bcp ) - kr_n( seff, bcp ) ) / h * dSedSw;
  if ( seff > 1.-h )
    return ( kr_n( seff, bcp ) - kr_n( seff - h, bcp ) ) / h * dSedSw;

  return ( kr_n( seff+h, bcp ) - kr_n( seff-h, bcp ) )/ (2. * h) * dSedSw;
}

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dkrnds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds_Numerical( TARGET_PLACEMENT& p, double64 h ) const
{
   const double64 sw(p.Obtain(User()->key_sH2O));
   const double64 swr(p.Obtain(User()->key_srH2O)), snr(p.Obtain(User()->key_srCO2));
   const double64 bcp(p.Obtain(User()->key_bcp));
   const double64 pd(p.Obtain(User()->key_pd));

  // pcBC( double64 sw, double64 swr, double64 snr, double64 bcp, double64 pd )
  if ( sw < 0.+h )
    return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw, swr, snr, bcp, pd ) ) / h;
      
  if ( sw > 1.-h )
    return ( pcBC( sw, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) ) / h;

  return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) )/ (2. * h);
}

template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds_Numerical( FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds_Numerical( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds_Numerical( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 BrooksCoreySaturationFunctions<3U,FlowFunctions>::dpcds_Numerical( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;





template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
void BrooksCoreySaturationFunctions<dim,USER>::Out( TARGET_PLACEMENT& p ) const
 {
    cout <<"\nBrooksCoreySaturationFunctions<"<< dim << ">::Out: Additional properties: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n       capillary entry pressure, pd: "<< p.Obtain(User()->key_pd);
    cout <<"\n      Brooks-Corey lambda parameter: "<< p.Obtain(User()->key_bcp);
    cout <<"\n             MAX_CAPILLARY_PRESSURE: "<< MaxCapillaryPressure() << endl << endl;
    cout <<"\n  MAX_CAPILLARY_PRESSURE_DERIVATIVE: "<< MaxCapillaryPressureDerivative() << endl << endl;

 } // end Out
 
 template void BrooksCoreySaturationFunctions<3U,FlowFunctions>::Out( FiniteElementPlacement<3U,ELEMENT>& ) const;
template void BrooksCoreySaturationFunctions<3U,FlowFunctions>::Out( FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template void BrooksCoreySaturationFunctions<3U,FlowFunctions>::Out( FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template void BrooksCoreySaturationFunctions<3U,FlowFunctions>::Out( FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;



template class BrooksCoreySaturationFunctions<1U,FlowFunctions>;
template class BrooksCoreySaturationFunctions<2U,FlowFunctions>;
template class BrooksCoreySaturationFunctions<3U,FlowFunctions>;

} // end namespace csmp


























