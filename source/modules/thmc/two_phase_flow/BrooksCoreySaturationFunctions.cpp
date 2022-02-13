#include "CSMP_mathUtilities.h"
#include "BrooksCoreySaturationFunctions.h"
#include "FlowFunctionsModule.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
/**
   The default constructor of Brooks Corey Saturation Functions class
*/
template<uint32_t dim, template<uint32_t> class USER> BrooksCoreySaturationFunctions<dim,USER>::BrooksCoreySaturationFunctions()
  {
  }




/// get seff at the element barycentre
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation( Element<dim>* const e ) const
  {
    const double sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert( sH2O >= 0. );
    assert( sH2O <= 1. );
    const double seff = (sH2O - e->Read(User()->key_srH2O)) /
                          (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return max( min( seff, 1. ), 0. );
  }





/// get seff from the supplied saturation value
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double seff = (sw - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    return max( min( seff, 1. ), 0. );
  }




template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::pc( Element<dim>* const e ) const
  {
    assert( !isnan(e->Read(User()->key_bcp) ) );
    const double  seff(EffectiveSaturation(e));
    const double  bcp(e->Read(User()->key_bcp));
    const double  entry_pressure(e->Read(User()->key_pd));

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure == 0. ) return 0.;

    // linear relperm model
    if ( bcp == 0. ) {
         if ( seff <= 0. ) return MaxCapillaryPressure();
         if ( seff >= 1. ) return entry_pressure;
         if ( entry_pressure == MaxCapillaryPressure() ) return entry_pressure;
         return entry_pressure + ( 1. - seff ) * ( MaxCapillaryPressure() - entry_pressure );
      }
 
   // Brooks-Corey model
   const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                           ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( MaxCapillaryPressure() / entry_pressure, -bcp ) )
     {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double Se_min =  pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * max_derivative_ / seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
  
} // end pc
  
 
 
 
 
 
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::pc_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );
    assert( !isnan(e->Read(User()->key_bcp) ) );
    assert( !isnan(e->Read(User()->key_pd) ) );

    const double  seff(EffectiveSaturation_at(e,sw));
    const double  bcp(e->Read(User()->key_bcp));
    const double  entry_pressure(e->Read(User()->key_pd));

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure == 0. ) return 0.;

    // linear relperm model
    if ( bcp == 0. ) {
         if ( seff <= 0. ) return MaxCapillaryPressure();
         if ( seff >= 1. ) return entry_pressure;
         if ( entry_pressure == MaxCapillaryPressure() ) return entry_pressure;
         return entry_pressure + ( 1. - seff ) * ( MaxCapillaryPressure() - entry_pressure );
      }
 
   // Brooks-Corey model
   const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                           ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( MaxCapillaryPressure() / entry_pressure, -bcp ) )
     {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double Se_min =  pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * max_derivative_ / seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
  
} // end pc

  
  
  
  
  
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dpcds( Element<dim>* const e ) const
  {
    assert( !isnan(e->Read(User()->key_bcp) ) );
    assert( !isnan(e->Read(User()->key_pd) ) );
    const double  seff(EffectiveSaturation(e));
    const double  bcp(e->Read(User()->key_bcp));
    const double  entry_pressure(e->Read(User()->key_pd));
    const double  seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    // linear relperm model
     if ( bcp == 0. ) {
          if ( ( seff < 0. ) || ( seff > 1. ) ) return 0.;
		      if ( entry_pressure == MaxCapillaryPressure() ) return 0.;

          return -( MaxCapillaryPressure() - entry_pressure ) * seff_mult;
       }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double Se_min = pow( ( entry_pressure / ( bcp * max_derivative_/seff_mult ) ), ( bcp / (1. + bcp) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( seff <= Se_min ) return -max_derivative_;

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;  
  
} // end dpcds





template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dpcds_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );
    assert( !isnan(e->Read(User()->key_bcp) ) );
    assert( !isnan(e->Read(User()->key_pd) ) );

    const double  seff(EffectiveSaturation_at(e,sw));
    const double  bcp(e->Read(User()->key_bcp));
    const double  entry_pressure(e->Read(User()->key_pd));
    const double  seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    // linear relperm model
     if ( bcp == 0. ) {
          if ( ( seff < 0. ) || ( seff > 1. ) ) return 0.;
          if ( entry_pressure == MaxCapillaryPressure() ) return 0.;

          return -( MaxCapillaryPressure() - entry_pressure ) * seff_mult;
       }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double Se_min = pow( ( entry_pressure / ( bcp * max_derivative_/seff_mult ) ), ( bcp / (1. + bcp) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( seff <= Se_min ) return -max_derivative_;

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;
  
} // end dpcds
 
 

  
  
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::krw( Element<dim>* const e ) const
 {
    const double bcp(e->Read(User()->key_bcp));
    assert( !isnan(bcp) );
   
    //  switch to linear relperm model if lambda = 0
    if ( bcp == static_cast<double>(0.) )
      return min( 1., max( 0., EffectiveSaturation(e) ) );

    // pm2 = lambda, the Brooks-Corey parameter
    const double seff = min( 1., max( 0., EffectiveSaturation(e) ) );
    return std::pow( seff, 2. / bcp + 3. );
}
  
  
   
 
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::krw_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );
    const double bcp(e->Read(User()->key_bcp));
    assert( !isnan(bcp) );

    //  switch to linear relperm model if lambda = 0
    if ( bcp == static_cast<double>(0.) )
      return min( 1., max( 0., EffectiveSaturation_at(e,sw) ) );
    
    // pm2 = lambda, the Brooks-Corey parameter
    const double seff = min( 1., max( 0., EffectiveSaturation_at(e,sw) ) );
    return std::pow( seff, 2. / bcp + 3. );
}



  
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::krn( Element<dim>* const e ) const
  {
    const double bcp(e->Read(User()->key_bcp));
    assert( !isnan(bcp) );

    //  switch to linear relperm model if lambda = 0
    if ( bcp == static_cast<double>(0.) )
      return min( 1., max( 0., 1. - EffectiveSaturation(e) ) );

    const double  seffn = min( 1., max( 0., 1. - EffectiveSaturation(e) ) );
    
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation(e), 2./ bcp + 1.) );
}




template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::krn_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );
    const double bcp(e->Read(User()->key_bcp));
    assert( !isnan(bcp) );

    //  switch to linear relperm model if lambda = 0
    if ( bcp == static_cast<double>(0.) )
      return min( 1., max( 0., 1. - EffectiveSaturation_at(e,sw) ) );

    const double  seffn = min( 1., max( 0., 1. - EffectiveSaturation_at(e,sw) ) );

    return (seffn * seffn) * (1. - pow(  EffectiveSaturation_at(e,sw), 2./ bcp + 1.) );
}




  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dkrwds( Element<dim>* const e ) const
  {
    assert( !isnan(e->Read(User()->key_bcp) ) );
    const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ e->Read(User()->key_bcp) + 3. ) *
           std::pow( EffectiveSaturation(e), 2./ e->Read(User()->key_bcp) + 2. ) * seff_mult;
}

  
  

  
  
/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dkrwds_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );
    assert( !isnan(e->Read(User()->key_bcp) ) );

    const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ e->Read(User()->key_bcp) + 3. ) *
           std::pow( EffectiveSaturation_at(e,sw), 2./ e->Read(User()->key_bcp) + 2. ) * seff_mult;
}

  
  
   
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dkrnds( Element<dim>* const e ) const
  {
    const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    //  switch to linear relperm model if lambda = 0
    assert( !isnan(e->Read(User()->key_bcp) ) );
    if ( e->Read(User()->key_bcp) == 0. ) return -seff_mult;

    const double  seff(EffectiveSaturation(e));
    const double  bcp(e->Read(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;
    
}
  
  
  

  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dkrnds_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    //  switch to linear relperm model if lambda = 0
    assert( !isnan(e->Read(User()->key_bcp) ) );
    if ( e->Read(User()->key_bcp) == 0. ) return -seff_mult;

    const double  seff(EffectiveSaturation_at(e,sw));
    const double  bcp(e->Read(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;
    
}

  
 
 
  
  // Numerical derivative added
  
  template<uint32_t dim, template<uint32_t> class USER>
  double BrooksCoreySaturationFunctions<dim,USER>::dkrwds_Numerical( Element<dim>* const e, double h ) const
  {
    const double dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
 
  
  template<uint32_t dim, template<uint32_t> class USER>
  double BrooksCoreySaturationFunctions<dim,USER>::dkrwds_at_Numerical( Element<dim>* const e, double sw, double h ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }



 
   
  
  
  template<uint32_t dim, template<uint32_t> class USER>
  double BrooksCoreySaturationFunctions<dim,USER>::dkrnds_Numerical( Element<dim>* const e, double h ) const
  {
    const double dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  template<uint32_t dim, template<uint32_t> class USER>
  double BrooksCoreySaturationFunctions<dim,USER>::dkrnds_at_Numerical( Element<dim>* const e, double sw, double h ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  // helper function for numerical differentiation
  double pcBC( double sw, double swr, double snr, double bcp, double pd )
  {
     const double pc_max(1.0e7), dpcds_max(1.0e6);
     const double seff( sw - swr / (1. - swr -snr) );
    
     if ( bcp == 0. ) {
          if ( seff <= 0. ) return pc_max;
          if ( seff >= 1. ) return pd;
          if ( pd == pc_max ) return pd;
          return pd + ( 1. - seff ) * ( pc_max - pd );
       }
   // for zero entry pressure capillary pressure always is zero
   if ( pd == 0. ) return 0.;

   const double seff_mult( 1./ (1. - swr - snr) );
   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double pcmax = pd * pow( ( pd / ( bcp * dpcds_max/seff_mult ) ),  ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( pc_max / pd, -bcp ) ) {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double Se_min =  pow( ( pd / ( bcp * dpcds_max/seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * dpcds_max / seff_mult;
     }
   return pd * std::pow( seff, -1. / bcp );
  }  



template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dpcds_Numerical( Element<dim>* const e, double h ) const
  {
   const double sw(e->PropertyValueAtBaryCenter(User()->key_sH2O));
   const double swr(e->Read(User()->key_srH2O)), snr(e->Read(User()->key_srCO2));
   const double bcp(e->Read(User()->key_bcp));
   const double pd(e->Read(User()->key_pd));

  // pcBC( double sw, double swr, double snr, double bcp, double pd )
  if ( sw < 0.+h )
    return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw, swr, snr, bcp, pd ) ) / h;
      
  if ( sw > 1.-h )
    return ( pcBC( sw, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) ) / h;

  return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) )/ (2. * h);
}



template<uint32_t dim, template<uint32_t> class USER>
double BrooksCoreySaturationFunctions<dim,USER>::dpcds_at_Numerical( Element<dim>* const e, double sw, double h ) const
 {
   assert( sw >= 0. );
   assert( sw <= 1. );

   const double swr(e->Read(User()->key_srH2O)), snr(e->Read(User()->key_srCO2));
   const double bcp(e->Read(User()->key_bcp));
   const double pd(e->Read(User()->key_pd));

  // pcBC( double sw, double swr, double snr, double bcp, double pd )
  if ( sw < 0.+h )
    return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw, swr, snr, bcp, pd ) ) / h;
  
  if ( sw > 1.-h )
    return ( pcBC( sw, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) ) / h;

  return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) )/ (2. * h);
}

  

template class BrooksCoreySaturationFunctions<1U,FlowFunctionsModule1>;
template class BrooksCoreySaturationFunctions<2U,FlowFunctionsModule1>;
template class BrooksCoreySaturationFunctions<3U,FlowFunctionsModule1>;

template class BrooksCoreySaturationFunctions<1U,FlowFunctionsModule4>;
template class BrooksCoreySaturationFunctions<2U,FlowFunctionsModule4>;
template class BrooksCoreySaturationFunctions<3U,FlowFunctionsModule4>;


} // csmp
