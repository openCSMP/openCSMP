#include "CSMP_mathUtilities.h"
#include "BrooksCoreySaturationFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"


using namespace std;

namespace csmp {
  
/**
   The default constructor of Brooks Corey Saturation Functions class
*/
template<size_t dim, template<size_t> class USER> BrooksCoreySaturationFunctions<dim,USER>::BrooksCoreySaturationFunctions()
  {
  }




/// get seff at the element barycentre
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation( const Element<dim>* const e ) const
  {
    const double64 sH2O = e->PropertyValueAtBaryCenter( User()->key_sH2O );
  
    return (sH2O - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
  }





/// get seff from the supplied saturation value
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::EffectiveSaturation_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

     return (sw - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
  }




template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::pc( const Element<dim>* const e ) const
  {
    const double64  seff(EffectiveSaturation(e));
    const double64  bcp(e->Read(User()->key_bcp));
    const double64  entry_pressure(e->Read(User()->key_pd));

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
   const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                           ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( MaxCapillaryPressure() / entry_pressure, -bcp ) )
     {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double64 Se_min =  pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * max_derivative_ / seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
  
} // end pc
  
  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::pc_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64  seff(EffectiveSaturation_at(e,sw));
    const double64  bcp(e->Read(User()->key_bcp));
    const double64  entry_pressure(e->Read(User()->key_pd));

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
   const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure * pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                           ( -1. / ( 1. + bcp ) ) );

   if ( seff <= std::pow( MaxCapillaryPressure() / entry_pressure, -bcp ) )
     {
        // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
         const double64 Se_min =  pow( ( entry_pressure / ( bcp * max_derivative_ / seff_mult ) ),
                                  ( bcp / ( 1. + bcp ) ) );
         // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
         return pcmax + ( Se_min - seff ) * max_derivative_ / seff_mult;
     }

   return entry_pressure * std::pow( seff, -1. / bcp );
  
} // end pc

  
  
  
  
  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds( const Element<dim>* const e ) const
  {
    const double64  seff(EffectiveSaturation(e));
    const double64  bcp(e->Read(User()->key_bcp));
    const double64  entry_pressure(e->Read(User()->key_pd));
    const double64  seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    // linear relperm model
     if ( bcp == 0. ) {
          if ( ( seff < 0. ) || ( seff > 1. ) ) return 0.;
		      if ( entry_pressure == MaxCapillaryPressure() ) return 0.;

          return -( MaxCapillaryPressure() - entry_pressure ) * seff_mult;
       }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double64 Se_min = pow( ( entry_pressure / ( bcp * max_derivative_/seff_mult ) ), ( bcp / (1. + bcp) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( seff <= Se_min ) return -max_derivative_;

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;  
  
} // end dpcds





template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64  seff(EffectiveSaturation_at(e,sw));
    const double64  bcp(e->Read(User()->key_bcp));
    const double64  entry_pressure(e->Read(User()->key_pd));
    const double64  seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    // linear relperm model
     if ( bcp == 0. ) {
          if ( ( seff < 0. ) || ( seff > 1. ) ) return 0.;
          if ( entry_pressure == MaxCapillaryPressure() ) return 0.;

          return -( MaxCapillaryPressure() - entry_pressure ) * seff_mult;
       }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure == 0. ) return 0.;

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double64 Se_min = pow( ( entry_pressure / ( bcp * max_derivative_/seff_mult ) ), ( bcp / (1. + bcp) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( seff <= Se_min ) return -max_derivative_;

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure * std::pow( seff, -1. - 1. / bcp ) / bcp * seff_mult;
  
} // end dpcds
 
 

  
  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::krw( const Element<dim>* const e ) const
  {
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == static_cast<double64>(0.) )
      return EffectiveSaturation(e);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( EffectiveSaturation(e), 2. / e->Read(User()->key_bcp) + 3. );  
  
}
  
  
   
 
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::krw_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == static_cast<double64>(0.) )
      return EffectiveSaturation_at(e,sw);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( EffectiveSaturation_at(e,sw), 2. / e->Read(User()->key_bcp) + 3. );
}



  
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::krn( const Element<dim>* const e ) const
  {
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == static_cast<double64>(0.) )
      return 1. - EffectiveSaturation(e);

    const double64  seffn(1. - EffectiveSaturation(e));
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation(e), 2./ e->Read(User()->key_bcp) + 1.) );
}




template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::krn_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == static_cast<double64>(0.) )
      return 1. - EffectiveSaturation_at(e,sw);

    const double64  seffn(1. - EffectiveSaturation_at(e,sw));
    return (seffn * seffn) * (1. - pow(  EffectiveSaturation_at(e,sw), 2./ e->Read(User()->key_bcp) + 1.) );
}




  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds( const Element<dim>* const e ) const
  {
    const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ e->Read(User()->key_bcp) + 3. ) *
           std::pow( EffectiveSaturation(e), 2./ e->Read(User()->key_bcp) + 2. ) * seff_mult;
}

  
  

  
  
/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );
    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return seff_mult;

    return ( 2./ e->Read(User()->key_bcp) + 3. ) *
           std::pow( EffectiveSaturation_at(e,sw), 2./ e->Read(User()->key_bcp) + 2. ) * seff_mult;
}

  
  
   
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds( const Element<dim>* const e ) const
  {
    const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return -seff_mult;

    const double64  seff(EffectiveSaturation(e));
    const double64  bcp(e->Read(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;
    
}
  
  
  

  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds_at( const Element<dim>* const e, double64 sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64 seff_mult( 1./ (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)) );

    //  switch to linear relperm model if lambda = 0
    if ( e->Read(User()->key_bcp) == 0. ) return -seff_mult;

    const double64  seff(EffectiveSaturation_at(e,sw));
    const double64  bcp(e->Read(User()->key_bcp));
    return 2.*(seff - 1.) * (1. + pow(seff, 2./bcp) * ( 1./bcp + 1./2. - seff * (1./bcp + 3./2.))) * seff_mult;
    
}

  
 
 
  
  // Numerical derivative added
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
 
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctions<dim,USER>::dkrwds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at(e,seff + h) - this->krw_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(e,seff) - this->krw_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(e,seff + h) - this->krw_at(e,seff - h) ) / (2. * h) * dSedSw;
  }



 
   
  
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 BrooksCoreySaturationFunctions<dim,USER>::dkrnds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(e,sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at(e,seff + h) - this->krn_at(e,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(e,seff) - this->krn_at(e,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(e,seff + h) - this->krn_at(e,seff - h) ) / (2. * h) * dSedSw;
  }
 
  
  
  // helper function for numerical differentiation
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



template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds_Numerical( const Element<dim>* const e, double64 h ) const
  {
   const double64 sw(e->PropertyValueAtBaryCenter(User()->key_sH2O));
   const double64 swr(e->Read(User()->key_srH2O)), snr(e->Read(User()->key_srCO2));
   const double64 bcp(e->Read(User()->key_bcp));
   const double64 pd(e->Read(User()->key_pd));

  // pcBC( double64 sw, double64 swr, double64 snr, double64 bcp, double64 pd )
  if ( sw < 0.+h )
    return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw, swr, snr, bcp, pd ) ) / h;
      
  if ( sw > 1.-h )
    return ( pcBC( sw, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) ) / h;

  return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) )/ (2. * h);
}



template<size_t dim, template<size_t> class USER>
double64 BrooksCoreySaturationFunctions<dim,USER>::dpcds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
 {
   assert( sw >= 0. );
   assert( sw <= 1. );

   const double64 swr(e->Read(User()->key_srH2O)), snr(e->Read(User()->key_srCO2));
   const double64 bcp(e->Read(User()->key_bcp));
   const double64 pd(e->Read(User()->key_pd));

  // pcBC( double64 sw, double64 swr, double64 snr, double64 bcp, double64 pd )
  if ( sw < 0.+h )
    return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw, swr, snr, bcp, pd ) ) / h;
  
  if ( sw > 1.-h )
    return ( pcBC( sw, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) ) / h;

  return ( pcBC( sw+h, swr, snr, bcp, pd ) - pcBC( sw-h, swr, snr, bcp, pd ) )/ (2. * h);
}

  

template class BrooksCoreySaturationFunctions<1U,CO2H2O_FunctionsModule0>;
template class BrooksCoreySaturationFunctions<2U,CO2H2O_FunctionsModule0>;
template class BrooksCoreySaturationFunctions<3U,CO2H2O_FunctionsModule0>;


} // csmp
