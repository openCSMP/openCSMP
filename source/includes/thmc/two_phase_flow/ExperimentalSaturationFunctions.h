#ifndef EXPERIMENTAL_SATURATION_FUNCTIONS_H
#define EXPERIMENTAL_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "CubicSpline.h"

namespace csmp {
  
  /// relperm- and pc-s relations from data files are interpolated with cubic splines
  /// @author SKM @date 14/1/2007
  /// @author MM  @date  7/9/2018

  
template<size_t dim, template<size_t> class USER>
class ExperimentalSaturationFunctions {
  public:
    /// TODO: reads capillary preessure and saturation functions from "rocktype" file
    explicit ExperimentalSaturationFunctions( const char* filename="rocktypes" );
  
    /// Effective Saturation function
    template<class TARGET_PLACEMENT>
    double64 EffectiveSaturation( const TARGET_PLACEMENT& ) const ;
  
    /// Effective Saturation function fro saturation S
    template<class TARGET_PLACEMENT>
    double64 EffectiveSaturation_at( const TARGET_PLACEMENT& , double64 ) const ;
  
    /// The Capillary pressure Eq. (2) from Skjaeveland et al. 2000
    template<class TARGET_PLACEMENT>
    double64 pc( const TARGET_PLACEMENT& );
  
    /// The first derivative of Capillary pressure.
    template<class TARGET_PLACEMENT>
    double64 dpcds( const TARGET_PLACEMENT& )  const ;
  
    /// The water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 krw( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 krw_at( const TARGET_PLACEMENT&, double64 sw ) const ;
  
    /// The CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 krn( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 krn_at( const TARGET_PLACEMENT&, double64 sw ) const ;
  
    /// The first relative of water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 dkrwds( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrwds_at( const TARGET_PLACEMENT&, double64 sw ) const ;
  
    /// The first relative of CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 dkrnds( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrnds_at( const TARGET_PLACEMENT&, double64 sw ) const ;

  
    /// Numerical derivatives of first derivatives of relative permeability of water and CO2
    template<class TARGET_PLACEMENT>
    double64 dkrwds_Numerical( const TARGET_PLACEMENT& p, double64 h ) const ;

    template<class TARGET_PLACEMENT>
    double64 dkrwds_at_Numerical( const TARGET_PLACEMENT& p, double64 sw, double64 h ) const ;

    template<class TARGET_PLACEMENT>
    double64 dkrnds_Numerical( const TARGET_PLACEMENT& p, double64 h ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 dkrnds_at_Numerical( const TARGET_PLACEMENT& p, double64 sw, double64 h ) const ;

    
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    void ConstructRTs( const char* rt_file_name );
  
    /// converts the floating-point rocktype identifier into a positive integer 0..254
    template<class TARGET_PLACEMENT>
    size_t RockType( const TARGET_PLACEMENT& ) const;
  
    std::vector<csmp::CubicSpline> kr1_, kr2_, pc_;
};
  
 
// INLINE FUNCTIONS
  
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline size_t ExperimentalSaturationFunctions<dim,USER>::RockType( const TARGET_PLACEMENT& p ) const
 {
    const size_t rocktype = static_cast<size_t>(p.Obtain( User()->key_RRT ));
    assert( rocktype < 254 );
    return rocktype;
 }



template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation( const TARGET_PLACEMENT& p ) const
{
  double64 seff = (p.Obtain(User()->key_sH2O) - p.Obtain(User()->key_srH2O)) /
                  (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
  
  return std::min( std::max( seff, 0. ), 1. );
}
  

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at( const TARGET_PLACEMENT& p, double64 sw ) const
  {
    double64 seff =  (sw - p.Obtain(User()->key_srH2O)) /
    (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }
  


/**
    TODO: capillary pressure also exists outside of the effective saturation range
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::pc( const TARGET_PLACEMENT& p )
{
    return pc_[ RockType(p) ].Value( EffectiveSaturation(p) );
}




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::dpcds( const TARGET_PLACEMENT& p) const
  {
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    return pc_[ RockType(p) ].Derivative( EffectiveSaturation(p) ) * seff_mult;
}
  


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::krw( const TARGET_PLACEMENT& p ) const
{
   return kr1_[ RockType(p) ].Value( EffectiveSaturation(p) );
}
  


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::krw_at( const TARGET_PLACEMENT& p , double64 S) const
{
   double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
   return kr1_[ RockType(p) ].Value( seff );
}


/**
   
   For calculating the oil relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
   The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
   the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
   
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::krn( const TARGET_PLACEMENT& p) const
  {
    return kr2_[ RockType(p) ].Value( EffectiveSaturation(p) );
  }

  
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::krn_at( const TARGET_PLACEMENT& p, double64 S) const
{
    double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    return kr2_[ RockType(p) ].Value( seff );
}
  


/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds( const TARGET_PLACEMENT& p ) const
  {
    const double64 srH2O  = p.Obtain(User()->key_srH2O);
    const double64 srCO2  = p.Obtain(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    return kr1_[ RockType(p) ].Derivative( EffectiveSaturation(p) )*seff_mult;
}


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_at( const TARGET_PLACEMENT& p , double64 S) const
  {
    const double64 srH2O  = p.Obtain(User()->key_srH2O);
    const double64 srCO2  = p.Obtain(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));

    return kr1_[ RockType(p) ].Derivative( seff )*seff_mult;
}




template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds( const TARGET_PLACEMENT& p) const
  {
    const double64 srH2O  = p.Obtain(User()->key_srH2O);
    const double64 srCO2  = p.Obtain(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    return kr2_[ RockType(p) ].Derivative( EffectiveSaturation(p) ) * seff_mult;
}



/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
inline double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_at( const TARGET_PLACEMENT& p, double64 S ) const
  {
    const double64 srH2O  = p.Obtain(User()->key_srH2O);
    const double64 srCO2  = p.Obtain(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    const double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return kr2_[ RockType(p) ].Derivative( seff ) * seff_mult;
    
}


} // end namespace csmp

#endif /* CSMP_EXPERIMENTAL_SATURATION_FUNCTIONS_H */
