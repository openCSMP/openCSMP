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
    double64 krw_at( const TARGET_PLACEMENT& , double64 S) const ;
  
    /// The CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 krn( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 krn_at( const TARGET_PLACEMENT& , double64 S) const ;
  
    /// The first relative of water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 dkrwds( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrwds_at( const TARGET_PLACEMENT& , double64 S) const ;
  
    /// The first relative of CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 dkrnds( const TARGET_PLACEMENT& ) const ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrnds_at( const TARGET_PLACEMENT& , double64 S) const ;

  
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
  
    void ConstructRTs(const char* rt_file_name);
  
    /// converts the floating-point rocktype identifier into a positive integer 0..254
    template<class TARGET_PLACEMENT>
    size_t RockType( const TARGET_PLACEMENT& ) const;
  
    std::vector<csmp::CubicSpline> kr1_, kr2_, pc_;
};
  
  
} // end namespace csmp

#endif /* CSMP_EXPERIMENTAL_SATURATION_FUNCTIONS_H */
