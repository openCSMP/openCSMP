#ifndef __H2OPropertiesLookUpTable_h__
#define __H2OPropertiesLookUpTable_h__

#include "CSMP_definitions.h"
#include "LookUpStorage.h"
#include "steam4.h"

namespace csmp {

template<typename fT>
class H2OPropertiesLookUpTable {

public:
    H2OPropertiesLookUpTable( bool use_lookup_tables=true );
    ~H2OPropertiesLookUpTable();
    H2OPropertiesLookUpTable( const H2OPropertiesLookUpTable<fT>& lookup );
    H2OPropertiesLookUpTable<fT>& operator=( const H2OPropertiesLookUpTable<fT>& lookup );
    // 2-phase properties
    fT LiquidSaturation( fT p, fT t );
    fT VaporSaturation( fT p, fT t );
    // for given p
    fT SaturationTemperatureFromP( fT p );
    fT LiquidDensityFromP( fT p );
    fT VaporDensityFromP( fT p );
    fT LiquidEnthalpyFromP( fT p );
    fT VaporEnthalpyFromP( fT p );
    fT LiquidHeatCapacityFromP( fT p );
    fT VaporHeatCapacityFromP( fT p );
    fT LiquidViscosityFromP( fT p );
    fT VaporViscosityFromP( fT p );
    fT LiquidCompressibilityFromP( fT p );
    fT VaporCompressibilityFromP( fT p );
    fT LiquidThermalExpansivityFromP( fT p );
    fT VaporThermalExpansivityFromP( fT p );
    fT LiquidPressureDensityDerivativeFromP( fT p );
    fT VaporPressureDensityDerivativeFromP( fT p );
    fT LiquidPressureTemperatureDerivativeFromP( fT p );
    fT VaporPressureTemperatureDerivativeFromP( fT p );
    // for given t
    fT SaturationPressureFromT( fT t );
    fT LiquidDensityFromT( fT t );
    fT VaporDensityFromT( fT t );
    fT LiquidEnthalpyFromT( fT t );
    fT VaporEnthalpyFromT( fT t );
    fT LiquidHeatCapacityFromT( fT t );
    fT VaporHeatCapacityFromT( fT t );
    fT LiquidViscosityFromT( fT t );
    fT VaporViscosityFromT( fT t );
    fT LiquidCompressibilityFromT( fT t );
    fT VaporCompressibilityFromT( fT t );
    fT LiquidThermalExpansivityFromT( fT t );
    fT VaporThermalExpansivityFromT( fT t );
    fT LiquidPressureDensityDerivativeFromT( fT t );
    fT VaporPressureDensityDerivativeFromT( fT t );
    fT LiquidPressureTemperatureDerivativeFromT( fT t );
    fT VaporPressureTemperatureDerivativeFromT( fT t );
    // for given p, t
    fT Density( fT t, fT p );
    fT HeatCapacity( fT t, fT p );
    fT Enthalpy( fT t, fT p );
    fT Viscosity( fT t, fT p );
    fT ThermalExpansivity( fT t, fT p );
    fT Compressibility( fT t, fT p );
    fT PressureDensityDerivative( fT t, fT p );
    fT PressureTemperatureDerivative( fT t, fT p );

    

private:
    const fT                 pcrit, tcrit, dp, kelvin, tatm, patm, p_cols, t_rows, dT, dP, 
                             dT_fine, dP_fine, p_sat_rows, t_sat_rows, sat_cols, p_sat_min,
                             tmin_c, tmax_c, pmin_c, pmax_c, tc_rows, pc_cols, tsat_atm;
    fT                       pmax, tmax, d, pcurrent, tcurrent;
    fT                       ak[4], bij[6][5], xy[2], prop[4];
    int32_t                    cols, rows, c_cols, c_rows, ij[4][2], it[2], ip[2];
    
    LookUpStorage <fT>       density, enthalpy, heat_capacity, viscosity, 
                             expansivity, compressibility, dp_dd_CT, dp_dT_Cd, t_sat, p_sat,
                             crit_density, crit_enthalpy, crit_heat_capacity, crit_viscosity,
                             crit_expansivity, crit_compressibility, crit_dp_dd_CT, crit_dp_dT_Cd;
    
    enum { regular, boundary, two_pt_liquid, two_pt_vapor, three_pt_liquid, three_pt_vapor } interpolation_type;  
    enum property { cp, rho, h, mu, a, b, dpddCT, dpdTCd };                             
                             

    bool LoadLookupTables();
    bool LoadLookupTablesTwoPhase();
    bool LoadLookupTablesCriticalPoint();
    void ComputeLookupTables();
    void ComputeLookupTablesTwoPhase();
    void ComputeLookupTablesCriticalPoint();
    fT   DynamicViscosity( fT t, fT rho ) const; 
    void CheckError( Prop* prop ) const;
    void GetLookupTableEntries( fT t, fT p );  
    void GetLookupTableEntriesCriticalPoint( fT t, fT p );  
    void GetDistances( fT t, fT p );
    void GetDistancesCriticalPoint( fT t, fT p );
    void GetPropertyValues( property p );
    void GetPropertyValuesCriticalPoint( property p );
    fT   Interpolate( fT t, fT p, property id );
    void GetLookupTableEntriesTwoPhaseCurveForT( fT t );                           
    void GetLookupTableEntriesTwoPhaseCurveForP( fT p ); 
    fT   DistancePTwoPhaseCurve( fT p );
    fT   DistanceTTwoPhaseCurve( fT t ); 
    fT   InterpolateForP( fT p, int32_t i );
    fT   InterpolateForT( fT t, int32_t i ); 
    fT   NumericalPressureDensityDerivative( fT r1, fT t, fT p, fT dpress );
    fT   NumericalCompressibility( fT r1, fT t, fT p, fT dpress );
    fT   NumericalHeatCapacity( fT h1, fT t, fT p, fT dtemp );                        

};

template<typename fT> 
inline fT H2OPropertiesLookUpTable<fT>::Interpolate( fT t, fT p, property id )
 {
   if ( t > tmin_c && t < tmax_c && p > pmin_c && p < pmax_c ) {
       GetLookupTableEntriesCriticalPoint( t, p );
       GetDistancesCriticalPoint( t, p );
       GetPropertyValuesCriticalPoint( id );
     }
   else {
       GetLookupTableEntries( t, p );
       GetDistances( t, p );
       GetPropertyValues( id );
     }
   
   return ( 1.0 - xy[0] ) * ( 1.0 - xy[1] ) * prop[0] + xy[0] * ( 1.0 - xy[1] ) * prop[1] + xy[0] * xy[1] * prop[2] + ( 1.0 - xy[0] ) * xy[1] * prop[3];
 }

template<typename fT> 
inline void H2OPropertiesLookUpTable<fT>::GetLookupTableEntriesTwoPhaseCurveForT( fT t )
 {
   if ( t > tcrit ) {
       it[0] = it[1] = -1; // no 2-phase properties above tcrit
     }
   else if ( t <= tatm ) {
       it[0] = it[1] = 0; // get 2-phase properties at tatm
     }
   else {
       it[0] = static_cast<int32_t>((t-tatm)/dT_fine);
       it[1] = it[0]+1;
     }
 }


template<typename fT> 
inline void H2OPropertiesLookUpTable<fT>::GetLookupTableEntriesTwoPhaseCurveForP( fT p )
 {
   fT press = p_sat_min + dP_fine; // pressure at 2nd entry in p_sat looup table
   
   if ( p > pcrit ) {
       ip[0] = ip[1] = -1; // no 2-phase properties above pcrit
     }
   else if ( p < press ) {
       ip[0] = 0;
       ip[1] = 1; // get 2-phase properties between patm (101325.0) and p_sat_min + dP_fine (105,000 Pa)
     }
   else {
       ip[0] = static_cast<int32_t>((p-press)/dP_fine)+1; // offset by 1 since entry 0 is patm, entry 1 os p_sat_min + dP_fine 
       ip[1] = ip[0]+1;
     }
 }

template<typename fT> 
inline fT H2OPropertiesLookUpTable<fT>::DistancePTwoPhaseCurve( fT p )
 {
   fT press = p_sat_min + dP_fine; // pressure at 2nd entry in p_sat looup table
   if ( p <= patm )                    return 0.0; // patm distance is zero
   else if ( p > patm && p < press )   return ( ( p - patm ) / ( press - patm ) ); // interpolate between 1st and 2nd entry in p_sat
   else {
       if ( ip[1] != p_sat.Rows()-1 )  return ( ( p - static_cast<fT>(ip[0])*dP_fine - p_sat_min)/dP_fine ); // regular interpolation below pcrit and above patm
       else                            return ( ( p - static_cast<fT>(ip[0])*dP_fine - p_sat_min)/( pcrit - static_cast<fT>(ip[0])*dP_fine - p_sat_min ) ); // rescale if last entry i pcrit
     }
 }

template<typename fT> 
inline fT H2OPropertiesLookUpTable<fT>::DistanceTTwoPhaseCurve( fT t )
 {
   if ( t <= tatm )                    return 0.0; // tatm distance is zero
   else {
       if ( it[1] != t_sat.Rows()-1 )  return ( ( t - static_cast<fT>(it[0])*dT_fine - tatm)/dT_fine ); // regular interpolation below tcrit and above tatm
       else                            return ( ( t - static_cast<fT>(it[0])*dT_fine - tatm)/( tcrit - static_cast<fT>(it[0])*dT_fine - tatm ) ); // rescale if last entry i tcrit
     }
 }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::InterpolateForT( fT t, int32_t i ) 
 {
   fT u;
   // obvious entries
   if ( t > tcrit )      return 0.0;
   else if ( t <= tatm ) return t_sat(0,i);
   else {
       GetLookupTableEntriesTwoPhaseCurveForT( t );
       u = DistanceTTwoPhaseCurve( t );
       return ( t_sat(it[0],i) - u * t_sat(it[0],i) + u * t_sat(it[1],i) );
     }
 }


template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::InterpolateForP( fT p, int32_t i ) 
 {
   fT u;
   // obvious entries
   if ( p > pcrit )      return 0.0;
   else if ( p <= patm ) return p_sat(0,i);
   else {
       GetLookupTableEntriesTwoPhaseCurveForP( p );
       u = DistancePTwoPhaseCurve( p );
       return ( p_sat(ip[0],i) - u * p_sat(ip[0],i) + u * p_sat(ip[1],i) );
     }
 }


// for given t and p

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::Density( fT t, fT p ) { return Interpolate( t, p, rho ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::Enthalpy( fT t, fT p ) { return Interpolate( t, p, h ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::HeatCapacity( fT t, fT p ) { return Interpolate( t, p, cp ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::Viscosity( fT t, fT p ) { return Interpolate( t, p, mu ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::Compressibility( fT t, fT p ) { return Interpolate( t, p, b ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::ThermalExpansivity( fT t, fT p ) { return Interpolate( t, p, a ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::PressureDensityDerivative( fT t, fT p ) { return Interpolate( t, p, dpddCT ); }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::PressureTemperatureDerivative( fT t, fT p ) { return Interpolate( t, p, dpdTCd ); }


// 2 -phase properties

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidSaturation( fT p, fT t ) {
    
    if ( p < pcrit && t < tcrit ) {
        if ( t <= InterpolateForP( p, 0 ) ) return 1.0; // t below saturation curve for given p -> liquid
        else                                return 0.0;
      }
    else return 1.0; // super critical fluid has liquid saturation set to 1
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporSaturation( fT p, fT t ) {
   
    return 1.0 - LiquidSaturation( p, t ); 
    
  }


// for given p
template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::SaturationTemperatureFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 0 ); 
    else              return 0.0; 
    
  }


template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidDensityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 1 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporDensityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 2 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidEnthalpyFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 3 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporEnthalpyFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 4 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidHeatCapacityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 5 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporHeatCapacityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 6 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidViscosityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 7 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporViscosityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 8 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidCompressibilityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 9 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporCompressibilityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 10 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidThermalExpansivityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 11 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporThermalExpansivityFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 12 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidPressureDensityDerivativeFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 13 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporPressureDensityDerivativeFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 14 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidPressureTemperatureDerivativeFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 15 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporPressureTemperatureDerivativeFromP( fT p ) {
    
    if ( p <= pcrit ) return InterpolateForP( p, 16 ); 
    else              return 0.0; 
    
  }



// for given t
template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::SaturationPressureFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 0 ); 
    else              return 0.0; 
    
  }


template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidDensityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 1 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporDensityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 2 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidEnthalpyFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 3 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporEnthalpyFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 4 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidHeatCapacityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 5 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporHeatCapacityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 6 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidViscosityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 7 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporViscosityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 8 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidCompressibilityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 9 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporCompressibilityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 10 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidThermalExpansivityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 11 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporThermalExpansivityFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 12 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidPressureDensityDerivativeFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 13 ); 
    else              return 0.0;   
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporPressureDensityDerivativeFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 14 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::LiquidPressureTemperatureDerivativeFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 15 ); 
    else              return 0.0; 
    
  }

template<typename fT>
inline fT H2OPropertiesLookUpTable<fT>::VaporPressureTemperatureDerivativeFromT( fT t ) {
    
    if ( t <= tcrit ) return InterpolateForT( t, 16 ); 
    else              return 0.0; 
    
  }


} // namespace csmp

#endif
