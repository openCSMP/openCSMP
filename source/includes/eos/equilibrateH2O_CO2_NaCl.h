//
//  equilibrateH2O_CO2_NaCl.h
//  CSMP_CO2GeoSequestrationSimulator
//
//  Created by Stephan Matthai on 21/10/18.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_EQUILIBRATE_H2O_CO2_NACL_H
#define CSMP_EQUILIBRATE_H2O_CO2_NACL_H

#include "VariableSet_CO2GeoSequestration.h"

namespace csmp {
  
  template<size_t> class Node;
  template<size_t> class Element;
  template<size_t> class Region;

  enum AQUEOUS_PHASE { XH2O=0, XCO2=1, XNACl_aq=2 }; ///< mass fractions stored 
  enum CARBONIC_PHASE { YCO2=0, YH2O=1 };
  enum THERMODYNAMIC_STATE { LIQUID, GASEOUS, SUPERCRITICAL };
  
  /// computes phase state of pure CO2
  THERMODYNAMIC_STATE  stateOfCO2( double64 pCO2, double64 TC );
  
  /// Pruess (2005), compute aqueous phase density from brine density and amount of CO2 dissolved
  double64 aqueousPhaseDensity_H2O_CO2_NaCl( double64 brine_density, double64 rhoCO2, double64 X_CO2 );
  
  /// compute kg/m3 of brine from mass fraction; NB: salinity is not the transport variable
  inline double64 salinity( double64 XNaCl_aq, double64 brine_density ) { return XNaCl_aq * brine_density; }

  /// computes new compositions of aqueous and carbonic phase, precipitates salt if any, updates and saturations, densities and viscosities
  template<size_t dim>
  void equilibrateH2O_CO2_NaCl( const variables::VariableSet_CO2GeoSequestration&, Node<dim>& );
  
  /// for given p,T and composition equilibrate phases with one-another
  template<size_t dim>
  void equilibrateFluid( const variables::VariableSet_CO2GeoSequestration&, Region<dim>& );  
  
  /// computes (barycentric) porosity from volume fraction of salt that is occupying the pore space
  template<size_t dim>
  double64 porosityWithSalt( const variables::VariableSet_CO2GeoSequestration&, Element<dim>& );
  
  /// accounts for pore-volume changes due to salt precipitation
  template<size_t dim>
  void updatePorosity( const variables::VariableSet_CO2GeoSequestration&, Region<dim>& );  

}

#endif /* CSMP_EQUILIBRATE_H2O_CO2_NACL_H */
