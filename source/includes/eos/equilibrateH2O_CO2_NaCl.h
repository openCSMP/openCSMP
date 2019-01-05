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
#include "PVTX_Calculator_H2O_CO2_NaCl.h"

namespace csmp {
  
  template<size_t> class Node;
  template<size_t> class Element;
  template<size_t> class Region;
  template<size_t> class PVTX_Calculator_H2O_CO2_NaCl;
  
  /// applying PVTX_Calculator_H2O_CO2_NaCl to all nodes of a region
  template<size_t dim>
  void equilibrateFluid( Region<dim>& model_subdomain, PVTX_Calculator_H2O_CO2_NaCl<dim>& pvtx_calculator, double64 delta_t );
  
  /// updating of fluid properties in all regions, not just where DES is invoked
  template<size_t dim>
  void updatePTXBasedFluidProperties(const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain );

  /// computes (barycentric) porosity from volume fraction of salt that is occupying the pore space
  template<size_t dim>
  double64 porosityWithSalt( const variables::VariableSet_CO2GeoSequestration&, Element<dim>& );
  
  /// accounts for porosity and permeability changes due to salt precipitation
  template<size_t dim>
  void updatePorosityAndPermeability( const variables::VariableSet_CO2GeoSequestration&, Region<dim>& );  
  
  /// accounts for pore-volume changes due to salt precipitation
  template<size_t dim>
  void updatePoreVolume( const variables::VariableSet_CO2GeoSequestration&, Region<dim>& );   
  
  /// update PTX based fluid properties
  template<size_t dim>
  void updateFluidProperties( const variables::VariableSet_CO2GeoSequestration&, Node<dim>&, double64 del_t );   
  
  template<size_t dim>
  void updatePTXBasedFluidProperties( const variables::VariableSet_CO2GeoSequestration&, Region<dim>&, double64 del_t );   

}

#endif /* CSMP_EQUILIBRATE_H2O_CO2_NACL_H */
