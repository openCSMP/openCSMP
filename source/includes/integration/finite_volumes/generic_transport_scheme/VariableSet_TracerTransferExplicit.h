//
//  VariableSet_TracerTransferExplicit
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#ifndef VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H
#define VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H

#include "Index.h"

namespace csmp {

template<size_t> class PropertyDatabase;

/**
    Variables used by the tracer-transfer scheme.
*/
template<size_t dim>
struct VariableSet_TracerTransferExplicit {
     explicit VariableSet_TracerTransferExplicit( const PropertyDatabase<dim>& );
     // finite volume variables
     const csmp::Index  fA_key;  ///< facet area
     const csmp::Index  fn_key;  ///< facet normal
     const csmp::Index  fnk_key; ///< facet normal permeability
     const csmp::Index  ff_key;  ///< facet flux
     const csmp::Index  ffC_key; ///< facet flux concentration product
     const csmp::Index  PV_key;  ///< FV pore volume (FV)
     const csmp::Index  spv_key; ///< sector pore volume
     const csmp::Index  fb_key;  ///< flux balance (FV)
     const csmp::Index  nsrc_key;  ///< nodal fluid volume source (absolute concentration*volume product assigned to FV)
     // general
     const csmp::Index  pf_key;  ///< fluid pressure
     const csmp::Index  k_key;   ///< permeability
     const csmp::Index  phi_key; ///< porosity
     const csmp::Index  thi_key; ///< thickness
     const csmp::Index  vD_key;  ///< Darcy velocity
     // fluid properties
     const csmp::Index  mu_key;   ///< fluid viscosity
     const csmp::Index  rhof_key; ///< fluid density
     const csmp::Index  C0_key;   ///< concentration
     const csmp::Index  C1_key;   ///< new concentration
 
   private:
     VariableSet_TracerTransferExplicit() {}; // = delete;
     VariableSet_TracerTransferExplicit( const VariableSet_TracerTransferExplicit& ) {}; // = delete;
     VariableSet_TracerTransferExplicit<dim>& operator=( const VariableSet_TracerTransferExplicit& ) { return *this; }; // = delete;

     /// verifies types and placements
     void CheckVariables() const;
};
 
} // end csmp

#endif /* defined(VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H) */
