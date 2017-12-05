//
//  VariableSet_TracerTransfer
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#ifndef VARIABLE_SET_TRACER_TRANSFER_H
#define VARIABLE_SET_TRACER_TRANSFER_H

#include "Index.h"

namespace csmp {

template<size_t> class PropertyDatabase;

/**
    Variables used by the tracer-transfer scheme.
*/
struct VariableSet_TracerTransfer {
     template<size_t dim>
     explicit VariableSet_TracerTransfer( const PropertyDatabase<dim>& );

     // finite volume variables
     const csmp::Index  ff_key;  ///< facet flux
     const csmp::Index  ffC_key; ///< facet flux concentration product
     const csmp::Index  PV_key;  ///< FV pore volume (FV)
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
     VariableSet_TracerTransfer() = delete;
     VariableSet_TracerTransfer( const VariableSet_TracerTransfer& ) = delete;
     VariableSet_TracerTransfer& operator=( const VariableSet_TracerTransfer& ) = delete;
  
  protected:
     /// verifies types and placements
     void CheckVariables() const;
};
 
} // end csmp

#endif /* defined(VARIABLE_SET_TRACER_TRANSFER_H) */
