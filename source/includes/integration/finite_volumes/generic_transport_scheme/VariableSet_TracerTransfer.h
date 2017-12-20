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
     const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>  ff_key;  ///< facet flux
     const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>  ffC_key; ///< facet flux concentration product
     const csmp::INDEX<SCALAR,NODE>  PV_key;  ///< FV pore volume (FV)
     const csmp::INDEX<SCALAR,NODE>  fb_key;  ///< flux balance (FV)
     const csmp::INDEX<SCALAR,NODE>  nsrc_key;  ///< nodal fluid volume source (absolute concentration*volume product assigned to FV)
     // general
     const csmp::INDEX<SCALAR,NODE>  pf_key;  ///< fluid pressure
     const csmp::INDEX<SCALAR,ELEMENT>  k_key;   ///< permeability
     const csmp::INDEX<SCALAR,ELEMENT>  phi_key; ///< porosity
     const csmp::INDEX<SCALAR,ELEMENT>  thi_key; ///< thickness
     const csmp::INDEX<VECTOR,ELEMENT>  vD_key;  ///< Darcy velocity
     // fluid properties
     const csmp::INDEX<SCALAR,MODEL>  mu_key;   ///< fluid viscosity
     const csmp::INDEX<SCALAR,MODEL>  rhof_key; ///< fluid density
     const csmp::INDEX<SCALAR,NODE>  C0_key;   ///< concentration
     const csmp::INDEX<SCALAR,NODE>  C1_key;   ///< new concentration
     // other
     const csmp::INDEX<SCALAR,ELEMENT>  diff_key; ///< diffusion coefficient (computed from fluid and material properties)

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
