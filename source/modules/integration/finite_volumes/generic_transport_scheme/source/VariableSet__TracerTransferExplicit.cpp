//
//  VariableSet_TracerTransferExplicit.cpp
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#include "VariableSet_TracerTransferExplicit.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp {

template<size_t dim>
VariableSet_TracerTransferExplicit<dim>::VariableSet_TracerTransferExplicit( const PropertyDatabase<dim>& p )
 : fA_key(p.StorageKey("facet area")),
   fn_key(p.StorageKey("facet normal")),
   fnk_key(p.StorageKey("facet normal permeability")),
   ff_key(p.StorageKey("facet flux")),
   ffC_key(p.StorageKey("facet flux concentration product")),
   PV_key(p.StorageKey("FV pore volume")),
   spv_key(p.StorageKey("sector pore volume")),
   fb_key(p.StorageKey("flux balance")),
   nsrc_key(p.StorageKey("nodal fluid volume source")),
   //
   pf_key(p.StorageKey("fluid pressure")),
   k_key(p.StorageKey("permeability")),
   phi_key(p.StorageKey("porosity")),
   thi_key(p.StorageKey("thickness")),
   vD_key(p.StorageKey("velocity")),
   //
   mu_key(p.StorageKey("fluid viscosity")),
   rhof_key(p.StorageKey("fluid density")),
   C0_key(p.StorageKey("concentration")),
   C1_key(p.StorageKey("new concentration"))
 {
   // potentially these properties could be created here from scratch to have them only
   // when there is a need for a transport calculation.
 }




 
template<size_t dim>
void VariableSet_TracerTransferExplicit<dim>::CheckVariables() const
 {
      // finite volume variables
    if ( PV_key.place != NODE || PV_key.type != SCALAR )
       throw csmp::Exception( CSMP_FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'pore volume' variable must be scalar placed on the nodes = node-centered finite volumes" );

    if ( fn_key.place != FACET_INTEGRATION_POINT || fn_key.type != VECTOR )
       throw csmp::Exception( CSMP_FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'facet normal' variable must be vector placed on the facet integration points" );
   
    if ( fA_key.place != FACET_INTEGRATION_POINT || fA_key.type != SCALAR )
       throw csmp::Exception( CSMP_FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'facet area' variable must be scalar placed on the facet integration points" );
   
    if ( fn_key.place != FACET_INTEGRATION_POINT || fn_key.type != VECTOR )
       throw csmp::Exception( CSMP_FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'facet normal' variable must be vector placed on the facet integration points" );
   
// TODO: complete missing checks

} // end CheckVariables

 
template class VariableSet_TracerTransferExplicit<1U>;
template class VariableSet_TracerTransferExplicit<2U>;
template class VariableSet_TracerTransferExplicit<3U>;
 
} // end csmp
