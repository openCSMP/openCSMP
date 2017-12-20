//
//  VariableSet_TracerTransfer.cpp
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#include "VariableSet_TracerTransfer.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp {

template<size_t dim>
VariableSet_TracerTransfer::VariableSet_TracerTransfer( const PropertyDatabase<dim>& p )

 : ff_key(INDEX<SCALAR,FACET_INTEGRATION_POINT>(p.StorageKey("facet flux"))),
   ffC_key(INDEX<SCALAR,FACET_INTEGRATION_POINT>(p.StorageKey("facet flux concentration"))),
   PV_key(INDEX<SCALAR,NODE>(p.StorageKey("FV pore volume"))),
   fb_key(INDEX<SCALAR,NODE>(p.StorageKey("flux balance"))),
   nsrc_key(INDEX<SCALAR,NODE>(p.StorageKey("nodal fluid volume source"))),
   //
   pf_key(INDEX<SCALAR,NODE>(p.StorageKey("fluid pressure"))),
   k_key(INDEX<SCALAR,NODE>(p.StorageKey("permeability"))),
   phi_key(INDEX<SCALAR,ELEMENT>(p.StorageKey("porosity"))),
   thi_key(INDEX<SCALAR,ELEMENT>(p.StorageKey("thickness"))),
   vD_key(INDEX<VECTOR,ELEMENT>(p.StorageKey("velocity"))),

   //
   mu_key(INDEX<SCALAR,MODEL>(p.StorageKey("fluid viscosity"))),
   rhof_key(INDEX<SCALAR,MODEL>(p.StorageKey("fluid density"))),
   C0_key(INDEX<SCALAR,NODE>(p.StorageKey("concentration"))),
   C1_key(INDEX<SCALAR,NODE>(p.StorageKey("new concentration"))),
   diff_key(INDEX<SCALAR,ELEMENT>(p.StorageKey("diffusivity")))

{
}


 
void VariableSet_TracerTransfer::CheckVariables() const
 {
      // finite volume variables
    if ( PV_key.place != NODE || PV_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'pore volume' variable must be scalar placed on the nodes = node-centered finite volumes" );


    if ( phi_key.place != ELEMENT || phi_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'porosity' variable must be scalar placed on the nodes = node-centered finite volumes" );

    if ( k_key.place != MODEL || k_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'permeability' variable must be scalar placed on the model (this will be extended in the future)" );

    if ( ff_key.place != FACET_INTEGRATION_POINT || ff_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'facet flux' variable must be scalar placed on the facet integration points" );

    if ( ffC_key.place != FACET_INTEGRATION_POINT || ffC_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "VariableSet_TracerTransferExplicit::CheckVariables:",
                      "The 'facet flux concentration' variable must be scalar placed on the facet integration points" );

} // end CheckVariables

 

template VariableSet_TracerTransfer::VariableSet_TracerTransfer( const PropertyDatabase<1U>& );
template VariableSet_TracerTransfer::VariableSet_TracerTransfer( const PropertyDatabase<2U>& );
template VariableSet_TracerTransfer::VariableSet_TracerTransfer( const PropertyDatabase<3U>& );
 
} // end csmp
