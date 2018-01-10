#ifndef CSMP_VARIABLES_TRACERTRANSFER_H
#define CSMP_VARIABLES_TRACERTRANSFER_H

/**
@file Variables_TracerTransfer.h
Automatically generated from /Users/andrew/Development/repository/csmp-api-library/source/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_TracerTransfer.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct Variables_TracerTransfer {
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_ff; // facet flux
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_ffC; // facet flux concentration
  csmp::INDEX<SCALAR,NODE> key_FVPV; // FV pore volume
  csmp::INDEX<SCALAR,NODE> key_FB; // flux balance
  csmp::INDEX<SCALAR,NODE> key_NQV; // nodal fluid volume source
  csmp::INDEX<SCALAR,NODE> key_PF; // fluid pressure
  csmp::INDEX<SCALAR,ELEMENT> key_k; // permeability
  csmp::INDEX<SCALAR,ELEMENT> key_PHI; // porosity
  csmp::INDEX<SCALAR,ELEMENT> key_THI; // thickness
  csmp::INDEX<VECTOR,ELEMENT> key_V; // velocity
  csmp::INDEX<SCALAR,MODEL> key_MU; // fluid viscosity
  csmp::INDEX<SCALAR,MODEL> key_RHO; // fluid density
  csmp::INDEX<SCALAR,NODE> key_C; // concentration
  csmp::INDEX<SCALAR,NODE> key_NC; // new concentration
  csmp::INDEX<SCALAR,ELEMENT> key_K; // conductivity
  csmp::INDEX<SCALAR,ELEMENT> key_D; // diffusivity

  template<size_t dim>
  explicit Variables_TracerTransfer( const PropertyDatabase<dim>& db )
    : key_ff( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux") ))
    , key_ffC( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux concentration") ))
    , key_FVPV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_FB( INDEX<SCALAR,NODE>( db.StorageKey("flux balance") ))
    , key_NQV( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid volume source") ))
    , key_PF( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_k( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability") ))
    , key_PHI( INDEX<SCALAR,ELEMENT>( db.StorageKey("porosity") ))
    , key_THI( INDEX<SCALAR,ELEMENT>( db.StorageKey("thickness") ))
    , key_V( INDEX<VECTOR,ELEMENT>( db.StorageKey("velocity") ))
    , key_MU( INDEX<SCALAR,MODEL>( db.StorageKey("fluid viscosity") ))
    , key_RHO( INDEX<SCALAR,MODEL>( db.StorageKey("fluid density") ))
    , key_C( INDEX<SCALAR,NODE>( db.StorageKey("concentration") ))
    , key_NC( INDEX<SCALAR,NODE>( db.StorageKey("new concentration") ))
    , key_K( INDEX<SCALAR,ELEMENT>( db.StorageKey("conductivity") ))
    , key_D( INDEX<SCALAR,ELEMENT>( db.StorageKey("diffusivity") ))
  {
    if ( key_ff.place != FACET_INTEGRATION_POINT || key_ff.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'facet flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_ffC.place != FACET_INTEGRATION_POINT || key_ffC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'facet flux concentration' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_FVPV.place != NODE || key_FVPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( key_FB.place != NODE || key_FB.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( key_NQV.place != NODE || key_NQV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( key_PF.place != NODE || key_PF.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_k.place != ELEMENT || key_k.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_PHI.place != ELEMENT || key_PHI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_THI.place != ELEMENT || key_THI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_V.place != ELEMENT || key_V.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_MU.place != MODEL || key_MU.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'fluid viscosity' variable must be SCALAR and placed on MODEL"  );
    if ( key_RHO.place != MODEL || key_RHO.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'fluid density' variable must be SCALAR and placed on MODEL"  );
    if ( key_C.place != NODE || key_C.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'concentration' variable must be SCALAR and placed on NODE"  );
    if ( key_NC.place != NODE || key_NC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'new concentration' variable must be SCALAR and placed on NODE"  );
    if ( key_K.place != ELEMENT || key_K.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_D.place != ELEMENT || key_D.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TracerTransfer::Variables_TracerTransfer:",
        "The 'diffusivity' variable must be SCALAR and placed on ELEMENT"  );
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLES_TRACERTRANSFER_H
