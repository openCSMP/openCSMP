// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_VARIABLESET_TRACERTRANSFER_H
#define CSMP_VARIABLESET_TRACERTRANSFER_H

/**
@file VariableSet_TracerTransfer.h
Automatically generated from VariableSet_TracerTransfer.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct VariableSet_TracerTransfer {
  // generic finite volume scheme
  csmp::INDEX<SCALAR,NODE> key_out_; // outflow
  csmp::INDEX<SCALAR,NODE> key_acc_; // accumulation
  // multiphase flow
  csmp::INDEX<SCALAR,ELEMENT> key_CT; // total system compressibility
  // numeric modelling
  csmp::INDEX<SCALAR,NODE> key_PF0; // previous fluid pressure
  // single phase flow
  csmp::INDEX<SCALAR,NODE> key_FVPV; // FV pore volume
  csmp::INDEX<SCALAR,NODE> key_FB; // flux balance
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_SPV; // sector pore volume
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_SV; // sector volume
  csmp::INDEX<SCALAR,NODE> key_NQV; // nodal fluid volume source
  csmp::INDEX<SCALAR,ELEMENT> key_QV; // fluid volume source
  csmp::INDEX<SCALAR,NODE> key_PF; // fluid pressure
  csmp::INDEX<TENSOR,ELEMENT> key_k; // permeability
  csmp::INDEX<TENSOR,ELEMENT> key_kk; // tensor permeability
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fAk; // facet area permeability
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; // facet normal
  csmp::INDEX<SCALAR,ELEMENT> key_THI; // thickness
  csmp::INDEX<VECTOR,ELEMENT> key_V; // velocity
  csmp::INDEX<SCALAR,MODEL> key_MU; // fluid viscosity
  csmp::INDEX<SCALAR,MODEL> key_RHO; // fluid density
  csmp::INDEX<SCALAR,NODE> key_C; // previous concentration
  csmp::INDEX<SCALAR,NODE> key_C0; // concentration
  csmp::INDEX<SCALAR,NODE> key_C1; // new concentration
  csmp::INDEX<TENSOR,ELEMENT> key_K; // conductivity
  csmp::INDEX<TENSOR,ELEMENT> key_KK; // tensor conductivity
  csmp::INDEX<SCALAR,ELEMENT> key_D; // diffusivity
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_ffC; // facet flux concentration
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_ff; // facet flux
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA; // facet area
  csmp::INDEX<SCALAR,ELEMENT> key_PHI; // porosity

  template<uint32_t dim>
  explicit VariableSet_TracerTransfer( const PropertyDatabase<dim>& db )
    : key_out_( INDEX<SCALAR,NODE>( db.StorageKey("outflow") ))
    , key_acc_( INDEX<SCALAR,NODE>( db.StorageKey("accumulation") ))
    , key_CT( INDEX<SCALAR,ELEMENT>( db.StorageKey("total system compressibility") ))
    , key_PF0( INDEX<SCALAR,NODE>( db.StorageKey("previous fluid pressure") ))
    , key_FVPV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_FB( INDEX<SCALAR,NODE>( db.StorageKey("flux balance") ))
    , key_SPV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector pore volume") ))
    , key_SV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector volume") ))
    , key_NQV( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid volume source") ))
    , key_QV( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid volume source") ))
    , key_PF( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_k( INDEX<TENSOR,ELEMENT>( db.StorageKey("permeability") ))
    , key_kk( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor permeability") ))
    , key_fAk( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area permeability") ))
    , key_fn( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal") ))
    , key_THI( INDEX<SCALAR,ELEMENT>( db.StorageKey("thickness") ))
    , key_V( INDEX<VECTOR,ELEMENT>( db.StorageKey("velocity") ))
    , key_MU( INDEX<SCALAR,MODEL>( db.StorageKey("fluid viscosity") ))
    , key_RHO( INDEX<SCALAR,MODEL>( db.StorageKey("fluid density") ))
    , key_C( INDEX<SCALAR,NODE>( db.StorageKey("previous concentration") ))
    , key_C0( INDEX<SCALAR,NODE>( db.StorageKey("concentration") ))
    , key_C1( INDEX<SCALAR,NODE>( db.StorageKey("new concentration") ))
    , key_K( INDEX<TENSOR,ELEMENT>( db.StorageKey("conductivity") ))
    , key_KK( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor conductivity") ))
    , key_D( INDEX<SCALAR,ELEMENT>( db.StorageKey("diffusivity") ))
    , key_ffC( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux concentration") ))
    , key_ff( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux") ))
    , key_fA( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area") ))
    , key_PHI( INDEX<SCALAR,ELEMENT>( db.StorageKey("porosity") ))
  {
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLESET_TRACERTRANSFER_H
