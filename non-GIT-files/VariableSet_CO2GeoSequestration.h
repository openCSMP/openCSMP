#ifndef CSMP_VARIABLESET_CO2GEOSEQUESTRATION_H
#define CSMP_VARIABLESET_CO2GEOSEQUESTRATION_H

/**
@file VariableSet_CO2GeoSequestration.h
Automatically generated from /home/unimelb.edu.au/qssha/CLionProjects/open-acgss-16Jan23/open-csmp/source/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_CO2GeoSequestration.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct VariableSet_CO2GeoSequestration {
  csmp::INDEX<ARRAY,ELEMENT> key_kri_param; // relative permeability model parameters
  csmp::INDEX<ARRAY,MODEL> key_t_out; // output times
  csmp::INDEX<SCALAR,NODE> key_Nn; // node number
  csmp::INDEX<VECTOR,FACE> key_dipf; // face dip vector
  csmp::INDEX<SCALAR,ELEMENT> key_Ne; // element number
  csmp::INDEX<SCALAR,MODEL> key_dt; // default time increment
  csmp::INDEX<SCALAR,NODE> key_fvPV; // FV pore volume
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; // facet normal
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA; // facet area
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_sPV; // sector pore volume
  csmp::INDEX<VECTOR,ELEMENT> key_dip; // dip vector
  csmp::INDEX<SCALAR,NODE> key_pf0; // previous fluid pressure
  csmp::INDEX<SCALAR,ELEMENT> key_bcp; // brooks corey parameter
  csmp::INDEX<SCALAR,MODEL> key_t; // model time
  csmp::INDEX<SCALAR,ELEMENT> key_ct; // total system compressibility
  csmp::INDEX<TENSOR,ELEMENT> key_LT; // tensor total mobility permeability product
  csmp::INDEX<SCALAR,ELEMENT> key_lambda_t; // total mobility permeability product
  csmp::INDEX<SCALAR,ELEMENT> key_lCO2; // mobility carbonic phase (lambda_CO2,sc * rho_CO2)
  csmp::INDEX<SCALAR,ELEMENT> key_lH2O; // mobility aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_ssH2O; // shock saturation aqueous phase
  csmp::INDEX<SCALAR,NODE> key_sH2O; // saturation aqueous phase
  csmp::INDEX<SCALAR,NODE> key_sCO2; // saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_UPC; // update count
  csmp::INDEX<SCALAR,MODEL> key_eoi; // end of injection
  csmp::INDEX<ARRAY,NODE> key_Xbulk; // bulk mass and fractions
  csmp::INDEX<SCALAR,ELEMENT> key_ePHS; // element system state
  csmp::INDEX<SCALAR,NODE> key_nPHS; // nodal system state
  csmp::INDEX<SCALAR,ELEMENT> key_Qh; // heat source
  csmp::INDEX<ARRAY,NODE> key_NTA; // non-wetting phase timing array
  csmp::INDEX<SCALAR,NODE> key_CFL; // cfl multiplier
  csmp::INDEX<SCALAR,NODE> key_SYC; // synchronize count
  csmp::INDEX<SCALAR,NODE> key_SCC; // schedule count
  csmp::INDEX<SCALAR,NODE> key_RAC; // rate count
  csmp::INDEX<ARRAY,NODE> key_H2O_comp; // composition aqueous phase
  csmp::INDEX<SCALAR,NODE> key_EVI; // event index
  csmp::INDEX<SCALAR,ELEMENT> key_RRT; // rocktype
  csmp::INDEX<SCALAR,ELEMENT> key_phi; // porosity
  csmp::INDEX<SCALAR,ELEMENT> key_k_cal; // permeability calibration factor
  csmp::INDEX<SCALAR,ELEMENT> key_kV; // vertical permeability
  csmp::INDEX<SCALAR,ELEMENT> key_k; // permeability
  csmp::INDEX<TENSOR,ELEMENT> key_kk; // tensor permeability
  csmp::INDEX<SCALAR,MODEL> key_g; // acceleration gravity
  csmp::INDEX<SCALAR,ELEMENT> key_thi; // thickness
  csmp::INDEX<SCALAR,MODEL> key_TMAX; // run duration
  csmp::INDEX<SCALAR,NODE> key_cCO2; // compressibility carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_Kr; // thermal conductivity
  csmp::INDEX<SCALAR,NODE> key_KCO2; // thermal conductivity carbonic phase
  csmp::INDEX<SCALAR,NODE> key_KH2O; // thermal conductivity aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_eC; // electric conductivity
  csmp::INDEX<SCALAR,NODE> key_muH2O; // viscosity aqueous phase
  csmp::INDEX<SCALAR,NODE> key_muCO2; // viscosity carbonic phase
  csmp::INDEX<SCALAR,NODE> key_IFT; // interfacial tension
  csmp::INDEX<SCALAR,NODE> key_rhoH2O; // density aqueous phase
  csmp::INDEX<SCALAR,NODE> key_rhoCO2; // density carbonic phase
  csmp::INDEX<SCALAR,NODE> key_cH2O; // compressibility aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_srCO2; // residual saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_mTot; // total mass
  csmp::INDEX<SCALAR,NODE> key_mNaCl; // mass NaCl
  csmp::INDEX<SCALAR,NODE> key_mH2O; // mass H2O
  csmp::INDEX<SCALAR,NODE> key_mCO2; // mass CO2
  csmp::INDEX<SCALAR,NODE> key_NaClaq; // salinity
  csmp::INDEX<SCALAR,NODE> key_NaCl; // salt
  csmp::INDEX<SCALAR,NODE> key_H2Og; // evaporated water
  csmp::INDEX<SCALAR,NODE> key_CO2aq; // dissolved CO2
  csmp::INDEX<ARRAY,NODE> key_CO2_comp; // composition carbonic phase
  csmp::INDEX<SCALAR,NODE> key_T; // temperature
  csmp::INDEX<SCALAR,ELEMENT> key_srH2O; // residual saturation aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_dsCO2; // change of saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_sCO2_0; // old saturation carbonic phase
  csmp::INDEX<VECTOR,ELEMENT> key_gt; // gravity term
  csmp::INDEX<SCALAR,ELEMENT> key_QM0; // previous fluid mass source
  csmp::INDEX<SCALAR,NODE> key_nQM0; // previous nodal fluid mass source
  csmp::INDEX<SCALAR,ELEMENT> key_QM; // fluid mass source
  csmp::INDEX<SCALAR,NODE> key_nQM; // nodal fluid mass source
  csmp::INDEX<SCALAR,NODE> key_rpf; // reduced fluid pressure
  csmp::INDEX<SCALAR,NODE> key_rpf0; // previous reduced fluid pressure
  csmp::INDEX<SCALAR,NODE> key_pf; // fluid pressure
  csmp::INDEX<VECTOR,FACE> key_fgt; // face gravity term
  csmp::INDEX<SCALAR,ELEMENT> key_Qe; // energy source
  csmp::INDEX<SCALAR,ELEMENT> key_cR; // compressibility rock
  csmp::INDEX<SCALAR,NODE> key_fv_ctr; // FV total rock compressibility
  csmp::INDEX<SCALAR,NODE> key_fv_ctm; // FV total system compressibility
  csmp::INDEX<SCALAR,ELEMENT> key_rhod; // dry rock density
  csmp::INDEX<SCALAR,ELEMENT> key_diff; // CO2 diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_diffpc; // capillary diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_pd; // entry pressure
  csmp::INDEX<VECTOR,FACE> key_flt; // face total mobility permeability product
  // fluid properties
  csmp::INDEX<SCALAR,NODE> key_rhom; // mixture density
  // multiphase flow
  csmp::INDEX<VECTOR,ELEMENT> key_vt; // total velocity
  csmp::INDEX<SCALAR,NODE> key_SwDrToImb; // previous drainage endpoint
  csmp::INDEX<SCALAR,NODE> key_SwImbToDr; // previous imbibition endpoint
  csmp::INDEX<SCALAR,ELEMENT> key_psrCO2; // pseudo residual saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_psrH2O; // pseudo residual saturation aqueous phase
  // output
  csmp::INDEX<SCALAR,NODE> key_mNaClaq; // mass NaClaq

  template<uint32_t dim>
  explicit VariableSet_CO2GeoSequestration( const PropertyDatabase<dim>& db )
    : key_kri_param( INDEX<ARRAY,ELEMENT>( db.StorageKey("relative permeability model parameters") ))
    , key_t_out( INDEX<ARRAY,MODEL>( db.StorageKey("output times") ))
    , key_Nn( INDEX<SCALAR,NODE>( db.StorageKey("node number") ))
    , key_dipf( INDEX<VECTOR,FACE>( db.StorageKey("face dip vector") ))
    , key_Ne( INDEX<SCALAR,ELEMENT>( db.StorageKey("element number") ))
    , key_dt( INDEX<SCALAR,MODEL>( db.StorageKey("default time increment") ))
    , key_fvPV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_fn( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal") ))
    , key_fA( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area") ))
    , key_sPV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector pore volume") ))
    , key_dip( INDEX<VECTOR,ELEMENT>( db.StorageKey("dip vector") ))
    , key_pf0( INDEX<SCALAR,NODE>( db.StorageKey("previous fluid pressure") ))
    , key_bcp( INDEX<SCALAR,ELEMENT>( db.StorageKey("brooks corey parameter") ))
    , key_t( INDEX<SCALAR,MODEL>( db.StorageKey("model time") ))
    , key_ct( INDEX<SCALAR,ELEMENT>( db.StorageKey("total system compressibility") ))
    , key_LT( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor total mobility permeability product") ))
    , key_lambda_t( INDEX<SCALAR,ELEMENT>( db.StorageKey("total mobility permeability product") ))
    , key_lCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility carbonic phase") ))
    , key_lH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility aqueous phase") ))
    , key_ssH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("shock saturation aqueous phase") ))
    , key_sH2O( INDEX<SCALAR,NODE>( db.StorageKey("saturation aqueous phase") ))
    , key_sCO2( INDEX<SCALAR,NODE>( db.StorageKey("saturation carbonic phase") ))
    , key_UPC( INDEX<SCALAR,NODE>( db.StorageKey("update count") ))
    , key_eoi( INDEX<SCALAR,MODEL>( db.StorageKey("end of injection") ))
    , key_Xbulk( INDEX<ARRAY,NODE>( db.StorageKey("bulk mass and fractions") ))
    , key_ePHS( INDEX<SCALAR,ELEMENT>( db.StorageKey("element system state") ))
    , key_nPHS( INDEX<SCALAR,NODE>( db.StorageKey("nodal system state") ))
    , key_Qh( INDEX<SCALAR,ELEMENT>( db.StorageKey("heat source") ))
    , key_NTA( INDEX<ARRAY,NODE>( db.StorageKey("non-wetting phase timing array") ))
    , key_CFL( INDEX<SCALAR,NODE>( db.StorageKey("cfl multiplier") ))
    , key_SYC( INDEX<SCALAR,NODE>( db.StorageKey("synchronize count") ))
    , key_SCC( INDEX<SCALAR,NODE>( db.StorageKey("schedule count") ))
    , key_RAC( INDEX<SCALAR,NODE>( db.StorageKey("rate count") ))
    , key_H2O_comp( INDEX<ARRAY,NODE>( db.StorageKey("composition aqueous phase") ))
    , key_EVI( INDEX<SCALAR,NODE>( db.StorageKey("event index") ))
    , key_RRT( INDEX<SCALAR,ELEMENT>( db.StorageKey("rocktype") ))
    , key_phi( INDEX<SCALAR,ELEMENT>( db.StorageKey("porosity") ))
    , key_k_cal( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability calibration factor") ))
    , key_kV( INDEX<SCALAR,ELEMENT>( db.StorageKey("vertical permeability") ))
    , key_k( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability") ))
    , key_kk( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor permeability") ))
    , key_g( INDEX<SCALAR,MODEL>( db.StorageKey("acceleration gravity") ))
    , key_thi( INDEX<SCALAR,ELEMENT>( db.StorageKey("thickness") ))
    , key_TMAX( INDEX<SCALAR,MODEL>( db.StorageKey("run duration") ))
    , key_cCO2( INDEX<SCALAR,NODE>( db.StorageKey("compressibility carbonic phase") ))
    , key_Kr( INDEX<SCALAR,ELEMENT>( db.StorageKey("thermal conductivity") ))
    , key_KCO2( INDEX<SCALAR,NODE>( db.StorageKey("thermal conductivity carbonic phase") ))
    , key_KH2O( INDEX<SCALAR,NODE>( db.StorageKey("thermal conductivity aqueous phase") ))
    , key_eC( INDEX<SCALAR,ELEMENT>( db.StorageKey("electric conductivity") ))
    , key_muH2O( INDEX<SCALAR,NODE>( db.StorageKey("viscosity aqueous phase") ))
    , key_muCO2( INDEX<SCALAR,NODE>( db.StorageKey("viscosity carbonic phase") ))
    , key_IFT( INDEX<SCALAR,NODE>( db.StorageKey("interfacial tension") ))
    , key_rhoH2O( INDEX<SCALAR,NODE>( db.StorageKey("density aqueous phase") ))
    , key_rhoCO2( INDEX<SCALAR,NODE>( db.StorageKey("density carbonic phase") ))
    , key_cH2O( INDEX<SCALAR,NODE>( db.StorageKey("compressibility aqueous phase") ))
    , key_srCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation carbonic phase") ))
    , key_mTot( INDEX<SCALAR,NODE>( db.StorageKey("total mass") ))
    , key_mNaCl( INDEX<SCALAR,NODE>( db.StorageKey("mass NaCl") ))
    , key_mH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass H2O") ))
    , key_mCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass CO2") ))
    , key_NaClaq( INDEX<SCALAR,NODE>( db.StorageKey("salinity") ))
    , key_NaCl( INDEX<SCALAR,NODE>( db.StorageKey("salt") ))
    , key_H2Og( INDEX<SCALAR,NODE>( db.StorageKey("evaporated water") ))
    , key_CO2aq( INDEX<SCALAR,NODE>( db.StorageKey("dissolved CO2") ))
    , key_CO2_comp( INDEX<ARRAY,NODE>( db.StorageKey("composition carbonic phase") ))
    , key_T( INDEX<SCALAR,NODE>( db.StorageKey("temperature") ))
    , key_srH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation aqueous phase") ))
    , key_dsCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("change of saturation carbonic phase") ))
    , key_sCO2_0( INDEX<SCALAR,NODE>( db.StorageKey("old saturation carbonic phase") ))
    , key_gt( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity term") ))
    , key_QM0( INDEX<SCALAR,ELEMENT>( db.StorageKey("previous fluid mass source") ))
    , key_nQM0( INDEX<SCALAR,NODE>( db.StorageKey("previous nodal fluid mass source") ))
    , key_QM( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid mass source") ))
    , key_nQM( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid mass source") ))
    , key_rpf( INDEX<SCALAR,NODE>( db.StorageKey("reduced fluid pressure") ))
    , key_rpf0( INDEX<SCALAR,NODE>( db.StorageKey("previous reduced fluid pressure") ))
    , key_pf( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_fgt( INDEX<VECTOR,FACE>( db.StorageKey("face gravity term") ))
    , key_Qe( INDEX<SCALAR,ELEMENT>( db.StorageKey("energy source") ))
    , key_cR( INDEX<SCALAR,ELEMENT>( db.StorageKey("compressibility rock") ))
    , key_fv_ctr( INDEX<SCALAR,NODE>( db.StorageKey("FV total rock compressibility") ))
    , key_fv_ctm( INDEX<SCALAR,NODE>( db.StorageKey("FV total system compressibility") ))
    , key_rhod( INDEX<SCALAR,ELEMENT>( db.StorageKey("dry rock density") ))
    , key_diff( INDEX<SCALAR,ELEMENT>( db.StorageKey("CO2 diffusivity") ))
    , key_diffpc( INDEX<SCALAR,ELEMENT>( db.StorageKey("capillary diffusivity") ))
    , key_pd( INDEX<SCALAR,ELEMENT>( db.StorageKey("entry pressure") ))
    , key_flt( INDEX<VECTOR,FACE>( db.StorageKey("face total mobility permeability product") ))
    , key_rhom( INDEX<SCALAR,NODE>( db.StorageKey("mixture density") ))
    , key_vt( INDEX<VECTOR,ELEMENT>( db.StorageKey("total velocity") ))
    , key_SwDrToImb( INDEX<SCALAR,NODE>( db.StorageKey("previous drainage endpoint") ))
    , key_SwImbToDr( INDEX<SCALAR,NODE>( db.StorageKey("previous imbibition endpoint") ))
    , key_psrCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("pseudo residual saturation carbonic phase") ))
    , key_psrH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("pseudo residual saturation aqueous phase") ))
    , key_mNaClaq( INDEX<SCALAR,NODE>( db.StorageKey("mass NaClaq") ))
  {
    if ( key_kri_param.place != ELEMENT || key_kri_param.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'relative permeability model parameters' variable must be ARRAY and placed on ELEMENT"  );
    if ( key_t_out.place != MODEL || key_t_out.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'output times' variable must be ARRAY and placed on MODEL"  );
    if ( key_Nn.place != NODE || key_Nn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'node number' variable must be SCALAR and placed on NODE"  );
    if ( key_dipf.place != FACE || key_dipf.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'face dip vector' variable must be VECTOR and placed on FACE"  );
    if ( key_Ne.place != ELEMENT || key_Ne.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'element number' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dt.place != MODEL || key_dt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'default time increment' variable must be SCALAR and placed on MODEL"  );
    if ( key_fvPV.place != NODE || key_fvPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_sPV.place != SECTOR_INTEGRATION_POINT || key_sPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'sector pore volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( key_dip.place != ELEMENT || key_dip.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'dip vector' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_pf0.place != NODE || key_pf0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_bcp.place != ELEMENT || key_bcp.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'brooks corey parameter' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_t.place != MODEL || key_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'model time' variable must be SCALAR and placed on MODEL"  );
    if ( key_ct.place != ELEMENT || key_ct.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'total system compressibility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_LT.place != ELEMENT || key_LT.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'tensor total mobility permeability product' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_lambda_t.place != ELEMENT || key_lambda_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'total mobility permeability product' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lCO2.place != ELEMENT || key_lCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mobility carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lH2O.place != ELEMENT || key_lH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mobility aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_ssH2O.place != ELEMENT || key_ssH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'shock saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sH2O.place != NODE || key_sH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'saturation aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_sCO2.place != NODE || key_sCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_UPC.place != NODE || key_UPC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'update count' variable must be SCALAR and placed on NODE"  );
    if ( key_eoi.place != MODEL || key_eoi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'end of injection' variable must be SCALAR and placed on MODEL"  );
    if ( key_Xbulk.place != NODE || key_Xbulk.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'bulk mass and fractions' variable must be ARRAY and placed on NODE"  );
    if ( key_ePHS.place != ELEMENT || key_ePHS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'element system state' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nPHS.place != NODE || key_nPHS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'nodal system state' variable must be SCALAR and placed on NODE"  );
    if ( key_Qh.place != ELEMENT || key_Qh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'heat source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_NTA.place != NODE || key_NTA.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'non-wetting phase timing array' variable must be ARRAY and placed on NODE"  );
    if ( key_CFL.place != NODE || key_CFL.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'cfl multiplier' variable must be SCALAR and placed on NODE"  );
    if ( key_SYC.place != NODE || key_SYC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'synchronize count' variable must be SCALAR and placed on NODE"  );
    if ( key_SCC.place != NODE || key_SCC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'schedule count' variable must be SCALAR and placed on NODE"  );
    if ( key_RAC.place != NODE || key_RAC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'rate count' variable must be SCALAR and placed on NODE"  );
    if ( key_H2O_comp.place != NODE || key_H2O_comp.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'composition aqueous phase' variable must be ARRAY and placed on NODE"  );
    if ( key_EVI.place != NODE || key_EVI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'event index' variable must be SCALAR and placed on NODE"  );
    if ( key_RRT.place != ELEMENT || key_RRT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'rocktype' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k_cal.place != ELEMENT || key_k_cal.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'permeability calibration factor' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kV.place != ELEMENT || key_kV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'vertical permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k.place != ELEMENT || key_k.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kk.place != ELEMENT || key_kk.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'tensor permeability' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_g.place != MODEL || key_g.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'acceleration gravity' variable must be SCALAR and placed on MODEL"  );
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_TMAX.place != MODEL || key_TMAX.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'run duration' variable must be SCALAR and placed on MODEL"  );
    if ( key_cCO2.place != NODE || key_cCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'compressibility carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_Kr.place != ELEMENT || key_Kr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'thermal conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KCO2.place != NODE || key_KCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'thermal conductivity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_KH2O.place != NODE || key_KH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'thermal conductivity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_eC.place != ELEMENT || key_eC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'electric conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_muH2O.place != NODE || key_muH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'viscosity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_muCO2.place != NODE || key_muCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'viscosity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_IFT.place != NODE || key_IFT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'interfacial tension' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoH2O.place != NODE || key_rhoH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'density aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoCO2.place != NODE || key_rhoCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'density carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_cH2O.place != NODE || key_cH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'compressibility aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_srCO2.place != ELEMENT || key_srCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'residual saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_mTot.place != NODE || key_mTot.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'total mass' variable must be SCALAR and placed on NODE"  );
    if ( key_mNaCl.place != NODE || key_mNaCl.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mass NaCl' variable must be SCALAR and placed on NODE"  );
    if ( key_mH2O.place != NODE || key_mH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mass H2O' variable must be SCALAR and placed on NODE"  );
    if ( key_mCO2.place != NODE || key_mCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mass CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_NaClaq.place != NODE || key_NaClaq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'salinity' variable must be SCALAR and placed on NODE"  );
    if ( key_NaCl.place != NODE || key_NaCl.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'salt' variable must be SCALAR and placed on NODE"  );
    if ( key_H2Og.place != NODE || key_H2Og.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'evaporated water' variable must be SCALAR and placed on NODE"  );
    if ( key_CO2aq.place != NODE || key_CO2aq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'dissolved CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_CO2_comp.place != NODE || key_CO2_comp.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'composition carbonic phase' variable must be ARRAY and placed on NODE"  );
    if ( key_T.place != NODE || key_T.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'temperature' variable must be SCALAR and placed on NODE"  );
    if ( key_srH2O.place != ELEMENT || key_srH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'residual saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dsCO2.place != ELEMENT || key_dsCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'change of saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sCO2_0.place != NODE || key_sCO2_0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'old saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_gt.place != ELEMENT || key_gt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'gravity term' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_QM0.place != ELEMENT || key_QM0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous fluid mass source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nQM0.place != NODE || key_nQM0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous nodal fluid mass source' variable must be SCALAR and placed on NODE"  );
    if ( key_QM.place != ELEMENT || key_QM.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'fluid mass source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nQM.place != NODE || key_nQM.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'nodal fluid mass source' variable must be SCALAR and placed on NODE"  );
    if ( key_rpf.place != NODE || key_rpf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'reduced fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_rpf0.place != NODE || key_rpf0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous reduced fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_pf.place != NODE || key_pf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_fgt.place != FACE || key_fgt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'face gravity term' variable must be VECTOR and placed on FACE"  );
    if ( key_Qe.place != ELEMENT || key_Qe.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'energy source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_cR.place != ELEMENT || key_cR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'compressibility rock' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fv_ctr.place != NODE || key_fv_ctr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'FV total rock compressibility' variable must be SCALAR and placed on NODE"  );
    if ( key_fv_ctm.place != NODE || key_fv_ctm.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'FV total system compressibility' variable must be SCALAR and placed on NODE"  );
    if ( key_rhod.place != ELEMENT || key_rhod.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'dry rock density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diff.place != ELEMENT || key_diff.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'CO2 diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diffpc.place != ELEMENT || key_diffpc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'capillary diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_pd.place != ELEMENT || key_pd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'entry pressure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_flt.place != FACE || key_flt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'face total mobility permeability product' variable must be VECTOR and placed on FACE"  );
    if ( key_rhom.place != NODE || key_rhom.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mixture density' variable must be SCALAR and placed on NODE"  );
    if ( key_vt.place != ELEMENT || key_vt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'total velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_SwDrToImb.place != NODE || key_SwDrToImb.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous drainage endpoint' variable must be SCALAR and placed on NODE"  );
    if ( key_SwImbToDr.place != NODE || key_SwImbToDr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'previous imbibition endpoint' variable must be SCALAR and placed on NODE"  );
    if ( key_psrCO2.place != ELEMENT || key_psrCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'pseudo residual saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_psrH2O.place != ELEMENT || key_psrH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'pseudo residual saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_mNaClaq.place != NODE || key_mNaClaq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration::VariableSet_CO2GeoSequestration:",
        "The 'mass NaClaq' variable must be SCALAR and placed on NODE"  );
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLESET_CO2GEOSEQUESTRATION_H
