#ifndef CSMP_VARIABLESET_CO2GEOSEQUESTRATION-FULL_MULTIPHYSICS_SET_H
#define CSMP_VARIABLESET_CO2GEOSEQUESTRATION-FULL_MULTIPHYSICS_SET_H

/**
@file VariableSet_CO2GeoSequestration-full_multiphysics_set.h
Automatically generated from VariableSet_CO2GeoSequestration-full_multiphysics_set.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct VariableSet_CO2GeoSequestration-full_multiphysics_set {
  // chemistry
  csmp::INDEX<ARRAY,NODE> key_CO2_comp; // composition carbonic phase
  csmp::INDEX<ARRAY,NODE> key_H2O_comp; // composition aqueous phase
  csmp::INDEX<ARRAY,NODE> key_Xbulk; // bulk mass and fractions
  csmp::INDEX<SCALAR,NODE> key_NaCl; // salt
  // fluid properties
  csmp::INDEX<SCALAR,NODE> key_muCO2; // viscosity carbonic phase
  csmp::INDEX<SCALAR,NODE> key_rhom; // mixture density
  csmp::INDEX<SCALAR,NODE> key_muH2O; // viscosity aqueous phase
  csmp::INDEX<SCALAR,NODE> key_cH2O; // compressibility aqueous phase
  csmp::INDEX<SCALAR,NODE> key_IFT; // interfacial tension
  csmp::INDEX<SCALAR,NODE> key_rhoCO2; // density carbonic phase
  csmp::INDEX<SCALAR,NODE> key_rhoH2O; // density aqueous phase
  csmp::INDEX<SCALAR,NODE> key_cCO2; // compressibility carbonic phase
  // geophysics
  csmp::INDEX<SCALAR,ELEMENT> key_eC; // electric conductivity
  // heat flow
  csmp::INDEX<SCALAR,NODE> key_KH2O; // thermal conductivity aqueous phase
  // heat transfer
  csmp::INDEX<SCALAR,ELEMENT> key_Kr; // thermal conductivity
  csmp::INDEX<SCALAR,NODE> key_T; // temperature
  csmp::INDEX<SCALAR,ELEMENT> key_Qe; // energy source
  csmp::INDEX<SCALAR,NODE> key_KCO2; // thermal conductivity carbonic phase
  // mechanics
  csmp::INDEX<TENSOR,ELEMENT> key_S; // stress
  csmp::INDEX<SCALAR,ELEMENT> key_Sy; // stress-y
  csmp::INDEX<SCALAR,ELEMENT> key_dV; // volume change
  csmp::INDEX<TENSOR,ELEMENT> key_e; // strain
  csmp::INDEX<VECTOR,ELEMENT> key_sig1; // sigma1
  csmp::INDEX<VECTOR,ELEMENT> key_sig3; // sigma3
  csmp::INDEX<VECTOR,ELEMENT> key_sig2; // sigma2
  csmp::INDEX<SCALAR,ELEMENT> key_Bio; // Biot alpha
  csmp::INDEX<SCALAR,ELEMENT> key_BioT; // Biot term
  csmp::INDEX<SCALAR,ELEMENT> key_nu; // Poissons ratio
  csmp::INDEX<SCALAR,ELEMENT> key_E; // Youngs modulus
  csmp::INDEX<SCALAR,ELEMENT> key_rhob; // bulk density
  csmp::INDEX<SCALAR,ELEMENT> key_coh; // cohesion
  csmp::INDEX<SCALAR,ELEMENT> key_cR; // compressibility rock
  csmp::INDEX<SCALAR,ELEMENT> key_dil; // dilatation
  csmp::INDEX<VECTOR,NODE> key_displ; // displacement
  csmp::INDEX<SCALAR,ELEMENT> key_rhod; // dry rock density
  csmp::INDEX<SCALAR,REGION> key_Se; // effective stress
  csmp::INDEX<SCALAR,ELEMENT_INTEGRATION_POINT> key_fail; // failure
  csmp::INDEX<VECTOR,NODE> key_F; // force
  csmp::INDEX<SCALAR,ELEMENT> key_frica; // friction angle
  csmp::INDEX<SCALAR,ELEMENT> key_SSmax; // max shear stress
  csmp::INDEX<SCALAR,ELEMENT> key_Smean; // mean stress
  csmp::INDEX<SCALAR,FACE> key_Sn; // normal stress
  csmp::INDEX<SCALAR,NODE> key_op; // overburden pressure
  csmp::INDEX<SCALAR,MODEL> key_Sv; // overburden stress
  csmp::INDEX<VECTOR,FACE> key_Ss; // shear stress
  // multiphase flow
  csmp::INDEX<SCALAR,NODE> key_fb; // flux balance
  csmp::INDEX<SCALAR,ELEMENT> key_psrH2O; // pseudo residual saturation aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_diff; // CO2 diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_ePHS; // element system state
  csmp::INDEX<SCALAR,NODE> key_nPHS; // nodal system state
  csmp::INDEX<SCALAR,ELEMENT> key_muH2Obc; // barcycenter aqueous phase viscosity
  csmp::INDEX<SCALAR,ELEMENT> key_muCO2bc; // barycenter carbonic phase viscosity
  csmp::INDEX<SCALAR,ELEMENT> key_kV; // vertical permeability
  csmp::INDEX<SCALAR,ELEMENT> key_diffpc; // capillary diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_pd; // entry pressure
  csmp::INDEX<VECTOR,FACE> key_fgt; // face gravity term
  csmp::INDEX<VECTOR,FACE> key_flt; // face total mobility permeability product
  csmp::INDEX<SCALAR,NODE> key_pf; // fluid pressure
  csmp::INDEX<SCALAR,NODE> key_rpf; // reduced fluid pressure
  csmp::INDEX<SCALAR,NODE> key_mum; // mixture viscosity
  csmp::INDEX<SCALAR,NODE> key_nQV; // nodal fluid volume source
  csmp::INDEX<SCALAR,NODE> key_nQM; // nodal fluid mass source
  csmp::INDEX<SCALAR,ELEMENT> key_QV; // fluid volume source
  csmp::INDEX<SCALAR,ELEMENT> key_QM; // fluid mass source
  csmp::INDEX<SCALAR,ELEMENT> key_rhoH2Obc; // barycenter aqueous phase density
  csmp::INDEX<VECTOR,ELEMENT> key_gt; // gravity term
  csmp::INDEX<SCALAR,NODE> key_sCO2_0; // old saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_srCO2; // residual saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_srH2O; // residual saturation aqueous phase
  csmp::INDEX<SCALAR,NODE> key_sCO2; // saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_sH2O; // saturation aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_ssH2O; // shock saturation aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_dsCO2; // change of saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_TDS; // total dissolved solids
  csmp::INDEX<SCALAR,ELEMENT> key_lH2O; // mobility aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_lCO2; // mobility carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_lambda_t; // total mobility permeability product
  csmp::INDEX<TENSOR,ELEMENT> key_LT; // tensor total mobility permeability product
  csmp::INDEX<SCALAR,ELEMENT> key_ct; // total system compressibility
  csmp::INDEX<VECTOR,ELEMENT> key_vt; // total velocity
  csmp::INDEX<SCALAR,ELEMENT> key_qf; // volume flux
  csmp::INDEX<SCALAR,MODEL> key_t; // model time
  csmp::INDEX<SCALAR,ELEMENT> key_bcp; // brooks corey parameter
  csmp::INDEX<VECTOR,ELEMENT> key_Fg; // gravity force
  csmp::INDEX<SCALAR,NODE> key_SwDrToImb; // previous drainage endpoint
  csmp::INDEX<SCALAR,NODE> key_SwImbToDr; // previous imbibition endpoint
  csmp::INDEX<SCALAR,ELEMENT> key_psrCO2; // pseudo residual saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_rhoCO2bc; // barycenter carbonic phase density
  csmp::INDEX<ARRAY,ELEMENT> key_kri_param; // relative permeability model parameters
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_FCAP; // capillary flux CO2 density area product
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_FGRAV; // gravitational CO2 flux density area product
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_FVISC; // viscous flux density area product
  // numeric modelling
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; // facet normal
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qft; // facet total flux
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qfCO2; // facet flux carbonic phase
  csmp::INDEX<SCALAR,NODE> key_fvPV; // FV pore volume
  csmp::INDEX<SCALAR,MODEL> key_dt; // default time increment
  csmp::INDEX<SCALAR,ELEMENT> key_Ne; // element number
  csmp::INDEX<VECTOR,FACE> key_dipf; // face dip vector
  csmp::INDEX<SCALAR,NODE> key_Nn; // node number
  csmp::INDEX<ARRAY,MODEL> key_t_out; // output times
  csmp::INDEX<SCALAR,NODE> key_pf0; // previous fluid pressure
  csmp::INDEX<SCALAR,ELEMENT> key_QV0; // previous fluid volume source
  csmp::INDEX<SCALAR,NODE> key_nQV0; // previous nodal fluid volume source
  csmp::INDEX<SCALAR,MODEL> key_TMAX; // run duration
  csmp::INDEX<SCALAR,ELEMENT> key_thi; // thickness
  csmp::INDEX<SCALAR,ELEMENT> key_V; // volume
  csmp::INDEX<SCALAR,MODEL> key_eoi; // end of injection
  csmp::INDEX<SCALAR,ELEMENT> key_QM0; // previous fluid mass source
  csmp::INDEX<VECTOR,ELEMENT> key_dip; // dip vector
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_sPV; // sector pore volume
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA; // facet area
  // output
  csmp::INDEX<SCALAR,NODE> key_mNaClaq; // mass NaClaq
  csmp::INDEX<SCALAR,NODE> key_mNaClsd; // mass NaClsd
  csmp::INDEX<SCALAR,NODE> key_CO2aq; // dissolved CO2
  csmp::INDEX<SCALAR,NODE> key_mCO2; // mass CO2
  csmp::INDEX<SCALAR,NODE> key_mH2O; // mass H2O
  csmp::INDEX<SCALAR,NODE> key_NaClaq; // salinity
  csmp::INDEX<SCALAR,NODE> key_H2Og; // evaporated water
  // physical constants
  csmp::INDEX<SCALAR,MODEL> key_g; // acceleration gravity
  // property modelling
  csmp::INDEX<SCALAR,ELEMENT> key_kFbr; // fault breccia permeability
  csmp::INDEX<SCALAR,ELEMENT> key_k_cal; // permeability calibration factor
  csmp::INDEX<SCALAR,REGION> key_Fs; // fault size
  csmp::INDEX<SCALAR,ELEMENT> key_kf; // fracture permeability
  csmp::INDEX<SCALAR,ELEMENT> key_k; // permeability
  csmp::INDEX<SCALAR,REGION> key_AS; // slip patch area
  csmp::INDEX<SCALAR,ELEMENT> key_RRT; // rocktype
  csmp::INDEX<SCALAR,ELEMENT> key_phi; // porosity
  csmp::INDEX<TENSOR,ELEMENT> key_kk; // tensor permeability

  template<uint32_t dim>
  explicit VariableSet_CO2GeoSequestration-full_multiphysics_set( const PropertyDatabase<dim>& db )
    : key_CO2_comp( INDEX<ARRAY,NODE>( db.StorageKey("composition carbonic phase") ))
    , key_H2O_comp( INDEX<ARRAY,NODE>( db.StorageKey("composition aqueous phase") ))
    , key_Xbulk( INDEX<ARRAY,NODE>( db.StorageKey("bulk mass and fractions") ))
    , key_NaCl( INDEX<SCALAR,NODE>( db.StorageKey("salt") ))
    , key_muCO2( INDEX<SCALAR,NODE>( db.StorageKey("viscosity carbonic phase") ))
    , key_rhom( INDEX<SCALAR,NODE>( db.StorageKey("mixture density") ))
    , key_muH2O( INDEX<SCALAR,NODE>( db.StorageKey("viscosity aqueous phase") ))
    , key_cH2O( INDEX<SCALAR,NODE>( db.StorageKey("compressibility aqueous phase") ))
    , key_IFT( INDEX<SCALAR,NODE>( db.StorageKey("interfacial tension") ))
    , key_rhoCO2( INDEX<SCALAR,NODE>( db.StorageKey("density carbonic phase") ))
    , key_rhoH2O( INDEX<SCALAR,NODE>( db.StorageKey("density aqueous phase") ))
    , key_cCO2( INDEX<SCALAR,NODE>( db.StorageKey("compressibility carbonic phase") ))
    , key_eC( INDEX<SCALAR,ELEMENT>( db.StorageKey("electric conductivity") ))
    , key_KH2O( INDEX<SCALAR,NODE>( db.StorageKey("thermal conductivity aqueous phase") ))
    , key_Kr( INDEX<SCALAR,ELEMENT>( db.StorageKey("thermal conductivity") ))
    , key_T( INDEX<SCALAR,NODE>( db.StorageKey("temperature") ))
    , key_Qe( INDEX<SCALAR,ELEMENT>( db.StorageKey("energy source") ))
    , key_KCO2( INDEX<SCALAR,NODE>( db.StorageKey("thermal conductivity carbonic phase") ))
    , key_S( INDEX<TENSOR,ELEMENT>( db.StorageKey("stress") ))
    , key_Sy( INDEX<SCALAR,ELEMENT>( db.StorageKey("stress-y") ))
    , key_dV( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume change") ))
    , key_e( INDEX<TENSOR,ELEMENT>( db.StorageKey("strain") ))
    , key_sig1( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma1") ))
    , key_sig3( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma3") ))
    , key_sig2( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma2") ))
    , key_Bio( INDEX<SCALAR,ELEMENT>( db.StorageKey("Biot alpha") ))
    , key_BioT( INDEX<SCALAR,ELEMENT>( db.StorageKey("Biot term") ))
    , key_nu( INDEX<SCALAR,ELEMENT>( db.StorageKey("Poissons ratio") ))
    , key_E( INDEX<SCALAR,ELEMENT>( db.StorageKey("Youngs modulus") ))
    , key_rhob( INDEX<SCALAR,ELEMENT>( db.StorageKey("bulk density") ))
    , key_coh( INDEX<SCALAR,ELEMENT>( db.StorageKey("cohesion") ))
    , key_cR( INDEX<SCALAR,ELEMENT>( db.StorageKey("compressibility rock") ))
    , key_dil( INDEX<SCALAR,ELEMENT>( db.StorageKey("dilatation") ))
    , key_displ( INDEX<VECTOR,NODE>( db.StorageKey("displacement") ))
    , key_rhod( INDEX<SCALAR,ELEMENT>( db.StorageKey("dry rock density") ))
    , key_Se( INDEX<SCALAR,REGION>( db.StorageKey("effective stress") ))
    , key_fail( INDEX<SCALAR,ELEMENT_INTEGRATION_POINT>( db.StorageKey("failure") ))
    , key_F( INDEX<VECTOR,NODE>( db.StorageKey("force") ))
    , key_frica( INDEX<SCALAR,ELEMENT>( db.StorageKey("friction angle") ))
    , key_SSmax( INDEX<SCALAR,ELEMENT>( db.StorageKey("max shear stress") ))
    , key_Smean( INDEX<SCALAR,ELEMENT>( db.StorageKey("mean stress") ))
    , key_Sn( INDEX<SCALAR,FACE>( db.StorageKey("normal stress") ))
    , key_op( INDEX<SCALAR,NODE>( db.StorageKey("overburden pressure") ))
    , key_Sv( INDEX<SCALAR,MODEL>( db.StorageKey("overburden stress") ))
    , key_Ss( INDEX<VECTOR,FACE>( db.StorageKey("shear stress") ))
    , key_fb( INDEX<SCALAR,NODE>( db.StorageKey("flux balance") ))
    , key_psrH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("pseudo residual saturation aqueous phase") ))
    , key_diff( INDEX<SCALAR,ELEMENT>( db.StorageKey("CO2 diffusivity") ))
    , key_ePHS( INDEX<SCALAR,ELEMENT>( db.StorageKey("element system state") ))
    , key_nPHS( INDEX<SCALAR,NODE>( db.StorageKey("nodal system state") ))
    , key_muH2Obc( INDEX<SCALAR,ELEMENT>( db.StorageKey("barcycenter aqueous phase viscosity") ))
    , key_muCO2bc( INDEX<SCALAR,ELEMENT>( db.StorageKey("barycenter carbonic phase viscosity") ))
    , key_kV( INDEX<SCALAR,ELEMENT>( db.StorageKey("vertical permeability") ))
    , key_diffpc( INDEX<SCALAR,ELEMENT>( db.StorageKey("capillary diffusivity") ))
    , key_pd( INDEX<SCALAR,ELEMENT>( db.StorageKey("entry pressure") ))
    , key_fgt( INDEX<VECTOR,FACE>( db.StorageKey("face gravity term") ))
    , key_flt( INDEX<VECTOR,FACE>( db.StorageKey("face total mobility permeability product") ))
    , key_pf( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_rpf( INDEX<SCALAR,NODE>( db.StorageKey("reduced fluid pressure") ))
    , key_mum( INDEX<SCALAR,NODE>( db.StorageKey("mixture viscosity") ))
    , key_nQV( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid volume source") ))
    , key_nQM( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid mass source") ))
    , key_QV( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid volume source") ))
    , key_QM( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid mass source") ))
    , key_rhoH2Obc( INDEX<SCALAR,ELEMENT>( db.StorageKey("barycenter aqueous phase density") ))
    , key_gt( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity term") ))
    , key_sCO2_0( INDEX<SCALAR,NODE>( db.StorageKey("old saturation carbonic phase") ))
    , key_srCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation carbonic phase") ))
    , key_srH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation aqueous phase") ))
    , key_sCO2( INDEX<SCALAR,NODE>( db.StorageKey("saturation carbonic phase") ))
    , key_sH2O( INDEX<SCALAR,NODE>( db.StorageKey("saturation aqueous phase") ))
    , key_ssH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("shock saturation aqueous phase") ))
    , key_dsCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("change of saturation carbonic phase") ))
    , key_TDS( INDEX<SCALAR,NODE>( db.StorageKey("total dissolved solids") ))
    , key_lH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility aqueous phase") ))
    , key_lCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility carbonic phase") ))
    , key_lambda_t( INDEX<SCALAR,ELEMENT>( db.StorageKey("total mobility permeability product") ))
    , key_LT( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor total mobility permeability product") ))
    , key_ct( INDEX<SCALAR,ELEMENT>( db.StorageKey("total system compressibility") ))
    , key_vt( INDEX<VECTOR,ELEMENT>( db.StorageKey("total velocity") ))
    , key_qf( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume flux") ))
    , key_t( INDEX<SCALAR,MODEL>( db.StorageKey("model time") ))
    , key_bcp( INDEX<SCALAR,ELEMENT>( db.StorageKey("brooks corey parameter") ))
    , key_Fg( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity force") ))
    , key_SwDrToImb( INDEX<SCALAR,NODE>( db.StorageKey("previous drainage endpoint") ))
    , key_SwImbToDr( INDEX<SCALAR,NODE>( db.StorageKey("previous imbibition endpoint") ))
    , key_psrCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("pseudo residual saturation carbonic phase") ))
    , key_rhoCO2bc( INDEX<SCALAR,ELEMENT>( db.StorageKey("barycenter carbonic phase density") ))
    , key_kri_param( INDEX<ARRAY,ELEMENT>( db.StorageKey("relative permeability model parameters") ))
    , key_FCAP( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("capillary flux CO2 density area product") ))
    , key_FGRAV( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("gravitational CO2 flux density area product") ))
    , key_FVISC( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("viscous flux density area product") ))
    , key_fn( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal") ))
    , key_qft( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet total flux") ))
    , key_qfCO2( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux carbonic phase") ))
    , key_fvPV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_dt( INDEX<SCALAR,MODEL>( db.StorageKey("default time increment") ))
    , key_Ne( INDEX<SCALAR,ELEMENT>( db.StorageKey("element number") ))
    , key_dipf( INDEX<VECTOR,FACE>( db.StorageKey("face dip vector") ))
    , key_Nn( INDEX<SCALAR,NODE>( db.StorageKey("node number") ))
    , key_t_out( INDEX<ARRAY,MODEL>( db.StorageKey("output times") ))
    , key_pf0( INDEX<SCALAR,NODE>( db.StorageKey("previous fluid pressure") ))
    , key_QV0( INDEX<SCALAR,ELEMENT>( db.StorageKey("previous fluid volume source") ))
    , key_nQV0( INDEX<SCALAR,NODE>( db.StorageKey("previous nodal fluid volume source") ))
    , key_TMAX( INDEX<SCALAR,MODEL>( db.StorageKey("run duration") ))
    , key_thi( INDEX<SCALAR,ELEMENT>( db.StorageKey("thickness") ))
    , key_V( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume") ))
    , key_eoi( INDEX<SCALAR,MODEL>( db.StorageKey("end of injection") ))
    , key_QM0( INDEX<SCALAR,ELEMENT>( db.StorageKey("previous fluid mass source") ))
    , key_dip( INDEX<VECTOR,ELEMENT>( db.StorageKey("dip vector") ))
    , key_sPV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector pore volume") ))
    , key_fA( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area") ))
    , key_mNaClaq( INDEX<SCALAR,NODE>( db.StorageKey("mass NaClaq") ))
    , key_mNaClsd( INDEX<SCALAR,NODE>( db.StorageKey("mass NaClsd") ))
    , key_CO2aq( INDEX<SCALAR,NODE>( db.StorageKey("dissolved CO2") ))
    , key_mCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass CO2") ))
    , key_mH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass H2O") ))
    , key_NaClaq( INDEX<SCALAR,NODE>( db.StorageKey("salinity") ))
    , key_H2Og( INDEX<SCALAR,NODE>( db.StorageKey("evaporated water") ))
    , key_g( INDEX<SCALAR,MODEL>( db.StorageKey("acceleration gravity") ))
    , key_kFbr( INDEX<SCALAR,ELEMENT>( db.StorageKey("fault breccia permeability") ))
    , key_k_cal( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability calibration factor") ))
    , key_Fs( INDEX<SCALAR,REGION>( db.StorageKey("fault size") ))
    , key_kf( INDEX<SCALAR,ELEMENT>( db.StorageKey("fracture permeability") ))
    , key_k( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability") ))
    , key_AS( INDEX<SCALAR,REGION>( db.StorageKey("slip patch area") ))
    , key_RRT( INDEX<SCALAR,ELEMENT>( db.StorageKey("rocktype") ))
    , key_phi( INDEX<SCALAR,ELEMENT>( db.StorageKey("porosity") ))
    , key_kk( INDEX<TENSOR,ELEMENT>( db.StorageKey("tensor permeability") ))
  {
    if ( key_CO2_comp.place != NODE || key_CO2_comp.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'composition carbonic phase' variable must be ARRAY and placed on NODE"  );
    if ( key_H2O_comp.place != NODE || key_H2O_comp.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'composition aqueous phase' variable must be ARRAY and placed on NODE"  );
    if ( key_Xbulk.place != NODE || key_Xbulk.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'bulk mass and fractions' variable must be ARRAY and placed on NODE"  );
    if ( key_NaCl.place != NODE || key_NaCl.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'salt' variable must be SCALAR and placed on NODE"  );
    if ( key_muCO2.place != NODE || key_muCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'viscosity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhom.place != NODE || key_rhom.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mixture density' variable must be SCALAR and placed on NODE"  );
    if ( key_muH2O.place != NODE || key_muH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'viscosity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_cH2O.place != NODE || key_cH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'compressibility aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_IFT.place != NODE || key_IFT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'interfacial tension' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoCO2.place != NODE || key_rhoCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'density carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoH2O.place != NODE || key_rhoH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'density aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_cCO2.place != NODE || key_cCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'compressibility carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_eC.place != ELEMENT || key_eC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'electric conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KH2O.place != NODE || key_KH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'thermal conductivity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_Kr.place != ELEMENT || key_Kr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'thermal conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_T.place != NODE || key_T.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'temperature' variable must be SCALAR and placed on NODE"  );
    if ( key_Qe.place != ELEMENT || key_Qe.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'energy source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KCO2.place != NODE || key_KCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'thermal conductivity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_S.place != ELEMENT || key_S.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'stress' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_Sy.place != ELEMENT || key_Sy.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'stress-y' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dV.place != ELEMENT || key_dV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'volume change' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_e.place != ELEMENT || key_e.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'strain' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_sig1.place != ELEMENT || key_sig1.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'sigma1' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sig3.place != ELEMENT || key_sig3.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'sigma3' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sig2.place != ELEMENT || key_sig2.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'sigma2' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_Bio.place != ELEMENT || key_Bio.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'Biot alpha' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_BioT.place != ELEMENT || key_BioT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'Biot term' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nu.place != ELEMENT || key_nu.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'Poissons ratio' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_E.place != ELEMENT || key_E.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'Youngs modulus' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_rhob.place != ELEMENT || key_rhob.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'bulk density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_coh.place != ELEMENT || key_coh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'cohesion' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_cR.place != ELEMENT || key_cR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'compressibility rock' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dil.place != ELEMENT || key_dil.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'dilatation' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_displ.place != NODE || key_displ.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'displacement' variable must be VECTOR and placed on NODE"  );
    if ( key_rhod.place != ELEMENT || key_rhod.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'dry rock density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Se.place != REGION || key_Se.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'effective stress' variable must be SCALAR and placed on REGION"  );
    if ( key_fail.place != ELEMENT_INTEGRATION_POINT || key_fail.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'failure' variable must be SCALAR and placed on ELEMENT_INTEGRATION_POINT"  );
    if ( key_F.place != NODE || key_F.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'force' variable must be VECTOR and placed on NODE"  );
    if ( key_frica.place != ELEMENT || key_frica.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'friction angle' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_SSmax.place != ELEMENT || key_SSmax.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'max shear stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Smean.place != ELEMENT || key_Smean.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mean stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Sn.place != FACE || key_Sn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'normal stress' variable must be SCALAR and placed on FACE"  );
    if ( key_op.place != NODE || key_op.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'overburden pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_Sv.place != MODEL || key_Sv.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'overburden stress' variable must be SCALAR and placed on MODEL"  );
    if ( key_Ss.place != FACE || key_Ss.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'shear stress' variable must be VECTOR and placed on FACE"  );
    if ( key_fb.place != NODE || key_fb.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( key_psrH2O.place != ELEMENT || key_psrH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'pseudo residual saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diff.place != ELEMENT || key_diff.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'CO2 diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_ePHS.place != ELEMENT || key_ePHS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'element system state' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nPHS.place != NODE || key_nPHS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'nodal system state' variable must be SCALAR and placed on NODE"  );
    if ( key_muH2Obc.place != ELEMENT || key_muH2Obc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'barcycenter aqueous phase viscosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_muCO2bc.place != ELEMENT || key_muCO2bc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'barycenter carbonic phase viscosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kV.place != ELEMENT || key_kV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'vertical permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diffpc.place != ELEMENT || key_diffpc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'capillary diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_pd.place != ELEMENT || key_pd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'entry pressure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fgt.place != FACE || key_fgt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'face gravity term' variable must be VECTOR and placed on FACE"  );
    if ( key_flt.place != FACE || key_flt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'face total mobility permeability product' variable must be VECTOR and placed on FACE"  );
    if ( key_pf.place != NODE || key_pf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_rpf.place != NODE || key_rpf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'reduced fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_mum.place != NODE || key_mum.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mixture viscosity' variable must be SCALAR and placed on NODE"  );
    if ( key_nQV.place != NODE || key_nQV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( key_nQM.place != NODE || key_nQM.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'nodal fluid mass source' variable must be SCALAR and placed on NODE"  );
    if ( key_QV.place != ELEMENT || key_QV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_QM.place != ELEMENT || key_QM.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fluid mass source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_rhoH2Obc.place != ELEMENT || key_rhoH2Obc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'barycenter aqueous phase density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_gt.place != ELEMENT || key_gt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'gravity term' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sCO2_0.place != NODE || key_sCO2_0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'old saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_srCO2.place != ELEMENT || key_srCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'residual saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_srH2O.place != ELEMENT || key_srH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'residual saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sCO2.place != NODE || key_sCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_sH2O.place != NODE || key_sH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'saturation aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_ssH2O.place != ELEMENT || key_ssH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'shock saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dsCO2.place != ELEMENT || key_dsCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'change of saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_TDS.place != NODE || key_TDS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'total dissolved solids' variable must be SCALAR and placed on NODE"  );
    if ( key_lH2O.place != ELEMENT || key_lH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mobility aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lCO2.place != ELEMENT || key_lCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mobility carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lambda_t.place != ELEMENT || key_lambda_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'total mobility permeability product' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_LT.place != ELEMENT || key_LT.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'tensor total mobility permeability product' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_ct.place != ELEMENT || key_ct.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'total system compressibility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_vt.place != ELEMENT || key_vt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'total velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_qf.place != ELEMENT || key_qf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'volume flux' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_t.place != MODEL || key_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'model time' variable must be SCALAR and placed on MODEL"  );
    if ( key_bcp.place != ELEMENT || key_bcp.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'brooks corey parameter' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Fg.place != ELEMENT || key_Fg.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'gravity force' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_SwDrToImb.place != NODE || key_SwDrToImb.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous drainage endpoint' variable must be SCALAR and placed on NODE"  );
    if ( key_SwImbToDr.place != NODE || key_SwImbToDr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous imbibition endpoint' variable must be SCALAR and placed on NODE"  );
    if ( key_psrCO2.place != ELEMENT || key_psrCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'pseudo residual saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_rhoCO2bc.place != ELEMENT || key_rhoCO2bc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'barycenter carbonic phase density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kri_param.place != ELEMENT || key_kri_param.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'relative permeability model parameters' variable must be ARRAY and placed on ELEMENT"  );
    if ( key_FCAP.place != FACET_INTEGRATION_POINT || key_FCAP.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'capillary flux CO2 density area product' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_FGRAV.place != FACET_INTEGRATION_POINT || key_FGRAV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'gravitational CO2 flux density area product' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_FVISC.place != FACET_INTEGRATION_POINT || key_FVISC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'viscous flux density area product' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qft.place != FACET_INTEGRATION_POINT || key_qft.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'facet total flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qfCO2.place != FACET_INTEGRATION_POINT || key_qfCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'facet flux carbonic phase' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fvPV.place != NODE || key_fvPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( key_dt.place != MODEL || key_dt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'default time increment' variable must be SCALAR and placed on MODEL"  );
    if ( key_Ne.place != ELEMENT || key_Ne.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'element number' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dipf.place != FACE || key_dipf.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'face dip vector' variable must be VECTOR and placed on FACE"  );
    if ( key_Nn.place != NODE || key_Nn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'node number' variable must be SCALAR and placed on NODE"  );
    if ( key_t_out.place != MODEL || key_t_out.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'output times' variable must be ARRAY and placed on MODEL"  );
    if ( key_pf0.place != NODE || key_pf0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_QV0.place != ELEMENT || key_QV0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nQV0.place != NODE || key_nQV0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( key_TMAX.place != MODEL || key_TMAX.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'run duration' variable must be SCALAR and placed on MODEL"  );
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_V.place != ELEMENT || key_V.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'volume' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_eoi.place != MODEL || key_eoi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'end of injection' variable must be SCALAR and placed on MODEL"  );
    if ( key_QM0.place != ELEMENT || key_QM0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'previous fluid mass source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dip.place != ELEMENT || key_dip.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'dip vector' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sPV.place != SECTOR_INTEGRATION_POINT || key_sPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'sector pore volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_mNaClaq.place != NODE || key_mNaClaq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mass NaClaq' variable must be SCALAR and placed on NODE"  );
    if ( key_mNaClsd.place != NODE || key_mNaClsd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mass NaClsd' variable must be SCALAR and placed on NODE"  );
    if ( key_CO2aq.place != NODE || key_CO2aq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'dissolved CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_mCO2.place != NODE || key_mCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mass CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_mH2O.place != NODE || key_mH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'mass H2O' variable must be SCALAR and placed on NODE"  );
    if ( key_NaClaq.place != NODE || key_NaClaq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'salinity' variable must be SCALAR and placed on NODE"  );
    if ( key_H2Og.place != NODE || key_H2Og.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'evaporated water' variable must be SCALAR and placed on NODE"  );
    if ( key_g.place != MODEL || key_g.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'acceleration gravity' variable must be SCALAR and placed on MODEL"  );
    if ( key_kFbr.place != ELEMENT || key_kFbr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fault breccia permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k_cal.place != ELEMENT || key_k_cal.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'permeability calibration factor' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Fs.place != REGION || key_Fs.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fault size' variable must be SCALAR and placed on REGION"  );
    if ( key_kf.place != ELEMENT || key_kf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'fracture permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k.place != ELEMENT || key_k.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_AS.place != REGION || key_AS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'slip patch area' variable must be SCALAR and placed on REGION"  );
    if ( key_RRT.place != ELEMENT || key_RRT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'rocktype' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kk.place != ELEMENT || key_kk.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_CO2GeoSequestration-full_multiphysics_set::VariableSet_CO2GeoSequestration-full_multiphysics_set:",
        "The 'tensor permeability' variable must be TENSOR and placed on ELEMENT"  );
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLESET_CO2GEOSEQUESTRATION-FULL_MULTIPHYSICS_SET_H
