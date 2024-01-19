#ifndef CSMP_VARIABLESET_TWOPHASEFLOW_H
#define CSMP_VARIABLESET_TWOPHASEFLOW_H

/**
@file VariableSet_TwoPhaseFlow.h
Automatically generated from VariableSet_TwoPhaseFlow.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct VariableSet_TwoPhaseFlow {
  // chemistry
  csmp::INDEX<SCALAR,NODE> key_XCO2; // mass fraction CO2 aqueous phase
  csmp::INDEX<SCALAR,NODE> key_XH2O; // mass fraction H2O aqueous phase
  csmp::INDEX<SCALAR,NODE> key_YCO2; // mass fraction CO2 carbonic phase
  csmp::INDEX<SCALAR,NODE> key_YH2O; // mass fraction H2O carbonic phase
  csmp::INDEX<SCALAR,NODE> key_xCO2; // mole fraction CO2 aqueous phase
  csmp::INDEX<SCALAR,NODE> key_xH2O; // mole fraction H2O aqueous phase
  csmp::INDEX<SCALAR,NODE> key_xSalt; // mole fraction NaCl aqueous phase
  csmp::INDEX<SCALAR,NODE> key_yCO2; // mole fraction CO2 carbonic phase
  csmp::INDEX<SCALAR,NODE> key_yH2O; // mole fraction H2O carbonic phase
  csmp::INDEX<SCALAR,NODE> key_CO2aq; // dissolved CO2
  csmp::INDEX<SCALAR,NODE> key_H2Og; // evaporated water
  csmp::INDEX<SCALAR,NODE> key_msalt; // molality salt
  // fluid properties
  csmp::INDEX<SCALAR,NODE> key_cCO2; // compressibility carbonic phase
  csmp::INDEX<SCALAR,NODE> key_cH2O; // compressibility aqueous phase
  csmp::INDEX<SCALAR,NODE> key_rhoCO2; // density carbonic phase
  csmp::INDEX<SCALAR,NODE> key_rhoH2O; // density aqueous phase
  csmp::INDEX<SCALAR,NODE> key_diffCO2; // diffusivity coefficient carbonic phase
  csmp::INDEX<SCALAR,NODE> key_diffH2O; // diffusivity coefficient aqueous phase
  csmp::INDEX<SCALAR,NODE> key_rhom; // density mixture
  csmp::INDEX<SCALAR,MODEL> key_IFT; // interfacial tension
  csmp::INDEX<SCALAR,NODE> key_muCO2; // viscosity carbonic phase
  csmp::INDEX<SCALAR,NODE> key_muH2O; // viscosity aqueous phase
  // geophysics
  csmp::INDEX<SCALAR,ELEMENT> key_eC; // electric conductivity
  // heat flow
  csmp::INDEX<SCALAR,MODEL> key_KH2O; // thermal conductivity aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_Qh; // heat source
  // heat transfer
  csmp::INDEX<SCALAR,MODEL> key_KCO2; // thermal conductivity carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_Kr; // thermal conductivity
  csmp::INDEX<SCALAR,ELEMENT> key_Qe; // energy source
  csmp::INDEX<SCALAR,NODE> key_T; // temperature
  // mechanics
  csmp::INDEX<SCALAR,ELEMENT> key_Bio; // Biot alpha
  csmp::INDEX<SCALAR,ELEMENT> key_BioT; // Biot term
  csmp::INDEX<SCALAR,ELEMENT> key_nu; // Poisson's ratio
  csmp::INDEX<SCALAR,ELEMENT> key_E; // Young's modulus
  csmp::INDEX<SCALAR,ELEMENT> key_rhob; // bulk density
  csmp::INDEX<SCALAR,ELEMENT> key_coh; // cohesion
  csmp::INDEX<SCALAR,ELEMENT> key_cR; // compressibility rock
  csmp::INDEX<SCALAR,ELEMENT> key_Vdef; // deformed volume
  csmp::INDEX<SCALAR,ELEMENT> key_dil; // dilatation
  csmp::INDEX<VECTOR,NODE> key_displ; // displacement
  csmp::INDEX<SCALAR,ELEMENT> key_rhod; // dry rock density
  csmp::INDEX<SCALAR,ELEMENT> key_Se; // effective stress
  csmp::INDEX<SCALAR,ELEMENT> key_fail; // failure
  csmp::INDEX<VECTOR,NODE> key_F; // force
  csmp::INDEX<SCALAR,ELEMENT> key_frica; // friction angle
  csmp::INDEX<SCALAR,ELEMENT> key_SSmax; // max shear stress
  csmp::INDEX<SCALAR,ELEMENT> key_Smean; // mean stress
  csmp::INDEX<SCALAR,FACE> key_Sn; // normal stress
  csmp::INDEX<SCALAR,NODE> key_op; // overburden pressure
  csmp::INDEX<SCALAR,MODEL> key_Sv; // overburden stress
  csmp::INDEX<SCALAR,FACE> key_Ss; // shear stress
  csmp::INDEX<VECTOR,ELEMENT> key_sig1; // sigma1
  csmp::INDEX<VECTOR,ELEMENT> key_sig2; // sigma2
  csmp::INDEX<VECTOR,ELEMENT> key_sig3; // sigma3
  csmp::INDEX<TENSOR,ELEMENT> key_e; // strain
  csmp::INDEX<VECTOR,ELEMENT> key_e1; // strain1
  csmp::INDEX<VECTOR,ELEMENT> key_e2; // strain2
  csmp::INDEX<VECTOR,ELEMENT> key_e3; // strain3
  csmp::INDEX<TENSOR,ELEMENT> key_S; // stress
  csmp::INDEX<SCALAR,ELEMENT> key_Sy; // stress-y
  csmp::INDEX<SCALAR,ELEMENT> key_fail1; // tensile failure
  csmp::INDEX<SCALAR,ELEMENT> key_dV; // volume change
  // multiphase flow
  csmp::INDEX<SCALAR,ELEMENT> key_kV; // vertical permeability
  csmp::INDEX<SCALAR,ELEMENT> key_diffpc; // capillary diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_K; // hydraulic conductivity
  csmp::INDEX<SCALAR,ELEMENT> key_pd; // entry pressure
  csmp::INDEX<VECTOR,FACE> key_fgt; // face gravity term
  csmp::INDEX<SCALAR,FACE> key_flt; // face total mobility
  csmp::INDEX<SCALAR,NODE> key_pf; // fluid pressure
  csmp::INDEX<SCALAR,NODE> key_mum; // mixture viscosity
  csmp::INDEX<SCALAR,NODE> key_nQV; // nodal fluid volume source
  csmp::INDEX<SCALAR,ELEMENT> key_QV; // fluid volume source
  csmp::INDEX<SCALAR,NODE> key_fb; // flux balance
  csmp::INDEX<VECTOR,ELEMENT> key_gt; // gravity term
  csmp::INDEX<SCALAR,ELEMENT> key_siH2O; // initial saturation aqueous phase
  csmp::INDEX<SCALAR,ELEMENT> key_laCO2; // mobility carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_laH2O; // mobility aqueous phase
  csmp::INDEX<SCALAR,NODE> key_sCO2_1; // new saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_srCO2; // residual saturation carbonic phase
  csmp::INDEX<SCALAR,ELEMENT> key_srH2O; // residual saturation aqueous phase
  csmp::INDEX<SCALAR,NODE> key_sCO2; // saturation carbonic phase
  csmp::INDEX<SCALAR,NODE> key_sH2O; // saturation aqueous phase
  csmp::INDEX<SCALAR,NODE> key_TDS; // total dissolved solids
  csmp::INDEX<SCALAR,ELEMENT> key_lat; // total mobility
  csmp::INDEX<SCALAR,ELEMENT> key_ct; // total system compressibility
  csmp::INDEX<VECTOR,ELEMENT> key_vt; // total velocity
  csmp::INDEX<SCALAR,ELEMENT> key_kh; // horizontal permeability
  csmp::INDEX<SCALAR,ELEMENT> key_qf; // volume flux
  csmp::INDEX<SCALAR,MODEL> key_t; // model time
  csmp::INDEX<SCALAR,ELEMENT> key_bcp; // brooks corey parameter
  csmp::INDEX<VECTOR,ELEMENT> key_Fg; // gravity force
  // numeric modelling
  csmp::INDEX<VECTOR,ELEMENT> key_dip; // dip vector
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_sPV; // sector pore volume
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_kfn; // facet normal permeability
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fAk; // facet area permeability
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA; // facet area
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; // facet normal
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qft; // facet total flux
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qfCO2; // facet flux carbonic phase
  csmp::INDEX<SCALAR,NODE> key_fvPV; // FV pore volume
  csmp::INDEX<SCALAR,NODE> key_cfl; // cfl
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
  // physical constants
  csmp::INDEX<SCALAR,MODEL> key_g; // acceleration gravity
  // property modelling
  csmp::INDEX<SCALAR,ELEMENT> key_kFbr; // fault breccia permeability
  csmp::INDEX<SCALAR,REGION> key_Fs; // fault size
  csmp::INDEX<SCALAR,ELEMENT> key_kf; // fracture permeability
  csmp::INDEX<TENSOR,ELEMENT> key_k; // permeability
  csmp::INDEX<SCALAR,ELEMENT> key_k_cal; // permeability calibration factor
  csmp::INDEX<SCALAR,ELEMENT> key_phi; // porosity
  csmp::INDEX<SCALAR,ELEMENT> key_RRT; // rocktype
  csmp::INDEX<SCALAR,REGION> key_AS; // slip patch area

  template<uint32_t dim>
  explicit VariableSet_TwoPhaseFlow( const PropertyDatabase<dim>& db )
    : key_XCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction CO2 aqueous phase") ))
    , key_XH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction H2O aqueous phase") ))
    , key_YCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction CO2 carbonic phase") ))
    , key_YH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction H2O carbonic phase") ))
    , key_xCO2( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction CO2 aqueous phase") ))
    , key_xH2O( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction H2O aqueous phase") ))
    , key_xSalt( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction NaCl aqueous phase") ))
    , key_yCO2( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction CO2 carbonic phase") ))
    , key_yH2O( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction H2O carbonic phase") ))
    , key_CO2aq( INDEX<SCALAR,NODE>( db.StorageKey("dissolved CO2") ))
    , key_H2Og( INDEX<SCALAR,NODE>( db.StorageKey("evaporated water") ))
    , key_msalt( INDEX<SCALAR,NODE>( db.StorageKey("molality salt") ))
    , key_cCO2( INDEX<SCALAR,NODE>( db.StorageKey("compressibility carbonic phase") ))
    , key_cH2O( INDEX<SCALAR,NODE>( db.StorageKey("compressibility aqueous phase") ))
    , key_rhoCO2( INDEX<SCALAR,NODE>( db.StorageKey("density carbonic phase") ))
    , key_rhoH2O( INDEX<SCALAR,NODE>( db.StorageKey("density aqueous phase") ))
    , key_diffCO2( INDEX<SCALAR,NODE>( db.StorageKey("diffusivity coefficient carbonic phase") ))
    , key_diffH2O( INDEX<SCALAR,NODE>( db.StorageKey("diffusivity coefficient aqueous phase") ))
    , key_rhom( INDEX<SCALAR,NODE>( db.StorageKey("density mixture") ))
    , key_IFT( INDEX<SCALAR,MODEL>( db.StorageKey("interfacial tension") ))
    , key_muCO2( INDEX<SCALAR,NODE>( db.StorageKey("viscosity carbonic phase") ))
    , key_muH2O( INDEX<SCALAR,NODE>( db.StorageKey("viscosity aqueous phase") ))
    , key_eC( INDEX<SCALAR,ELEMENT>( db.StorageKey("electric conductivity") ))
    , key_KH2O( INDEX<SCALAR,MODEL>( db.StorageKey("thermal conductivity aqueous phase") ))
    , key_Qh( INDEX<SCALAR,ELEMENT>( db.StorageKey("heat source") ))
    , key_KCO2( INDEX<SCALAR,MODEL>( db.StorageKey("thermal conductivity carbonic phase") ))
    , key_Kr( INDEX<SCALAR,ELEMENT>( db.StorageKey("thermal conductivity") ))
    , key_Qe( INDEX<SCALAR,ELEMENT>( db.StorageKey("energy source") ))
    , key_T( INDEX<SCALAR,NODE>( db.StorageKey("temperature") ))
    , key_Bio( INDEX<SCALAR,ELEMENT>( db.StorageKey("Biot alpha") ))
    , key_BioT( INDEX<SCALAR,ELEMENT>( db.StorageKey("Biot term") ))
    , key_nu( INDEX<SCALAR,ELEMENT>( db.StorageKey("Poisson's ratio") ))
    , key_E( INDEX<SCALAR,ELEMENT>( db.StorageKey("Young's modulus") ))
    , key_rhob( INDEX<SCALAR,ELEMENT>( db.StorageKey("bulk density") ))
    , key_coh( INDEX<SCALAR,ELEMENT>( db.StorageKey("cohesion") ))
    , key_cR( INDEX<SCALAR,ELEMENT>( db.StorageKey("compressibility rock") ))
    , key_Vdef( INDEX<SCALAR,ELEMENT>( db.StorageKey("deformed volume") ))
    , key_dil( INDEX<SCALAR,ELEMENT>( db.StorageKey("dilatation") ))
    , key_displ( INDEX<VECTOR,NODE>( db.StorageKey("displacement") ))
    , key_rhod( INDEX<SCALAR,ELEMENT>( db.StorageKey("dry rock density") ))
    , key_Se( INDEX<SCALAR,ELEMENT>( db.StorageKey("effective stress") ))
    , key_fail( INDEX<SCALAR,ELEMENT>( db.StorageKey("failure") ))
    , key_F( INDEX<VECTOR,NODE>( db.StorageKey("force") ))
    , key_frica( INDEX<SCALAR,ELEMENT>( db.StorageKey("friction angle") ))
    , key_SSmax( INDEX<SCALAR,ELEMENT>( db.StorageKey("max shear stress") ))
    , key_Smean( INDEX<SCALAR,ELEMENT>( db.StorageKey("mean stress") ))
    , key_Sn( INDEX<SCALAR,FACE>( db.StorageKey("normal stress") ))
    , key_op( INDEX<SCALAR,NODE>( db.StorageKey("overburden pressure") ))
    , key_Sv( INDEX<SCALAR,MODEL>( db.StorageKey("overburden stress") ))
    , key_Ss( INDEX<SCALAR,FACE>( db.StorageKey("shear stress") ))
    , key_sig1( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma1") ))
    , key_sig2( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma2") ))
    , key_sig3( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma3") ))
    , key_e( INDEX<TENSOR,ELEMENT>( db.StorageKey("strain") ))
    , key_e1( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain1") ))
    , key_e2( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain2") ))
    , key_e3( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain3") ))
    , key_S( INDEX<TENSOR,ELEMENT>( db.StorageKey("stress") ))
    , key_Sy( INDEX<SCALAR,ELEMENT>( db.StorageKey("stress-y") ))
    , key_fail1( INDEX<SCALAR,ELEMENT>( db.StorageKey("tensile failure") ))
    , key_dV( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume change") ))
    , key_kV( INDEX<SCALAR,ELEMENT>( db.StorageKey("vertical permeability") ))
    , key_diffpc( INDEX<SCALAR,ELEMENT>( db.StorageKey("capillary diffusivity") ))
    , key_K( INDEX<SCALAR,ELEMENT>( db.StorageKey("hydraulic conductivity") ))
    , key_pd( INDEX<SCALAR,ELEMENT>( db.StorageKey("entry pressure") ))
    , key_fgt( INDEX<VECTOR,FACE>( db.StorageKey("face gravity term") ))
    , key_flt( INDEX<SCALAR,FACE>( db.StorageKey("face total mobility") ))
    , key_pf( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_mum( INDEX<SCALAR,NODE>( db.StorageKey("mixture viscosity") ))
    , key_nQV( INDEX<SCALAR,NODE>( db.StorageKey("nodal fluid volume source") ))
    , key_QV( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid volume source") ))
    , key_fb( INDEX<SCALAR,NODE>( db.StorageKey("flux balance") ))
    , key_gt( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity term") ))
    , key_siH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("initial saturation aqueous phase") ))
    , key_laCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility carbonic phase") ))
    , key_laH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility aqueous phase") ))
    , key_sCO2_1( INDEX<SCALAR,NODE>( db.StorageKey("new saturation carbonic phase") ))
    , key_srCO2( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation carbonic phase") ))
    , key_srH2O( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation aqueous phase") ))
    , key_sCO2( INDEX<SCALAR,NODE>( db.StorageKey("saturation carbonic phase") ))
    , key_sH2O( INDEX<SCALAR,NODE>( db.StorageKey("saturation aqueous phase") ))
    , key_TDS( INDEX<SCALAR,NODE>( db.StorageKey("total dissolved solids") ))
    , key_lat( INDEX<SCALAR,ELEMENT>( db.StorageKey("total mobility") ))
    , key_ct( INDEX<SCALAR,ELEMENT>( db.StorageKey("total system compressibility") ))
    , key_vt( INDEX<VECTOR,ELEMENT>( db.StorageKey("total velocity") ))
    , key_kh( INDEX<SCALAR,ELEMENT>( db.StorageKey("horizontal permeability") ))
    , key_qf( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume flux") ))
    , key_t( INDEX<SCALAR,MODEL>( db.StorageKey("model time") ))
    , key_bcp( INDEX<SCALAR,ELEMENT>( db.StorageKey("brooks corey parameter") ))
    , key_Fg( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity force") ))
    , key_dip( INDEX<VECTOR,ELEMENT>( db.StorageKey("dip vector") ))
    , key_sPV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector pore volume") ))
    , key_kfn( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal permeability") ))
    , key_fAk( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area permeability") ))
    , key_fA( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area") ))
    , key_fn( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal") ))
    , key_qft( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet total flux") ))
    , key_qfCO2( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux carbonic phase") ))
    , key_fvPV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_cfl( INDEX<SCALAR,NODE>( db.StorageKey("cfl") ))
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
    , key_g( INDEX<SCALAR,MODEL>( db.StorageKey("acceleration gravity") ))
    , key_kFbr( INDEX<SCALAR,ELEMENT>( db.StorageKey("fault breccia permeability") ))
    , key_Fs( INDEX<SCALAR,REGION>( db.StorageKey("fault size") ))
    , key_kf( INDEX<SCALAR,ELEMENT>( db.StorageKey("fracture permeability") ))
    , key_k( INDEX<TENSOR,ELEMENT>( db.StorageKey("permeability") ))
    , key_k_cal( INDEX<SCALAR,ELEMENT>( db.StorageKey("permeability calibration factor") ))
    , key_phi( INDEX<SCALAR,ELEMENT>( db.StorageKey("porosity") ))
    , key_RRT( INDEX<SCALAR,ELEMENT>( db.StorageKey("rocktype") ))
    , key_AS( INDEX<SCALAR,REGION>( db.StorageKey("slip patch area") ))
  {
    if ( key_XCO2.place != NODE || key_XCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mass fraction CO2 aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_XH2O.place != NODE || key_XH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mass fraction H2O aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_YCO2.place != NODE || key_YCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mass fraction CO2 carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_YH2O.place != NODE || key_YH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mass fraction H2O carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xCO2.place != NODE || key_xCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mole fraction CO2 aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xH2O.place != NODE || key_xH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mole fraction H2O aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xSalt.place != NODE || key_xSalt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mole fraction NaCl aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_yCO2.place != NODE || key_yCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mole fraction CO2 carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_yH2O.place != NODE || key_yH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mole fraction H2O carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_CO2aq.place != NODE || key_CO2aq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'dissolved CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_H2Og.place != NODE || key_H2Og.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'evaporated water' variable must be SCALAR and placed on NODE"  );
    if ( key_msalt.place != NODE || key_msalt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'molality salt' variable must be SCALAR and placed on NODE"  );
    if ( key_cCO2.place != NODE || key_cCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'compressibility carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_cH2O.place != NODE || key_cH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'compressibility aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoCO2.place != NODE || key_rhoCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'density carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhoH2O.place != NODE || key_rhoH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'density aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_diffCO2.place != NODE || key_diffCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'diffusivity coefficient carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_diffH2O.place != NODE || key_diffH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'diffusivity coefficient aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhom.place != NODE || key_rhom.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'density mixture' variable must be SCALAR and placed on NODE"  );
    if ( key_IFT.place != MODEL || key_IFT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'interfacial tension' variable must be SCALAR and placed on MODEL"  );
    if ( key_muCO2.place != NODE || key_muCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'viscosity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_muH2O.place != NODE || key_muH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'viscosity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_eC.place != ELEMENT || key_eC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'electric conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KH2O.place != MODEL || key_KH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'thermal conductivity aqueous phase' variable must be SCALAR and placed on MODEL"  );
    if ( key_Qh.place != ELEMENT || key_Qh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'heat source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KCO2.place != MODEL || key_KCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'thermal conductivity carbonic phase' variable must be SCALAR and placed on MODEL"  );
    if ( key_Kr.place != ELEMENT || key_Kr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'thermal conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Qe.place != ELEMENT || key_Qe.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'energy source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_T.place != NODE || key_T.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'temperature' variable must be SCALAR and placed on NODE"  );
    if ( key_Bio.place != ELEMENT || key_Bio.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'Biot alpha' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_BioT.place != ELEMENT || key_BioT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'Biot term' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nu.place != ELEMENT || key_nu.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'Poisson's ratio' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_E.place != ELEMENT || key_E.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'Young's modulus' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_rhob.place != ELEMENT || key_rhob.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'bulk density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_coh.place != ELEMENT || key_coh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'cohesion' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_cR.place != ELEMENT || key_cR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'compressibility rock' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Vdef.place != ELEMENT || key_Vdef.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'deformed volume' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dil.place != ELEMENT || key_dil.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'dilatation' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_displ.place != NODE || key_displ.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'displacement' variable must be VECTOR and placed on NODE"  );
    if ( key_rhod.place != ELEMENT || key_rhod.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'dry rock density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Se.place != ELEMENT || key_Se.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'effective stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fail.place != ELEMENT || key_fail.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'failure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_F.place != NODE || key_F.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'force' variable must be VECTOR and placed on NODE"  );
    if ( key_frica.place != ELEMENT || key_frica.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'friction angle' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_SSmax.place != ELEMENT || key_SSmax.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'max shear stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Smean.place != ELEMENT || key_Smean.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mean stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Sn.place != FACE || key_Sn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'normal stress' variable must be SCALAR and placed on FACE"  );
    if ( key_op.place != NODE || key_op.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'overburden pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_Sv.place != MODEL || key_Sv.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'overburden stress' variable must be SCALAR and placed on MODEL"  );
    if ( key_Ss.place != FACE || key_Ss.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'shear stress' variable must be SCALAR and placed on FACE"  );
    if ( key_sig1.place != ELEMENT || key_sig1.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'sigma1' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sig2.place != ELEMENT || key_sig2.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'sigma2' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sig3.place != ELEMENT || key_sig3.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'sigma3' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_e.place != ELEMENT || key_e.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'strain' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_e1.place != ELEMENT || key_e1.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'strain1' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_e2.place != ELEMENT || key_e2.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'strain2' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_e3.place != ELEMENT || key_e3.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'strain3' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_S.place != ELEMENT || key_S.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'stress' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_Sy.place != ELEMENT || key_Sy.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'stress-y' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fail1.place != ELEMENT || key_fail1.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'tensile failure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dV.place != ELEMENT || key_dV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'volume change' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kV.place != ELEMENT || key_kV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'vertical permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diffpc.place != ELEMENT || key_diffpc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'capillary diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_K.place != ELEMENT || key_K.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'hydraulic conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_pd.place != ELEMENT || key_pd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'entry pressure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fgt.place != FACE || key_fgt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'face gravity term' variable must be VECTOR and placed on FACE"  );
    if ( key_flt.place != FACE || key_flt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'face total mobility' variable must be SCALAR and placed on FACE"  );
    if ( key_pf.place != NODE || key_pf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_mum.place != NODE || key_mum.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mixture viscosity' variable must be SCALAR and placed on NODE"  );
    if ( key_nQV.place != NODE || key_nQV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( key_QV.place != ELEMENT || key_QV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fb.place != NODE || key_fb.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( key_gt.place != ELEMENT || key_gt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'gravity term' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_siH2O.place != ELEMENT || key_siH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'initial saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_laCO2.place != ELEMENT || key_laCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mobility carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_laH2O.place != ELEMENT || key_laH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'mobility aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sCO2_1.place != NODE || key_sCO2_1.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'new saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_srCO2.place != ELEMENT || key_srCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'residual saturation carbonic phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_srH2O.place != ELEMENT || key_srH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'residual saturation aqueous phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sCO2.place != NODE || key_sCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'saturation carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_sH2O.place != NODE || key_sH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'saturation aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_TDS.place != NODE || key_TDS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'total dissolved solids' variable must be SCALAR and placed on NODE"  );
    if ( key_lat.place != ELEMENT || key_lat.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'total mobility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_ct.place != ELEMENT || key_ct.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'total system compressibility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_vt.place != ELEMENT || key_vt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'total velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_kh.place != ELEMENT || key_kh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'horizontal permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_qf.place != ELEMENT || key_qf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'volume flux' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_t.place != MODEL || key_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'model time' variable must be SCALAR and placed on MODEL"  );
    if ( key_bcp.place != ELEMENT || key_bcp.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'brooks corey parameter' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Fg.place != ELEMENT || key_Fg.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'gravity force' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_dip.place != ELEMENT || key_dip.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'dip vector' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_sPV.place != SECTOR_INTEGRATION_POINT || key_sPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'sector pore volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( key_kfn.place != FACET_INTEGRATION_POINT || key_kfn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet normal permeability' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fAk.place != FACET_INTEGRATION_POINT || key_fAk.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet area permeability' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qft.place != FACET_INTEGRATION_POINT || key_qft.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet total flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qfCO2.place != FACET_INTEGRATION_POINT || key_qfCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'facet flux carbonic phase' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fvPV.place != NODE || key_fvPV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( key_cfl.place != NODE || key_cfl.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'cfl' variable must be SCALAR and placed on NODE"  );
    if ( key_dt.place != MODEL || key_dt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'default time increment' variable must be SCALAR and placed on MODEL"  );
    if ( key_Ne.place != ELEMENT || key_Ne.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'element number' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dipf.place != FACE || key_dipf.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'face dip vector' variable must be VECTOR and placed on FACE"  );
    if ( key_Nn.place != NODE || key_Nn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'node number' variable must be SCALAR and placed on NODE"  );
    if ( key_t_out.place != MODEL || key_t_out.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'output times' variable must be ARRAY and placed on MODEL"  );
    if ( key_pf0.place != NODE || key_pf0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'previous fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_QV0.place != ELEMENT || key_QV0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'previous fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nQV0.place != NODE || key_nQV0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'previous nodal fluid volume source' variable must be SCALAR and placed on NODE"  );
    if ( key_TMAX.place != MODEL || key_TMAX.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'run duration' variable must be SCALAR and placed on MODEL"  );
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_V.place != ELEMENT || key_V.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'volume' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_g.place != MODEL || key_g.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'acceleration gravity' variable must be SCALAR and placed on MODEL"  );
    if ( key_kFbr.place != ELEMENT || key_kFbr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'fault breccia permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Fs.place != REGION || key_Fs.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'fault size' variable must be SCALAR and placed on REGION"  );
    if ( key_kf.place != ELEMENT || key_kf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'fracture permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k.place != ELEMENT || key_k.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'permeability' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_k_cal.place != ELEMENT || key_k_cal.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'permeability calibration factor' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_RRT.place != ELEMENT || key_RRT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'rocktype' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_AS.place != REGION || key_AS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "VariableSet_TwoPhaseFlow::VariableSet_TwoPhaseFlow:",
        "The 'slip patch area' variable must be SCALAR and placed on REGION"  );
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLESET_TWOPHASEFLOW_H
