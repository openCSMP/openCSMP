#ifndef CSMP_VARIABLES_TWOPHASEFLOW_H
#define CSMP_VARIABLES_TWOPHASEFLOW_H

/**
@file Variables_TwoPhaseFlow.h
Automatically generated from /Users/andrew/Development/repository/csmp-api-library/source/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_TwoPhaseFlow.csv
DO NOT EDIT!
*/

#include "Index.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp { namespace variables {

struct Variables_TwoPhaseFlow {
  // chemistry
  csmp::INDEX<SCALAR,NODE> key_XCO2; // mass fraction CO2 aqueous phase
  csmp::INDEX<SCALAR,NODE> key_XH2O; // mass fraction H2O aqueous phase
  csmp::INDEX<SCALAR,NODE> key_YCO2; // mass fraction CO2 carbonic phase
  csmp::INDEX<SCALAR,NODE> key_YH2O; // mass fraction H2O carbonic phase
  csmp::INDEX<SCALAR,NODE> key_xCO2; // mole fraction CO2 aqueous phase
  csmp::INDEX<SCALAR,NODE> key_xH2O; // mole fraction H2O aqueous phase
  csmp::INDEX<SCALAR,NODE> key_xSalt; // mole fraction NaCl aqueous phase
  csmp::INDEX<SCALAR,NODE> key_yCO2; // mole fraction CO2 carbonic phase
  csmp::INDEX<SCALAR,NODE> key_yH2O; // mole fraction H2O aqueous phase
  csmp::INDEX<SCALAR,NODE> key_CO2aq; // dissolved CO2
  csmp::INDEX<SCALAR,NODE> key_H2Og; // evaporated water
  csmp::INDEX<SCALAR,NODE> key_msalt; // molality salt
  // fluid properties
  csmp::INDEX<SCALAR,NODE> key_cnw; // compressibility carbonic phase
  csmp::INDEX<SCALAR,NODE> key_cw; // compressibility aqueous phase
  csmp::INDEX<SCALAR,NODE> key_rhonw; // density carbonic phase
  csmp::INDEX<SCALAR,NODE> key_rhow; // density aqueous phase
  csmp::INDEX<SCALAR,NODE> key_diffnw; // diffusivity coefficient carbonic phase
  csmp::INDEX<SCALAR,NODE> key_diffw; // diffusivity coefficient aqueous phase
  csmp::INDEX<SCALAR,NODE> key_rhom; // density mixture
  csmp::INDEX<SCALAR,MODEL> key_IFT; // interfacial tension
  csmp::INDEX<SCALAR,NODE> key_munw; // viscosity carbonic phase
  csmp::INDEX<SCALAR,NODE> key_muw; // viscosity aqueous phase
  // geophysics
  csmp::INDEX<SCALAR,ELEMENT> key_eC; // electric conductivity
  // heat flow
  csmp::INDEX<SCALAR,MODEL> key_KH2O; // thermal conductivity water
  csmp::INDEX<SCALAR,ELEMENT> key_Qh; // heat source
  // heat transfer
  csmp::INDEX<SCALAR,MODEL> key_KCO2; // thermal conductivity CO2
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
  csmp::INDEX<SCALAR,ELEMENT> key_VD; // deformed volume
  csmp::INDEX<SCALAR,ELEMENT> key_dil; // dilatation
  csmp::INDEX<VECTOR,NODE> key_displ; // displacement
  csmp::INDEX<SCALAR,ELEMENT> key_rhod; // dry rock density
  csmp::INDEX<SCALAR,ELEMENT> key_SSE; // effective stress
  csmp::INDEX<SCALAR,ELEMENT> key_F; // failure
  csmp::INDEX<VECTOR,NODE> key_f; // force
  csmp::INDEX<SCALAR,ELEMENT> key_frica; // friction angle
  csmp::INDEX<SCALAR,ELEMENT> key_SSmax; // max shear stress
  csmp::INDEX<SCALAR,ELEMENT> key_Smean; // mean stress
  csmp::INDEX<SCALAR,FACE> key_Sn; // normal stress
  csmp::INDEX<SCALAR,NODE> key_op; // overburden pressure
  csmp::INDEX<SCALAR,MODEL> key_ll; // overburden stress
  csmp::INDEX<SCALAR,FACE> key_SS; // shear stress
  csmp::INDEX<VECTOR,ELEMENT> key_SI1; // sigma1
  csmp::INDEX<VECTOR,ELEMENT> key_SI2; // sigma2
  csmp::INDEX<VECTOR,ELEMENT> key_SI3; // sigma3
  csmp::INDEX<TENSOR,ELEMENT> key_EPS; // strain
  csmp::INDEX<VECTOR,ELEMENT> key_e1; // strain1
  csmp::INDEX<VECTOR,ELEMENT> key_e2; // strain2
  csmp::INDEX<VECTOR,ELEMENT> key_e3; // strain3
  csmp::INDEX<TENSOR,ELEMENT> key_S; // stress
  csmp::INDEX<SCALAR,ELEMENT> key_Sy; // stress-y
  csmp::INDEX<SCALAR,ELEMENT> key_F1; // tensile failure
  csmp::INDEX<SCALAR,ELEMENT> key_dV; // volume change
  // multiphase flow
  csmp::INDEX<SCALAR,ELEMENT> key_kV; // vertical permeability
  csmp::INDEX<SCALAR,ELEMENT> key_diffpc; // capillary diffusivity
  csmp::INDEX<SCALAR,ELEMENT> key_K; // hydraulic conductivity
  csmp::INDEX<SCALAR,ELEMENT> key_pd; // entry pressure
  csmp::INDEX<VECTOR,FACE> key_fgt; // face gravity term
  csmp::INDEX<SCALAR,FACE> key_flt; // face total mobility
  csmp::INDEX<SCALAR,NODE> key_pf; // fluid pressure
  csmp::INDEX<SCALAR,NODE> key_muM; // mixture viscosity
  csmp::INDEX<SCALAR,ELEMENT> key_QV; // fluid volume source
  csmp::INDEX<SCALAR,NODE> key_fb; // flux balance
  csmp::INDEX<VECTOR,ELEMENT> key_gt; // gravity term
  csmp::INDEX<SCALAR,ELEMENT> key_SWI; // initial water saturation
  csmp::INDEX<SCALAR,ELEMENT> key_lnw; // mobility CO2
  csmp::INDEX<SCALAR,ELEMENT> key_lw; // mobility water
  csmp::INDEX<SCALAR,NODE> key_sCO2_1; // new saturation CO2
  csmp::INDEX<SCALAR,ELEMENT> key_SNR; // residual saturation nonwetting phase
  csmp::INDEX<SCALAR,ELEMENT> key_SWR; // residual saturation wetting phase
  csmp::INDEX<SCALAR,NODE> key_snw; // saturation nonwetting phase
  csmp::INDEX<SCALAR,NODE> key_sw; // saturation wetting phase
  csmp::INDEX<SCALAR,NODE> key_TDS; // total dissolved solids
  csmp::INDEX<SCALAR,ELEMENT> key_lt; // total mobility
  csmp::INDEX<SCALAR,ELEMENT> key_ct; // total system compressibility
  csmp::INDEX<VECTOR,ELEMENT> key_vt; // total velocity
  csmp::INDEX<SCALAR,ELEMENT> key_kH; // horizontal permeability
  csmp::INDEX<SCALAR,ELEMENT> key_qf; // volume flux
  csmp::INDEX<SCALAR,MODEL> key_t; // model time
  csmp::INDEX<SCALAR,ELEMENT> key_bcp; // brooks corey parameter
  // multiphase fow
  csmp::INDEX<VECTOR,ELEMENT> key_gf; // gravity force
  // numeric modelling
  csmp::INDEX<VECTOR,ELEMENT> key_vd; // dip vector
  csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> key_poreSFV; // sector pore volume
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_kfn; // facet normal permeability
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_fA; // facet area
  csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> key_fn; // facet normal
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qft; // facet total flux
  csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> key_qfCO2; // facet flux CO2
  csmp::INDEX<SCALAR,NODE> key_poreFV; // FV pore volume
  csmp::INDEX<SCALAR,NODE> key_cfl; // cfl
  csmp::INDEX<SCALAR,MODEL> key_DT; // default time increment
  csmp::INDEX<SCALAR,ELEMENT> key_Ne; // element number
  csmp::INDEX<VECTOR,FACE> key_fdvec; // face dip vector
  csmp::INDEX<SCALAR,NODE> key_Nn; // node number
  csmp::INDEX<ARRAY,MODEL> key_t_out; // output times
  csmp::INDEX<SCALAR,NODE> key_pf0; // previous fluid pressure
  csmp::INDEX<SCALAR,ELEMENT> key_QV0; // previous fluid volume source
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
  csmp::INDEX<SCALAR,REGION> key_RRT; // rocktype
  csmp::INDEX<SCALAR,REGION> key_AS; // slip patch area

  template<size_t dim>
  explicit Variables_TwoPhaseFlow( const PropertyDatabase<dim>& db )
    : key_XCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction CO2 aqueous phase") ))
    , key_XH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction H2O aqueous phase") ))
    , key_YCO2( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction CO2 carbonic phase") ))
    , key_YH2O( INDEX<SCALAR,NODE>( db.StorageKey("mass fraction H2O carbonic phase") ))
    , key_xCO2( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction CO2 aqueous phase") ))
    , key_xH2O( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction H2O aqueous phase") ))
    , key_xSalt( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction NaCl aqueous phase") ))
    , key_yCO2( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction CO2 carbonic phase") ))
    , key_yH2O( INDEX<SCALAR,NODE>( db.StorageKey("mole fraction H2O aqueous phase") ))
    , key_CO2aq( INDEX<SCALAR,NODE>( db.StorageKey("dissolved CO2") ))
    , key_H2Og( INDEX<SCALAR,NODE>( db.StorageKey("evaporated water") ))
    , key_msalt( INDEX<SCALAR,NODE>( db.StorageKey("molality salt") ))
    , key_cnw( INDEX<SCALAR,NODE>( db.StorageKey("compressibility carbonic phase") ))
    , key_cw( INDEX<SCALAR,NODE>( db.StorageKey("compressibility aqueous phase") ))
    , key_rhonw( INDEX<SCALAR,NODE>( db.StorageKey("density carbonic phase") ))
    , key_rhow( INDEX<SCALAR,NODE>( db.StorageKey("density aqueous phase") ))
    , key_diffnw( INDEX<SCALAR,NODE>( db.StorageKey("diffusivity coefficient carbonic phase") ))
    , key_diffw( INDEX<SCALAR,NODE>( db.StorageKey("diffusivity coefficient aqueous phase") ))
    , key_rhom( INDEX<SCALAR,NODE>( db.StorageKey("density mixture") ))
    , key_IFT( INDEX<SCALAR,MODEL>( db.StorageKey("interfacial tension") ))
    , key_munw( INDEX<SCALAR,NODE>( db.StorageKey("viscosity carbonic phase") ))
    , key_muw( INDEX<SCALAR,NODE>( db.StorageKey("viscosity aqueous phase") ))
    , key_eC( INDEX<SCALAR,ELEMENT>( db.StorageKey("electric conductivity") ))
    , key_KH2O( INDEX<SCALAR,MODEL>( db.StorageKey("thermal conductivity water") ))
    , key_Qh( INDEX<SCALAR,ELEMENT>( db.StorageKey("heat source") ))
    , key_KCO2( INDEX<SCALAR,MODEL>( db.StorageKey("thermal conductivity CO2") ))
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
    , key_VD( INDEX<SCALAR,ELEMENT>( db.StorageKey("deformed volume") ))
    , key_dil( INDEX<SCALAR,ELEMENT>( db.StorageKey("dilatation") ))
    , key_displ( INDEX<VECTOR,NODE>( db.StorageKey("displacement") ))
    , key_rhod( INDEX<SCALAR,ELEMENT>( db.StorageKey("dry rock density") ))
    , key_SSE( INDEX<SCALAR,ELEMENT>( db.StorageKey("effective stress") ))
    , key_F( INDEX<SCALAR,ELEMENT>( db.StorageKey("failure") ))
    , key_f( INDEX<VECTOR,NODE>( db.StorageKey("force") ))
    , key_frica( INDEX<SCALAR,ELEMENT>( db.StorageKey("friction angle") ))
    , key_SSmax( INDEX<SCALAR,ELEMENT>( db.StorageKey("max shear stress") ))
    , key_Smean( INDEX<SCALAR,ELEMENT>( db.StorageKey("mean stress") ))
    , key_Sn( INDEX<SCALAR,FACE>( db.StorageKey("normal stress") ))
    , key_op( INDEX<SCALAR,NODE>( db.StorageKey("overburden pressure") ))
    , key_ll( INDEX<SCALAR,MODEL>( db.StorageKey("overburden stress") ))
    , key_SS( INDEX<SCALAR,FACE>( db.StorageKey("shear stress") ))
    , key_SI1( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma1") ))
    , key_SI2( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma2") ))
    , key_SI3( INDEX<VECTOR,ELEMENT>( db.StorageKey("sigma3") ))
    , key_EPS( INDEX<TENSOR,ELEMENT>( db.StorageKey("strain") ))
    , key_e1( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain1") ))
    , key_e2( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain2") ))
    , key_e3( INDEX<VECTOR,ELEMENT>( db.StorageKey("strain3") ))
    , key_S( INDEX<TENSOR,ELEMENT>( db.StorageKey("stress") ))
    , key_Sy( INDEX<SCALAR,ELEMENT>( db.StorageKey("stress-y") ))
    , key_F1( INDEX<SCALAR,ELEMENT>( db.StorageKey("tensile failure") ))
    , key_dV( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume change") ))
    , key_kV( INDEX<SCALAR,ELEMENT>( db.StorageKey("vertical permeability") ))
    , key_diffpc( INDEX<SCALAR,ELEMENT>( db.StorageKey("capillary diffusivity") ))
    , key_K( INDEX<SCALAR,ELEMENT>( db.StorageKey("hydraulic conductivity") ))
    , key_pd( INDEX<SCALAR,ELEMENT>( db.StorageKey("entry pressure") ))
    , key_fgt( INDEX<VECTOR,FACE>( db.StorageKey("face gravity term") ))
    , key_flt( INDEX<SCALAR,FACE>( db.StorageKey("face total mobility") ))
    , key_pf( INDEX<SCALAR,NODE>( db.StorageKey("fluid pressure") ))
    , key_muM( INDEX<SCALAR,NODE>( db.StorageKey("mixture viscosity") ))
    , key_QV( INDEX<SCALAR,ELEMENT>( db.StorageKey("fluid volume source") ))
    , key_fb( INDEX<SCALAR,NODE>( db.StorageKey("flux balance") ))
    , key_gt( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity term") ))
    , key_SWI( INDEX<SCALAR,ELEMENT>( db.StorageKey("initial water saturation") ))
    , key_lnw( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility CO2") ))
    , key_lw( INDEX<SCALAR,ELEMENT>( db.StorageKey("mobility water") ))
    , key_sCO2_1( INDEX<SCALAR,NODE>( db.StorageKey("new saturation CO2") ))
    , key_SNR( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation nonwetting phase") ))
    , key_SWR( INDEX<SCALAR,ELEMENT>( db.StorageKey("residual saturation wetting phase") ))
    , key_snw( INDEX<SCALAR,NODE>( db.StorageKey("saturation nonwetting phase") ))
    , key_sw( INDEX<SCALAR,NODE>( db.StorageKey("saturation wetting phase") ))
    , key_TDS( INDEX<SCALAR,NODE>( db.StorageKey("total dissolved solids") ))
    , key_lt( INDEX<SCALAR,ELEMENT>( db.StorageKey("total mobility") ))
    , key_ct( INDEX<SCALAR,ELEMENT>( db.StorageKey("total system compressibility") ))
    , key_vt( INDEX<VECTOR,ELEMENT>( db.StorageKey("total velocity") ))
    , key_kH( INDEX<SCALAR,ELEMENT>( db.StorageKey("horizontal permeability") ))
    , key_qf( INDEX<SCALAR,ELEMENT>( db.StorageKey("volume flux") ))
    , key_t( INDEX<SCALAR,MODEL>( db.StorageKey("model time") ))
    , key_bcp( INDEX<SCALAR,ELEMENT>( db.StorageKey("brooks corey parameter") ))
    , key_gf( INDEX<VECTOR,ELEMENT>( db.StorageKey("gravity force") ))
    , key_vd( INDEX<VECTOR,ELEMENT>( db.StorageKey("dip vector") ))
    , key_poreSFV( INDEX<SCALAR,SECTOR_INTEGRATION_POINT>( db.StorageKey("sector pore volume") ))
    , key_kfn( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal permeability") ))
    , key_fA( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet area") ))
    , key_fn( INDEX<VECTOR,FACET_INTEGRATION_POINT>( db.StorageKey("facet normal") ))
    , key_qft( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet total flux") ))
    , key_qfCO2( INDEX<SCALAR,FACET_INTEGRATION_POINT>( db.StorageKey("facet flux CO2") ))
    , key_poreFV( INDEX<SCALAR,NODE>( db.StorageKey("FV pore volume") ))
    , key_cfl( INDEX<SCALAR,NODE>( db.StorageKey("cfl") ))
    , key_DT( INDEX<SCALAR,MODEL>( db.StorageKey("default time increment") ))
    , key_Ne( INDEX<SCALAR,ELEMENT>( db.StorageKey("element number") ))
    , key_fdvec( INDEX<VECTOR,FACE>( db.StorageKey("face dip vector") ))
    , key_Nn( INDEX<SCALAR,NODE>( db.StorageKey("node number") ))
    , key_t_out( INDEX<ARRAY,MODEL>( db.StorageKey("output times") ))
    , key_pf0( INDEX<SCALAR,NODE>( db.StorageKey("previous fluid pressure") ))
    , key_QV0( INDEX<SCALAR,ELEMENT>( db.StorageKey("previous fluid volume source") ))
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
    , key_RRT( INDEX<SCALAR,REGION>( db.StorageKey("rocktype") ))
    , key_AS( INDEX<SCALAR,REGION>( db.StorageKey("slip patch area") ))
  {
    if ( key_XCO2.place != NODE || key_XCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mass fraction CO2 aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_XH2O.place != NODE || key_XH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mass fraction H2O aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_YCO2.place != NODE || key_YCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mass fraction CO2 carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_YH2O.place != NODE || key_YH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mass fraction H2O carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xCO2.place != NODE || key_xCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mole fraction CO2 aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xH2O.place != NODE || key_xH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mole fraction H2O aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_xSalt.place != NODE || key_xSalt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mole fraction NaCl aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_yCO2.place != NODE || key_yCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mole fraction CO2 carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_yH2O.place != NODE || key_yH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mole fraction H2O aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_CO2aq.place != NODE || key_CO2aq.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'dissolved CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_H2Og.place != NODE || key_H2Og.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'evaporated water' variable must be SCALAR and placed on NODE"  );
    if ( key_msalt.place != NODE || key_msalt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'molality salt' variable must be SCALAR and placed on NODE"  );
    if ( key_cnw.place != NODE || key_cnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'compressibility carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_cw.place != NODE || key_cw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'compressibility aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhonw.place != NODE || key_rhonw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'density carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhow.place != NODE || key_rhow.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'density aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_diffnw.place != NODE || key_diffnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'diffusivity coefficient carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_diffw.place != NODE || key_diffw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'diffusivity coefficient aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_rhom.place != NODE || key_rhom.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'density mixture' variable must be SCALAR and placed on NODE"  );
    if ( key_IFT.place != MODEL || key_IFT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'interfacial tension' variable must be SCALAR and placed on MODEL"  );
    if ( key_munw.place != NODE || key_munw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'viscosity carbonic phase' variable must be SCALAR and placed on NODE"  );
    if ( key_muw.place != NODE || key_muw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'viscosity aqueous phase' variable must be SCALAR and placed on NODE"  );
    if ( key_eC.place != ELEMENT || key_eC.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'electric conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KH2O.place != MODEL || key_KH2O.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'thermal conductivity water' variable must be SCALAR and placed on MODEL"  );
    if ( key_Qh.place != ELEMENT || key_Qh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'heat source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_KCO2.place != MODEL || key_KCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'thermal conductivity CO2' variable must be SCALAR and placed on MODEL"  );
    if ( key_Kr.place != ELEMENT || key_Kr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'thermal conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Qe.place != ELEMENT || key_Qe.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'energy source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_T.place != NODE || key_T.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'temperature' variable must be SCALAR and placed on NODE"  );
    if ( key_Bio.place != ELEMENT || key_Bio.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'Biot alpha' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_BioT.place != ELEMENT || key_BioT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'Biot term' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_nu.place != ELEMENT || key_nu.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'Poisson's ratio' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_E.place != ELEMENT || key_E.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'Young's modulus' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_rhob.place != ELEMENT || key_rhob.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'bulk density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_coh.place != ELEMENT || key_coh.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'cohesion' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_cR.place != ELEMENT || key_cR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'compressibility rock' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_VD.place != ELEMENT || key_VD.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'deformed volume' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dil.place != ELEMENT || key_dil.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'dilatation' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_displ.place != NODE || key_displ.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'displacement' variable must be VECTOR and placed on NODE"  );
    if ( key_rhod.place != ELEMENT || key_rhod.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'dry rock density' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_SSE.place != ELEMENT || key_SSE.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'effective stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_F.place != ELEMENT || key_F.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'failure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_f.place != NODE || key_f.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'force' variable must be VECTOR and placed on NODE"  );
    if ( key_frica.place != ELEMENT || key_frica.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'friction angle' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_SSmax.place != ELEMENT || key_SSmax.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'max shear stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Smean.place != ELEMENT || key_Smean.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mean stress' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Sn.place != FACE || key_Sn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'normal stress' variable must be SCALAR and placed on FACE"  );
    if ( key_op.place != NODE || key_op.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'overburden pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_ll.place != MODEL || key_ll.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'overburden stress' variable must be SCALAR and placed on MODEL"  );
    if ( key_SS.place != FACE || key_SS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'shear stress' variable must be SCALAR and placed on FACE"  );
    if ( key_SI1.place != ELEMENT || key_SI1.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'sigma1' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_SI2.place != ELEMENT || key_SI2.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'sigma2' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_SI3.place != ELEMENT || key_SI3.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'sigma3' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_EPS.place != ELEMENT || key_EPS.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'strain' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_e1.place != ELEMENT || key_e1.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'strain1' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_e2.place != ELEMENT || key_e2.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'strain2' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_e3.place != ELEMENT || key_e3.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'strain3' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_S.place != ELEMENT || key_S.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'stress' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_Sy.place != ELEMENT || key_Sy.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'stress-y' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_F1.place != ELEMENT || key_F1.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'tensile failure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_dV.place != ELEMENT || key_dV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'volume change' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_kV.place != ELEMENT || key_kV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'vertical permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_diffpc.place != ELEMENT || key_diffpc.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'capillary diffusivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_K.place != ELEMENT || key_K.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'hydraulic conductivity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_pd.place != ELEMENT || key_pd.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'entry pressure' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fgt.place != FACE || key_fgt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'face gravity term' variable must be VECTOR and placed on FACE"  );
    if ( key_flt.place != FACE || key_flt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'face total mobility' variable must be SCALAR and placed on FACE"  );
    if ( key_pf.place != NODE || key_pf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_muM.place != NODE || key_muM.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mixture viscosity' variable must be SCALAR and placed on NODE"  );
    if ( key_QV.place != ELEMENT || key_QV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fb.place != NODE || key_fb.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'flux balance' variable must be SCALAR and placed on NODE"  );
    if ( key_gt.place != ELEMENT || key_gt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'gravity term' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_SWI.place != ELEMENT || key_SWI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'initial water saturation' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lnw.place != ELEMENT || key_lnw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mobility CO2' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_lw.place != ELEMENT || key_lw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'mobility water' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_sCO2_1.place != NODE || key_sCO2_1.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'new saturation CO2' variable must be SCALAR and placed on NODE"  );
    if ( key_SNR.place != ELEMENT || key_SNR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'residual saturation nonwetting phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_SWR.place != ELEMENT || key_SWR.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'residual saturation wetting phase' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_snw.place != NODE || key_snw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'saturation nonwetting phase' variable must be SCALAR and placed on NODE"  );
    if ( key_sw.place != NODE || key_sw.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'saturation wetting phase' variable must be SCALAR and placed on NODE"  );
    if ( key_TDS.place != NODE || key_TDS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'total dissolved solids' variable must be SCALAR and placed on NODE"  );
    if ( key_lt.place != ELEMENT || key_lt.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'total mobility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_ct.place != ELEMENT || key_ct.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'total system compressibility' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_vt.place != ELEMENT || key_vt.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'total velocity' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_kH.place != ELEMENT || key_kH.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'horizontal permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_qf.place != ELEMENT || key_qf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'volume flux' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_t.place != MODEL || key_t.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'model time' variable must be SCALAR and placed on MODEL"  );
    if ( key_bcp.place != ELEMENT || key_bcp.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'brooks corey parameter' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_gf.place != ELEMENT || key_gf.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'gravity force' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_vd.place != ELEMENT || key_vd.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'dip vector' variable must be VECTOR and placed on ELEMENT"  );
    if ( key_poreSFV.place != SECTOR_INTEGRATION_POINT || key_poreSFV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'sector pore volume' variable must be SCALAR and placed on SECTOR_INTEGRATION_POINT"  );
    if ( key_kfn.place != FACET_INTEGRATION_POINT || key_kfn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'facet normal permeability' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fA.place != FACET_INTEGRATION_POINT || key_fA.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'facet area' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_fn.place != FACET_INTEGRATION_POINT || key_fn.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'facet normal' variable must be VECTOR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qft.place != FACET_INTEGRATION_POINT || key_qft.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'facet total flux' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_qfCO2.place != FACET_INTEGRATION_POINT || key_qfCO2.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'facet flux CO2' variable must be SCALAR and placed on FACET_INTEGRATION_POINT"  );
    if ( key_poreFV.place != NODE || key_poreFV.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'FV pore volume' variable must be SCALAR and placed on NODE"  );
    if ( key_cfl.place != NODE || key_cfl.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'cfl' variable must be SCALAR and placed on NODE"  );
    if ( key_DT.place != MODEL || key_DT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'default time increment' variable must be SCALAR and placed on MODEL"  );
    if ( key_Ne.place != ELEMENT || key_Ne.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'element number' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_fdvec.place != FACE || key_fdvec.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'face dip vector' variable must be VECTOR and placed on FACE"  );
    if ( key_Nn.place != NODE || key_Nn.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'node number' variable must be SCALAR and placed on NODE"  );
    if ( key_t_out.place != MODEL || key_t_out.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'output times' variable must be ARRAY and placed on MODEL"  );
    if ( key_pf0.place != NODE || key_pf0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'previous fluid pressure' variable must be SCALAR and placed on NODE"  );
    if ( key_QV0.place != ELEMENT || key_QV0.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'previous fluid volume source' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_TMAX.place != MODEL || key_TMAX.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'run duration' variable must be SCALAR and placed on MODEL"  );
    if ( key_thi.place != ELEMENT || key_thi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'thickness' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_V.place != ELEMENT || key_V.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'volume' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_g.place != MODEL || key_g.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'acceleration gravity' variable must be SCALAR and placed on MODEL"  );
    if ( key_kFbr.place != ELEMENT || key_kFbr.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'fault breccia permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_Fs.place != REGION || key_Fs.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'fault size' variable must be SCALAR and placed on REGION"  );
    if ( key_kf.place != ELEMENT || key_kf.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'fracture permeability' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_k.place != ELEMENT || key_k.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'permeability' variable must be TENSOR and placed on ELEMENT"  );
    if ( key_k_cal.place != ELEMENT || key_k_cal.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'permeability calibration factor' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_phi.place != ELEMENT || key_phi.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'porosity' variable must be SCALAR and placed on ELEMENT"  );
    if ( key_RRT.place != REGION || key_RRT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'rocktype' variable must be SCALAR and placed on REGION"  );
    if ( key_AS.place != REGION || key_AS.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Variables_TwoPhaseFlow::Variables_TwoPhaseFlow:",
        "The 'slip patch area' variable must be SCALAR and placed on REGION"  );
  }
};

} } // end namespace csmp

#endif // CSMP_VARIABLES_TWOPHASEFLOW_H
