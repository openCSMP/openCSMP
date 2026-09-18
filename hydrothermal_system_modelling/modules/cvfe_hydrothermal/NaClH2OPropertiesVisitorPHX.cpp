// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <cmath>
#include <algorithm>
#include <sstream>

#include "NaClH2OPropertiesVisitorPHX.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{

namespace
{
// Molar masses [kg/mol] — used in Pass 1 to convert advected mass concentrations
// (mv, ma in kg/m³ pore) to molar fractions.  M_VAPOR cancels in the volume-fraction
// identity used by Pass 2 (ComputeGasMixtureProperties), but NOT in Pass 1 where we
// work from mv() and ma() directly.
constexpr double M_VAPOR = 0.01802;   // H2O
constexpr double M_AIR   = 0.02897;   // dry air
}

template<uint32_t dim>
NaClH2OPropertiesVisitorPHX<dim>::NaClH2OPropertiesVisitorPHX( Model<dim>& model )
//
#include "NaClH2OPropertiesVisitorPHX_initializer_list.h"
//
{
    // ── Temperature / pressure ────────────────────────────────────────────────
    t_key  = pref.StorageKey("temperature");
    tp_key = pref.StorageKey("previous temperature");
    p_key  = pref.StorageKey("fluid pressure");
    pp_key = pref.StorageKey("previous fluid pressure");

    // ── Enthalpy contents (primary transported energy variables) ─────────────
    hCl_key  = pref.StorageKey("enthalpy content liquid");
    hCv_key  = pref.StorageKey("enthalpy content vapor");
    hCa_key  = pref.StorageKey("enthalpy content air");
    hClp_key = pref.StorageKey("previous enthalpy content liquid");
    hCvp_key = pref.StorageKey("previous enthalpy content vapor");
    hCap_key = pref.StorageKey("previous enthalpy content air");

    // ── Nodal scalars ─────────────────────────────────────────────────────────
    phi_key = pref.StorageKey("nodal porosity");
    cpr_key = pref.StorageKey("nodal heat capacity rock");
    ncp_key = pref.StorageKey("nodal heat capacity");
    rr_key  = pref.StorageKey("nodal density rock");
    cpv_key = pref.StorageKey("pore volume");

    // ── Fluid mass concentrations (primary transported mass variables) ────────
    ml_key  = pref.StorageKey("fluid mass liquid");
    mv_key  = pref.StorageKey("fluid mass vapor");
    ma_key  = pref.StorageKey("fluid mass air");
    mlp_key = pref.StorageKey("previous fluid mass liquid");
    mvp_key = pref.StorageKey("previous fluid mass vapor");
    map_key = pref.StorageKey("previous fluid mass air");

    mt_key       = pref.StorageKey("fluid density");
    mtp_key      = pref.StorageKey("previous fluid density");
    rho_bulk_key = pref.StorageKey("bulk fluid density");
    Htp_key      = pref.StorageKey("previous total enthalpy");

    // ── Saturations ───────────────────────────────────────────────────────────
    sl_key = pref.StorageKey("saturation liquid");
    sv_key = pref.StorageKey("saturation vapor");
    sa_key = pref.StorageKey("saturation air");

    // ── Phase densities ───────────────────────────────────────────────────────
    rl_key = pref.StorageKey("density liquid");
    rv_key = pref.StorageKey("density vapor");
    ra_key = pref.StorageKey("density air");

    // ── Viscosities ───────────────────────────────────────────────────────────
    mul_key = pref.StorageKey("viscosity liquid");
    muv_key = pref.StorageKey("viscosity vapor");
    mua_key = pref.StorageKey("viscosity air");

    // ── Specific enthalpies ───────────────────────────────────────────────────
    hf_key = pref.StorageKey("fluid enthalpy");
    hl_key = pref.StorageKey("enthalpy liquid");
    hv_key = pref.StorageKey("enthalpy vapor");
    ha_key = pref.StorageKey("enthalpy air");

    // ── Volumetric enthalpies (hV = h*rho, flux density; hC = hV*S, accumulation) ──
    hVl_key = pref.StorageKey("volumetric enthalpy liquid");
    hVv_key = pref.StorageKey("volumetric enthalpy vapor");
    hVa_key = pref.StorageKey("volumetric enthalpy air");

    // ── Compressibility / heat capacity ──────────────────────────────────────
    cpf_key       = pref.StorageKey("fluid heat capacity");
    beta_key      = pref.StorageKey("compressibility");
    CT_key        = pref.StorageKey("nodal total compressibility");
    beta_rock_key = pref.StorageKey("nodal compressibility rock");
    beta_air_key  = pref.StorageKey("nodal compressibility air");
    apc_key       = pref.StorageKey("after phasechange counter");
    dpc_key       = pref.StorageKey("dangerous phase change");

    // ── Transport densities (post volume-factor correction) ───────────────────
    rl_transport_key = pref.StorageKey("density liquid transport");
    rv_transport_key = pref.StorageKey("density vapor transport");
    ra_transport_key = pref.StorageKey("density air transport");

    // ── Pressure equation ─────────────────────────────────────────────────────
    nQ_key      = pref.StorageKey("nodal fluid volume source");
    vol_fac_key = pref.StorageKey("volume factor");

    // ── State / source ────────────────────────────────────────────────────────
    state_key    = pref.StorageKey("fluid state");
    state_p_key  = pref.StorageKey("previous fluid state");
    src_h_key    = pref.StorageKey("fluid source h");
    src_rate_key = pref.StorageKey("fluid source rate");

    // ── Phase mobilities: mm = kr*rho/mu,  em = kr*hV/mu,  rv = kr/mu ────────
    mml_key  = pref.StorageKey("liquid mass mobility");
    mmv_key  = pref.StorageKey("vapor mass mobility");
    mma_key  = pref.StorageKey("air mass mobility");
    mmld_key = pref.StorageKey("liquid mass mobility density");
    mmvd_key = pref.StorageKey("vapor mass mobility density");
    mmad_key = pref.StorageKey("air mass mobility density");

    eml_key  = pref.StorageKey("liquid enthalpy mobility");
    emv_key  = pref.StorageKey("vapor enthalpy mobility");
    ema_key  = pref.StorageKey("air enthalpy mobility");
    emld_key = pref.StorageKey("liquid enthalpy mobility density");
    emvd_key = pref.StorageKey("vapor enthalpy mobility density");
    emad_key = pref.StorageKey("air enthalpy mobility density");

    rvl_key = pref.StorageKey("relperm viscosity liquid");
    rvv_key = pref.StorageKey("relperm viscosity vapor");
    rva_key = pref.StorageKey("relperm viscosity air");

    // ── Boundary flow accumulators ────────────────────────────────────────────
    bfm_key       = pref.StorageKey("boundary flow mass");
    bfe_key       = pref.StorageKey("boundary flow enthalpy");
    bfs_key       = pref.StorageKey("boundary flow salt");
    bfa_key       = pref.StorageKey("boundary flow air");
    ref_enthalpy_liquid_top_key         = pref.StorageKey("specific enthalpy liquid top");
    ref_enthalpy_vapor_in_air_top_key   = pref.StorageKey("specific enthalpy vapor top");
    ref_enthalpy_air_top_key            = pref.StorageKey("specific enthalpy air top");

    ref_frac_vapor_in_air_top_key       = pref.StorageKey("mass fraction inflow vapor in air top");
    ref_frac_air_top_key                = pref.StorageKey("mass fraction inflow air top");

    rain_recharge_top_key               = pref.StorageKey("rain recharge top");

    // ── NaCl / halite ─────────────────────────────────────────────────────────
    wt_key   = pref.StorageKey("salinity");
    sh_key   = pref.StorageKey("saturation halite");
    rh_key   = pref.StorageKey("density halite");
    mh_key   = pref.StorageKey("solid mass halite");
    hh_key   = pref.StorageKey("enthalpy halite");
    hVh_key  = pref.StorageKey("volumetric enthalpy halite");
    hCh_key  = pref.StorageKey("enthalpy content halite");
    xf_key   = pref.StorageKey("salt fraction fluid");
    xl_key   = pref.StorageKey("salt fraction liquid");
    xv_key   = pref.StorageKey("salt fraction vapor");
    xh_key   = pref.StorageKey("salt fraction halite");
    xVl_key  = pref.StorageKey("volumetric salinity liquid");
    xVv_key  = pref.StorageKey("volumetric salinity vapor");
    xVh_key  = pref.StorageKey("volumetric salinity halite");
    xCl_key  = pref.StorageKey("salt content liquid");
    xCv_key  = pref.StorageKey("salt content vapor");
    xCh_key  = pref.StorageKey("salt content halite");
    xCf_key  = pref.StorageKey("salt content fluid");
    xClp_key = pref.StorageKey("previous salt content liquid");
    xCvp_key = pref.StorageKey("previous salt content vapor");
    xml_key  = pref.StorageKey("liquid salt mobility");
    xmv_key  = pref.StorageKey("vapor salt mobility");
    src_wt_key = pref.StorageKey("fluid source wt");
    msp_key    = pref.StorageKey("previous mass salt");

    // ── Vapor-air mixture quantities ──────────────────────────────────────────
    // Registered unconditionally (same pattern as air keys); only meaningful when
    // gas_mixture = true.  These are diagnostic output variables:
    // Pass 1 writes p_partial_vapor/air (from mv/ma), Pass 2 overwrites them with
    // the equilibrated values (from sv/sa).  Neither pass reads them back —
    // previous-step partial pressures are not used anywhere in the computation.

    xa_mix_key          = pref.StorageKey("gas mix molar fraction air");
    xv_mix_key          = pref.StorageKey("gas mix molar fraction vapor");
    p_partial_air_key   = pref.StorageKey("partial pressure air");
    p_partial_vapor_key = pref.StorageKey("partial pressure vapor");
    rho_gas_mix_key     = pref.StorageKey("gas mix density");
    mu_gas_mix_key      = pref.StorageKey("gas mix viscosity");
    h_gas_mix_key       = pref.StorageKey("gas mix enthalpy");
    beta_gas_mix_key    = pref.StorageKey("gas mix compressibility");

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    Liquid.InitToZero();
    Vapor.InitToZero();
    Bulk.InitToZero();
    Salt.InitToZero();
}

template<uint32_t dim>
NaClH2OPropertiesVisitorPHX<dim>::~NaClH2OPropertiesVisitorPHX()
{}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::Visit( Model<dim>* n )
{}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::Visit( Region<dim>* n )
{}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::Visit( Node<dim>* n )
{
    // ============================================================================
    // All variables are normalised to UNIT PORE VOLUME (i.e. "per m³ pore space")
    // unless noted otherwise.
    //
    // ── PRE-EQUILIBRATOR ──────────────────────────────────────────────────────
    //
    // STEP 1 — ReadAllVariables
    //   Reads the post-transport nodal state from the CSMP mesh:
    //     ml, mv, ma       : phase mass concentrations [kg/m³ pore]
    //     mlp, mvp, map    : previous-step values (for computing increments)
    //     hCl, hCv, hCa    : phase enthalpy contents [J/m³ pore]
    //     hClp, hCvp, hCap : previous-step enthalpy contents
    //     Htp              : total enthalpy of previous accepted step [J/m³ pore]
    //     t, p, wt, ...    : temperature, pressure, salinity, etc.
    //
    // STEP 2 — CalculateAbsoluteVariables
    //   Computes mass and enthalpy increments from transported variables:
    //     dml_ = ml - mlp,  dmv_ = mv - mvp,  dma_ = ma - map
    //     dhCl_ = phi*(hCl - hClp),  dhCv_ = phi*(hCv - hCvp),  dhCa_ = ...
    //     dh_rock_diff_ = rock enthalpy change from temperature shift
    //   Updates total mass: mt() = mtp() + dml_ + dmv_ [+ dma_]
    //   Updates total salt: ms() = msp() + dxCl_ + dxCv_
    //   Freezes interior DIRICH nodes (p+T fixed → zero all increments).
    //
    // STEP 3 — UpdateEquilibratorVariables  (PRE-Equilibrator)
    //   Assembles the five Equilibrator inputs from the increments above:
    //     H_current_ = Htp + dh_rock_diff_ + dh_fluid_diff_ + dhCl_ + dhCv_ [+ dhCa_]
    //     m_fluid_   = (mt - ma) * pore_volume    [air mass excluded]
    //     wt_        = ms / (mt - ma) * 100       [air excluded from denominator]
    //     phi_       = phi * (1 - sa)              [air volume shaved]
    //     p_current_ = pressure for the Equilibrator (see below)
    //
    //   3a. MIXTURE — Pass 1: partial pressures from advected masses:
    //     nv = mv/M_VAPOR,  na = ma/M_AIR,  xv = nv/(nv+na)
    //     p_partial_air() = p*xa,  p_partial_vapor() = p*xv  (stored for diagnostics).
    //     p_current_: if xv >= XV_MIN → p*xv (vapor partial pressure for Equilibrator).
    //                 if xv < XV_MIN  → p_total (liquid under total pressure; physically
    //                 correct when H2O vapor is negligible; avoids near-zero pressure).
    //   IMMISCIBLE / no air: p_current_ = p_total.
    //
    //   3b. sa() and phi_ shaving:
    //     sa() and ra() computed with ComputeAirSaturationAndDensity(ma, t, p, p_partial_air())  [sa() is clamped]
    //     phi_ = phi * (1 - sa)
    //     In mixture mode p_partial_air() is already set (3a) — exact.
    //     In immiscible mode p_partial_air() = p_total — also exact.
    //
    // ── EQUILIBRATOR ──────────────────────────────────────────────────────────
    //
    // STEP 4 — Equilibrate
    //   Calls H2O-NaCl Equilibrator on the air-shaved pore space at p_current_.
    //   Returns Liquid.s, Vapor.s, Salt.s (sum to 1 over shaved space) and
    //   Vapor.rho/mu/beta/h at p_current_ (= p*xv in mixture mode).
    //   At open boundaries, calls BoundaryIteration instead (root-finding for bfm).
    //
    // STEP 4a — CheckPhaseChange / CheckVolumeMismatchCompensation (optional)
    //
    // ── POST-EQUILIBRATOR ─────────────────────────────────────────────────────
    //
    // STEP 5 — UpdateCSMPVariables  (POST-Equilibrator)
    //   Translates Equilibrator output back to full CSMP variables:
    //     t() = Bulk.t, phase densities/viscosities/enthalpies, beta() = Bulk.beta
    //     ml() = rl()*sl(),  mv() = rv()*sv()   [re-derived from equilibrium]
    //     Volumetric/content enthalpies, salt fractions, ncp(), rho_bulk()
    //   5a. Recompute ra() and sa() at the equilibrated temperature t_ = Bulk.t.
    //       MIXTURE:    ra() = air.Density(t_, p_partial_air())
    //                   sa() = ComputeAirSaturation(ma, t_, p, p_partial_air())
    //       IMMISCIBLE: ra() = air.Density(t_, p_total)
    //                   sa() = ComputeAirSaturation(ma, t_, p, p_total)
    //   5b. Rescale Equilibrator saturations to full pore space:
    //         sl, sv, sh = Equilibrator.s * (1 - sa)
    //       Closure: sl+sv+sh+sa = 1 exactly by construction.
    //   5c. MIXTURE — Pass 2: ComputeGasMixtureProperties with rescaled sv/sa:
    //         xv_mix, xa_mix, p_partial_vapor/air, rho/mu/h/beta_gas_mix
    //         ra() updated to Pass-2 partial pressure; sa() NOT updated.
    //
    // STEP 6 — VolumeFactorComputations
    //   nQ = mt - rho_bulk_with_air  (mass mismatch driving pressure equation).
    //   volume_factor_LHS/RHS: H2O-NaCl density-mismatch corrections (air excluded).
    //   Scales transport variables (ml, mv, hCl, hCv, rl_transport, etc.) by vf.
    //   MIXTURE: rv_transport and ra_transport are rescaled by s_gas/sv and s_gas/sa
    //   respectively, so that mmv = krv*rv_transport/mu correctly gives
    //   kr_gas*rv_individual/mu (density-fraction allocation, not volume-fraction).
    //   hVv/hVa are updated to use the rescaled transport densities so that
    //   emv/mmv = hv and ema/mma = ha hold exactly.
    //
    // STEP 7 — PrepareVariablesForStorage
    //   CT (total compressibility), relative permeabilities, phase mobilities,
    //   gravity terms — all using mixture properties where applicable.
    //   Stores Htp = H_current_, mtp = mt(), msp = ms() for next timestep.
    //
    // STEP 8 — StorePropertiesAndFlags
    //   Writes all variables back to CSMP mesh for the next timestep.
    //
    // ============================================================================

    // ****************************************
    // 1. Read all nodal variables of interest
    // ****************************************

    ReadAllVariables( n );

    const int nid = n->Idx();
    const bool node_dbg = (nid == -1);//-1: none

    auto log_stage = [&](const char* stage) {
        if (!node_dbg) return;
        cerr << "[MASS " << stage
             << "]"
             << "  mt="  << mt()
             << "  ml="  << ml();
        if (with_air) cerr << "  ma=" << ma();
        cerr << "  sl="  << sl();

        if (with_air) cerr << "  sa=" << sa();
        cerr << "  rl="  << rl();
        if (with_air) cerr << "  ra=" << ra();

        cerr << "  mml=" << mml();
        if (with_air) cerr << "  mma=" << mma();
        {
            const double krl_now = RelativePermeabilityLiquid(sl(), sv(), with_air ? sa() : 0.);
            const double krg_now = RelativePermeabilityGas    (sl(), sv(), with_air ? sa() : 0.);
            cerr << "  krl=" << krl_now
                 << "  krg=" << krg_now;
        }
        cerr << "  p="   << p()
             << "  T="   << t()
             << "  state=" << Bulk.state
             << endl;

    };

    // One-time header showing the transport increments (these are fixed for this timestep)
    if (node_dbg && verbose_output)
    {
        cerr << "\n========== Visit node " << nid << " ==========" << endl
             << "[TRANSPORT id=" << nid << "]"
             << "  mlp=" << mlp()
             << "  mvp=" << mvp();
        if (with_air) cerr << "  map=" << map();
        cerr << "  mtp=" << mtp()
             << "\n             ml_in=" << ml()
             << "  mv_in=" << mv();
        if (with_air) cerr << "  ma_in=" << ma();
        cerr << "\n             dml=" << (ml() - mlp())
             << "  dmv=" << (mv() - mvp());
        if (with_air) cerr << "  dma=" << (ma() - map());
        cerr << endl;
    }

    if(verbose_output)
    {
        log_stage("read");                  // raw mesh values: transport solution
        CheckPhysicalState( n, "ReadAllVariables" );   // ← raw mesh inputs
    }

    old_state = int(state()+0.01); // should always give the correct result
    state_p() = double(old_state);

    // ****************************************
    // 2. Initialize some bools
    // ****************************************

    bogus_variables   = false;
    pure_halite       = false;
    CheckForOutOfRange( n );
    CheckBoundaryFlags( n );

    // ****************************************
    // 3. Prepare a few variables for thermal equilibration code
    // ****************************************

    CalculateAbsoluteVariables( n );
    if(verbose_output)
        log_stage("post_abs");              // mt assembled from increments

    UpdateEquilibratorVariables( n );

    if(verbose_output)
    {
        log_stage("pre_equil");             // m_fluid, phi_, p_current_ prepared
        CheckPhysicalState( n, "PreEquilibrate" );     // ← what we feed the EOS
    }

    t_diffusion_       = t();

    // ****************************************
    // 4. EQUILIBRATE, core task of the visitor
    // ****************************************

    Equilibrate( n );

    if(verbose_output)
        log_stage("post_equil");            // Bulk, Liquid, Vapor, Salt populated


    // ****************************************
    // 5. Prepare variables to be passed back to CSMP
    // ****************************************

    UpdateCSMPVariables( n );

    if(verbose_output)
        log_stage("post_update");           // sl rescaled, ml = rl*sl OVERWRITE HERE


    // ****************************************
    // 6. Volume Mismatch
    // ****************************************

    VolumeFactorComputations( n );

    if(verbose_output)
        log_stage("post_vf");               // ml *= vol_fac_LHS CORRECTION HERE

    PrepareVariablesForStorage();

    if(verbose_output)
    {
        log_stage("post_prep");
        CheckPhysicalState( n, "PostPrepare" );        // ← final state before mesh store
    }

    if(Bulk.state!=none)
    {
        StorePropertiesAndFlags( n );
    }

    else
    {
        cerr << endl << "[Visit] Fluid state undefined — not storing result from equilibration.";
        cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        DumpNodeState("Visit/UndefinedState");
        csmp_error.Note( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                         "\nFluid's state is undefined, not storing result from equilibration!");
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::ReadAllVariables( Node<dim>* n )
{

    if(with_rock_liquidus_solidus)
    {
        n->Read( tl_key, tl );
        n->Read( ts_key, ts );
        tl_ = tl();
        ts_ = ts();
        equilibrator.SetRockLiquidusSolidusTemperatures( tl_, ts_ );
    }

    if(with_air)
    {
        n->Read( hCa_key,      hCa  );
        n->Read( hCap_key,     hCap );
        n->Read( ma_key,       ma   );
        n->Read( map_key,      map  );
        n->Read( sa_key,       sap  );
        n->Read( ra_key,       ra   );
        n->Read( cpa_key,      cpa  );
        n->Read( beta_air_key, beta_air );

        // open-boundary previous values
        n->Read( mua_key, muap );
        n->Read( ha_key,  hap  );
        n->Read( ra_key,  rap  );

        n->Read( ref_enthalpy_vapor_in_air_top_key, ref_enthalpy_vapor_in_air_top );
        n->Read( ref_enthalpy_air_top_key, ref_enthalpy_air_top );

        n->Read( ref_frac_vapor_in_air_top_key, ref_frac_vapor_in_air_top );

        n->Read( rain_recharge_top_key, rain_recharge_top );
    }

    // ── Primary state ─────────────────────────────────────────────────────────
    n->Read( t_key,   t   );
    n->Read( tp_key,  tp  );
    n->Read( p_key,   p   );
    n->Read( pp_key,  pp  );
    n->Read( mt_key,  mt  );
    n->Read( phi_key, phi );
    n->Read( rr_key,  rr  );
    n->Read( cpr_key, cpr );
    n->Read( cpv_key, cpv );

    // ── Transported mass ──────────────────────────────────────────────────────
    n->Read( ml_key,  ml  );
    n->Read( mv_key,  mv  );
    n->Read( mlp_key, mlp );
    n->Read( mvp_key, mvp );
    n->Read( mtp_key, mtp );

    // ── Transported enthalpy ──────────────────────────────────────────────────
    n->Read( hCl_key,  hCl  );
    n->Read( hCv_key,  hCv  );
    n->Read( hClp_key, hClp );
    n->Read( hCvp_key, hCvp );
    n->Read( Htp_key,  Htp  );

    // ── Sources / state ───────────────────────────────────────────────────────
    n->Read( src_h_key,    src_h    );
    n->Read( src_rate_key, src_rate );
    n->Read( beta_rock_key, beta_rock );
    n->Read( state_key,     state    );
    n->Read( apc_key,       after_phasechange_counter );
    n->Read( dpc_key,       dangerous_phase_change    );

    // ── NaCl ──────────────────────────────────────────────────────────────────
    n->Read( xCl_key,  xCl  );
    n->Read( xCv_key,  xCv  );
    n->Read( xCf_key,  xCf  );
    n->Read( xClp_key, xClp );
    n->Read( xCvp_key, xCvp );
    n->Read( mh_key,   mh   );
    n->Read( sh_key,   sh   );
    n->Read( msp_key,  msp  );
    n->Read( src_wt_key, src_wt );
    n->Read( wt_key,     wt     );

    // ── Open-boundary / previous-step values ─────────────────────────────────
    // slp/svp: previous saturations — needed by BoundaryFlow mobility and by
    // the mixture molar-fraction computation in UpdateEquilibratorVariables.
    n->Read( sl_key,   slp  );
    n->Read( sv_key,   svp  );
    n->Read( rho_bulk_key, rho_bulk );
    n->Read( mul_key, mulp );
    n->Read( muv_key, muvp );
    n->Read( xl_key,  xlp  );
    n->Read( xv_key,  xvp  );
    n->Read( hl_key,  hlp  );
    n->Read( hv_key,  hvp  );
    n->Read( hh_key,  hhp  );
    n->Read( rl_key,  rlp  );
    n->Read( rv_key,  rvp  );

    n->Read( ref_enthalpy_liquid_top_key, ref_enthalpy_liquid_top );

}// end ReadAllVariables

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::CalculateAbsoluteVariables( Node<dim>* n )
{
    // Initialize boundary fluxes to 0
    bfm() = bfe() = bfs() = 0.0;
    if(with_air)
    {
        bfa() = 0.0;
        ref_frac_air_top() = 0.;
    }


    // Rock enthalpy change (accounts for optional crystallisation curve / latent heat)
    if(with_rock_liquidus_solidus)
    {
        hrock_prev = rock.Enthalpy(tp(), tl(), ts());
        hrock_curr = rock.Enthalpy(t(),  tl(), ts());
    }
    else
    {
        hrock_prev = rock.Enthalpy(tp());
        hrock_curr = rock.Enthalpy(t());
    }

    // Conductive enthalpy increments (fluid and air diffusion not yet implemented)
    dh_rock_diff_  = (1.0 - phi()) * rr() * (hrock_curr - hrock_prev);
    dh_fluid_diff_ = 0.;
    if(with_air) dh_air_diff_ = 0.;

    // Advective enthalpy increments
    dhCl_ = phi() * (hCl() - hClp());
    dhCv_ = phi() * (hCv() - hCvp());
    if(with_air)
        dhCa_ = phi() * (hCa() - hCap());

    dP = phi() * (p() - pp());  // informational; not wired into H_current_

    // Advective mass increments
    dml_ = ml() - mlp();
    dmv_ = mv() - mvp();
    if(with_air) dma_ = ma() - map();

    mt() = mtp() + dml_ + dmv_;
    if(with_air) mt() += dma_;

    // Round-off dust from exact depletion: when a node drains fully, the
    // cancellation mtp - mlp can land a machine-epsilon hair below zero
    // ── DUST/CONSISTENCY RULE 2 of 3 — mt increment residue ─────────────────
    // (Rule 1: SBL beat landings, in UpdateTransportVariables. Rule 3: the
    //  EMPTY-NODE override in UpdateCSMPVariables below. Transport has NO
    //  rule — clean by construction.) Relative tolerance: the increment
    //  mt = mtp + dml + dmv cancels large terms; residue scales with them.
    const double mt_dust_tol = 1.0e-12 * (fabs(mtp()) + fabs(mlp()) + fabs(mvp()));
    if ( mt() < 0. && mt() > -mt_dust_tol )
        mt() = 0.;

    if ( mt() < -mt_dust_tol )
        cerr << "\n[mt<0] node=" << n->Idx()
             << " ml=" << ml() << " mv=" << mv()
             << " mlp=" << mlp() << " mvp=" << mvp() << " mtp=" << mtp();

    // ── DIRICH interior nodes: freeze all state ───────────────────────────────
    // Both p and T fixed → node is thermodynamically pinned; zero all increments.
    if (p.Flag() == DIRICH && t.Flag() == DIRICH && n->AtBoundary() == NOT)
    {
        dh_rock_diff_ = dh_fluid_diff_ = dhCl_ = dhCv_ = dml_ = dmv_ = 0.;
        hCl() = hClp();    hCv() = hCvp();
        ml()  = mlp();     mv()  = mvp();
        mt()  = mtp();

        if(with_air)
        {
            dh_air_diff_ = dhCa_ = dma_ = 0.;
            hCa() = hCap();
            ma()  = map();
        }
    }

    // NaCl mass increments
    dxCl_    = xCl() - xClp();
    dxCv_    = xCv() - xCvp();
    dx_diff_ = 0.;//No salt diffusion

    ms() = msp() + dxCl_ + dxCv_ + dx_diff_;

    // Safety clamp: salt mass cannot go negative. CSMP transport of xCv/xCl
    // can overshoot msp at phase-transition nodes (e.g. large dmv, small msp).
    // The salt volume factor in VolumeFactorComputations prevents recurrence;
    // this clamp handles the first timestep after a large phase transition.
    // Commented for now since the salt volume factor should fix.
    // if (ms() < 0.) ms() = 0.;

    if (p.Flag() == DIRICH && t.Flag() == DIRICH && wt.Flag() == DIRICH && n->AtBoundary() == NOT)
    {
        dxCl_ = dxCv_ = 0.;
        xCl() = xClp();    xCv() = xCvp();
        ms()  = msp();
    }

}// end CalculateAbsoluteVariables

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::UpdateEquilibratorVariables( Node<dim>* n )
{
    // ====================================================================
    // PURPOSE: Prepare all inputs the Equilibrator needs.
    //   Called PRE-Equilibrator, after CalculateAbsoluteVariables.
    //   Leaves the following ready for Equilibrate():
    //     phi_       : porosity shaved by air (phi * (1-sa))
    //     m_fluid_   : H2O-NaCl mass only (air excluded)
    //     wt_        : salinity excluding air from denominator
    //     H_current_ : total enthalpy including rock, fluid AND air energy.
    //                  The Equilibrator strips air energy internally before
    //                  bisection — do NOT subtract it here.
    //     p_current_ : pressure the Equilibrator will use
    //                  MIXTURE: p * xv  (vapor partial pressure, Dalton's law)
    //                  IMMISCIBLE / no air: p_total
    //
    // AIR HANDLING:
    //   Pass 1 (mixture mode) runs FIRST so that p_partial_air() is known before
    //   ComputeAirSaturation is called.  This means sa() and phi_ are computed at
    //   the correct partial pressure in mixture mode — no approximation is needed.
    //   In immiscible mode p_partial_air() = p_total, which is also exact.
    //   UpdateCSMPVariables refines ra() and sa() only because the Equilibrator
    //   may shift xv slightly (by redistributing H2O between phases), which is a
    //   genuine second-order physical correction, not a numerical approximation.
    // ====================================================================

    t_  = t();
    tp_ = tp();
    p_  = p();

    // ── Rock ──────────────────────────────────────────────────────────────────
    m_rock_   = rr() * (1.0 - phi());
    rho_rock_ = rr();

    if(with_rock_liquidus_solidus)
        cp_rock_ = rock.HeatCapacity(t_, tl_, ts_);
    else
        cp_rock_ = rock.HeatCapacity(t_);

    // ── Air ───────────────────────────────────────────────────────────────────
    if(with_air)
    {
        m_air_  = ma() * phi();
        cp_air_ = air.HeatCapacity(t_);
    }

    // ── Pass 1: partial pressures from advected masses (MIXTURE) ─────────────
    // Runs before phi shaving so p_partial_air() is available for ComputeAirSaturation.
    // IMMISCIBLE / no air: p_partial_air() = p_current_ = p_total.
    //
    // xv is derived from the current advected masses mv() and ma():
    //   nv = mv / M_VAPOR,   na = ma / M_AIR   [mol/m³ pore]
    //   xv = nv / (nv + na)
    // Using advected masses is exact for the current timestep composition.
    // On a timestep reset mv→mvp and ma→map automatically, giving the correct
    // previous-step composition with no extra stored variable needed.
    //
    // Edge cases:
    //   ma()=0 (no air / pure vapor): na=0 → xv=1 → p_current_=p  ✓
    //   mv()=0 (air-dominated, no H2O vapor): xv=0 → see below
    //
    // When xv < XV_MIN (H2O vapor is negligible), p_current_ falls back to p_total.
    // Rationale: the H2O liquid is under the total pore pressure, not zero.  Passing
    // p*xv≈0 to the Equilibrator would be physically wrong (inventing zero vapor
    // pressure on a liquid that is actually compressed) and numerically unstable.
    // p_partial_vapor/air are still stored at their true Dalton values for diagnostics.

    constexpr double XV_MIN = 1.0e-6;  // below this xv → treat as liquid under total p
    if (with_air && full_gas_mixture)
    {
        const double nv   = (mv() > 0.) ? mv() / M_VAPOR : 0.;
        const double na   = (ma() > 0.) ? ma() / M_AIR   : 0.;
        const double ntot = nv + na;

        const double xv = (ntot > 0.) ? nv / ntot : 1.0;
        const double xa = 1.0 - xv;

        p_partial_vapor() = p_ * xv;
        p_partial_air()   = p_ * xa;

        if (xv < XV_MIN)
        {
            // No meaningful H2O vapor in the gas phase — system is liquid.
            // Liquid is in mechanical equilibrium at p_total; no vapor-liquid
            // equilibrium to resolve.  p_current_ = p_total is exact.
            p_current_ = p_;
        }
        else
        {
            p_current_ = p_partial_vapor();
        }
    }
    else
    {
        p_partial_vapor() = p_;
        p_partial_air()   = p_;
        p_current_        = p_;
    }

    // ── Porosity shaving ──────────────────────────────────────────────────────
    // p_partial_air() is now set (above), so ComputeAirSaturation uses the
    // correct pressure in both mixture and immiscible modes.
    phi_ = phi();
    if(with_air) ComputeAirSaturationAndDensity( ma(), t_, p_, p_partial_air() );
    phi_ *= (1.0 - sa());
    if (with_air && phi_ < 1.e-4 * phi())
        cerr << "\n[UpdateEquilibratorVariables] phi_ near-zero after air shaving:"
             << "  phi_=" << phi_ << "  phi=" << phi() << "  sa=" << sa()
             << "  sa_max_lim=" << (1.-sh()-min_h2o_nacl_fraction) << endl;

    // ── H2O-NaCl mass (air mass hidden from Equilibrator) ─────────────────────
    m_fluid_ = with_air ? (mt() - ma()) * phi()
                        :  mt()         * phi();

    // ── Salinity (air excluded from denominator) ──────────────────────────────
    wt_ = (ms() > 0. && mt() > ma())
            ? (with_air ? ms()/(mt()-ma())*100. : ms()/mt()*100.)
            : 0.;
    if(wt_ > 100.)
    {
        cerr << "\n[UpdateEquilibratorVariables] wt_ clamped " << wt_ << " -> 100 wt%"
             << "  ms=" << ms() << "  mt=" << mt();
        if(with_air) cerr << "  ma=" << ma();
        cerr << endl;
        wt_ = 100.;
    }

    // ── Total enthalpy ────────────────────────────────────────────────────────
    // Air volume/mass are shaved from phi_ and m_fluid_; air energy is stripped
    // internally by Equilibrator — do NOT subtract it here.
    H_previous_ = Htp();
    H_current_  = Htp()
            + dh_fluid_diff_              // H2O-NaCl conduction (currently zero)
            + dh_rock_diff_               // rock conduction
            + dhCl_ + dhCv_;             // liquid + vapor advection
    if(with_air) H_current_ += dh_air_diff_ + dhCa_;  // air conduction (zero) + advection

    //H_current_ += dP;  // work term — correct in principle, kept off for now

    if (p.Flag() == DIRICH && t.Flag() == DIRICH && n->AtBoundary() == NOT)
        H_current_ = H_previous_;  // interior DIRICH: pin enthalpy

} //end UpdateEquilibratorVariables

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::Equilibrate( Node<dim>* n )
{

    // ========================================================================
    // EQUILIBRATE — Compute thermodynamic equilibrium between all phases.
    //
    // AIR SHAVING:
    //
    //   The H2O-NaCl equilibrator (Equilibrator) has no knowledge of air:
    //
    //   Equilibrator strips air ENERGY internally (subtracts mass_air *
    //   air.Enthalpy(t) before bisection) — do NOT subtract it here or the
    //   energy will be double-stripped.
    //
    //   This visitor has already shaved air VOLUME and MASS from phi_ and m_fluid_,
    //   so Equilibrator sees a smaller pore space containing only H2O-NaCl-Halite.
    //   Returned sl, sv, sh sum to 1 relative to that reduced pore space.
    // ========================================================================

    // At boundary: go to BoundaryIteration() rather than the thermal equilibrator.
    // Nodes with Dirichlet salinity (e.g. ocean BC) skip BoundaryIteration and
    // go through normal equilibration — their salinity is externally controlled.
    //if (open_boundaries && t.Flag() == DIRICH && n->AtBoundary() != NOT && n->AtBoundary() != INTERNAL)
    //{

    if (open_boundaries && t.Flag() == DIRICH && wt.Flag() != DIRICH)
    {
        BoundaryIteration();
    }

    else
    {
        thermal_eq_succeeded = equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);

        if(!thermal_eq_succeeded)
        {
            cerr << "\n[Equilibrate] ThreePhaseProperties failed.";
            cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
            cerr << endl << "  p=" << p() << " Pa  t=" << t() << " oC  wt=" << wt_ << " wt%"
                 << "  mt=" << mt() << "  mtp=" << mtp();
            if (std::abs(p_current_ - p()) > 1.0)
                cerr << "  p_current=" << p_current_ << " Pa";
            cerr << "  H=" << H_current_ << "  m_fluid=" << m_fluid_;
            if (with_air) cerr << "  phi_shaved=" << phi_ << "  ma=" << ma() << "  sa=" << sa();
            cerr << "  rho_bulk=" << rho_bulk() << endl;
            //DumpNodeState("Equilibrate/ThreePhasePropertiesFailed");
        }
    }

    if(equilibrator.Fatal())
    {
        cerr << endl << "[Equilibrate] FATAL: received Fatal() from equilibrator. Conditions:";
        cerr << endl << "Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        DumpNodeState("Equilibrate/Fatal");
        csmp_error.Note( FATAL_ERROR,
                         "NaClH2OPropertiesVisitorPHX<dim>::Visit( Node<dim>* n ) -",
                         "received message Fatal() from equilibrator ... terminating!!!");
    }

    current_state = Bulk.state; // that should be type safe
    state()       = double(current_state);

}// end Equilibrate

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables( Node<dim>* n )
{
    // ====================================================================
    // PURPOSE: Translate Equilibrator output back to full CSMP variables.
    //
    // CONTEXT (what happened before this function):
    //   UpdateEquilibratorVariables (PRE-Equilibrator):
    //     - Pass 1 ran first: computed xv from mv/ma, set p_partial_air() = p*xa,
    //       p_partial_vapor() = p*xv (diagnostic).
    //       p_current_ = p*xv if xv >= XV_MIN, else p_total (liquid fallback).
    //     - sa() = ComputeAirSaturation(ma, t, p, p_partial_air()) — exact,
    //       because p_partial_air() was already set by Pass 1.
    //     - phi_ = phi*(1-sa) shaved so Equilibrator saw only H2O-NaCl pore space.
    //   Equilibrate() was then called by Visit() at p_current_ and returned:
    //     - Liquid.s, Vapor.s, Salt.s summing to 1 over the shaved pore space.
    //     - Vapor.rho/mu/beta/h evaluated at p_current_ (= p*xv in mixture mode).
    //
    // WHAT THIS FUNCTION DOES (in order):
    //   1. Update temperature and rock/air heat capacities.
    //   2. Update ra() and sa() at the new equilibrated temperature t_ = Bulk.t.
    //      The pre-Equil sa() used the pre-Equil temperature; air density is
    //      temperature-dependent so both must be recomputed before rescaling.
    //      Partial pressures are still Pass-1 values at this point.
    //      MIXTURE:    ra() = air.Density(t_, p_partial_air());
    //      IMMISCIBLE: ra() = air.Density(t_, p_total);
    //   3. Rescale Equilibrator saturations to full pore space:
    //        sl, sv, sh = Equilibrator.s * (1 - sa)
    //      Closure: sl+sv+sh+sa = 1 exactly by construction.
    //   4. MIXTURE Pass 2: call ComputeGasMixtureProperties with rescaled
    //      sv/sa and equilibrated Vapor.rho/mu/beta/h to get:
    //        xv_mix, xa_mix         : final molar fractions from equilibrated state
    //        p_partial_vapor/air    : refined partial pressures (p*xv_pass2, p*xa_pass2)
    //        rho/mu/h/beta_gas_mix  : mixture bulk properties
    //        ra()                   : air density at p*xa_pass2 (updated; sa NOT)
    //      sa() is not updated in Pass 2 — changing sa() would break the
    //      sl+sv+sh+sa=1 closure established by the rescaling above.
    //
    // ma() vs ra()*sa() — definitive rule (unchanged):
    //   sa() = sa_clamped = min(ma()/ra(), sa_max),  sa_max=(1-sh)*(1-res_sl)
    //   Use ma()       for: mass balance, transport (dma_, bfa, m_fluid_), hCa=ha*ma
    //   Use ra()*sa()  for: rho_bulk, nQ, CT pore fractions, ncp, hVa=ha*ra
    // ====================================================================

    t() = t_ = Bulk.t;   // new equilibrated temperature

    if(with_rock_liquidus_solidus)
        cpr() = cp_rock_ = rock.HeatCapacity(t_, tl_, ts_);
    else
        cpr() = cp_rock_ = rock.HeatCapacity(t_);

    if(with_air) cp_air_ = air.HeatCapacity(t_);

    ncp() = cpr() * rr() * (1. - phi());   // rock contribution; fluid/air added below

    // ── Saturations ───────────────────────────────────────────────────────────
    if(with_air)
    {
        if (full_gas_mixture)
        {
            // ── MIXTURE: update sa() and ra() at equilibrated temperature ─────────
            // The pre-Equil sa() was computed at the pre-Equil temperature.
            // t_ is now Bulk.t (equilibrated); air density changes with temperature,
            // so ra() and sa() must be recomputed before the saturation rescaling.
            // p_partial_air() is still the Pass-1 value here (Pass 2 hasn't run yet).
            // Both calls use previous-step sh() for sa_max — see ComputeAirSaturation.
            // Guard: only recompute when air is present; zero sa when ma=0 so that
            // saturation rescaling below (sl/sv/sh *= (1-sa)) is unaffected.
            if (ma() > 0.)
            {
                ComputeAirSaturationAndDensity(ma(), t_, p_, p_partial_air());
            }
            else
            {
                ra() = 0.;
                sa() = 0.;
            }
        }
        else
        {
            // ── IMMISCIBLE: update ra(), sa() at equilibrated temperature ────────────
            // sa() is recomputed here (not retained from pre-Equilibrator) because
            // t_ is now Bulk.t which may differ from the pre-Equil temperature, and
            // air density is temperature-dependent.
            // sa() is used immediately below for saturation rescaling:
            //   sl+sv+sh = (Liquid.s+Vapor.s+Salt.s)*(1-sa) = 1-sa  →  closure ✓
            // Properties are zeroed when ma=0 so visualization shows clean zeros.
            if (ma() > 0.)
            {
                ComputeAirSaturationAndDensity(ma(), t_, p_, p_);
                mua() = air.Viscosity(t_);
                ha()  = air.Enthalpy(t_);
            }
            else
            {
                ra()  = 0.;
                mua() = 0.;
                ha()  = 0.;
                sa()  = 0.;
            }
        }

        // ── Rescale Equilibrator saturations to full pore space ───────────────
        // Equilibrator worked on phi*(1-sa), so its sl/sv/sh sum to 1 over that
        // reduced space.  Multiply by (1-sa) to recover fractions of total pore space.
        // sa() is now finalised for this rescaling; closure holds exactly.
        sl() = Liquid.s * (1.0 - sa());
        sv() = Vapor.s  * (1.0 - sa());
        sh() = Salt.s   * (1.0 - sa());
        // Closure: sl + sv + sh + sa = (Liquid.s+Vapor.s+Salt.s)*(1-sa) + sa = 1 ✓

        if (full_gas_mixture)
        {
            if (ma() > 0.)
            {
                // ── MIXTURE Pass 2 (POST-Equilibrator) ───────────────────────────────
                // The Equilibrator may have redistributed H2O between liquid and vapor,
                // shifting sv() and sa() slightly from their Pass-1 values. The Dalton
                // partial pressures must therefore be recomputed from the equilibrated
                // saturations sv/sa before calling ComputeGasMixtureProperties, so that
                // ra_new is evaluated at the correct Pass-2 partial pressure.
                //
                // This is the authoritative partial-pressure update: Pass-1 used
                // molar-mass-derived fractions from transported mv/ma (the best estimate
                // before the EOS ran); Pass-2 uses volume fractions from equilibrated
                // sv/sa (exact for ideal gases, consistent with the EOS state).
                //
                // sa() is deliberately NOT updated after this call: the change in
                // p_partial_air between Pass 1 and Pass 2 is second-order, and
                // modifying sa() would break the sl+sv+sh+sa=1 closure established
                // by the saturation rescaling above.
                const double s_gas_eq = sv() + sa();  //FIX 30-03-26
                if (s_gas_eq > 0.)
                {
                    p_partial_vapor() = p_ * sv() / s_gas_eq;  //FIX 30-03-26
                    p_partial_air()   = p_ * sa() / s_gas_eq;  //FIX 30-03-26
                }

                double ra_new;
                ComputeGasMixtureProperties(
                            t_, p_,
                            sv(), sa(),
                            Vapor.rho, Vapor.mu, Vapor.beta, Vapor.h,
                            p_partial_air(), p_partial_vapor(),   // Pass-2 Dalton pressures (just updated)
                            xa_mix(), xv_mix(),
                            ra_new,
                            mu_gas_mix(), rho_gas_mix(), h_gas_mix(), beta_gas_mix());

                ra() = ra_new;   // refined air density at Pass-2 p*xa; sa() left unchanged

                // Each phase carries its own specific enthalpy — h_gas_mix is diagnostic
                // only and must NOT replace ha() or hv() (would cross-contaminate budgets).
                // Both phases share mu_gas_mix so they travel at the same Darcy velocity.
                mua()   = mu_gas_mix();
                ha()    = air.Enthalpy(t_);
                cp_air_ = air.HeatCapacity(t_);
            }
            else
            {
                // No air at this node — zero all air-derived quantities for clean
                // visualization.  sa() was already zeroed in the guard above.
                ra()  = 0.;   mua() = 0.;   ha()  = 0.;
                xa_mix() = 0.; xv_mix() = 1.;
                p_partial_air()   = 0.;
                p_partial_vapor() = p_;

                //pure vapor properies for mix:
                rho_gas_mix()  = Vapor.rho;
                mu_gas_mix()   = Vapor.mu;
                h_gas_mix()    = Vapor.h;
                beta_gas_mix() = Vapor.beta;
            }
        }
    }
    else
    {
        sl() = Liquid.s;
        sv() = Vapor.s;
        sh() = Salt.s;
    }

    // ── Saturation floor: zero numerically empty phases ───────────────────
    // Prevents near-zero saturations from causing division overflow in
    // mixture rescaling (s_gas/sv, s_gas/sa) and mobility computation.
    // A phase below this threshold is physically absent — zeroing it is
    // correct and avoids inf/nan propagating into transport variables.
    if (sl() < sat_min) sl() = 0.;
    if (sv() < sat_min) sv() = 0.;
    if (sh() < sat_min) sh() = 0.;
    if (with_air && sa() < sat_min) sa() = 0.;

    // ── Nodal heat capacity (fluid + air, sa() now correct) ───────────────────
    if(add_fluid_contribution_to_heat_capacity)
    {
        cpf() = Bulk.cp;
        if(with_air)
        {
            ncp() += cpf()   * Bulk.rho * phi() * (1.0 - sa());
            ncp() += cp_air_ * ra()     * phi() * sa();
        }
        else
        {
            ncp() += cpf() * Bulk.rho * phi();
        }
    }

    // ── Thermodynamic equilibrium density ─────────────────────────────────────
    // Uses ra()*sa_clamped (not ma()) — see ma() vs ra()*sa() rule above.
    rho_bulk() = with_air ? Bulk.rho * (1. - sa()) + ra() * sa()
                          : Bulk.rho;

    // ── Pure-halite edge cases ────────────────────────────────────────────────
    {
        x_ = Weight2XNaCl(wt_);
        if(x_ > 1.) x_ = 1.;

        if( (abs(sh() - 1.0)) < 1.0e-10 )   // sh≈1: avoid NaN from Equilibrator
        {
            cerr << "\n[UpdateCSMPVariables] Pure halite (sh≈1).";
            cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
            cerr << endl << "  state=" << Bulk.state
                 << "  sl=" << sl() << "  sv=" << sv() << "  sh=" << sh()
                 << " -> resetting sl=sv=0 sh=1" << endl;

            Vapor.InitToZero();  Liquid.InitToZero();
            Salt.state = H;      Salt.s = 1.0;      Salt.mf = 1.0;
            Bulk = Salt;
            sl() = Liquid.s;     sv() = Vapor.s;    sh() = Salt.s;
            if(with_air) sa() = 0.;

            pure_halite = true;
        }

        if( t() > 799.99 && x_ > 0.999 && Bulk.state == V )   // pure halite in VH field
        {
            cerr << "\n[UpdateCSMPVariables] Pure halite in VAPOR field.";
            cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
            cerr << endl << "  t=" << t() << " oC  x=" << x_
                 << "  state=" << Bulk.state
                 << " -> resetting sl=sv=0 sh=1" << endl;

            Vapor.InitToZero();  Liquid.InitToZero();
            Salt.state = H;      Salt.s = 1.0;      Salt.mf = 1.0;
            Bulk = Salt;
            sl() = Liquid.s;     sv() = Vapor.s;    sh() = Salt.s;

            pure_halite = true;
        }
    }

    // ── Phase densities, viscosities, enthalpies ──────────────────────────────
    rl() = Liquid.rho;    rv() = Vapor.rho;    rh() = Salt.rho;
    mul() = Liquid.mu;    muv() = Vapor.mu;    // overridden below in mixture mode
    hf() = Bulk.h;        hl() = Liquid.h;     hv() = Vapor.h;    hh() = Salt.h;

    // In mixture mode, vapor and air share one gas-phase Darcy velocity, so they
    // must use the same effective viscosity.  muv() is stored to the mesh so that
    // muvp() = mu_gas_mix_prev is available in BoundaryFlow next timestep.
    if (with_air && full_gas_mixture) muv() = mu_gas_mix();

    // In hydraulic mixture mode, compute all gas mixture properties via the
    // shared ComputeGasMixtureProperties function. Both phases are at p_total
    // (no Dalton partial pressures), so p_partial_air = p_partial_vapor = p_.
    // ra_out is ignored here — ra() is already correct at p_total from
    // ComputeAirSaturationAndDensity called above; sa() must not change.
    if (with_air && hydraulic_gas_mixture)  //FIX 30-03-26
    {
        const double s_gas = sv() + sa();
        if (s_gas > 0.)
        {
            double ra_ignored;  //FIX 30-03-26
            ComputeGasMixtureProperties(
                        t_, p_,
                        sv(), sa(),
                        Vapor.rho, Vapor.mu, Vapor.beta, Vapor.h,
                        p_, p_,   // both at p_total — no Dalton partial pressures
                        xa_mix(), xv_mix(),
                        ra_ignored,
                        mu_gas_mix(), rho_gas_mix(), h_gas_mix(), beta_gas_mix());

            p_partial_vapor() = p_;
            p_partial_air()   = p_;

            // Both phases share mu_gas_mix — update phase viscosities for
            // BoundaryFlow (reads muvp/muap as mu_gas_mix_prev next timestep).
            if (sa() > 0.) mua() = mu_gas_mix(); else mua() = 0.;
            if (sv() > 0.) muv() = mu_gas_mix(); else muv() = 0.;

            ha()    = air.Enthalpy(t_);
            cp_air_ = air.HeatCapacity(t_);
        }
        else
        {
            // No gas phase — zero all mixture diagnostics
            xv_mix() = 0.;   xa_mix() = 0.;
            p_partial_vapor() = 0.;   p_partial_air() = 0.;
            rho_gas_mix() = 0.;   mu_gas_mix() = 0.;
            h_gas_mix()   = 0.;   beta_gas_mix() = 0.;
            ha() = 0.;
        }
    }

    if(!add_fluid_contribution_to_heat_capacity)
        cpf() = Bulk.cp;   // visualisation only

    // ── Salt fractions ────────────────────────────────────────────────────────
    wt() = Bulk.wt;
    xf() = Bulk.smf;    xl() = Liquid.smf;    xv() = Vapor.smf;    xh() = Salt.smf;

    // Transport compressibility (H2O-NaCl, so including halite; air term added separately in CT)
    beta() = Bulk.beta;

    // // ── DIAGNOSTIC: does Bulk.beta include the solid halite phase? ────────────
    // //   Bulk.beta ≈ Liquid.beta, flat in sh      → FLUID-ONLY → halite→rock exact
    // //   Bulk.beta drifts toward Salt.beta as sh↑  → INCLUDES halite (vol-weighted)
    // if (Bulk.state == LH || Bulk.state == VH || Bulk.state == VLH)
    // {
    //     const double sf    = sl() + sv();                          // fluid fraction
    //     const double blend = (sf + sh() > 0.)
    //         ? (sf*Liquid.beta + sh()*Salt.beta) / (sf + sh()) : Liquid.beta;
    //     cerr << "[beta-check] state=" << Bulk.state
    //          << "  Bulk="   << Bulk.beta
    //          << "  Liquid=" << Liquid.beta
    //          << "  Salt="   << Salt.beta
    //          << "  blend="  << blend
    //          << "  sh="     << sh() << "\n";
    // }

    // ── Mass concentrations [kg/m³ pore space] ────────────────────────────────
    // ml and mv are re-derived from the equilibrated density × saturation.
    // ma() is NOT touched: it is the conserved transported air mass and is the
    // ground truth for air.  ra() and sa() are derived from it, not the other
    // way around.  Re-syncing ma() from ra()*sa() would silently destroy the
    // excess air mass that was clamped out of sa(), breaking the pressure source
    // term (see nQ in VolumeFactorComputations).
    ml() = rl() * sl();
    mv() = rv() * sv();

    // ── Volumetric and content enthalpies ─────────────────────────────────────
    // hV [J/m³ phase volume]  = h * rho         → flux density for mobility
    // hC [J/m³ total pore]    = hV * S          → accumulation variable
    // Air exception: hCa = ha * ma()  (not hVa*sa) to preserve energy of clamped excess air

    hVl() = hl() * rl();    hCl() = hVl() * sl();
    hVv() = hv() * rv();    hCv() = hVv() * sv();
    hVa() = ha() * ra();    hCa() = ha()  * ma();   // hCa uses ma(), not hVa*sa()
    hVh() = hh() * rh();    hCh() = hVh() * sh();

    xVl() = xl() * rl();    xCl() = xVl() * sl();
    xVv() = xv() * rv();    xCv() = xVv() * sv();
    xVh() = xh() * rh();    xCh() = xVh() * sh();
    xCf() = Bulk.rho * xf();

    // ── Sanity checks ─────────────────────────────────────────────────────────
    if(Bulk.state != 8)
    {
        if(mul() < 0.0) mul() = fabs(mul());
        if(muv() < 0.0) muv() = fabs(muv());

        const char* bogus_reason = nullptr;
        if      (mul() < 1.e-10 && muv() < 1.e-10)           bogus_reason = "viscosities near-zero";
        else if (sl() == 0.0 && sv() == 0.0 && sh() == 0.0)  bogus_reason = "all saturations zero";
        else if (sl() < 0.0 || sv() < 0.0 || sh() < 0.0)     bogus_reason = "negative saturation";
        else if (ml() == 0.0 && mv() == 0.0 && mh() == 0.0)  bogus_reason = "all masses zero";

        if (bogus_reason)
        {
            bogus_variables = true;
            cerr << endl <<"[UpdateCSMPVariables] Bogus EOS output (" << bogus_reason << ")";
            cerr << endl <<"Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
            DumpNodeState("UpdateCSMPVariables/BogusEOS");
        }
    }

    if( definitelyGreaterThan(sl()+sv()+sa(), 1.))
    {
        cerr << endl <<"[UpdateCSMPVariables] Saturation closure violated";
        cerr << endl <<"Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        cerr << endl <<"sl+sv+sa=" << (sl()+sv()+sa()) << " > 1-sh=" << (1.-sh());
        DumpNodeState("UpdateCSMPVariables/SaturationClosure");
    }

    // ── DUST/CONSISTENCY RULE 3 of 3 — EMPTY-NODE OVERRIDE ──────────────────
    // (Rule 1: SBL beat landings. Rule 2: the mt-increment snap above.)
    // mt is the CONSERVED transported total; ml/mv above are re-derived from
    // the equilibrated state (rl*sl, rv*sv) — the intentional desync that
    // feeds nQ. But when mt is dust-empty, the state (e.g. sv=1 at a fully
    // drained node) is stale bookkeeping, and the re-derivation RESURRECTS
    // phase mass from nothing: observed node restarting every stage at
    // mv = rv*sv = 276.52191 exactly after draining to ~30 — a perpetual
    // phantom mass fountain whose exports corrupt neighbours, whose mt
    // increment goes hugely negative ([mt<0] -247), and whose inconsistent
    // energy row eventually hands the thermal equilibrator an unbracketable
    // enthalpy (the tmin>tmax fatal). If the conserved total says EMPTY, the
    // node IS empty: no phase mass, no fluid energy/salt content to
    // accumulate or advect. (Placed after the sanity checks so the deliberate
    // zeroing is not reported as bogus EOS output. hV/xV flux densities stay —
    // intensive EOS quantities; with zero upwind mass nothing is exported.)
    constexpr double MT_EMPTY_TOL = 1.e-9;   // per-PV; transport dry convention
    if ( fabs(mt()) < MT_EMPTY_TOL )
    {
        ml()  = 0.;   mv()  = 0.;
        hCl() = 0.;   hCv() = 0.;
        xCl() = 0.;   xCv() = 0.;
    }

    if(beta() < 0.)
    {
        if(Bulk.state == L  || Bulk.state == F  || Bulk.state == LH)  beta() = Liquid.beta;
        if(Bulk.state == V  || Bulk.state == VH)                       beta() = Vapor.beta;
        if(Bulk.state == VL || Bulk.state == VLH)
            beta() = TwoPhasePureWaterCompressibility(Liquid.cp, Vapor.cp);
        if(Bulk.state == H || Bulk.state == 9)                         beta() = Salt.beta;
        cerr << "\n[UpdateCSMPVariables] Negative beta reset to " << beta()
             << "  state=" << Bulk.state;
        cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
    }
    if(muv() < 0.)
    {
        cerr << "\n[UpdateCSMPVariables] Negative muv=" << muv() << " reset to 1e-5"
             << "  state=" << Bulk.state;
        cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        muv() = 1.0e-5;
    }

}// end UpdateCSMPVariables

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::VolumeFactorComputations( Node<dim>* n )
{
    // ============================================================================
    // PURPOSE
    //   The H2O-NaCl EOS (SOWAT) returns an equilibrium density Bulk.rho that may
    //   differ from the transported mass mt() due to numerical discretisation.
    //   This function reconciles the two by computing correction factors (volume
    //   factors) and applying them to all transport variables before they are
    //   stored to the mesh.
    //
    // VARIABLE TYPES
    //   LHS (accumulation): ml, mv, hCl, hCv, xCl, xCv
    //     — storage terms in the time derivative; represent what is IN the cell.
    //   RHS (flux density): rl_transport, rv_transport, hVl, hVv, xVl, xVv
    //     — integrated over cell faces by CSMP; represent what FLOWS between cells.
    //
    // SCALING PIPELINE (applied in order):
    //
    //   1. H2O-NaCl volume factors (vf_LHS, vf_RHS)
    //      Correct for the SOWAT density mismatch.
    //      vf_RHS = mvl / (ml+mv)        applied to RHS variables
    //      vf_LHS = m_transp / m_eq      applied to LHS variables
    //      vf_LHS ≠ vf_RHS because halite is immobile: it enters m_transp (LHS
    //      numerator) but not mvl (RHS numerator).
    //      Air excluded: ma is conserved, ra = air.Density(T,p) — no SOWAT mismatch.
    //
    //   2. Salt volume factor (salt_vf)
    //      Correct for the drift between the two independent salt tracking paths:
    //        ms      : updated via ms = msp + d(xCl) + d(xCv)
    //        xCl/xCv : independently transported by CSMP, then overwritten by EOS
    //      salt_vf = ms / (xCl+xCv+xCh) — forces xCl+xCv+xCh = ms exactly,
    //      preventing ms < 0 at the next step.
    //      Same factor applied to both LHS (xCl, xCv) and RHS (xVl, xVv):
    //      no LHS/RHS split needed because xVh = 0 (halite has no salt flux).
    //      Applied AFTER vf scaling so it corrects the already-vf-corrected values.
    //
    //   3. Mixture mode rescaling (s_gas/sv, s_gas/sa)
    //      In mixture mode, vapor and air share one Darcy velocity u_gas, so each
    //      component's mass flux = individual density × u_gas.
    //      PrepareVariablesForStorage splits kr_gas by volume fraction:
    //        krv = kr_gas * sv/s_gas,  kra = kr_gas * sa/s_gas
    //      This gives mmv = krv * rv_transport / mu = kr_gas*(sv/s_gas)*rv/mu — wrong.
    //      Rescaling by s_gas/sv absorbs the volume-fraction split back out:
    //        rv_transport *= s_gas/sv  → mmv = kr_gas * rv / mu  ✓
    //        ra_transport *= s_gas/sa  → mma = kr_gas * ra / mu  ✓
    //      The same factor applies to hVv, xVv (vapor) and hVa (air) for consistency:
    //        emv = krv * hVv / mu  and  xmv = krv * xVv / mu get the same correction.
    //      Total mobility mmv+mma = kr_gas*(rv+ra)/mu is unchanged.
    //      Gravity uses rv()+ra() from EOS directly in PrepareVariablesForStorage
    //      — rv_transport/ra_transport do not carry gravity information.
    //      Applied LAST, on top of vf and salt_vf corrections.
    // ============================================================================

    // Initialise transport densities from EOS values — corrected below.
    rl_transport() = rl();
    rv_transport() = rv();
    if(with_air) ra_transport() = ra();

    // ── Equilibrium bulk density ──────────────────────────────────────────────
    // Air via ra()*sa() (clamped volume), NOT ma() (conserved mass).
    // When clamped: ra()*sa < ma() → rho_bulk < mt() → nQ > 0 (overpressure).
    const double rho_bulk_with_air = with_air ? Bulk.rho * (1.0 - sa()) + sa() * ra()
                                              : Bulk.rho;

    // ── Mass quantities for vf_LHS and vf_RHS ────────────────────────────────
    // mh       : immobile halite              [kg/m³ pore]
    // mvl      : mobile H2O-NaCl (no halite, no air) — denominator of vf_RHS
    // m_transp : transported H2O-NaCl+halite (no air) — numerator of vf_LHS
    // m_eq     : equilibrated H2O-NaCl+halite (no air) — denominator of vf_LHS
    double m_eq = 0., m_transp = 0., mvl = 0.;

    if(with_air)
    {
        mh()     = Salt.mf * (mt() - ma());
        mvl      = mt() - mh() - ma();
        m_eq     = Bulk.rho * (1. - sa());
        m_transp = mt() - ma();
    }
    else
    {
        mh()     = Salt.mf * mt();
        mvl      = mt() - mh();
        m_eq     = Bulk.rho;
        m_transp = mt();
    }

    // ── Volume factors ────────────────────────────────────────────────────────
    volume_factor_RHS = (ml() + mv() > 0.) ? mvl      / (ml() + mv()) : 0.;
    volume_factor_LHS = (m_eq > 0.)        ? m_transp / m_eq          : 1.;

    vol_fac() = volume_factor_LHS;   // diagnostic output

    if(mvl < 1.0e-10)
    {
        cerr << "\n[VolumeFactorComputations] mvl -> 0";
        cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        cerr << endl << "  mvl=" << mvl << "  mt=" << mt() << "  mh=" << mh()
             << "  rl=" << rl() << "  rv=" << rv() << "  sl=" << sl() << "  sv=" << sv()
             << "  phi=" << phi() << "  state=" << Bulk.state << endl;
    }

    if(volume_factor_RHS < 0.0)
    {
        cerr << endl <<"[VolumeFactorComputations] Negative volume_factor_RHS: " << volume_factor_RHS;
        cerr << endl <<"Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        DumpNodeState("VolumeFactorComputations/NegativeVF");
        csmp_error.Note(FATAL_ERROR, "NaClH2OPropertiesVisitorPHX::VolumeFactorComputations()",
                        "\nNegative Volume Factor, terminating!");
    }

    // ── nQ: mass mismatch driving the pressure equation ───────────────────────
    // nQ = mt() - rho_bulk_with_air
    //   Normal (unclamped):     nQ ≈ 0 — no pressure signal from air.
    //   Air excess (clamped):   nQ > 0 — overpressure expels excess air next step.
    //   H2O-NaCl phase change:  nQ ≠ 0 — normal pressure signal from density change.
    nQ() = mt() - rho_bulk_with_air;


    // // ── nQ magnitude smoothing (lagged-flash oscillation control) ─────────────
    // // tanh on the RELATIVE mismatch: ~no scaling when small (tanh→0 ⇒ s→1),
    // // saturating reduction when large (tanh→1 ⇒ s→1-cap). Phase- and
    // // location-agnostic by design — discrimination is by mismatch magnitude.
    // //   s = 1 - cap · tanh( |nQ| / (knee · rho_bulk_with_air) )
    // // STOPGAP: scaling under-corrects a real volume error, so it also damps
    // // legitimate large-nQ events (e.g. matrix boiling) — validate fault- AND
    // // matrix-boiling pressure vs HT on/off. FaultPressureIteration is the
    // // converged (zero-residual) version of this.
    {
        constexpr bool   enable_nQ_smoothing = false;
        constexpr double nQ_smooth_cap       = 0.5;    // retain >= 50% of nQ at saturation
        constexpr double nQ_smooth_knee      = 0.02;   // 2% relative mismatch = tanh knee

        if (enable_nQ_smoothing && rho_bulk_with_air > 0. && nQ() != 0.)
        {
            const double rel = std::fabs(nQ()) / (nQ_smooth_knee * rho_bulk_with_air);
            //const double s   = 1.0 - nQ_smooth_cap * std::tanh(rel);

            double s   = 1.;
            if(p()>2.e5) {s   = 0.2;}

            // if (/*verbose_output && */s < 0.9)
            //     cerr << "\n[nQ-smooth] nQ=" << nQ()
            //          << "  rel_mismatch=" << (rel * nQ_smooth_knee)
            //          << "  scale=" << s;
            nQ() *= s;
        }
    }
    // ── expected_dp: the pressure change this node's nQ will actually produce ──
    //
    // History — the original predictor was
    //     expected_dp = (mt() - Bulk.rho) / (mt() * Bulk.beta);   // Δρ / (mt·β)
    // mass-mismatch ÷ (transported density × fluid compressibility). Two defects:
    //
    //   (1) Density was mt() (transported) instead of the EOS equilibrium density.
    //       dρ_eq/dp must be the EOS slope where the fluid actually sits (Bulk.rho),
    //       not at the disequilibrium mass; at off-equilibrium nodes mt() ≠ Bulk.rho
    //       so the slope is misjudged. (Second-order near equilibrium.)
    //
    //   (2) Bulk.beta is FLUID compressibility ONLY — it omits the rock-frame
    //       storage beta_rock·(1-phi) and (after the halite fix) the halite volume
    //       on the rock term. The realised Δp is set by the FULL storage; at stiff
    //       single-phase liquid nodes beta is tiny and the rock term dominates it,
    //       so the fluid-only denominator was a factor ~(1 + beta_rock·(1-phi)/
    //       (beta·phi)) too small — i.e. expected_dp was over-predicted, badly so
    //       for nearly incompressible liquid.
    //
    // Fix — use the SAME CT that becomes the implicit pressure diagonal. The lumped
    // pressure equation, source-only, balances
    //     CT·V_bulk·Δp = nQ·cpv   (cpv = phi·V_bulk)  ⇒  Δp = nQ·phi/CT = nQ/(CT/phi),
    // so rho_beta := CT/phi is the correct effective storage density. It carries the
    // rock frame, the halite-on-rock volume, the per-phase split, and the air cross
    // term — everything the diagonal carries — so the predictor equals the Δp the
    // solve produces by construction, and the floor regulator lands exactly on
    // p_floor instead of overshooting.

    if(p.Flag() == ANY && p.Flag() != DIRICH) //Benoit: for some reason p.Flag()=PLAIN does not seem to work...
    {
        const double CT_local = ComputeTotalCompressibility();
        const double rho_beta = (phi() > 0.) ? CT_local / phi() : 0.;

        // VLH saturation-line oscillation guard (unchanged behaviour).
        if(p() > 38.0e6 && p() < 40.e6 &&
                t() > 590.   && t() < 605.  && wt() > 0.0 &&
                avoid_pressure_oscillations_at_VLH)
        {
            const double factor = std::fabs(p() - 39.e6) / 1.0e6;
            if (nQ() > 0. && nQ() > factor*factor*mt())
                nQ() = factor*factor*mt();
        }

        expected_dp = (rho_beta > 0.) ? nQ() / rho_beta : 0.;

        // Pre-solve floor regulator: rewrite nQ so the predicted pressure lands on
        // p_floor instead of going sub-floor. Same rho_beta as the diagonal, so the
        // predicted landing is exact.
        constexpr bool   avoid_negative_pressures = true;
        constexpr double p_floor = 1.0e4;

        if (avoid_negative_pressures && rho_beta > 0. && p() + expected_dp < p_floor)
        {
            nQ()        = (p_floor - p()) * rho_beta;
            expected_dp = p_floor - p();
        }
    }

    nQ() *= cpv();   // [kg/m³ pore] → [kg]. Full pore volume, NO *(1-sh).

    // ── STEP 1: H2O-NaCl volume factor scaling ────────────────────────────────
    // Air variables (ra_transport, hVa) not scaled — self-consistent, no SOWAT mismatch.
    if(Bulk.rho > 0. && mt() > 0.)
    {
        // LHS accumulation
        ml()  *= volume_factor_LHS;    mv()  *= volume_factor_LHS;
        hCl() *= volume_factor_LHS;    hCv() *= volume_factor_LHS;
        xCl() *= volume_factor_LHS;    xCv() *= volume_factor_LHS;

        // RHS flux densities
        rl_transport() *= volume_factor_RHS;    rv_transport() *= volume_factor_RHS;
        hVl()          *= volume_factor_RHS;    hVv()          *= volume_factor_RHS;
        xVl()          *= volume_factor_RHS;    xVv()          *= volume_factor_RHS;
    }

    // ── STEP 2: Salt volume factor ────────────────────────────────────────────
    // salt_vf = ms / (xCl+xCv+xCh): forces xCl+xCv+xCh = ms, preventing drift.
    // Same factor for LHS (xCl, xCv) and RHS (xVl, xVv) — no split needed (xVh=0).
    // xCh excluded: immobile, not a CSMP transport variable.
    {
        const double xC_eos = xCl() + xCv() + xCh();
        if (xC_eos > 0. && ms() >= 0.)
        {
            const double salt_vf = ms() / xC_eos;
            xCl() *= salt_vf;    xCv() *= salt_vf;    // LHS accumulation
            xVl() *= salt_vf;    xVv() *= salt_vf;    // RHS flux density
        }
    }

    // ── STEP 3: Full Mixture mode rescaling (Dalton Mode ONLY) ────────────────
    // Corrects the volume-fraction kr split so each phase gets density-fraction
    // mass flux. All RHS flux variables for vapor (rv_transport, hVv, xVv) and
    // air (ra_transport, hVa) get the same factor s_gas/sv or s_gas/sa.
    //
    // IMPORTANT: This is ONLY applied in full_gas_mixture. In hydraulic_gas_mixture,
    // rv and ra are pure phase densities (evaluated at p_total), so applying this
    // inverse-volume-fraction scaling would cause mobilities to blow up to infinity
    // as a phase disappears.
    //
    // Applied on top of the vf_RHS and salt_vf corrections already in place.
    if(with_air && full_gas_mixture && (sv() + sa()) > 0.)
    {
        const double s_gas = sv() + sa();

        // Vapor: rv_transport, hVv, xVv all scale by s_gas/sv
        if (sv() > sat_min) { rv_transport() *= s_gas/sv(); hVv() *= s_gas/sv(); xVv() *= s_gas/sv(); }
        else                { rv_transport() = 0.; hVv() = 0.; xVv() = 0.; }

        // Air: ra_transport, hVa scale by s_gas/sa (no salt in air — no xVa)
        if (sa() > sat_min) { ra_transport() *= s_gas/sa(); hVa() *= s_gas/sa(); }
        else                { ra_transport() = 0.; hVa() = 0.; }
    }

    // ── Diagnostic: DIRICH interior nodes ────────────────────────────────────
    // Both factors should be 1.0 when p and T are pinned (no mass mismatch).
    if(p.Flag() == DIRICH && t.Flag() == DIRICH && n->AtBoundary() == NOT)
    {
        if(!essentiallyEqual(volume_factor_RHS, 1., 1.e-6))
        {
            cerr << "\n[VolumeFactorComputations] DIRICH interior vf_RHS != 1: " << volume_factor_RHS;
            cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        }
        if(!essentiallyEqual(volume_factor_LHS, 1., 1.e-6))
        {
            cerr << "\n[VolumeFactorComputations] DIRICH interior vf_LHS != 1: " << volume_factor_LHS;
            cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        }
    }
}// end VolumeFactorComputations

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::BoundaryIteration()
{
    // ============================================================================
    // FUNCTION: BoundaryIteration
    // ROLE: The Numerical "Strategist"
    // PURPOSE: Finds the exact mass flow (bfm) such that BoundaryFlow() returns
    //          a volume ratio of exactly 1.0 — i.e. the fluid mass precisely fills
    //          the available pore space.
    //
    // HIGH-LEVEL LOGIC:
    //   BoundaryFlow() is treated as a black-box physical evaluator:
    //     for any trial bfm, it applies phase partitioning, updates state,
    //     and returns vol = mt() / rho_bulk_with_air.
    //
    //   BoundaryIteration() solves the scalar root-finding problem:
    //       f(bfm) = BoundaryFlow(bfm) - 1.0 = 0
    //
    // NUMERICAL STRATEGY:
    //   1. Save original conserved state so trials can be safely undone.
    //   2. Evaluate bfm = 0 to determine flow direction (outflow vs inflow).
    //      If air is present and ma > ra*sa_max, BoundaryFlow removes the excess
    //      at this step. The backup is then updated to the cleaned state so all
    //      subsequent trials start from a physically consistent baseline.
    //   3. Generate an initial guess from the density mismatch.
    //   4. Bracket the root so vol crosses 1.0 within (bfm_min, bfm_max).
    //   5. Refine with a secant/bisection hybrid until converged.
    //   6. Post-convergence air cap: clamp ma <= ma_max to prevent the next
    //      step from starting with a physically impossible air saturation.
    //      Needed because the inflow bracketing loop deliberately skips excess
    //      air removal (to allow vol > 1.0), so the converged inflow solution
    //      can accumulate ma >> ma_max.
    //
    // IMPORTANT: BoundaryFlow() is stateful — it modifies mt(), ms(), H_current_,
    //   ma() (when with_air), and mv() (when full_gas_mixture). The saved state
    //   MUST be restored before every trial evaluation via restore_backup.
    //
    // AIR EXCESS REMOVAL DESIGN:
    //   BoundaryFlow removes excess air (ma > ra*sa_max) only when bfm >= 0
    //   (outflow or zero-flow direction check). This keeps removal out of the
    //   inflow bracketing loop, which doubles bfm negatively and would otherwise
    //   strip inflow air and prevent vol from ever reaching > 1.0.
    //
    //   In practice, removal only fires at the bfm=0 trial: outflow trials start
    //   from ma_backup=ma_max (set by the backup refresh below) and can only
    //   reduce ma(), so they structurally cannot re-introduce excess. The wider
    //   bfm >= 0 guard is kept as a defensive catch against future code changes.
    //
    //   Pre-existing excess is handled by the backup refresh after the bfm=0 trial:
    //   excess is removed once, and all subsequent trials restore to the cleaned
    //   baseline where ma <= ra*sa_max.
    //
    // NOTE: outflow mobilities are based on previous-step saturations (slp, svp, sap).
    //   A phase that just arrived this timestep (map=0, dma>0) has zero previous-step
    //   saturation and therefore zero outflow mobility — it cannot leave until the
    //   next timestep.
    // ============================================================================

    auto restore_backup = [&](double mt_backup, double ms_backup,
            double H_backup,  double ma_backup, double mv_backup)
    {
        mt()       = mt_backup;
        ms()       = ms_backup;
        H_current_ = H_backup;
        if (with_air)                     ma() = ma_backup;
        if (with_air && full_gas_mixture) mv() = mv_backup;
    };

    auto max_outflow = [&](double trial_bfm, double mt_backup, double ma_backup)
    {
        // With air: ensure at least a minimal H2O-NaCl mass remains for the
        // Equilibrator. This is a numerical floor only — physical protection of
        // H2O-NaCl pore fraction is handled by sa_max in ComputeAirSaturationAndDensity.
        if (with_air)
        {
            const double h2o_nacl  = mt_backup - ma_backup;
            const double protected_h2o_nacl = std::max(h2o_nacl * 1.e-6, 1.e-12);
            return std::min(trial_bfm, mt_backup - protected_h2o_nacl);
        }
        return std::min(trial_bfm, mt_backup * (1. - 1.e-6));
    };

    // eval_trial restores the backup, sets bfm, calls BoundaryFlow, and returns vol.
    // It leaves the system in the post-BoundaryFlow state. The root-finding loop
    // relies on this: the final eval_trial call at convergence (or failure) leaves
    // the accepted state in place, so no separate "accept" step is needed.
    auto eval_trial = [&](double mt_backup, double ms_backup,
            double H_backup,  double ma_backup, double mv_backup,
            double trial_bfm)
    {
        restore_backup(mt_backup, ms_backup, H_backup, ma_backup, mv_backup);
        bfm() = max_outflow(trial_bfm, mt_backup, ma_backup);
        const double trial_vol = BoundaryFlow();
        boundary_flow_trial++;
        return trial_vol;
    };

    // ── State variables ───────────────────────────────────────────────────────
    double mt_backup (0.), ms_backup (0.), H_backup (0.), ma_backup (0.), mv_backup (0.);
    double bfm_min (0.), bfm_max (0.);
    double vol (0.), vol_min (0.), vol_max (0.);
    double crit (1.);
    bool outflow(false);

    // Save initial conserved state
    mt_backup = mt();
    ms_backup = ms();
    H_backup  = H_current_;
    ma_backup = with_air ? ma() : 0.;
    mv_backup = mv();

    // Update inflow reference enthalpies from nodal values set by
    // InitialPropertiesFromPTX. These are per-node and reflect the local
    // equilibrated T/P at each boundary node.
    ref_spec_h_liquid = ref_enthalpy_liquid_top();
    if (with_air)
    {
        ref_spec_h_air          = ref_enthalpy_air_top();
        ref_spec_h_vapor_in_air = ref_enthalpy_vapor_in_air_top();
        ref_frac_vapor_in_air   = ref_frac_vapor_in_air_top();
    }

    // ════════════════════════════════════════════════════════════════════════
    // RAIN RECHARGE — UNCONDITIONAL, ONCE PER TIMESTEP
    // ════════════════════════════════════════════════════════════════════════
    //
    // Rain is a fixed atmospheric mass source. It is NOT conditioned on the
    // sign of bfm: rain falls on the mountain whether the underlying cell
    // is exfiltrating or accepting inflow. Previously the rain logic lived
    // inside BoundaryFlow's inflow branch (bfm < 0) and was clamped by bfm,
    // with two pathological consequences for heterogeneous drainage:
    //
    //   (a) Top cells in net outflow (bfm >= 0, common on wet flanks during
    //       drainage) received zero rain — the branch was skipped entirely.
    //   (b) Top cells with |bfm| < rain*dt had their rain truncated to bfm,
    //       silently discarding the excess. Mass was lost.
    //
    // Both effects produced a precipitation field that depended on the local
    // boundary flow direction rather than on the weather, producing the
    // fingered/patchy drainage structure typical of this bug.
    //
    // The fix is to deposit the full rain mass into the conserved state
    // ONCE, before the root-finder starts. rain_recharge_top is a per-node
    // volumetric source [kg/m³-pore/s] set by SetPrecipitationRecharge,
    // which has already applied the node-area, horizontal-projection, and
    // control-pore-volume normalisations, so `rain_recharge_top * dt` is
    // the correct mass increment per unit pore volume for this timestep.
    //
    // We modify the BACKUP, not the live state, because eval_trial restores
    // the backup on every trial. Depositing into the backup guarantees that
    // every trial sees the post-rain starting state, and the iteration is
    // free to converge on any sign of bfm:
    //
    //   - Cell initially undersaturated even after rain: bfm < 0 (humid air
    //     sucked in to fill the remainder). Handled by the inflow branch in
    //     BoundaryFlow, which — after this fix — handles only humid air.
    //
    //   - Cell exactly filled by rain:                    bfm ≈ 0.
    //
    //   - Cell oversaturated by rain (runoff):            bfm > 0 (excess
    //     expelled by the normal outflow branch, which is physically what
    //     surface runoff is). Nothing is lost; mass balance is exact.
    //
    // Rain is fresh liquid at local atmospheric conditions:
    //   - Salt content:  zero (ms_backup unchanged)
    //   - Air content:   zero (ma_backup unchanged) — the air already in the
    //                    pore space is compressed as mt grows, which is the
    //                    physically correct response.
    //   - Enthalpy:      rain_mass * ref_spec_h_liquid * phi. See Issue 5
    //                    note: ref_spec_h_liquid is the reference-state
    //                    liquid enthalpy; if the atmospheric temperature
    //                    differs substantially from the reference, replace
    //                    with an explicit rain temperature later.
    //
    // ════════════════════════════════════════════════════════════════════════
    if (with_air && rain_recharge_top() > 0. && dt_ > 0.)
    {
        const double rain_mass = rain_recharge_top() * dt_;   // [kg/m³-pore]

        // Deposit into the BACKUP so every restore_backup call sees post-rain
        // state. The live mt()/H_current_ are overwritten by the first
        // restore_backup inside eval_trial, so writing to them here would
        // be lost — the backup is the source of truth for the root-finder.
        mt_backup += rain_mass;
        H_backup  += rain_mass * ref_spec_h_liquid * phi();

        // Also update the live state so any pre-iteration diagnostics
        // (e.g. a debug print before the first eval_trial) see the correct
        // post-rain state. Not strictly necessary for correctness.
        mt()       = mt_backup;
        H_current_ = H_backup;
    }

    boundary_flow_trial = 0;

    // ====================================================================
    // STEP 1: DIRECTION CHECK  (bfm = 0)
    //
    // Evaluate at bfm = 0 to see whether the system is over- or under-filled.
    //   vol > 1.0  =>  too much mass  =>  OUTFLOW required
    //   vol < 1.0  =>  too little mass =>  INFLOW required
    //
    // If air is present and ma > ra*sa_max, BoundaryFlow removes the excess
    // here (bfm=0 satisfies the bfm>=0 guard). We then refresh the backup to
    // the post-removal state so all subsequent trials restore to the physically
    // clean baseline. Without this refresh, pre-existing excess (e.g. map >>
    // ma_max from a previous timestep) would fire removal on every inflow trial,
    // preventing vol from ever reaching > 1.0 and causing the inflow doubling
    // loop to overflow.
    // ====================================================================
    bfm() = 0;
    vol = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());

    // Refresh backup to post-removal clean state (no-op if no excess was present).
    if (with_air)
    {
        mt_backup = mt();
        ms_backup = ms();
        H_backup  = H_current_;
        ma_backup = ma();
        mv_backup = mv();
    }

    // Early exit: node is already at equilibrium with zero flow.REMOVED
    // State and Equilibrator output are valid from the bfm=0 trial above.
    // if (essentiallyEqual(vol, 1., 1.e-4))
    //     return;
    // if (essentiallyEqual(mt(), rho_bulk(), 1.e-4))//stale rho bulk
    //     return;

    if (definitelyGreaterThan(vol, 1., numeric_limits<double>::epsilon()))
        outflow = true;

    if (outflow) { vol_min = vol; bfm_min = bfm(); }
    else         { vol_max = vol; bfm_max = bfm(); }

    // ====================================================================
    // STEP 2: INITIAL GUESS
    // Estimate bfm from the density mismatch; enforce correct sign.
    // ====================================================================
    bfm() = (mt() - rho_bulk()) / 2.;

    if ( outflow && definitelyLessThan   (bfm(), 0., numeric_limits<double>::epsilon()))
        bfm() *= -1.;
    if (!outflow && definitelyGreaterThan(bfm(), 0., numeric_limits<double>::epsilon()))
        bfm() *= -1.;

    const double tiny = 1.e-12 * mt();   // minimal flow to start bracketing
    // If the guess is near zero, seed a tiny flow.
    if (essentiallyEqual(bfm(), 0., 1.e-12))
    {
        bfm() = outflow ? tiny : -tiny;
    }

    vol = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());

    // ====================================================================
    // STEP 3: BRACKETING
    //
    // Expand bfm magnitude until we have a valid bracket:
    //     vol_min > 1.0   and   vol_max < 1.0
    // State is restored before every BoundaryFlow() call via eval_trial.
    // ====================================================================
    if (outflow)
    {
        if (definitelyLessThan(vol, 1., numeric_limits<double>::epsilon()))
        {
            // Initial guess already crossed the target — bracket complete.
            vol_max = vol;
            bfm_max = bfm();
        }
        else
        {
            // Double bfm until vol falls below 1.0.
            bool find_vol_max(true);
            int count_vol_max(0);

            while (find_vol_max)
            {
                count_vol_max++;
                bfm() *= 2.;

                vol = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());

                if (definitelyLessThan(vol, 1., numeric_limits<double>::epsilon()) ||
                        count_vol_max > 100)
                {
                    if (count_vol_max > 100)
                    {
                        cerr << "\n[BoundaryIteration] Outflow bracket: vol<1 not reached after 100 doublings"
                             << "  t=" << t_ << " oC  p=" << p_ << " Pa"
                             << "  bfm=" << bfm() << "  vol=" << vol << endl;
                        PAUSE();
                    }

                    vol_max = vol;
                    bfm_max = bfm();
                    find_vol_max = false;
                }
            }
        }
    }
    else // inflow
    {
        if (definitelyGreaterThan(vol, 1., numeric_limits<double>::epsilon()))
        {
            vol_min = vol;
            bfm_min = bfm();
        }
        else
        {
            // Double bfm (negatively) until vol rises above 1.0.
            // Excess air removal is deliberately skipped inside BoundaryFlow
            // during inflow trials (bfm < 0), allowing vol to exceed 1.0.
            bool find_vol_min(true);
            int count_vol_min(0);

            while (find_vol_min)
            {
                count_vol_min++;
                bfm() *= 2.;

                vol = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());

                if (definitelyGreaterThan(vol, 1., numeric_limits<double>::epsilon()) ||
                        count_vol_min > 100)
                {
                    if (count_vol_min > 100)
                    {
                        cerr << "\n[BoundaryIteration] Inflow bracket: vol>1 not reached after 100 doublings"
                             << "  t=" << t_ << " oC  p=" << p_ << " Pa"
                             << "  bfm=" << bfm() << "  vol=" << vol << endl;
                        PAUSE();
                    }

                    vol_min = vol;
                    bfm_min = bfm();
                    find_vol_min = false;
                }
            }
        }
    }

    // Bracket validation
    if (definitelyGreaterThan(bfm_min, bfm_max, numeric_limits<double>::epsilon()) ||
            definitelyLessThan   (vol_min, 1., numeric_limits<double>::epsilon()) ||
            definitelyGreaterThan(vol_max, 1., numeric_limits<double>::epsilon()))
    {
        cerr << "\n[BoundaryIteration] Invalid bracket!"
             << "  outflow=" << outflow
             << "  t=" << t_ << " oC  p=" << p_ << " Pa  wt=" << wt_ << " wt%"
             << "  mt_backup=" << mt_backup << "  H_backup=" << H_backup
             << "  m_fluid=" << m_fluid_ << "  mv_backup=" << mv_backup;
        if (with_air) cerr << "  ma_backup=" << ma_backup << "  sa=" << sa();
        cerr << "\n  bfm=[" << bfm_min << ", " << bfm_max
             << "] span=" << (bfm_max - bfm_min)
             << "  vol=[" << vol_max << ", " << vol_min << "]" << endl;
        DumpNodeState("BoundaryIteration/InvalidBracket");
    }

    // ====================================================================
    // STEP 4: ROOT-FINDING  (secant / bisection hybrid)
    //
    // Every 3rd iteration forces bisection to guarantee monotone bracket
    // shrinkage. All other iterations use the secant method for speed.
    // If the secant prediction escapes the bracket, fall back to bisection.
    // Convergence criterion: |1 - vol| < 1e-4.
    //
    // eval_trial leaves the system in the post-BoundaryFlow state, so the
    // final call in the failure path also serves as the accept step.
    // ====================================================================
    crit = std::abs(1. - vol);
    int count(0);

    while (crit > 1.e-4)
    {
        count++;

        if (count % 3 == 0)
        {
            // Bisection step — guaranteed bracket reduction
            bfm() = 0.5 * (bfm_min + bfm_max);
        }
        else
        {
            // Secant step — fast convergence near the root
            bfm()  = bfm_min * (1. - vol_max);
            bfm() += bfm_max * (vol_min - 1.);
            bfm() /= (vol_min - vol_max);

            // Safety: if secant escapes the bracket, fall back to bisection
            if (bfm() <= std::min(bfm_min, bfm_max) ||
                    bfm() >= std::max(bfm_min, bfm_max))
                bfm() = 0.5 * (bfm_min + bfm_max);
        }

        if (count > 100)
        {
            cerr << "\n[BoundaryIteration] Failed to converge after 100 iterations."
                 << "  outflow=" << outflow
                 << "  t=" << t_ << " oC  p=" << p_ << " Pa" << endl;
            cerr << "  bracket: bfm=[" << bfm_min << ", " << bfm_max
                 << "] span=" << (bfm_max - bfm_min)
                 << "  vol=[" << vol_max << ", " << vol_min
                 << "]  crit=" << crit << endl;
            cerr << "  state backup: mt=" << mt_backup << "  ms=" << ms_backup
                 << "  H=" << H_backup << "  mv=" << mv_backup;
            if (with_air) cerr << "  ma=" << ma_backup;
            cerr << endl;
            // dma_ is the original pre-BoundaryIteration increment from
            // CalculateAbsoluteVariables — never modified here, so these
            // diagnostics correctly reflect the initial timestep state.
            cerr << "  phase inputs: mlp=" << mlp() << "  mvp=" << mvp();
            if (with_air) cerr << "  map=" << map();
            cerr << "  dml=" << dml_ << "  dmv=" << dmv_;
            if (with_air) cerr << "  dma=" << dma_;
            cerr << "\n  old_l=" << (mlp() + std::min(0., dml_))
                 << "  old_v=" << (mvp() + std::min(0., dmv_));
            if (with_air) cerr << "  old_a=" << (map() + std::min(0., dma_));
            cerr << "  new_l=" << std::max(0., dml_)
                 << "  new_v=" << std::max(0., dmv_);
            if (with_air) cerr << "  new_a=" << std::max(0., dma_);
            cerr << endl;

            // Fresh evals at bracket endpoints for diagnostics.
            // The final eval_trial also serves as the accept step.
            bfm() = bfm_min;
            double v_min = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());
            cerr << "  [bfm_min=" << bfm_min << "] vol=" << v_min
                 << "  state=" << Bulk.state << "  Bulk.rho=" << Bulk.rho;
            if (with_air) cerr << "  ra=" << ra() << "  sa=" << sa()
                               << "  sa_unc=" << (ra() > 0. ? ma() / ra() : 0.)
                               << "  rho_with_air=" << ((1. - sa()) * Bulk.rho + sa() * ra());
            cerr << "  Liquid.s=" << Liquid.s << "  Vapor.s=" << Vapor.s
                 << "  mt=" << mt() << "  ms=" << ms() << "  wt=" << wt_
                 << "  H=" << H_current_ << "  phi_=" << phi_ << endl;

            bfm() = bfm_max;
            double v_max = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());
            cerr << "  [bfm_max=" << bfm_max << "] vol=" << v_max
                 << "  state=" << Bulk.state << "  Bulk.rho=" << Bulk.rho;
            if (with_air) cerr << "  ra=" << ra() << "  sa=" << sa()
                               << "  sa_unc=" << (ra() > 0. ? ma() / ra() : 0.)
                               << "  rho_with_air=" << ((1. - sa()) * Bulk.rho + sa() * ra());
            cerr << "  Liquid.s=" << Liquid.s << "  Vapor.s=" << Vapor.s
                 << "  mt=" << mt() << "  ms=" << ms() << "  wt=" << wt_
                 << "  H=" << H_current_ << "  phi_=" << phi_ << endl;

            // Accept the bracket endpoint closest to vol = 1.
            if (std::abs(v_min - 1.) <= std::abs(v_max - 1.))
                eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm_min);
            else
                eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm_max);

            crit = 0.;
            break;
        }

        vol  = eval_trial(mt_backup, ms_backup, H_backup, ma_backup, mv_backup, bfm());
        crit = std::abs(1. - vol);

        // Tighten bracket
        if (vol > 1.) { bfm_min = bfm(); vol_min = vol; }
        if (vol < 1.) { bfm_max = bfm(); vol_max = vol; }
    }

    // ====================================================================
    // STEP 5: POST-CONVERGENCE AIR CAP
    //
    // After convergence (or failure exit), clamp ma() to ra*sa_max.
    //
    // The inflow bracketing loop (Step 3) doubles bfm negatively and
    // deliberately skips excess air removal to allow vol > 1. The converged
    // inflow solution can therefore store ma >> ma_max. If unchecked, the
    // next timestep's bfm=0 trial cannot bracket vol < 1 (dilute steam at
    // low pressure has very low density, so vol stays > 1 even after draining
    // all H2O-NaCl mass).
    //
    // ENERGY CONSISTENCY: excess air enthalpy is charged to bfe() and removed
    // from H_current_, matching the convention in BoundaryFlow's in-trial
    // excess removal block.
    //
    // NOTE: m_fluid_ = (mt-ma)*phi is invariant under this cap — both mt and
    // ma decrease by the same excess. The last Equilibrator output (Bulk.rho,
    // phase saturations, temperature) therefore remains valid; no re-equilibration
    // is needed.
    // ====================================================================
    if (with_air && ra() > 0.)
    {
        const double sa_max_physical = std::max(0., 1.0 - sh());
        const double sa_max          = std::min(sa_max_physical, 1.0 - min_h2o_nacl_fraction);
        const double ma_max          = ra() * sa_max;

        if (ma() > ma_max)
        {
            const double ma_excess = ma() - ma_max;
            const double H_excess = ma_excess * phi() * air.Enthalpy(t_);

            H_current_ -= H_excess;
            bfe()      += H_excess;
            bfa()      += ma_excess;
            mt()       -= ma_excess;
            ma()        = ma_max;

            cerr << "\n[BoundaryIteration] Post-convergence ma cap:"
                 << "  ma=" << (ma_max + ma_excess) << " -> " << ma_max
                 << "  excess=" << ma_excess << endl;
        }
    }

} // end BoundaryIteration


template<uint32_t dim>
double NaClH2OPropertiesVisitorPHX<dim>::BoundaryFlow()
{
    // ============================================================================
    // FUNCTION: BoundaryFlow
    // ROLE: The Physical "Accountant"
    // PURPOSE: Calculates the consequences of a specific trial mass flow (bfm).
    //          It partitions mobile mass into phases, computes energy and salt
    //          transported across the boundary, updates system state, and returns
    //          the resulting volume ratio (mt / rho_bulk_with_air) for the
    //          root-finder in BoundaryIteration.
    //
    // HIGH-LEVEL LOGIC:
    //   1. Partition mobile mass into old/new liquid, vapor, and air.
    //   2. Based on sign of bfm(), handle outflow or inflow.
    //   3. For outflow, compute phase fractional flow using relative mobilities.
    //   4. Compute energy (bfe), salt (bfs), and air mass (bfa) transported.
    //   5. Update system totals: mt, H_current_, ms, ma (and mv if full_gas_mixture).
    //   6. Remove excess air if bfm >= 0. Structurally only fires at bfm=0 —
    //      see AIR EXCESS REMOVAL note below.
    //   7. Re-equilibrate thermodynamic state via ThreePhaseProperties.
    //   8. Return volume ratio for BoundaryIteration's numerical controller.
    //
    // CONSERVED-STATE TRACKING:
    //   BoundaryIteration restores {mt, ms, H_current_, ma, mv} before each call.
    //   This function updates all five consistently:
    //     mt, ms, H_current_ : always (outflow and inflow).
    //     ma                 : always when with_air.
    //     mv                 : only when full_gas_mixture. Tracked so that the
    //                          per-trial partial pressure p*xa can be computed from
    //                          the current gas-phase composition. See air density block.
    //
    // AIR EXCESS REMOVAL (bfm >= 0 only):
    //   After state updates, if ma > ra*sa_max the excess is removed directly.
    //   In practice this only fires at the bfm=0 direction-check trial: outflow
    //   trials start from ma_backup=ma_max (set by BoundaryIteration's backup
    //   refresh) and can only decrease ma(), so excess cannot reappear. The wider
    //   bfm >= 0 guard is kept as a defensive catch. During inflow trials (bfm < 0),
    //   removal is deliberately skipped — stripping inflow air would prevent vol
    //   from exceeding 1.0 and break the inflow bracketing loop.
    //   All four conserved quantities updated consistently on removal:
    //   ma, mt, bfa, H_current_.
    //
    // AIR PHASE PARTITIONING (dma_local):
    //   old_air and new_air are computed from dma_local = ma() - map(), not from
    //   the member dma_. This makes the split immune to the bfm=0 excess removal
    //   shifting the effective baseline: after BoundaryIteration's backup refresh,
    //   ma_backup = ma_max, so dma_local = ma_max - map() is always consistent
    //   with the restored ma(). The member dma_ (original pre-BoundaryIteration
    //   increment from CalculateAbsoluteVariables) is preserved unmodified and
    //   used only in diagnostics.
    //
    // NEW-MASS OUTFLOW ENERGY:
    //   When outflow demand exceeds old-phase inventory, BoundaryFlow draws from
    //   newly-arrived mass (dml_ > 0, dmv_ > 0, dma_ > 0). Their specific energy
    //   uses the transport-derived value (dhCl_/new_liquid, etc.) rather than the
    //   local equilibrated Liquid.h. Although dhCl_/new_liquid may underestimate
    //   enthalpy at hot nodes, it is a compile-time constant during root-finding
    //   (fixed in CalculateAbsoluteVariables, never changes between trials). This
    //   keeps bfe() exactly linear in bfm, ensuring smooth and well-conditioned
    //   convergence. Using Liquid.h (updated each trial by the Equilibrator) would
    //   introduce non-linear coupling between bfe() and bfm and risk oscillations.
    //   Salt fractions (dxCl_/new_liquid) are composition, not enthalpy — kept as-is.
    //
    // AIR DENSITY IN FULL GAS MIXTURE MODE:
    //   ComputeAirSaturationAndDensity is called with a per-trial partial pressure
    //   p*xa_trial derived from {mv() (tracked per-trial), ma() (trial state)}.
    //   This is more accurate than using p_total when the gas composition shifts
    //   significantly across trials. p_current_ (the vapor EOS pressure for the
    //   Equilibrator) is NOT updated per trial — changing it would make Bulk.rho
    //   jump non-smoothly between trials and break root-finder monotonicity. The
    //   residual pressure non-closure (p_partial_air_trial + p_partial_vapor_fixed
    //   != p_total) is second-order in Δxa and smaller than the O(xa) error from
    //   naively using p_total. Pass 2 in UpdateCSMPVariables refines both partial
    //   pressures after full convergence.
    //
    // Called repeatedly by BoundaryIteration() during root-finding.
    // ============================================================================

    double bfm_l(0.), bfm_v(0.), bfm_a(0.), bfm_h(0.);
    double old_liquid(0.), old_vapor(0.), old_air(0.);
    double new_liquid(0.), new_vapor(0.), new_air(0.);

    // =========================================================================
    // 1. PHASE PARTITIONING
    //
    // Split each mobile phase into:
    //   old = pre-existing inventory (previous-step mass + any net decrease)
    //   new = mass that arrived this timestep (positive increment only)
    //
    // For air, dma_local = ma() - map() is recomputed from the restored ma()
    // rather than from the member dma_. See AIR PHASE PARTITIONING note above.
    // Solid halite is tracked for the sanity check only; it does not flow.
    // =========================================================================
    old_liquid = mlp() + std::min(0., dml_);
    old_vapor  = mvp() + std::min(0., dmv_);
    new_liquid = std::max(0., dml_);
    new_vapor  = std::max(0., dmv_);

    if (with_air)
    {
        const double dma_local = ma() - map();
        old_air = map() + std::min(0., dma_local);
        new_air = std::max(0., dma_local);
    }

    // ── Rain recharge: classify as new liquid ────────────────────────────────
    // Rain was deposited into mt_backup (and mt()) once at the start of
    // BoundaryIteration, BEFORE the root-finder started. That deposit did not
    // touch dml_ / mlp() — both are properties of the transport step, set by
    // CalculateAbsoluteVariables before rain is applied. The old_/new_
    // partitioning above therefore sums to the pre-rain mt, and the sanity
    // check below fires with diff = rain_mass every trial.
    //
    // Rain IS new liquid arriving this timestep — physically indistinguishable
    // from a transported liquid increment as far as the outflow priority
    // scheme is concerned. Classifying it as new_liquid:
    //   (a) restores the mass-balance sanity check (sum matches mt());
    //   (b) ensures that on runoff (bfm > 0), the outflow branch at line
    //       ~2175 expels rain before pre-existing liquid — correct upwind-
    //       in-time behaviour, since rain is the most recently arrived mass.
    //
    // rain_recharge_top() and dt_ are the same values used in
    // BoundaryIteration to populate mt_backup, so the classification here is
    // guaranteed consistent with the deposit there. BoundaryFlow is only
    // called from BoundaryIteration (single call site), so there is no risk
    // of double-counting from another caller.
    if (with_air && rain_recharge_top() > 0. && dt_ > 0.)
    {
        const double rain_mass_this_step = rain_recharge_top() * dt_;
        new_liquid += rain_mass_this_step;
    }

    // // Sanity check: partitioned mobile mass + solid halite must equal mt().
    // {
    //     double old_salt = 0.;
    //     if (with_air)
    //     {
    //         if (definitelyGreaterThan(mtp(), mlp() + mvp() + map(),
    //                                   numeric_limits<double>::epsilon()))
    //             old_salt = mtp() - mlp() - mvp() - map();
    //     }
    //     else
    //     {
    //         if (definitelyGreaterThan(mtp(), mlp() + mvp(),
    //                                   numeric_limits<double>::epsilon()))
    //             old_salt = mtp() - mlp() - mvp();
    //     }

    //     double total = old_salt + old_liquid + old_vapor + new_liquid + new_vapor;
    //     if (with_air) total += old_air + new_air;

    //     if (fabs(mt() - total) > 1.0e-10 * mt())
    //     {
    //         cerr << "\n[BoundaryFlow] Mass balance mismatch! trial=" << boundary_flow_trial
    //              << "  t=" << t_ << " oC  p=" << p_ << " Pa" << endl;
    //         cerr << "  old_salt=" << old_salt
    //              << "  old_l=" << old_liquid << "  old_v=" << old_vapor;
    //         if (with_air) cerr << "  old_a=" << old_air;
    //         cerr << "  new_l=" << new_liquid << "  new_v=" << new_vapor;
    //         if (with_air) cerr << "  new_a=" << new_air;
    //         cerr << "  sum=" << total << "  mt=" << mt()
    //              << "  diff=" << (mt() - total)
    //              << "  rel=" << (mt() > 0. ? (mt() - total) / mt() : 0.) << endl;
    //         PAUSE();
    //     }
    // }

    bfe() = 0.;
    bfs() = 0.;
    if (with_air) bfa() = 0.;

    // =========================================================================
    // CASE 1: OUTFLOW  (bfm > 0, system is not fully solid halite)
    // =========================================================================
    if (definitelyGreaterThan(bfm(), 0., numeric_limits<double>::epsilon()) &&
            definitelyLessThan(sh(), 1., numeric_limits<double>::epsilon()))
    {
        // ── Relative permeabilities and fractional mobilities ─────────────────
        // Vapor and air form a combined non-wetting gas phase for kr purposes.
        // Total gas kr is split between vapor and air by previous-step volume
        // fraction. All quantities use previous-step values (slp, svp, sap,
        // rlp, rvp, rap, mulp, muvp, muap) — current-step equilibrium is unknown.
        double krl = with_air
                ? RelativePermeabilityLiquid(slp(), svp(), sap())
                : RelativePermeabilityLiquid(slp(), svp(), 0.);

        double mob_l = 0., mob_v = 0., mob_a = 0.;

        if (with_air)
        {
            const double gas_space = svp() + sap();
            if (definitelyGreaterThan(gas_space, 0., numeric_limits<double>::epsilon()))
            {
                // Total gas kr via linear or Brooks-Corey. Previous-step
                // saturations (slp, svp, sap) used consistently with krl above.
                // The split between vapor and air by volume fraction reflects
                // pore-network sharing and is independent of the kr model.
                const double krg = RelativePermeabilityGas(slp(), svp(), sap());
                const double krv = (svp() / gas_space) * krg;
                const double kra = (sap() / gas_space) * krg;

                if (full_gas_mixture || hydraulic_gas_mixture)
                {
                    // ── COUPLED GAS PHASE ─────────────────────────────────────────
                    // Vapor and air share one Darcy velocity → one effective viscosity.
                    //
                    // FULL GAS MIXTURE (Dalton): rvp() and rap() are partial densities;
                    //   volume-fraction information is already embedded in them.
                    //   Using krg_total avoids double-penalising by volume fraction.
                    //
                    // HYDRAULIC MIXTURE: rvp() and rap() are pure densities at p_total;
                    //   volume-fraction allocation must come from kr, so fractional
                    //   krv and kra are used.
                    const double krg_total   = krv + kra;
                    const double mu_gas_prev = (svp() > 0.) ? muvp() : muap();

                    if (full_gas_mixture)
                    {
                        if (svp() > 0.) mob_v = krg_total * rvp() / mu_gas_prev;
                        if (sap() > 0.) mob_a = krg_total * rap() / mu_gas_prev;
                    }
                    else // hydraulic_gas_mixture
                    {
                        if (svp() > 0.) mob_v = krv * rvp() / mu_gas_prev;
                        if (sap() > 0.) mob_a = kra * rap() / mu_gas_prev;
                    }
                }
                else
                {
                    // ── IMMISCIBLE MODE ───────────────────────────────────────────
                    if (svp() > 0.) mob_v = krv * rvp() / muvp();
                    if (sap() > 0.) mob_a = kra * rap() / muap();
                }
            }
            if (slp() > 0.) mob_l = krl * rlp() / mulp();
        }
        else
        {

            // No air: vapor is the only gas phase, so kr_gas = kr_vapor.
            // Call RelativePermeabilityGas with sa=0 for a consistent curve.
            const double krv = std::max(0., RelativePermeabilityGas(slp(), svp(), 0.));
            if (slp() > 0.) mob_l = krl * rlp() / mulp();
            if (svp() > 0.) mob_v = krv * rvp() / muvp();

        }

        // Normalise to fractional mobilities
        const double mob = mob_l + mob_v + (with_air ? mob_a : 0.);
        if (mob > 0.)
        {
            mob_l /= mob;
            mob_v /= mob;
            if (with_air) mob_a /= mob;
        }

        // ── Outflow allocation: two-stage (old mass, then new mass) ───────────
        //
        // Forward cascade (Air → Vapor → Liquid) redistributes unmet demand to
        // the next phase when CASCADE=true. Each phase is capped by its inventory;
        // cascade never creates mass or allows overdraw.
        const bool CASCADE_OLD_MASS = true;
        const bool CASCADE_NEW_MASS = true;

        const double old_mobile = old_liquid + old_vapor + (with_air ? old_air : 0.);

        if (!definitelyGreaterThan(bfm(), old_mobile, numeric_limits<double>::epsilon()))
        {
            bfm_l = bfm() * mob_l;
            bfm_v = bfm() * mob_v;
            if (with_air) bfm_a = bfm() * mob_a;

            if (CASCADE_OLD_MASS)
            {
                if (with_air &&
                        definitelyGreaterThan(bfm_a, old_air, numeric_limits<double>::epsilon()))
                {
                    double excess = bfm_a - old_air;
                    bfm_a = old_air;
                    bfm_v += excess;
                }
                if (definitelyGreaterThan(bfm_v, old_vapor, numeric_limits<double>::epsilon()))
                {
                    double excess = bfm_v - old_vapor;
                    bfm_v = old_vapor;
                    bfm_l += excess;
                }
                if (definitelyGreaterThan(bfm_l, old_liquid, numeric_limits<double>::epsilon()))
                {
                    double excess = bfm_l - old_liquid;
                    bfm_l = old_liquid;
                    bfm_h = std::min(excess, mh());  // halite is terminal sink at boundary
                }
            }
            else
            {
                if (with_air) bfm_a = std::min(bfm_a, old_air);
                bfm_v = std::min(bfm_v, old_vapor);
                bfm_l = std::min(bfm_l, old_liquid);
                // Remaining unmet demand goes to halite regardless of cascade mode —
                // halite draining is independent of fluid-phase redistribution.
                const double remaining = bfm() - bfm_l - bfm_v - (with_air ? bfm_a : 0.);
                bfm_h = std::min(mh(), remaining);

            }
        }
        else
        {
            // bfm exceeds all available old mobile mass: drain everything
            bfm_l = old_liquid;
            bfm_v = old_vapor;
            if (with_air) bfm_a = old_air;
            bfm_h = std::min(mh(), bfm() - old_liquid - old_vapor - old_air);
        }

        // Energy and salt carried out by old phase mass.
        // Uses previous-step specific enthalpies and salt fractions.
        bfe() = phi() * (bfm_l * hlp() + bfm_v * hvp() + bfm_h * hhp());
        if (with_air) bfe() += phi() * bfm_a * hap();
        bfs() = bfm_l * xlp() + bfm_v * xvp() + bfm_h;

        // ── Additional outflow from new mass ──────────────────────────────────
        const double flowed_so_far = bfm_l + bfm_v + bfm_h + (with_air ? bfm_a : 0.);
        const double extra_flow    = bfm() - flowed_so_far;

        if (definitelyGreaterThan(extra_flow, 0., numeric_limits<double>::epsilon()))
        {
            double extra_l = extra_flow * mob_l;
            double extra_v = extra_flow * mob_v;
            double extra_a = with_air ? extra_flow * mob_a : 0.;

            if (CASCADE_NEW_MASS)
            {
                // Forward: Air → Vapor → Liquid
                if (with_air &&
                        definitelyGreaterThan(extra_a, new_air, numeric_limits<double>::epsilon()))
                {
                    double excess = extra_a - new_air;
                    extra_a = new_air;
                    extra_v += excess;
                }
                if (definitelyGreaterThan(extra_v, new_vapor, numeric_limits<double>::epsilon()))
                {
                    double excess = extra_v - new_vapor;
                    extra_v = new_vapor;
                    extra_l += excess;
                }
                // Reverse: Liquid → Vapor → Air
                if (definitelyGreaterThan(extra_l, new_liquid, numeric_limits<double>::epsilon()))
                {
                    double excess = extra_l - new_liquid;
                    extra_l = new_liquid;
                    extra_v += excess;
                }
                if (definitelyGreaterThan(extra_v, new_vapor, numeric_limits<double>::epsilon()))
                {
                    double excess = extra_v - new_vapor;
                    extra_v = new_vapor;
                    if (with_air) extra_a += excess;
                }
                // Final cap on air
                if (with_air && definitelyGreaterThan(extra_a, new_air, numeric_limits<double>::epsilon()))
                    extra_a = new_air;
            }
            else
            {
                if (with_air) extra_a = std::min(extra_a, new_air);
                extra_v = std::min(extra_v, new_vapor);
                extra_l = std::min(extra_l, new_liquid);
            }

            bfm_l += extra_l;
            bfm_v += extra_v;
            if (with_air) bfm_a += extra_a;

            // Energy and salt from new mass.
            // Uses transport-derived specific enthalpy (dhCx_/new_x). Although this
            // may underestimate enthalpy at hot nodes, it is a compile-time constant
            // during root-finding — keeping bfe() linear in bfm for smooth convergence.
            // See NEW-MASS OUTFLOW ENERGY note in function header.
            if (new_liquid > 0. && extra_l > 0.)
            {
                bfe() += extra_l * (dhCl_ / new_liquid);
                bfs() += extra_l * (dxCl_ / new_liquid);
            }
            if (new_vapor > 0. && extra_v > 0.)
            {
                bfe() += extra_v * (dhCv_ / new_vapor);
                bfs() += extra_v * (dxCv_ / new_vapor);
            }
            if (with_air && new_air > 0. && extra_a > 0.)
                bfe() += extra_a * (dhCa_ / new_air);
        }

        // Sync bfm() to actual total flow (phase caps may have reduced it).
        // mt() -= bfm() below must match what was actually allocated.
        bfm() = bfm_l + bfm_v + bfm_h + (with_air ? bfm_a : 0.);
        if (with_air) bfa() = bfm_a;

        // Track mv for per-trial partial pressure in full gas mixture mode.
        // See AIR DENSITY IN FULL GAS MIXTURE MODE in function header.
        if (with_air && full_gas_mixture)
            mv() = std::max(0., mv() - bfm_v);
    }

    // =========================================================================
    // CASE 2: INFLOW  (bfm < 0)
    // =========================================================================
    //
    // Rain has already been deposited into mt_backup by BoundaryIteration
    // BEFORE the root-finder started. This branch now handles ONLY the
    // residual gas inflow required to fill whatever pore space the rain
    // did not already fill. That residual comes in as humid air — a blend
    // of vapor and dry air in ratio ref_frac_vapor_in_air, with no salt
    // and at atmospheric reference enthalpies.
    //
    // If the cell was already full (or overfilled by rain alone), the
    // root-finder lands on bfm >= 0 and this branch is never entered —
    // the outflow branch handles the runoff case naturally.
    //
    if (definitelyLessThan(bfm(), 0., numeric_limits<double>::epsilon()))
    {
        if (with_air)
        {
            // All inflow is humid air. No rain here anymore.
            const double gas_flow     = bfm();                                // < 0
            const double vap_flow     = gas_flow * ref_frac_vapor_in_air;
            const double dry_air_flow = gas_flow * (1.0 - ref_frac_vapor_in_air);

            bfa() = bfm_a = dry_air_flow;
            bfs() = 0.;                                                       // humid air is salt-free
            bfe() = phi() * (vap_flow     * ref_spec_h_vapor_in_air
                             + dry_air_flow * ref_spec_h_air);

            // Diagnostic: converged dry-air fraction of the gas inflow.
            // Now equals (1 - ref_frac_vapor_in_air) exactly when bfm < 0,
            // since rain is no longer part of the inflow mix. Kept for
            // continuity of output format.
            ref_frac_air_top() = (bfm() != 0.) ? bfa() / bfm() : 0.;

            // Track mv for per-trial partial pressure in full gas mixture mode.
            if (full_gas_mixture) mv() += vap_flow;
        }
        else
        {
            // No-air case unchanged: inflow is pure liquid at ref composition.
            bfe() = bfm() * ref_spec_h_liquid * phi();
            bfs() = bfm() * ref_sal;
        }
    }

    // =========================================================================
    // UPDATE SYSTEM STATE
    // Apply boundary flow to all conserved quantities.
    // =========================================================================
    mt()       -= bfm();
    H_current_ -= bfe();

    // Salt: guard against numerical overshoot to negative values
    if (bfs() < ms())
        ms() -= bfs();
    else
    {
        bfs() = ms();
        ms()  = 0.;
    }

    if (with_air)
    {
        // Guard against numerical overshoot to negative values
        if (bfa() < ma())
            ma() -= bfa();
        else
        {
            bfa() = ma();
            ma()  = 0.;
        }
    }

    // ── Direct excess air removal (bfm >= 0 only) ────────────────────────────
    // Remove air mass exceeding ra*sa_max. In practice only fires at the bfm=0
    // direction-check trial — see AIR EXCESS REMOVAL note in function header.
    // During inflow trials (bfm < 0), removal is deliberately skipped.
    // ra() is constant throughout BoundaryIteration (T and p are Dirichlet).
    if (with_air && bfm() >= 0. && ra() > 0.)
    {
        const double sa_max_physical = std::max(0., 1.0 - sh());
        const double sa_max          = std::min(sa_max_physical, 1.0 - min_h2o_nacl_fraction);
        const double ma_max          = ra() * sa_max;

        if (ma() > ma_max)
        {
            const double ma_excess = ma() - ma_max;
            const double H_excess  = ma_excess * phi() * air.Enthalpy(t_);

            cerr << "\n[BoundaryFlow] Direct excess air removal"
                 << "  trial=" << boundary_flow_trial
                 << "  ma=" << ma() << "  ma_max=" << ma_max
                 << "  excess=" << ma_excess
                 << "  sa_unclamped=" << ma() / ra()
                 << "  sap=" << sap() << "  map=" << map()
                 << "  bfm=" << bfm() << "  mt=" << mt()
                 << "  mt-ma=" << mt() - ma()
                 << "  ra=" << ra()
                 << "  t=" << t_ << "  p=" << p_ << endl;

            H_current_ -= H_excess;
            bfe()      += H_excess;
            bfa()      += ma_excess;
            ma()       -= ma_excess;
            mt()       -= ma_excess;
        }
    }

    // ── Re-prepare Equilibrator inputs ───────────────────────────────────────
    m_fluid_ = with_air
            ? std::max(1.e-12, (mt() - ma())) * phi()
            : std::max(1.e-12, mt()) * phi();

    wt() = (ms() > 0. && mt() > ma())
            ? (with_air ? ms() / (mt() - ma()) * 100. : ms() / mt() * 100.)
            : 0.;
    if (wt() > 100.) wt() = 100.;
    wt_ = wt();

    // Diagnostics
    if (mt() < 0. || H_current_ < 0. || ms() < 0.)
    {
        cerr << "\n[BoundaryFlow] Nonphysical negative state! trial=" << boundary_flow_trial
             << "  t=" << t_ << " oC  p=" << p_ << " Pa"
             << "  H=" << H_current_ << "  mt=" << mt() << "  ms=" << ms();
        if (with_air) cerr << "  ma=" << ma();
        cerr << "  bfm=" << bfm() << "  bfe=" << bfe() << "  bfs=" << bfs();
        if (with_air) cerr << "  bfa=" << bfa();
        cerr << endl;
    }
    if (with_air && (ma() < 0. || ma() > mt() * (1. + 1.e-10)))
    {
        cerr << "\n[BoundaryFlow] Nonphysical air state! trial=" << boundary_flow_trial
             << "  t=" << t_ << " oC  p=" << p_ << " Pa"
             << "  mt=" << mt() << "  ma=" << ma() << "  mt-ma=" << mt() - ma()
             << "  m_fluid=" << m_fluid_ << "  ms=" << ms() << "  H=" << H_current_ << endl;
        cerr << "  bfm=" << bfm() << "  bfe=" << bfe()
             << "  bfs=" << bfs() << "  bfa=" << bfa() << endl;
        cerr << "  sap=" << sap()
             << "  (mob_a=0 when sap=0 => air locked until next step)" << endl;
    }

    // ── Air saturation and porosity shaving ──────────────────────────────────
    // phi_ is recomputed from phi() each call — NOT multiplied onto the incoming
    // value. BoundaryFlow is called repeatedly and phi_ arrives already shaved;
    // using *= would give phi*(1-sa)^N instead of phi*(1-sa).
    //
    // FULL GAS MIXTURE: air density uses per-trial p*xa_trial derived from
    //   {mv() (tracked per-trial), ma() (trial state)}. p_current_ is unchanged.
    //   See AIR DENSITY IN FULL GAS MIXTURE MODE in function header.
    // ALL OTHER MODES: air density uses p_total (exact for immiscible/hydraulic).
    if (with_air)
    {
        if (ma() > 0.)
        {
            if (full_gas_mixture)
            {
                const double nv_t = (mv() > 0.) ? mv() / M_VAPOR : 0.;
                const double na_t = (ma() > 0.) ? ma() / M_AIR   : 0.;
                const double ntot = nv_t + na_t;
                const double xa_t = (ntot > 0.) ? na_t / ntot : 0.;
                ComputeAirSaturationAndDensity(ma(), t_, p_, p_ * xa_t);
            }
            else
            {
                ComputeAirSaturationAndDensity(ma(), t_, p_, p_);
            }
            phi_   = phi() * (1.0 - sa());
            m_air_ = ma() * phi();
        }
        else
        {
            ra()   = 0.;
            sa()   = 0.;
            phi_   = phi();
            m_air_ = 0.;
        }
    }

    // ── Thermodynamic equilibration ───────────────────────────────────────────
    thermal_eq_succeeded = equilibrator.ThreePhaseProperties(Bulk, Liquid, Vapor, Salt);

    if (!thermal_eq_succeeded)
    {
        cerr << "\n[BoundaryFlow] ThreePhaseProperties failed at boundary"
             << "  p=" << p() << " Pa  t=" << t() << " oC  wt=" << wt_ << " wt%"
             << "  mt=" << mt() << "  mtp=" << mtp();
        if (std::abs(p_current_ - p()) > 1.0)
            cerr << "  p_current=" << p_current_ << " Pa";
        cerr << "  H=" << H_current_ << "  m_fluid=" << m_fluid_;
        if (with_air) cerr << "  phi_shaved=" << phi_ << "  ma=" << ma() << "  sa=" << sa();
        cerr << "  rho_bulk=" << rho_bulk() << endl;
    }

    // ── Volume ratio ──────────────────────────────────────────────────────────
    // Returns mt / rho_bulk_with_air. Convergence target is 1.0.
    //
    // Uses CLAMPED sa() = min(ma/ra, sa_max), not unclamped ma/ra. This ensures
    // vol is monotone in bfm, which the root-finder requires:
    //   Clamped:   when ma > ma_max, sa = sa_max (fixed), rho_bulk denominator is
    //              constant, and vol grows linearly with mt. Well-conditioned.
    //   Unclamped: as ma grows, sa = ma/ra increases, (1-sa)*Bulk.rho drops sharply
    //              (since Bulk.rho >> ra), rho_bulk falls, and vol becomes non-monotone.
    // This formula is consistent with rho_bulk_with_air in VolumeFactorComputations
    // and PrepareVariablesForStorage.
    if (with_air)
    {
        const double rho_bulk_with_air = (1.0 - sa()) * Bulk.rho + sa() * ra();
        return (rho_bulk_with_air > 0.) ? mt() / rho_bulk_with_air : 1.0;
    }
    else
    {
        return mt() / Bulk.rho;
    }

} // end BoundaryFlow

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::PrepareVariablesForStorage()
{
    state() = static_cast<double>(Bulk.state);
    Htp()   = H_current_;
    msp()   = ms();
    mtp()   = mt();

    // rho_bulk stored for transport — same formula as VolumeFactorComputations
    rho_bulk() = with_air ? Bulk.rho * (1.0 - sa()) + sa() * ra()
                          : Bulk.rho;

    // Air compressibility — stored as a DIAGNOSTIC only. It is not a CT term; the
    // air storage enters CT through the regime-dependent cross term inside
    // ComputeTotalCompressibility().
    if (with_air)
        beta_air() = (p() > 0.) ? (full_gas_mixture ? beta_gas_mix()
                                                    : air.Compressibility(p()))
                                : 0.;

    // Total compressibility — single source of truth (see ComputeTotalCompressibility).
    CT() = ComputeTotalCompressibility();

    if (CT() <= 0.)
        cerr << "\n[PrepareVariablesForStorage] Non-positive CT=" << CT()
             << "  beta=" << beta() << "  beta_air=" << beta_air()
             << "  phi=" << phi() << "  sl=" << sl() << "  sv=" << sv()
             << "  sa=" << sa() << "  sh=" << sh()
             << "  p=" << p_ << " Pa  state=" << Bulk.state << endl;

    // ── Mobility accumulators — zero before accumulation loop ─────────────────
    mml() = mmv() = mmld() = mmvd() = 0.0;
    eml() = emv() = emld() = emvd() = 0.0;
    rvl() = rvv() = 0.0;
    xml() = xmv() = 0.0;
    if(with_air) { mma() = mmad() = ema() = emad() = rva() = 0.0; }

    // ── Relative permeabilities ───────────────────────────────────────────────
    // Vapor and air are a combined non-wetting gas phase for kr; total gas kr
    // is then split between vapor and air by volume fraction.
    if(with_air)
    {
        krl = RelativePermeabilityLiquid(sl(), sv(), sa());
        // Total gas kr via linear or Brooks-Corey. Split into krv and kra
        // below by volume fraction; kr_gas itself represents the combined
        // vapor+air mobility through the gas pore network.
        const double kr_gas = RelativePermeabilityGas(sl(), sv(), sa());
        const double s_gas  = sv() + sa();

        if(s_gas > 0.)
        {
            // Simple fractional flow — air gets its share of total gas kr
            krv = kr_gas * sv() / s_gas;
            kra = kr_gas * sa() / s_gas;

        }

        else krv = kra = 0.;
    }

    else
    {
        krl = RelativePermeabilityLiquid(sl(), sv(), 0.);
        krv = RelativePermeabilityGas(sl(), sv(), 0.);
    }

    // ── Phase mobilities: mm = kr*rho/mu,  em = kr*hV/mu,  rv = kr/mu ────────
    // In mixture mode, vapor and air share one Darcy velocity u_gas.
    // VolumeFactorComputations has rescaled rv_transport and ra_transport so that
    //   mmv = krv * rv_transport / muv  correctly gives  kr_gas * rv_individual / mu_gas
    //   mma = kra * ra_transport / mua  correctly gives  kr_gas * ra_individual / mu_gas
    // i.e. each phase gets the mass flux proportional to its own density, not to
    // the combined gas density. Total gas mobility mmv+mma = kr_gas*(rv+ra)/mu_gas
    // is preserved. muv == mua == mu_gas_mix in both mixture modes.

    if(sl() > 0.)
    {
        if (mul() <= 0.) cerr << "\n[PrepareVariablesForStorage] mul=" << mul()
                              << " non-positive — liquid mobility undefined  sl=" << sl()
                              << "  state=" << Bulk.state << endl;
        else { mml()=krl*rl_transport()/mul(); eml()=krl*hVl()/mul(); rvl()=krl/mul(); xml()=krl*xVl()/mul(); }
    }
    if(sv() > 0.)
    {
        if (muv() <= 0.) cerr << "\n[PrepareVariablesForStorage] muv=" << muv()
                              << " non-positive — vapor mobility undefined  sv=" << sv()
                              << "  state=" << Bulk.state << endl;
        else { mmv()=krv*rv_transport()/muv(); emv()=krv*hVv()/muv(); rvv()=krv/muv(); xmv()=krv*xVv()/muv(); }
    }
    if(with_air && sa() > 0.)
    {
        if (mua() <= 0.) cerr << "\n[PrepareVariablesForStorage] mua=" << mua()
                              << " non-positive — air mobility undefined  sa=" << sa() << endl;
        else { mma()=kra*ra_transport()/mua(); ema()=kra*hVa()/mua(); rva()=kra/mua(); }
    }

    // Old comments ── Gravity-term densities ────────────────────────────────────────────────
    // The gravity term uses the density that drives buoyancy.
    // In mixture mode, both vapor and air travel at u_gas and experience the
    // combined gas buoyancy rho_gas = rv + ra. We use rv() + ra() from the EOS
    // directly — rv_transport/ra_transport have been rescaled in VolumeFactorComputations
    // for correct pressure-driven flux allocation and no longer equal rho_gas.
    //
    // For full_gas_mixture:    rv() + ra() = total gas density at Dalton partial pressures.
    // For hydraulic_gas_mixture: rv() + ra() = sum at p_total.

    //FIX 01-04-2026 ── Gravity-term densities ────────────────────────────────────────────────
    // The gravity term uses the density that drives buoyancy.
    // In mixture mode, both vapor and air travel at u_gas and experience the
    // combined gas buoyancy rho_gas.

    double rho_gas_grav = 0.0;
    const bool is_mixture = with_air && (full_gas_mixture || hydraulic_gas_mixture);

    if (is_mixture)
    {
        if (full_gas_mixture)
        {
            // Dalton mode: rv() and ra() are partial densities. Sum = true mixture density.
            rho_gas_grav = rv() + ra();
        }
        else if (hydraulic_gas_mixture)
        {
            // Hydraulic mode: rv() and ra() are pure densities.
            // Summing them creates an artificially dense gas. We MUST volume-weight them.
            const double s_gas = sv() + sa();
            rho_gas_grav = (s_gas > 0.) ? (rv() * (sv() / s_gas) + ra() * (sa() / s_gas)) : 0.0;
        }
    }

    const bool thermo = thermodymamic_density_in_gravity_term;
    const double rv_grav = is_mixture ? rho_gas_grav : rv();
    const double ra_grav = is_mixture ? rho_gas_grav : ra();

    if(sl() > 0.0)             { mmld() = mml() * (thermo ? rl()     : rl_transport());
        emld() = eml() * (thermo ? rl()     : rl_transport()); }
    if(sv() > 0.0)             { mmvd() = mmv() * (thermo ? rv_grav  : rv_transport());
        emvd() = emv() * (thermo ? rv_grav  : rv_transport()); }
    if(with_air && sa() > 0.0) { mmad() = mma() * (thermo ? ra_grav  : ra_transport());
        emad() = ema() * (thermo ? ra_grav  : ra_transport()); }

}// end PrepareVariablesForStorage

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetTimeIncrement( double time_increment )
{
    dt_ = time_increment;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeTotalCompressibility
//
// Total compressibility CT [kg·m⁻³·Pa⁻¹, per BULK volume] — the diagonal
// coefficient of the pressure equation
//
//     CT · V_bulk · Δp  =  nQ · cpv  +  (flux and source terms)        (cpv = phi·V_bulk)
//
// SINGLE SOURCE OF TRUTH. Called from:
//   • VolumeFactorComputations   → predictor rho_beta = CT/phi, so the
//                                   expected_dp prediction matches the solve.
//   • PrepareVariablesForStorage → stored as the implicit pressure diagonal.
// Keeping both on this one routine guarantees predictor and solver can't diverge.
//
// STORAGE PARTITION (pore saturations sum to 1: sl + sv + sa + sh = 1):
//   • H2O-NaCl (liquid+vapor+halite): beta() covers all phases that occupy the
//       non‑air pore space (1 - sa). For states with halite (LH, VH, VLH, H) it
//       already includes the halite contribution (dissolution storage or mechanical
//       compressibility). The fallback to Liquid.beta in UpdateCSMPVariables is
//       only a safety net for negative beta values, not the typical case.
//   • Solid halite (phi·sh): a stiff solid — compressed at beta_rock, folded into
//       the rock-frame term in EVERY branch. Compressing it at beta() instead
//       would be ~10× too soft and is the bug this routine fixes.
//   • Rock frame ((1-phi)): beta_rock.
//   • Air: regime-dependent cross term d(rho_a·sa)/dp — see block comment below.
//       beta_air is NOT a standalone CT term (it would double-count); it is only
//       stored as a diagnostic, outside this routine.
//
// The whole expression is scaled by the chosen reference density rho_scale
// (rho_bulk for thermodynamic consistency, else mt). The air cross term is added
// AFTER scaling, written so rho_scale cancels.
// ─────────────────────────────────────────────────────────────────────────────
template<uint32_t dim>
double NaClH2OPropertiesVisitorPHX<dim>::ComputeTotalCompressibility()
{
    const double rho_scale = thermodynamic_density_in_compressibility_term
            ? rho_bulk() : mt();
    double ct = 0.;

    if (with_air && full_gas_mixture)
    {
        // ── FULL GAS MIXTURE (Dalton) ─────────────────────────────────────────
        // Liquid on sl; gas-mixture compressibility on (sv+sa); halite on rock.
        // No Grant & Sorey phase-conversion term here (Liquid.beta is mechanical),
        // but beta_gas_mix >> Liquid.beta so the gas dominates CT anyway. The
        // immiscible air cross term is intentionally absent: in mixture mode
        // dsa/dp also redistributes vapor, so a clean three-phase Grant & Sorey
        // extension would be required. Revisit if mixture runs oscillate.
        ct  = Liquid.beta    * phi() * sl();
        ct += beta_gas_mix() * phi() * (sv() + sa());
        ct += beta_rock()    * ((1.0 - phi()) + phi() * sh());
        ct *= rho_scale;
    }
    else if (with_air)
    {
        // ── IMMISCIBLE / HYDRAULIC GAS MIXTURE ────────────────────────────────
        // beta() (= Bulk.beta) already accounts for halite: for LH/VH/VLH it is the
        // Grant & Sorey dissolution storage built from salt.rho/h/cp/s, and for pure
        // H it is halite.Compressibility(). It governs the full non-air pore fraction
        // (1 - sa) — do NOT carve out phi*sh. The rock term is the mechanical frame
        // only; the thermal rock buffering of the phase change is already inside
        // TwophaseCompressibility() (its (1-phi)·cpr·rr term), so adding beta_rock·
        // (1-phi) here is the separate mechanical compression, not a double count.
        ct  = beta()      * phi() * (1.0 - sa());
        ct += beta_rock() * (1.0 - phi());
        ct *= rho_scale;

        // ── AIR CROSS TERM  d(rho_a·sa)/dp ─────────────────────────────────────
        // Regime-dependent. beta_air is NOT added as a standalone term (double
        // count). For ideal-gas air, beta_air = 1/p.
        //
        //   UNCLAMPED (sa = ma/ra, ma conserved):
        //       rho_a·sa = ma = const  ⇒  air self-storage d(rho_a·sa)/dp = 0.
        //       The real storage is liquid backfilling the pore volume the air
        //       vacates as it compresses → rho_cross = Bulk.rho, term Bulk.rho·sa/p.
        //       Dominant stabilising feedback at air+liquid nodes (Bulk.rho >> ra).
        //
        //   CLAMPED (sa pinned at sa_max, fixed-volume pocket):
        //       only the trapped air compresses → sa·ra/p, i.e. rho_cross = ra
        //       (≡ beta_air·ra·sa). ~800× smaller near-atmospheric; using Bulk.rho
        //       here would over-damp CT and produce fingered drainage fronts.
        //
        // Regime detection reconstructs the unclamped sa = ma/ra and compares to
        // the SAME sa_max used by BoundaryFlow / the post-convergence cap; sat_min
        // guards against per-step flip-flop at ma = ra·sa_max.
        //
        // Added post-scale (rho_scale cancels: the immiscible block above already
        // carries rho_scale on its beta terms).
        const double sa_max_physical = std::max(0., 1.0 - sh());
        const double sa_max_local    = std::min(sa_max_physical,
                                                1.0 - min_h2o_nacl_fraction);
        const double sa_unclamped    = (ra() > 0.) ? ma() / ra() : 0.;
        const bool   air_is_clamped  = (sa_unclamped > sa_max_local + sat_min);
        const double rho_cross       = air_is_clamped ? ra() : Bulk.rho;

        if (p() > 0.)
            ct += rho_cross * sa() * phi() / p();

        // Optional regime diagnostic — enable once to confirm, then silence.
        // if (air_is_clamped && sa() > sat_min)
        //     cerr << "[CT] clamped air: sa_uncl=" << sa_unclamped
        //          << "  sa_max=" << sa_max_local << "  rho_cross=" << rho_cross
        //          << "  (Bulk.rho=" << Bulk.rho << ", ra=" << ra() << ")\n";
    }
    else
    {
        // ── NO AIR ────────────────────────────────────────────────────────────
        ct  = beta()      * phi();
        ct += beta_rock() * (1.0 - phi());
        ct *= rho_scale;
    }

    return ct;
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::ComputeAirSaturationAndDensity(
        double ma_,
        double t_,
        double p_,
        double p_air_partial)
{
    // =========================================================================
    // PURPOSE: convert transported air mass ma_ [kg/m³ pore] into the volume
    //   fraction sa [-] that air actually occupies, subject to a physical cap.
    //
    // INPUTS:
    //   ma_           : transported air mass [kg/m³ pore] — conserved, never modified here
    //   t_            : temperature [°C]
    //   p_            : total pore pressure [Pa] — used only for immiscible mode
    //   p_air_partial : air partial pressure [Pa] — used in mixture mode
    //
    // CALLED FROM TWO PLACES with different pressures:
    //
    //   1. UpdateEquilibratorVariables (PRE-Equilibrator):
    //        ComputeAirSaturation(ma(), t_, p_, p_partial_air())
    //        Pass 1 runs before this call, so p_partial_air() is already set
    //        from the advected masses mv/ma.  In mixture mode the air density is
    //        therefore at the correct partial pressure — no approximation.
    //        In immiscible mode p_partial_air() = p_total, also exact.
    //        sh() is the PREVIOUS-STEP halite saturation (read from mesh at the
    //        start of Visit()), so sa_max is one step stale.
    //
    //   2. UpdateCSMPVariables (POST-Equilibrator, BEFORE saturation rescaling):
    //        ComputeAirSaturation(ma(), t_, p_, p_partial_air())
    //        t_ = Bulk.t is the newly equilibrated temperature; the pre-Equil
    //        sa() used the pre-Equil temperature so must be recomputed here.
    //        p_partial_air() is still the Pass-1 value (Pass 2 hasn't run yet).
    //        sh() is STILL the PREVIOUS-STEP halite saturation — the rescaling
    //        sh() = Salt.s*(1-sa) that would update it happens AFTER this call.
    //        This is the authoritative sa() used for saturation rescaling and all
    //        downstream variables.
    //
    // CLAMPING AND PRESSURE SOURCE TERM:
    //   sa_unclamped = ma_ / ra(t, p_air)               [volume air would need]
    //   sa_max       = min(1-sh, 1-min_h2o_nacl_fraction) [available pore volume]
    //   sa           = min(sa_unclamped, sa_max)         [actual occupied volume]
    //
    //   res_sl is NOT applied here — residual liquid saturation is a mobility
    //   concept (kr_l → 0 at sl=res_sl) and does not restrict the physical volume
    //   that air can occupy. It belongs in RelativePermeabilityLiquid only.
    //
    //   min_h2o_nacl_fraction ensures phi_ = phi*(1-sa) > 0 so the Equilibrator
    //   always receives a well-posed problem. Excess air beyond sa_max stays in
    //   ma() and drives nQ > 0 (overpressure) which expels it via transport.
    //
    //   ma_ is NEVER modified: if sa_unclamped > sa_max the excess air mass is
    //   retained in ma().  This causes:
    //     rho_bulk_with_air = Bulk.rho*(1-sa) + ra()*sa  <  mt()
    //     nQ = mt() - rho_bulk_with_air  >  0
    //   nQ > 0 is a positive pressure source in the pressure equation, which
    //   raises the nodal pressure and expels the excess air mass via transport
    //   in the next timestep.  This is the correct physical mechanism: excess
    //   air creates overpressure until it can escape.
    // =========================================================================

    const double p_air_use = full_gas_mixture ? p_air_partial : p_;

    ra() = air.Density(t_, p_air_use);

    const double sa_unclamped = (ra() > 0.) ? ma_ / ra() : 0.0;

    // sa_max: air can occupy any non-halite pore space.
    // res_sl does NOT appear here — residual liquid saturation affects kr only,
    // not the physical volume air can occupy. sh() is previous-step (one-step lag,
    // acceptable and consistent with all other previous-step usage in the code).
    const double sa_max_physical = std::max(0., 1.0 - sh());

    // Minimum H2O-NaCl pore fraction: always reserve a small volume for H2O-NaCl
    // fluid so that phi_ = phi*(1-sa) > 0 and the Equilibrator always receives a
    // well-posed problem. Without this, when sa → 1, phi_ → 0 and the Equilibrator
    // gets fluid mass squeezed into zero volume → returns bogus sl=sv=sh=0.
    // The excess air beyond sa_limit stays in ma() and drives nQ overpressure
    // which expels it via transport in the next timestep — correct physical mechanism.

    const double sa_max = std::min(sa_max_physical, 1.0 - min_h2o_nacl_fraction); // reserve 2% pore for H2O-NaCl

    const double sa_clamped = std::min(sa_unclamped, sa_max);

    if (verbose_output && sa_unclamped > sa_max + 1.e-10)
        cerr << "\n[ComputeAirSaturationAndDensity] sa clamped:"
             << "  sa_unclamped=" << sa_unclamped << " -> sa=" << sa_clamped
             << "  sa_max=" << sa_max
             << "  excess_ma=" << (sa_unclamped - sa_max) * ra() << " kg/m3"
             << "  t=" << t_ << " oC  p_air=" << p_air_use << " Pa" << endl;

    sa() = sa_clamped;

}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::ComputeGasMixtureProperties(
        double t,               // equilibrated temperature [°C]
        double p,               // total pore pressure [Pa]
        double sv_in,           // vapor volume fraction [-]
        double sa_in,           // air volume fraction [-]
        double rv_in,           // vapor density [kg/m³]
        double muv_in,          // vapor viscosity [Pa·s]: 0 → air-only proxy (pre-Equil)
        double beta_v_in,       // vapor compressibility [1/Pa]: 0 → ideal-gas proxy (pre-Equil)
        double hv_in,           // vapor specific enthalpy [J/kg]: 0 → air-only (pre-Equil)
        double p_partial_air_in,    // air partial pressure [Pa]: p*xa (Dalton) or p (hydraulic)
        double p_partial_vapor_in,  // vapor partial pressure [Pa]: p*xv (Dalton) or p (hydraulic)
        double& xa_out,             double& xv_out,
        double& ra_out,             double& mu_gas_mix_out,
        double& rho_gas_mix_out,    double& h_gas_mix_out,
        double& beta_gas_mix_out)
{
    // ============================================================================
    // GAS PHASE MIXTURE PROPERTIES  (air + water vapor)  — MODE-AGNOSTIC
    //
    // Computes mixture properties given the partial pressures supplied by the caller.
    // This makes the function usable in two modes without internal mode logic:
    //
    //   full_gas_mixture   : caller passes Dalton partial pressures (p*xv, p*xa)
    //   hydraulic_gas_mixture: caller passes p_total for both
    //
    // Called in three contexts:
    //
    //   Pass 1  (UpdateEquilibratorVariables, PRE-Equilibrator, full_gas only):
    //     sv_in = svp(), sa_in = sap(), rv_in = rvp()
    //     muv_in = beta_v_in = hv_in = 0 (not yet available)
    //     p_partial_air_in / p_partial_vapor_in from mass-based Dalton estimate
    //     Purpose: estimate xv so Equilibrator receives correct partial pressure
    //
    //   Pass 2  (UpdateCSMPVariables, POST-Equilibrator, full_gas only):
    //     sv_in = sv(), sa_in = sa(), rv_in = Vapor.rho
    //     muv_in = Vapor.mu, beta_v_in = Vapor.beta, hv_in = Vapor.h
    //     p_partial_air_in / p_partial_vapor_in: Pass-1 Dalton values
    //     Purpose: compute final mixture properties; ra_out refined at p*xa
    //
    //   hydraulic_gas_mixture (UpdateCSMPVariables, POST-Equilibrator):
    //     sv_in = sv(), sa_in = sa(), rv_in = Vapor.rho (at p_total)
    //     muv_in = Vapor.mu, beta_v_in = Vapor.beta, hv_in = Vapor.h
    //     p_partial_air_in = p_partial_vapor_in = p_total
    //     ra_out: air density at p_total (no Dalton refinement)
    //
    // Volume fractions = molar fractions for ideal gases (Avogadro's law).
    // ============================================================================

    const double s_gas = sv_in + sa_in;

    // ── No gas phase: safe pure-air defaults ─────────────────────────────────
    if (s_gas <= 0.)
    {
        xa_out          = 1.;
        xv_out          = 0.;
        ra_out          = air.Density(t, p_partial_air_in > 0. ? p_partial_air_in : p);
        rho_gas_mix_out = ra_out;
        mu_gas_mix_out  = air.Viscosity(t);
        h_gas_mix_out   = air.Enthalpy(t);
        beta_gas_mix_out= (p_partial_air_in > 0.) ? 1.0 / p_partial_air_in : 0.;
        return;
    }

    // ── Volume / molar fractions ─────────────────────────────────────────────
    xv_out = sv_in / s_gas;
    xa_out = sa_in / s_gas;

    // ── Air density at the supplied partial pressure ──────────────────────────
    // full_gas: p_partial_air_in = p*xa (Dalton) → ra at partial pressure
    // hydraulic: p_partial_air_in = p_total       → ra at total pressure
    ra_out = air.Density(t, p_partial_air_in);

    //FIX 01-04 ── Mixture density ───────────────────────────────────────────────────────
    // If Dalton: rv_in and ra_out are partial densities (already scaled by x). Summing is correct.
    // If Hydraulic: rv_in and ra_out are pure-phase densities. Summing is wrong.
    // Using volume-weighting (xv_out*rv + xa_out*ra) works correctly for both.

    if (full_gas_mixture) {
        // In Dalton mode, the 'input' rv_in is already a partial density.
        // Summing them is the standard definition.
        rho_gas_mix_out = rv_in + ra_out;
    } else {
        // In Hydraulic mode, we must volume-weight the pure-phase densities.
        rho_gas_mix_out = xv_out * rv_in + xa_out * ra_out;
    }

    // ── Mass fractions ────────────────────────────────────────────────────────
    if (rho_gas_mix_out <= 0.)
        cerr << "\n[ComputeGasMixtureProperties] rho_gas_mix=0 — falling back to molar fractions"
             << "  rv_in=" << rv_in << "  ra_out=" << ra_out
             << "  sv=" << sv_in << "  sa=" << sa_in
             << "  t=" << t << " oC  p=" << p << " Pa" << endl;
    const double omega_v = (rho_gas_mix_out > 0.) ? rv_in  / rho_gas_mix_out : xv_out;
    const double omega_a = (rho_gas_mix_out > 0.) ? ra_out / rho_gas_mix_out : xa_out;

    // ── Mixture viscosity: molar-fraction weighted ─────────────────────────────
    const double mua_now = air.Viscosity(t);
    mu_gas_mix_out = (muv_in > 0.) ? xv_out * muv_in + xa_out * mua_now
                                   : mua_now;

    // ── Mixture specific enthalpy: mass-fraction weighted ─────────────────────
    h_gas_mix_out = (hv_in != 0.) ? omega_v * hv_in + omega_a * air.Enthalpy(t)
                                  : air.Enthalpy(t);

    // ── Mixture compressibility: mass-fraction weighted ───────────────────────
    // beta_a = 1/p_a; for hydraulic mode p_a = p_total
    const double beta_a = (p_partial_air_in    > 0.) ? 1.0 / p_partial_air_in    : 0.;
    const double beta_v = (beta_v_in > 0.)
            ? beta_v_in
            : (p_partial_vapor_in > 0.) ? 1.0 / p_partial_vapor_in : 0.;
    beta_gas_mix_out = omega_v * beta_v + omega_a * beta_a;
}


template<uint32_t dim>
double NaClH2OPropertiesVisitorPHX<dim>::TwoPhasePureWaterCompressibility(double cpl, double cpv)
{
    // compressibility of liquid/vapor mixture after Grant and Sorey, WRR 15(3) p. 684-686
    // compressibility modified as to add the vapor contribution cpv * rv * sv * phi as well

    // called when compressibility fluid is negative!! and if(Bulk.state == VL|| Bulk.state == VLH)

    const double denom_2ph = (hv() - hl()) * rl() * rv();
    if (std::abs(denom_2ph) < 1.e-12)
    {
        cerr << "\n[TwoPhasePureWaterCompressibility] Degenerate denominator (hv-hl)*rl*rv=" << denom_2ph
             << "  hv=" << hv() << "  hl=" << hl() << "  rl=" << rl() << "  rv=" << rv()
             << "  t=" << t_ << " oC  p=" << p_ << " Pa  state=" << Bulk.state << endl;
        return Liquid.beta;   // safe fallback
    }
    double product, b;
    product = ( rl() - rv() ) / denom_2ph;
    b  = ( 1.0 - phi() ) * cpr() * rr();
    b += phi() * rl() * cpl * sl();
    b += phi() * rv() * cpv * sv();
    b *= product * product;
    b *= ( t() + kelvin ) * 1.0 / phi();

    return b;
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::StorePropertiesAndFlags( Node<dim>* n )
{
    n->Store( mt_key, mt );
    n->Store( mtp_key,mtp );
    n->Store( Htp_key, Htp );
    n->Store( msp_key, msp );

    if( bogus_variables )
    {
        // Bogus variables already reported in UpdateCSMPVariables — skip mesh store silently.
    }
    else
    {
        // temperature
        n->Store( t_key,  t );
        n->Store( sl_key, sl );
        n->Store( sv_key, sv );
        n->Store( sh_key, sh );

        n->Store( state_key, state ); // really new?
        n->Store( apc_key, after_phasechange_counter );
        n->Store( dpc_key, dangerous_phase_change );

        // densities
        n->Store( rho_bulk_key, rho_bulk );
        n->Store( rl_key, rl );
        n->Store( rv_key, rv );

        n->Store( rh_key, rh );
        n->Store( ml_key, ml );
        n->Store( mv_key, mv );

        n->Store( mh_key, mh );

        //viscosities
        n->Store( mul_key,mul );
        n->Store( muv_key,muv );

        //enthalpy variables
        n->Store( hf_key, hf );
        n->Store( hl_key, hl );
        n->Store( hVl_key,hVl );
        n->Store( hCl_key,hCl );
        n->Store( hv_key, hv );
        n->Store( hVv_key,hVv );
        n->Store( hCv_key,hCv );

        n->Store( hh_key, hh );
        n->Store( hVh_key,hVh );
        n->Store( hCh_key,hCh );

        n->Store( cpf_key,cpf );
        n->Store( cpr_key,cpr );
        n->Store( ncp_key,ncp );

        // expansivities
        n->Store( beta_key,beta );         // bulk compressibility (includes halite if present), air excluded
        n->Store( CT_key,  CT );

        // extras
        n->Store( rl_transport_key, rl_transport );
        n->Store( rv_transport_key, rv_transport );

        n->Store( vol_fac_key, vol_fac );
        n->Store( nQ_key,    nQ );

        n->Store( phi_key,   phi );

        // mobilities
        n->Store( mml_key,   mml );
        n->Store( mmv_key,   mmv );
        n->Store( mmld_key,   mmld );
        n->Store( mmvd_key,   mmvd );
        n->Store( eml_key,   eml );
        n->Store( emv_key,   emv );
        n->Store( emld_key,   emld );
        n->Store( emvd_key,   emvd );

        n->Store( xml_key,   xml );
        n->Store( xmv_key,   xmv );
        n->Store( rvl_key,   rvl );
        n->Store( rvv_key,   rvv );

        n->Store( bfm_key,   bfm );
        n->Store( bfe_key,   bfe );
        n->Store( bfs_key,   bfs );

        n->Store( wt_key, wt );
        n->Store( xf_key, xf );
        n->Store( xl_key, xl );
        n->Store( xVl_key,xVl );
        n->Store( xCl_key,xCl );
        n->Store( xv_key, xv );
        n->Store( xVv_key,xVv );
        n->Store( xCv_key,xCv );
        n->Store( xh_key, xh );
        n->Store( xVh_key,xVh );
        n->Store( xCh_key,xCh );
        n->Store( xCf_key,xCf );

        if(with_air)
        {
            n->Store( sa_key, sa );

            n->Store( ra_key, ra );

            n->Store( ma_key, ma );

            n->Store( mua_key, mua );

            n->Store( ha_key, ha );
            n->Store( hVa_key,hVa );
            n->Store( hCa_key,hCa );

            n->Store( beta_air_key, beta_air );
            n->Store( ra_transport_key, ra_transport );

            n->Store( mma_key,   mma );
            n->Store( mmad_key,   mmad );
            n->Store( ema_key,   ema );
            n->Store( emad_key,   emad );

            n->Store( rva_key,   rva );

            n->Store( bfa_key,   bfa );

            // Only store boundary reference values at open boundary nodes —
            // these are only read by BoundaryIteration and meaningless elsewhere.
            if (open_boundaries && t.Flag() == DIRICH
                    && n->AtBoundary() != NOT && n->AtBoundary() != INTERNAL)
            {
                n->Store( rain_recharge_top_key,              rain_recharge_top              );
                n->Store( ref_frac_air_top_key,               ref_frac_air_top               );
                n->Store( ref_frac_vapor_in_air_top_key,      ref_frac_vapor_in_air_top      );
                n->Store( ref_enthalpy_vapor_in_air_top_key,  ref_enthalpy_vapor_in_air_top  );
                n->Store( ref_enthalpy_air_top_key,           ref_enthalpy_air_top           );
            }

            // ── Mixture quantities ────────────────────────────────────────────
            // Stored for diagnostics and visualisation only.  None of these are
            // read back as computational inputs: Pass 1 derives xv fresh from
            // mv()/ma(), and Pass 2 derives all mixture properties from the
            // equilibrated saturations.  See ReadAllVariables comment.
            if (full_gas_mixture || hydraulic_gas_mixture)
            {
                n->Store( xa_mix_key,          xa_mix          );
                n->Store( xv_mix_key,          xv_mix          );
                n->Store( p_partial_air_key,   p_partial_air   );
                n->Store( p_partial_vapor_key, p_partial_vapor );
                n->Store( rho_gas_mix_key,     rho_gas_mix     );
                n->Store( mu_gas_mix_key,      mu_gas_mix      );
                n->Store( h_gas_mix_key,       h_gas_mix       );
                n->Store( beta_gas_mix_key,    beta_gas_mix    );
            }
        }
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::StoreInitialPropertiesAndFlags( Node<dim>* n )
{
    // ── Saturations ───────────────────────────────────────────────────────────
    n->Store( t_key,  t  );
    n->Store( sl_key, sl );
    n->Store( sv_key, sv );

    n->Store( sa_key, sa );

    n->Store( sh_key, sh );

    n->Store( state_key,               state                     );
    n->Store( state_p_key,             state_p                   );
    n->Store( apc_key,                 after_phasechange_counter );
    n->Store( dpc_key,                 dangerous_phase_change    );
    n->Store( cpr_key,                 cpr                       );
    n->Store( ncp_key,                 ncp                       );

    // ── Mass concentrations [kg/m³ pore space] ────────────────────────────────
    n->Store( mt_key,       mt       );
    n->Store( mtp_key,      mtp      );
    n->Store( rho_bulk_key, rho_bulk );
    n->Store( rl_key,       rl       );
    n->Store( rv_key,       rv       );

    n->Store( ra_key, ra );

    n->Store( rh_key,       rh       );
    n->Store( ml_key,       ml       );
    n->Store( mv_key,       mv       );

    n->Store( ma_key, ma );

    n->Store( mh_key,       mh       );

    // ── Viscosities ───────────────────────────────────────────────────────────
    n->Store( mul_key, mul );
    n->Store( muv_key, muv );

    n->Store( mua_key, mua );

    // ── Enthalpy variables ────────────────────────────────────────────────────
    n->Store( hf_key,  hf  );   // does not include air
    n->Store( hl_key,  hl  );
    n->Store( hVl_key, hVl );
    n->Store( hCl_key, hCl );
    n->Store( hv_key,  hv  );
    n->Store( hVv_key, hVv );
    n->Store( hCv_key, hCv );

    n->Store( ha_key,  ha  );
    n->Store( hVa_key, hVa );
    n->Store( hCa_key, hCa );

    n->Store( hh_key,  hh  );
    n->Store( hVh_key, hVh );
    n->Store( hCh_key, hCh );
    n->Store( cpf_key, cpf );   // does not include air

    // ── Salt bookkeeping ──────────────────────────────────────────────────────
    n->Store( msp_key,  msp  );
    n->Store( wt_key,   wt   );
    n->Store( xf_key,   xf   );
    n->Store( xl_key,   xl   );
    n->Store( xVl_key,  xVl  );
    n->Store( xCl_key,  xCl  );
    n->Store( xv_key,   xv   );
    n->Store( xVv_key,  xVv  );
    n->Store( xCv_key,  xCv  );
    n->Store( xh_key,   xh   );
    n->Store( xVh_key,  xVh  );
    n->Store( xCh_key,  xCh  );
    n->Store( xCf_key,  xCf  );

    // ── Total enthalpy ────────────────────────────────────────────────────────
    n->Store( Htp_key, Htp );

    // ── Compressibilities ─────────────────────────────────────────────────────
    n->Store( beta_key,   beta   );   // H2O-NaCl fluid compressibility only
    n->Store( CT_key,     CT     );   // total compressibility (rock + fluid + air)

    n->Store( beta_air_key, beta_air );

    // ── Transport densities ───────────────────────────────────────────────────
    n->Store( rl_transport_key, rl_transport );
    n->Store( rv_transport_key, rv_transport );

    n->Store( ra_transport_key, ra_transport );

    n->Store( vol_fac_key, vol_fac );
    n->Store( nQ_key,      nQ      );
    n->Store( phi_key,     phi     );

    // ── Phase mobilities ──────────────────────────────────────────────────────
    n->Store( mml_key,  mml  );
    n->Store( mmv_key,  mmv  );
    n->Store( mmld_key, mmld );
    n->Store( mmvd_key, mmvd );
    n->Store( eml_key,  eml  );
    n->Store( emv_key,  emv  );
    n->Store( emld_key, emld );
    n->Store( emvd_key, emvd );

    n->Store( mma_key,  mma  );
    n->Store( mmad_key, mmad );
    n->Store( ema_key,  ema  );
    n->Store( emad_key, emad );

    n->Store( xml_key,  xml  );
    n->Store( xmv_key,  xmv  );
    n->Store( rvl_key,  rvl  );
    n->Store( rvv_key,  rvv  );

    n->Store( rva_key, rva );

    n->Store( bfm_key,  bfm  );
    n->Store( bfe_key,  bfe  );
    n->Store( bfs_key,  bfs  );

    n->Store( bfa_key, bfa );

    // ── Boundary ──────────────────────────────────────────────────────
    // Stored on every nodes for visualization purpose: store 0 instead of NaN
    n->Store( ref_enthalpy_liquid_top_key, ref_enthalpy_liquid_top );
    n->Store( ref_enthalpy_vapor_in_air_top_key, ref_enthalpy_vapor_in_air_top );
    n->Store( ref_enthalpy_air_top_key, ref_enthalpy_air_top );

    n->Store( ref_frac_vapor_in_air_top_key, ref_frac_vapor_in_air_top );
    n->Store( ref_frac_air_top_key, ref_frac_air_top );
    n->Store( rain_recharge_top_key, rain_recharge_top );

    // ── Mixture quantities (populated by UpdateCSMPVariables Pass 2 above) ────
    // For the moment, at init ma=0/sa=0 so Pass 2 produces pure-vapor defaults (xv=1, xa=0)
    // or pure-air safe defaults if sv=0 as well.  Values are never zero here.
    if (full_gas_mixture || hydraulic_gas_mixture)
    {
        n->Store( xa_mix_key,          xa_mix          );
        n->Store( xv_mix_key,          xv_mix          );
        n->Store( p_partial_air_key,   p_partial_air   );
        n->Store( p_partial_vapor_key, p_partial_vapor );
        n->Store( rho_gas_mix_key,     rho_gas_mix     );
        n->Store( mu_gas_mix_key,      mu_gas_mix      );
        n->Store( h_gas_mix_key,       h_gas_mix       );
        n->Store( beta_gas_mix_key,    beta_gas_mix    );
    }
    // Mass, density, and volumetric-enthalpy transport variables inherit the
    // temperature flag; salt transport variables inherit the salinity flag.
    if(!open_boundaries)
    {
        n->Status( ml_key,           t.Flag()  );
        n->Status( mv_key,           t.Flag()  );

        n->Status( ma_key, t.Flag() );

        n->Status( mh_key,           t.Flag()  );
        n->Status( rl_key,           t.Flag()  );
        n->Status( rv_key,           t.Flag()  );

        n->Status( ra_key, t.Flag() );

        n->Status( rh_key,           t.Flag()  );
        n->Status( rl_transport_key, t.Flag()  );
        n->Status( rv_transport_key, t.Flag()  );

        n->Status( ra_transport_key, t.Flag() );

        n->Status( hVl_key,          t.Flag()  );
        n->Status( hVv_key,          t.Flag()  );
        n->Status( hVh_key,          t.Flag()  );
        n->Status( hCl_key,          t.Flag()  );
        n->Status( hCv_key,          t.Flag()  );
        n->Status( hCh_key,          t.Flag()  );

        n->Status( hVa_key, t.Flag() );
        n->Status( hCa_key, t.Flag() );

        n->Status( xVl_key,          wt.Flag() );
        n->Status( xVv_key,          wt.Flag() );
        n->Status( xVh_key,          wt.Flag() );
        n->Status( xCl_key,          wt.Flag() );
        n->Status( xCv_key,          wt.Flag() );
        n->Status( xCh_key,          wt.Flag() );
        n->Status( xCf_key,          wt.Flag() );
    }

}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX(Model<dim>& model)
{
    // Only use for calculating properties when calculating static pressure.
    // NOTE: the local fluid object does not need to know about air; air energy
    // is handled by stripping/restoring inside the Equilibrator.
    //H2ONaClFluidProperties fluid(t_, p_, x_, h_fluid_, cp_rock_, rho_rock_, phi_, false);

    cout << "\nNaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX(): Initialising fluid properties";
    if(with_air) cout << " (air phase enabled)";
    cout << endl;

    // ── Diagnostic accumulators for boundary summary ─────
    double dbg_t_min     =  1e30, dbg_t_max     = -1e30, dbg_t_sum     = 0.;
    double dbg_p_min     =  1e30, dbg_p_max     = -1e30, dbg_p_sum     = 0.;
    double dbg_psat_min  =  1e30, dbg_psat_max  = -1e30, dbg_psat_sum  = 0.;
    double dbg_vfrac_min =  1e30, dbg_vfrac_max = -1e30, dbg_vfrac_sum = 0.;
    double dbg_hl_min    =  1e30, dbg_hl_max    = -1e30, dbg_hl_sum    = 0.;
    double dbg_hv_min    =  1e30, dbg_hv_max    = -1e30, dbg_hv_sum    = 0.;
    double dbg_ha_min    =  1e30, dbg_ha_max    = -1e30, dbg_ha_sum    = 0.;
    int    dbg_n_boundary = 0;

    for( auto it = model.Region("Model").NodesBegin(); it != model.Region("Model").NodesEnd(); it++ )
    {
        ReadAllVariables( *it );

        // ── Zero-initialise phase masses and rates ────────────────────────────
        tp()  = t();
        ml() = mlp() = mv() = mvp() = src_rate() = 0.0;

        // All air variables are zeroed here unconditionally for now
        ma()           = 0.0;   // air mass concentration [kg/m³ pore space]
        map()          = 0.0;   // previous-step air mass  (prevents dma_ = NaN)
        hCa()          = 0.0;   // volumetric air enthalpy content [J/m³]
        hCap()         = 0.0;   // previous-step hCa       (prevents dhCa_ = NaN)
        hVa()          = 0.0;   // specific volumetric air enthalpy
        ha()           = 0.0;   // specific air enthalpy [J/kg]
        sa()           = 0.0;   // air saturation
        ra()           = 0.0;   // air density (recomputed on first Visit)
        ra_transport() = 0.0;
        beta_air()     = 0.0;
        mua()          = 0.0;
        rva()          = 0.0;
        bfa()          = 0.0;
        mma()  = mmad() = 0.0;
        ema()  = emad() = 0.0;
        kra = 0.0;

        rain_recharge_top()              = 0.;
        ref_frac_air_top()               = 0.;
        ref_frac_vapor_in_air_top()      = 0.;   // ← add
        ref_enthalpy_liquid_top()        = 0.;
        ref_enthalpy_air_top()           = 0.;
        ref_enthalpy_vapor_in_air_top()  = 0.;

        // ── Rock heat capacity ────────────────────────────────────────────────
        cpr() = with_rock_liquidus_solidus ? rock.HeatCapacity(t(), tl(), ts())
                                           : rock.HeatCapacity(t());
        ncp() = cpr() * rr() * (1. - phi());

        if(add_fluid_contribution_to_heat_capacity)
            ncp() += cpf() * rho_bulk() * phi();
        // Air contribution to ncp is zero at initialization (ma=sa=0).

        // ── Fixed-temperature flag ────────────────────────────────────────────
        fixed_temperature = (t.Flag() == DIRICH);
        if(fixed_temperature) t_fixed = t();



        // ── Range checks ─────────────────────────────────────────────────────
        // Note: ScreenOutputEquilibratorVariables() is not called here because
        // most state variables are not yet initialised at this point.
        if (p() < 0.0 || t() < 0.0)
        {
            cerr << "\n[InitialPropertiesFromPTX] Negative p or t — cannot initialise"
                 << "  node=(" << (*it)->x() << ", " << (*it)->y() << ")"
                 << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
            csmp_error.Note(FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX",
                            "Negative pressure or temperature — terminating.");
        }
        if (p() > 5000.0e5 || t() > 1000.0)
        {
            cerr << "\n[InitialPropertiesFromPTX] p or t above EOS maximum"
                 << "  node=(" << (*it)->x() << ", " << (*it)->y() << ")"
                 << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
            csmp_error.Note(FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX",
                            "Pressure or temperature above maximum of lookup table — terminating.");
        }
        if (/*p() < 100000.0 ||*/ t() < 5.0)
        {
            cerr << "\n[InitialPropertiesFromPTX] p or t below EOS minimum"
                 << "  node=(" << (*it)->x() << ", " << (*it)->y() << ")"
                 << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
            csmp_error.Note(ERROR, "NaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX",
                            "Pressure or temperature below minimum of lookup table — erroneous results are possible.");
        }

        // ── H2O-NaCl initial equilibration ────────────────────────────────────────────

        CalculateAbsoluteVariables( *it );
        UpdateEquilibratorVariables( *it );

        wt_                     = wt();
        x_                      = Weight2XNaCl(wt_);
        if(x_ > 1.) x_          = 1.;
        h_fluid_                = 2086000.; // starting guess: critical enthalpy of water [J/kg]

        Bulk   = fluid.BulkProperties();
        Liquid = fluid.ReportLiquidProperties();
        Vapor  = fluid.ReportVaporProperties();
        Salt   = fluid.ReportSaltProperties();

        t()      = Bulk.t;
        h_fluid_ = Bulk.h;

        UpdateCSMPVariables( *it );

        after_phasechange_counter() = 0.0;
        dangerous_phase_change()    = 0.0;

        // ── Fluid masses ──────────────────────────────────────────────

        mt() = mtp() = Bulk.rho;

        if(with_air)
        {
            m_fluid_ = (mt() - ma()) * phi();   // H2O-NaCl mass only
            m_air_   =  ma()         * phi();   // = 0 at init
        }
        else
            m_fluid_ = mt() * phi();            // includes halite

        // ── Total enthalpy ────────────────────────────────────────────────────
        // Air contributes zero here because ma = 0 at init.
        H_current_ = Bulk.h * m_fluid_
                + rr() * (with_rock_liquidus_solidus ? rock.Enthalpy(t(), tl(), ts())
                                                     : rock.Enthalpy(t())) * (1.0 - phi());
        Htp() = H_current_;

        // ── Salt and fluid fraction bookkeeping ───────────────────────────────
        ms() = msp() = xCl() + xCv() + xCh();
        mh() = Salt.mf * mt();

        //mf() = mt() - mh();
        //if(with_air) mf() -= ma();

        // ── Transport densities ───────────────────────────────────────────────
        rl_transport() = rl();
        rv_transport() = rv();
        if(with_air) ra_transport() = ra();   // = 0 at init

        rho_bulk() = Bulk.rho;

        nQ()      = 0.;
        vol_fac() = 1.0;

        if (with_air)
            beta_air() = (p() > 0.) ? (full_gas_mixture ? beta_gas_mix()
                                                        : air.Compressibility(p()))
                                    : 0.;

        // Total compressibility — single source of truth (see ComputeTotalCompressibility).
        CT() = ComputeTotalCompressibility();

        // ── Phase mobilities ──────────────────────────────────────────────────
        mml() = mmv() = mmld() = mmvd() = 0.0;
        eml() = emv() = rvl()  = rvv()  = 0.0;
        emld() = emvd() = 0.0;
        bfm()  = bfe()  = bfs() = xml() = xmv() = 0.0;
        // Air mobilities already zeroed above

        // ── Relative permeabilities ───────────────────────────────────────────
        if(with_air)
        {
            krl  = RelativePermeabilityLiquid(sl(), sv(), sa());
            double kr_gas = RelativePermeabilityGas(sl(), sv(), sa());
            double s_gas  = sv() + sa();

            if(s_gas > 0.)
            {
                // Simple fractional flow — air gets its share of total gas kr
                krv = kr_gas * sv() / s_gas;
                kra = kr_gas * sa() / s_gas;

            }

            else krv = kra = 0.;
        }

        else
        {
            krl = RelativePermeabilityLiquid(sl(), sv(), 0.);
            krv = RelativePermeabilityGas(sl(), sv(), 0.);
        }

        // ── Liquid and vapor mobilities ───────────────────────────────────────
        if(sl() > 0.) mml() = krl * rl_transport() / mul();
        if(sv() > 0.) mmv() = krv * rv_transport() / muv();

        if(with_air && sa() > 0.) mma() = kra * ra_transport() / mua();
        // Air mobility: sa = 0 at init, so mma remains zero

        if(thermodymamic_density_in_gravity_term)
        {
            if(sl() > 0.) mmld() = mml() * rl();
            if(sv() > 0.) mmvd() = mmv() * rv();
            if(with_air && sa() > 0.) mmad() = mma() * ra();
        }
        else
        {
            if(sl() > 0.) mmld() = mml() * rl_transport();
            if(sv() > 0.) mmvd() = mmv() * rv_transport();
            if(with_air && sa() > 0.) mmad() = mma() * ra_transport();
        }

        if(sl() > 0.) eml() = krl * hVl() / mul();
        if(sv() > 0.) emv() = krv * hVv() / muv();
        if(with_air && sa() > 0.) ema() = kra * hVa() / mua();
        // Air enthalpy mobility: sa = 0 at init, so ema remains zero

        if(thermodymamic_density_in_gravity_term)
        {
            if(sl() > 0.) emld() = eml() * rl();
            if(sv() > 0.) emvd() = emv() * rv();
            if(with_air && sa() > 0.) emad() = ema() * ra();
        }
        else
        {
            if(sl() > 0.) emld() = eml() * rl_transport();
            if(sv() > 0.) emvd() = emv() * rv_transport();
            if(with_air && sa() > 0.) emad() = ema() * ra_transport();
        }

        if(sl() > 0.) { xml() = krl * xVl() / mul();  rvl() = krl / mul(); }
        if(sv() > 0.) { xmv() = krv * xVv() / muv();  rvv() = krv / muv(); }
        // no salt in air

        // ── Reference enthalpies for open-boundary inflow ─────────────────────
        // Stored as nodal ScalarVariables so BoundaryIteration can read them
        // back per-node (boundary T/P can vary spatially).

        ref_enthalpy_liquid_top() = Bulk.h;

        if(with_air)
        {
            ref_enthalpy_air_top()    = air.Enthalpy(t());
        }

        double p_sat = p();
        ref_enthalpy_vapor_in_air_top() = 0.;
        ref_frac_vapor_in_air_top() = 0.;

        if (open_boundaries && with_air)
        {
            // Pure-water saturation pressure at this node's equilibrated temperature.
            // Must use EOS lookup — Liquid.p is total pressure for subcooled liquid,
            // NOT p_sat (they are only equal when the fluid is in VL state).


            if( ref_humidity > 0.)
            {
                p_sat = fluid.SaturationPressureFromT(t());

                ref_enthalpy_vapor_in_air_top() = fluid.SaturationVaporEnthalpyFromT(t());

                const double xv_mol = ref_humidity * p_sat / p();
                ref_frac_vapor_in_air_top() = std::min(1., xv_mol * M_VAPOR
                                                       / (xv_mol * M_VAPOR + (1. - xv_mol) * M_AIR));
            }

            if((*it)->AtBoundary()==TOP)
            {
                // accumulate Top boundary for summary
                dbg_t_min     = std::min(dbg_t_min,    t());
                dbg_t_max     = std::max(dbg_t_max,    t());
                dbg_t_sum    += t();

                dbg_p_min     = std::min(dbg_p_min,    p());
                dbg_p_max     = std::max(dbg_p_max,    p());
                dbg_p_sum    += p();

                dbg_psat_min  = std::min(dbg_psat_min, p_sat);
                dbg_psat_max  = std::max(dbg_psat_max, p_sat);
                dbg_psat_sum += p_sat;

                dbg_vfrac_min  = std::min(dbg_vfrac_min, ref_frac_vapor_in_air_top());
                dbg_vfrac_max  = std::max(dbg_vfrac_max, ref_frac_vapor_in_air_top());
                dbg_vfrac_sum += ref_frac_vapor_in_air_top();

                dbg_hl_min    = std::min(dbg_hl_min,   ref_enthalpy_liquid_top());
                dbg_hl_max    = std::max(dbg_hl_max,   ref_enthalpy_liquid_top());
                dbg_hl_sum   += ref_enthalpy_liquid_top();

                dbg_hv_min    = std::min(dbg_hv_min,   ref_enthalpy_vapor_in_air_top());
                dbg_hv_max    = std::max(dbg_hv_max,   ref_enthalpy_vapor_in_air_top());
                dbg_hv_sum   += ref_enthalpy_vapor_in_air_top();

                dbg_ha_min    = std::min(dbg_ha_min,   ref_enthalpy_air_top());
                dbg_ha_max    = std::max(dbg_ha_max,   ref_enthalpy_air_top());
                dbg_ha_sum   += ref_enthalpy_air_top();

                ++dbg_n_boundary;
            }
        }

        state()   = Bulk.state;
        state_p() = state();

        StoreInitialPropertiesAndFlags( *it );
    }

    if (open_boundaries && dbg_n_boundary > 0)
    {
        const double n = static_cast<double>(dbg_n_boundary);
        cerr << "\n[InitialPropertiesFromPTX] Variable-pressure open boundary summary:\n"
             << "  ref_humidity        = " << ref_humidity            << " [-]\n"
             << "  ref_sal             = " << ref_sal * 100.           << " wt%\n"
             << "  boundary nodes      = " << dbg_n_boundary           << "\n"
             << "  T        [°C]  min/mean/max = "
             << dbg_t_min     << " / " << dbg_t_sum/n     << " / " << dbg_t_max     << "\n"
             << "  p        [Pa]  min/mean/max = "
             << dbg_p_min     << " / " << dbg_p_sum/n     << " / " << dbg_p_max     << "\n"
             << "  p_sat    [Pa]  min/mean/max = "
             << dbg_psat_min  << " / " << dbg_psat_sum/n  << " / " << dbg_psat_max  << "\n"
             << "  vap_frac [kg_vap/kg_gas]   min/mean/max = "
             << dbg_vfrac_min << " / " << dbg_vfrac_sum/n << " / " << dbg_vfrac_max << "\n"
             << "  h_liq  [J/kg]  min/mean/max = "
             << dbg_hl_min    << " / " << dbg_hl_sum/n    << " / " << dbg_hl_max    << "\n"
             << "  h_vap  [J/kg]  min/mean/max = "
             << dbg_hv_min    << " / " << dbg_hv_sum/n    << " / " << dbg_hv_max    << "\n"
             << "  h_air  [J/kg]  min/mean/max = "
             << dbg_ha_min    << " / " << dbg_ha_sum/n    << " / " << dbg_ha_max    << "\n\n";
    }

    // ── Rain recharge: convert mm/year to per-node mass rates ─────────────────
    // Must come after the main loop because it needs cpv (pore volume) which was
    // just stored by StoreInitialPropertiesAndFlags.  Overwrites rain_recharge_top
    // on TOP boundary nodes; interior nodes keep the zero set in the loop above.
    if (with_air && rain_mm_per_year_top > 0.)
    {
        SetPrecipitationRecharge(model);
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::CheckPhysicalState( Node<dim>* n, const char* stage )
{
    // =========================================================================
    // PHYSICAL STATE VALIDATOR
    //
    // Called at three pipeline checkpoints inside Visit():
    //
    //   "ReadAllVariables" — raw values just read from the CSMP mesh.
    //     Checks: finiteness of mesh inputs, thermodynamic bounds (p, t, phi,
    //     rr, mt, wt, ms), air mass invariants.
    //     NOTE: ms() is STALE at this stage (not yet recomputed). Only msp()
    //     (freshly read from mesh) is checked for salt invariants.
    //
    //   "PreEquilibrate"   — after UpdateEquilibratorVariables, just before
    //     Equilibrate(). Checks everything in ReadAllVariables PLUS the derived
    //     inputs we are about to feed the EOS: H_current_, m_fluid_, phi_
    //     (shaved), p_current_, wt_, x_, and air shaving correctness. Catching
    //     problems here means we know the Equilibrator received bad input.
    //     NOTE: ms() is now valid (computed in CalculateAbsoluteVariables).
    //     Both msp() and ms() are checked.
    //
    //   "PostPrepare"      — after PrepareVariablesForStorage, just before
    //     StorePropertiesAndFlags(). Checks everything: EOS outputs (rl, rv,
    //     mul, muv, beta, saturations, closure), volume factors, nQ, mobilities,
    //     CT, relative permeabilities.  This is the final gate before mesh store.
    //     NOTE: ms() is the authoritative current salt mass.
    //
    // Every violation is written to cerr with a consistent tagged prefix:
    //   [STAGE] (x, y[, z]) | VARIABLE = value | reason
    //
    // If any violation is found the full DumpNodeState() dump goes to cout
    // (details), with a brief summary header to cerr.
    // =========================================================================

    bool any_violation = false;

    // ── Location tag ─────────────────────────────────────────────────────────
    std::ostringstream loc_ss;
    loc_ss << "[" << stage << "] (";
    if (dim >= 2) loc_ss << n->x() << ", " << n->y();
    if (dim == 3) loc_ss << ", " << n->z();
    loc_ss << ")";
    const std::string loc = loc_ss.str();

    // ── Helpers ───────────────────────────────────────────────────────────────
    auto flag = [&](const std::string& var, double val, const std::string& reason)
    {
        cerr << loc << " | " << var << " = " << val << " | " << reason << "\n";
        any_violation = true;
    };

    auto check_finite = [&](const std::string& name, double v)
    {
        if (!std::isfinite(v))
            flag(name, v, "NaN or Inf — upstream computation blew up");
    };

    const bool is_read    = (stage == std::string("ReadAllVariables"));
    const bool is_pre_eq  = (stage == std::string("PreEquilibrate"));
    const bool is_post    = (stage == std::string("PostPrepare"));

    // =========================================================================
    // STAGE 1 — ReadAllVariables  AND  PreEquilibrate
    // Basic mesh inputs and thermodynamic bounds.  Checked at both early stages
    // so problems are caught as soon as possible.
    // =========================================================================
    if (is_read || is_pre_eq)
    {
        // Finiteness of raw mesh variables
        check_finite("t",           t());
        check_finite("p",           p());
        check_finite("phi",         phi());
        check_finite("rr",          rr());
        check_finite("mt",          mt());
        check_finite("msp",         msp());   // ms() is stale at ReadAllVariables; msp() is freshly read
        if (is_pre_eq) check_finite("ms",  ms());   // valid after CalculateAbsoluteVariables
        check_finite("ml",          ml());
        check_finite("mv",          mv());
        check_finite("hCl",         hCl());
        check_finite("hCv",         hCv());

        if (with_air)
        {
            check_finite("ma",  ma());
            check_finite("hCa", hCa());
        }

        // Thermodynamic bounds
        if (std::isfinite(p()) && p() <= 0.)
            flag("p [Pa]", p(), "non-positive pressure — EOS undefined");
        // if (std::isfinite(p()) && p() < 1.0e5)
        //     flag("p [Pa]", p(), "below EOS minimum (~1 bar)");
        if (std::isfinite(t()) && t() < 5.)
            flag("t [°C]", t(), "below EOS minimum (5 °C)");
        if (std::isfinite(t()) && t() > 1500.)
            flag("t [°C]", t(), "above EOS validity range (1500 °C)");
        if (std::isfinite(phi()) && (phi() <= 0. || phi() >= 1.))
            flag("phi", phi(), "porosity outside (0, 1)");
        if (std::isfinite(rr()) && rr() <= 0.)
            flag("rr [kg/m³]", rr(), "non-positive rock density");
        if (std::isfinite(mt()) && mt() < 0.)
            flag("mt [kg/m³]", mt(), "negative total mass — mass balance violated");

        // ── Salt mass checks ─────────────────────────────────────────────────
        // ReadAllVariables: ms() is STALE (not recomputed until CalculateAbsoluteVariables).
        //   Check msp() only — the freshly read stored value from the mesh.
        // PreEquilibrate: ms() is now valid (just computed in CalculateAbsoluteVariables).
        //   Check both msp() (what was stored last step) and ms() (current value).
        if (is_read)
        {
            if (std::isfinite(msp()) && msp() < 0.)
                flag("msp [kg/m³]", msp(), "negative stored salt mass");
            if (std::isfinite(msp()) && std::isfinite(mt()) && msp() > mt() + 1.e-10)
                flag("msp/mt", msp()/mt(), "stored salt mass exceeds total mass");
        }
        if (is_pre_eq)
        {
            if (std::isfinite(msp()) && msp() < 0.)
                flag("msp [kg/m³]", msp(), "negative stored salt mass (previous step)");
            if (std::isfinite(msp()) && std::isfinite(mt()) && msp() > mt() + 1.e-10)
                flag("msp/mt", msp()/mt(), "stored salt mass exceeds total mass (previous step)");
            if (std::isfinite(ms()) && ms() < 0.)
                flag("ms [kg/m³]", ms(), "negative salt mass after CalculateAbsoluteVariables");
            if (std::isfinite(ms()) && std::isfinite(mt()) && ms() > mt() + 1.e-10)
                flag("ms/mt", ms()/mt(), "salt mass exceeds total mass after CalculateAbsoluteVariables");
        }

        if (std::isfinite(wt()) && (wt() < 0. || wt() > 100.))
            flag("wt [wt%]", wt(), "salinity out of [0, 100]");
    }

    // =========================================================================
    // STAGE 2 (PreEquilibrate only) — Equilibrator inputs
    // =========================================================================
    if (is_pre_eq)
    {
        check_finite("H_current_", H_current_);
        check_finite("m_fluid_",   m_fluid_);
        check_finite("phi_",       phi_);
        check_finite("p_current_", p_current_);
        check_finite("wt_",        wt_);
        check_finite("x_",         x_);

        if (std::isfinite(m_fluid_) && m_fluid_ < 0.)
            flag("m_fluid_ [kg]", m_fluid_, "negative fluid mass fed to Equilibrator — will produce bogus EOS output");
        if (std::isfinite(phi_) && phi_ <= 0.)
            flag("phi_ (shaved)", phi_, "non-positive shaved porosity fed to Equilibrator");
        if (std::isfinite(p_current_) && p_current_ <= 0.)
            flag("p_current_ [Pa]", p_current_, "non-positive pressure fed to Equilibrator");
        if (std::isfinite(wt_) && (wt_ < 0. || wt_ > 100.))
            flag("wt_ [wt%]", wt_, "salinity out of [0, 100] before EOS call");

        if (with_air)
        {
            check_finite("sa",  sa());
            check_finite("ra",  ra());
            check_finite("phi_*(1-sa)", phi() * (1.0 - sa()));

            if (std::isfinite(ma()) && ma() < 0.)
                flag("ma [kg/m³]", ma(), "negative air mass");
            if (std::isfinite(ma()) && std::isfinite(mt()) && ma() > mt() + 1.e-10)
                flag("ma/mt", ma()/mt(), "air mass exceeds total mass");

            if (std::isfinite(sa()) && (sa() < 0. || sa() > 1. + 1.e-8))
                flag("sa", sa(), "air saturation out of [0, 1] before EOS call");
            const double phi_shaved = phi() * (1.0 - sa());
            if (std::isfinite(phi_shaved) && phi_shaved < 0.)
                flag("phi*(1-sa)", phi_shaved, "negative shaved porosity — Equilibrator ill-posed");
            if (std::isfinite(ra()) && ra() <= 0. && ma() > 1.e-12)
                flag("ra [kg/m³]", ra(), "non-positive air density with non-zero air mass");
            if (full_gas_mixture)
            {
                check_finite("p_partial_air",   p_partial_air());
                check_finite("p_partial_vapor", p_partial_vapor());
                if (std::isfinite(p_partial_air()) && p_partial_air() < 0.)
                    flag("p_partial_air [Pa]", p_partial_air(), "negative air partial pressure");
                if (std::isfinite(p_partial_vapor()) && p_partial_vapor() < 0.)
                    flag("p_partial_vapor [Pa]", p_partial_vapor(), "negative vapor partial pressure");
            }
        }
    }

    // =========================================================================
    // STAGE 3 (PostPrepare only) — Full pipeline output
    // =========================================================================
    if (is_post)
    {
        // ── Finiteness ───────────────────────────────────────────────────────
        check_finite("t",               t());
        check_finite("p",               p());
        check_finite("mt",              mt());
        check_finite("H_current_",      H_current_);
        check_finite("rl",              rl());
        check_finite("rv",              rv());
        check_finite("mul",             mul());
        check_finite("muv",             muv());
        check_finite("Bulk.rho",        Bulk.rho);
        check_finite("ncp",             ncp());
        check_finite("beta",            beta());
        check_finite("nQ",              nQ());
        check_finite("volume_factor_LHS", volume_factor_LHS);
        check_finite("volume_factor_RHS", volume_factor_RHS);
        check_finite("CT",              CT());
        check_finite("mml",             mml());
        check_finite("mmv",             mmv());
        check_finite("eml",             eml());
        check_finite("emv",             emv());
        if (with_air)
        {
            check_finite("ma",  ma());
            check_finite("sa",  sa());
            check_finite("ra",  ra());
            check_finite("ha",  ha());
            check_finite("mma", mma());
            check_finite("ema", ema());
        }

        // ── Thermodynamic bounds ─────────────────────────────────────────────
        if (std::isfinite(mt()) && mt() < 0.)
            flag("mt [kg/m³]", mt(), "negative total mass");
        if (std::isfinite(beta()) && beta() < 0.)
            flag("beta [1/Pa]", beta(), "negative fluid compressibility");

        // ── Salt mass ────────────────────────────────────────────────────────
        // ms() is authoritative here — PrepareVariablesForStorage has run.
        if (std::isfinite(ms()) && ms() < 0.)
            flag("ms [kg/m³]", ms(), "negative salt mass");
        if (std::isfinite(ms()) && std::isfinite(mt()) && ms() > mt() + 1.e-10)
            flag("ms/mt", ms()/mt(), "salt mass exceeds total mass");

        // ── Saturations ──────────────────────────────────────────────────────
        if (std::isfinite(sl()) && sl() < -1.e-10)
            flag("sl", sl(), "negative liquid saturation");
        if (std::isfinite(sv()) && sv() < -1.e-10)
            flag("sv", sv(), "negative vapor saturation");
        if (std::isfinite(sh()) && sh() < -1.e-10)
            flag("sh", sh(), "negative halite saturation");
        if (std::isfinite(sl()) && sl() > 1. + 1.e-8)
            flag("sl", sl(), "liquid saturation > 1");
        if (std::isfinite(sv()) && sv() > 1. + 1.e-8)
            flag("sv", sv(), "vapor saturation > 1");
        if (std::isfinite(sh()) && sh() > 1. + 1.e-8)
            flag("sh", sh(), "halite saturation > 1");
        {
            double sat_sum = sl() + sv() + sh();
            if (with_air && std::isfinite(sa())) sat_sum += sa();
            if (std::isfinite(sat_sum) && std::abs(sat_sum - 1.) > 1.e-6)
                flag("sl+sv+sh+sa", sat_sum, "saturation closure violated (should be 1.0)");
        }

        // ── EOS outputs ──────────────────────────────────────────────────────
        if (std::isfinite(rl()) && sl() > 1.e-10 && rl() <= 0.)
            flag("rl [kg/m³]", rl(), "non-positive liquid density while sl > 0");
        if (std::isfinite(rv()) && sv() > 1.e-10 && rv() <= 0.)
            flag("rv [kg/m³]", rv(), "non-positive vapor density while sv > 0");
        if (std::isfinite(mul()) && sl() > 1.e-10 && mul() <= 0.)
            flag("mul [Pa.s]", mul(), "non-positive liquid viscosity while sl > 0");
        if (std::isfinite(muv()) && sv() > 1.e-10 && muv() <= 0.)
            flag("muv [Pa.s]", muv(), "non-positive vapor viscosity while sv > 0");

        // ── ENERGY-WITHOUT-MASS INVARIANT ────────────────────────────────────
        // hCv = hv * rv * sv, so hCv/mv = hv (specific vapor enthalpy). A cell
        // with vapor present (sv>0, mv>0) but ZERO enthalpy content means the
        // specific enthalpy hv itself came out zero — a degenerate EOS/datum
        // state, not a transport or leakage fault. It propagates through the SB
        // relay (enthalpy leg = M*hCv/mv = 0) and shows up downstream as
        // "mv>0, hCv=0". Catching it HERE (equilibration output, PostPrepare)
        // pins the origin to the EOS. Liquid mirror included.
        if (std::isfinite(hv()) && sv() > 1.e-10 && rv() > 0. && hv() == 0.)
            flag("hv [J/kg]", hv(), "ZERO specific vapor enthalpy while vapor present (sv>0) -> hCv=0 with mv>0");
        if (std::isfinite(hCv()) && mv() > 1.e-10 && hCv() == 0.)
            flag("hCv [J/m3]", hCv(), "ZERO vapor enthalpy content while mv>0 -> energy/mass desync source");
        if (std::isfinite(hl()) && sl() > 1.e-10 && rl() > 0. && hl() == 0.)
            flag("hl [J/kg]", hl(), "ZERO specific liquid enthalpy while liquid present (sl>0)");

        {
            const double h_threshold = -1.e9;
            if (std::isfinite(hCl()) && hCl() < h_threshold)
                flag("hCl [J/m³]", hCl(), "extremely negative liquid enthalpy content");
            if (std::isfinite(hCv()) && hCv() < h_threshold)
                flag("hCv [J/m³]", hCv(), "extremely negative vapor enthalpy content");
            if (with_air && std::isfinite(hCa()) && hCa() < h_threshold)
                flag("hCa [J/m³]", hCa(), "extremely negative air enthalpy content");
        }

        // ── Air invariants ───────────────────────────────────────────────────
        if (with_air)
        {
            if (std::isfinite(ma()) && ma() < 0.)
                flag("ma [kg/m³]", ma(), "negative air mass");
            if (std::isfinite(ma()) && std::isfinite(mt()) && ma() > mt() + 1.e-10)
                flag("ma/mt", ma()/mt(), "air mass exceeds total mass");
            if (std::isfinite(sa()) && sa() < 0.)
                flag("sa", sa(), "negative air saturation");
            if (std::isfinite(sa()) && sa() > 1. + 1.e-8)
                flag("sa", sa(), "air saturation > 1");
            if (std::isfinite(ra()) && ra() <= 0. && ma() > 1.e-12)
                flag("ra [kg/m³]", ra(), "non-positive air density with non-zero air mass");
            if (full_gas_mixture)
            {
                if (std::isfinite(p_partial_air()) && p_partial_air() < 0.)
                    flag("p_partial_air [Pa]", p_partial_air(), "negative air partial pressure");
                if (std::isfinite(p_partial_vapor()) && p_partial_vapor() < 0.)
                    flag("p_partial_vapor [Pa]", p_partial_vapor(), "negative vapor partial pressure");
                if (std::isfinite(p_partial_air()) && std::isfinite(p_partial_vapor()) &&
                        std::isfinite(p()) &&
                        std::abs(p_partial_air() + p_partial_vapor() - p()) > 1.e-3 * p())
                    flag("p_air+p_vap [Pa]", p_partial_air()+p_partial_vapor(),
                         "Dalton sum does not match total pressure");
            }
        }

        // ── Volume factors ───────────────────────────────────────────────────
        if (std::isfinite(volume_factor_LHS) && volume_factor_LHS < 0.)
            flag("volume_factor_LHS", volume_factor_LHS, "negative LHS volume factor");
        if (std::isfinite(volume_factor_RHS) && volume_factor_RHS < 0.)
            flag("volume_factor_RHS", volume_factor_RHS, "negative RHS volume factor");

        // ── Mobilities, CT, relative permeabilities ──────────────────────────
        if (std::isfinite(CT()) && CT() <= 0.)
            flag("CT [1/Pa·m³]", CT(), "non-positive total compressibility — pressure equation ill-conditioned");
        if (krl < -1.e-10 || krl > 1. + 1.e-8)
            flag("krl", krl, "liquid relative permeability outside [0, 1]");
        if (krv < -1.e-10 || krv > 1. + 1.e-8)
            flag("krv", krv, "vapor relative permeability outside [0, 1]");
        if (with_air && (kra < -1.e-10 || kra > 1. + 1.e-8))
            flag("kra", kra, "air relative permeability outside [0, 1]");
        {
            const double kr_sum = krl + krv + (with_air ? kra : 0.);
            const double kr_max = 1.0 - sh();
            if (std::isfinite(kr_sum) && kr_sum > kr_max + 1.e-8)
                flag("krl+krv+kra", kr_sum, "relative permeability sum exceeds (1-sh)");
        }
        if (std::isfinite(mml()) && mml() < 0.)
            flag("mml", mml(), "negative liquid mass mobility");
        if (std::isfinite(mmv()) && mmv() < 0.)
            flag("mmv", mmv(), "negative vapor mass mobility");
        if (with_air && std::isfinite(mma()) && mma() < 0.)
            flag("mma", mma(), "negative air mass mobility");
    }

    // =========================================================================
    // FULL STATE DUMP — cerr summary header, cout full state
    // =========================================================================
    if (any_violation)
    {
        cerr << loc << " | *** " << stage
             << ": physical-state violation(s) detected — see stdout for full state dump ***\n";
        const char* dump_site = is_post    ? "CheckPhysicalState/PostPrepare"    :
                                             is_pre_eq  ? "CheckPhysicalState/PreEquilibrate" :
                                                          "CheckPhysicalState/ReadAllVariables";
        DumpNodeState(dump_site);
        PAUSE();PAUSE();
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::DumpNodeState( const char* call_site )
{
    // =========================================================================
    // DumpNodeState — full diagnostic dump of all visitor state variables.
    //
    // call_site identifies the pipeline position and drives the stale-variable
    // warning printed in the header. Call sites and their positions:
    //
    //   PRE-EOS (EOS not yet run — rho/sat/visc are previous-timestep values):
    //     CheckForOutOfRange               after ReadAllVariables, before CalcAbsVars
    //     CheckPhysicalState/ReadAllVariables   same
    //     CheckPhysicalState/PreEquilibrate     after UpdateEquilibratorVariables
    //                                           H_current_, m_fluid_, phi_ are set
    //
    //   POST-EOS, PRE-VF (EOS done; volume factors/CT/mobilities not computed):
    //     Equilibrate/Fatal                after Equilibrate(); rl/rv/sl/sv still
    //                                      stale — UpdateCSMPVariables not run yet
    //     UpdateCSMPVariables/BogusEOS     rl/rv/sl/sv/mul/muv current; vf not yet
    //     UpdateCSMPVariables/SaturationClosure   same
    //
    //   POST-VF, PRE-PREP (volume factors current; CT/mobilities not yet):
    //     VolumeFactorComputations/NegativeVF
    //
    //   MID-ITERATION (trial state — conserved vars modified, EOS from last trial):
    //     BoundaryIteration/InvalidBracket
    //
    //   FULLY CURRENT (PrepareVariablesForStorage done — everything valid):
    //     CheckPhysicalState/PostPrepare
    //     Visit/UndefinedState
    // =========================================================================

    const std::string cs(call_site);

    const bool pre_eos    = (cs == "CheckForOutOfRange"                  ||
                             cs == "CheckPhysicalState/ReadAllVariables"  ||
                             cs == "CheckPhysicalState/PreEquilibrate");
    const bool post_eos_pre_csmp = (cs == "Equilibrate/Fatal");
    const bool post_csmp  = (cs == "UpdateCSMPVariables/BogusEOS"        ||
                             cs == "UpdateCSMPVariables/SaturationClosure");
    const bool post_vf    = (cs == "VolumeFactorComputations/NegativeVF");
    const bool mid_iter   = (cs == "BoundaryIteration/InvalidBracket");
    const bool post_prep  = (cs == "CheckPhysicalState/PostPrepare"      ||
                             cs == "Visit/UndefinedState");

    cout << "\n========================================================================\n";
    cout << "  NODE STATE DUMP  [called from: " << call_site << "]\n";
    if (pre_eos)
        cout << "  *** PRE-EOS: densities, saturations, viscosities, volume factors,\n"
             << "  *** CT and mobilities are FROM PREVIOUS TIMESTEP — not yet updated.\n"
             << "  *** H_current_, m_fluid_, phi_"
             << (cs == "CheckPhysicalState/PreEquilibrate" ? " are current (Equilibrator inputs)."
                                                           : " not yet computed.")
             << "\n";
    else if (post_eos_pre_csmp)
        cout << "  *** POST-EQUILIBRATE: EOS ran on Bulk/Liquid/Vapor/Salt structs but\n"
             << "  *** UpdateCSMPVariables has not run — rl/rv/sl/sv/mul/muv on mesh\n"
             << "  *** are still FROM PREVIOUS TIMESTEP. Volume factors, CT, mobilities\n"
             << "  *** also not yet updated.\n";
    else if (post_csmp)
        cout << "  *** POST-UpdateCSMPVariables: EOS outputs (rl/rv/sl/sv/mul/muv)\n"
             << "  *** are current. Volume factors, CT and mobilities not yet computed.\n";
    else if (post_vf)
        cout << "  *** POST-VolumeFactorComputations: EOS outputs and volume factors\n"
             << "  *** are current. CT and mobilities not yet computed.\n";
    else if (mid_iter)
        cout << "  *** MID ROOT-FINDING TRIAL: EOS outputs reflect last BoundaryFlow\n"
             << "  *** trial, NOT the converged state. Conserved vars (mt, H, ma)\n"
             << "  *** reflect the current trial modification.\n";
    // post_prep: no warning — everything is current
    cout << "========================================================================\n";

    std::ios_base::fmtflags old_flags = cout.flags();
    cout << std::fixed << std::setprecision(8);

    // ── General node state ────────────────────────────────────────────────────
    cout << "--- General ---\n";
    cout << "Bulk.state = " << Bulk.state << "  (previous = " << previous_state << ")\n";
    cout << "t          = " << t_ << " oC  (tp = " << tp_ << ")\n";
    cout << "p          = " << p_ << " Pa";
    if (std::abs(p_current_ - p_) > 1.0)
        cout << "  p_current = " << p_current_ << " Pa  (mixture mode: vapor partial pressure)";
    cout << "\n";
    cout << "x          = " << x_ << "  wt = " << wt_ << " wt%\n\n";

    // ── Mass & volume ─────────────────────────────────────────────────────────
    cout << "--- Mass & Volume ---\n";
    cout << "phi        = " << phi() << "  phi_ (shaved) = " << phi_ << "\n";
    cout << "mt         = " << mt() << " kg/m3  (mtp = " << mtp() << ")\n";
    cout << "ml         = " << ml() << "  mv = " << mv() << "  mh = " << mh();
    if (with_air) cout << "  ma = " << ma();
    cout << "\n";
    cout << "m_fluid_   = " << m_fluid_ << " kg  (passed to Equilibrator)\n";
    if (with_air) cout << "m_air_     = " << m_air_ << " kg\n";
    cout << "m_rock_    = " << m_rock_ << " kg\n";
    cout << "H_current_ = " << H_current_ << " J  (H_previous = " << H_previous_ << ")\n";
    cout << "dml        = " << dml_ << "  dmv = " << dmv_;
    if (with_air) cout << "  dma = " << dma_;
    cout << "\n\n";

    // ── Salt bookkeeping ──────────────────────────────────────────────────────
    cout << "--- Salt ---\n";
    cout << "ms         = " << ms() << " kg/m3  (msp = " << msp() << ")\n";
    cout << "wt         = " << wt_ << " wt%  x = " << x_ << "\n";
    cout << "dxCl       = " << dxCl_ << "  dxCv = " << dxCv_ << "\n";
    cout << "xCl        = " << xCl() << "  xCv = " << xCv() << "\n\n";

    // ── EOS outputs — may be stale, see header warning ────────────────────────
    cout << "--- EOS outputs" << (pre_eos ? " [PREVIOUS TIMESTEP]" : mid_iter ? " [TRIAL STATE]" : "") << " ---\n";
    cout << "Bulk.rho   = " << Bulk.rho << " kg/m3  rho_bulk = " << rho_bulk() << "\n";
    cout << "rl         = " << rl() << "  rv = " << rv() << "  rh = " << rh();
    if (with_air) cout << "  ra = " << ra();
    cout << "\n";
    cout << "sl         = " << sl() << "  sv = " << sv() << "  sh = " << sh();
    if (with_air) cout << "  sa = " << sa();
    {
        double sum_sat = sl() + sv() + sh();
        if (with_air) sum_sat += sa();
        cout << "  sum = " << sum_sat << " (should be 1.0)";
    }
    cout << "\n";
    cout << "mul        = " << mul() << "  muv = " << muv();
    if (with_air) cout << "  mua = " << mua();
    cout << "\n";
    cout << "beta       = " << beta() << " 1/Pa\n";
    cout << "h_fluid_   = " << h_fluid_ << " J/kg  hf = " << hf() << "\n";
    cout << "hl         = " << hl() << "  hv = " << hv() << "  hh = " << hh();
    if (with_air) cout << "  ha = " << ha();
    cout << "\n";
    cout << "hCl        = " << hCl() << "  hCv = " << hCv();
    if (with_air) cout << "  hCa = " << hCa();
    cout << "\n\n";

    // ── Gas mixture (full_gas_mixture or hydraulic_gas_mixture) ─────────────
    if (with_air && (full_gas_mixture || hydraulic_gas_mixture))  //FIX 30-03-26
    {
        cout << "--- Gas Mixture [" << (full_gas_mixture ? "full/Dalton" : "hydraulic/p_total") << "] ---\n";  //FIX 30-03-26
        cout << "  p_partial_vapor = " << p_partial_vapor() << " Pa"
             << "  p_partial_air = " << p_partial_air() << " Pa\n";
        cout << "  xv_mix = " << xv_mix() << "  xa_mix = " << xa_mix() << "\n";
        cout << "  rho_gas_mix = " << rho_gas_mix() << " kg/m3"
             << "  mu_gas_mix = " << mu_gas_mix() << " Pa.s\n";
        cout << "  h_gas_mix = " << h_gas_mix() << " J/kg"
             << "  beta_gas_mix = " << beta_gas_mix() << " 1/Pa\n\n";
    }

    // ── Volume factors, nQ, CT, mobilities — may be stale, see header warning ─
    cout << "--- Volume factors / CT / Mobilities" << (pre_eos || mid_iter ? " [POTENTIALLY STALE]" : !post_prep ? " [MAY NOT BE UPDATED]" : "") << " ---\n";
    cout << "nQ             = " << nQ() << "\n";
    cout << "vol_factor_LHS = " << volume_factor_LHS << "\n";
    cout << "vol_factor_RHS = " << volume_factor_RHS << "\n";
    cout << "CT             = " << CT() << "\n";
    cout << "krl = " << krl << "  krv = " << krv;
    if (with_air) cout << "  kra = " << kra;
    cout << "  sum = " << (krl + krv + (with_air ? kra : 0.))
         << "  (1-sh = " << (1. - sh()) << ")\n";
    cout << "mml = " << mml() << "  mmv = " << mmv();
    if (with_air) cout << "  mma = " << mma();
    cout << "\n";
    cout << "========================================================================\n\n";

    cout.flags(old_flags);
}

template <uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::PAUSE()
{
    string dummy;
    cerr << " press ENTER to continue... " << endl;
    getline(cin, dummy);
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::CheckForOutOfRange( Node<dim>* n )
{
    // ── Transport enthalpy and salt fraction sanity ───────────────────────────
    // These fire before EOS — cerr summary only, no full state dump available.
    if (hCl() < 0.)
    {
        cerr << "\n[CheckForOutOfRange] Negative hCl=" << hCl()
             << "  node=(" << n->x() << ", " << n->y();
        if (dim == 3) cerr << ", "<< n->z();
        cerr<<")"
           << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
    }

    if (hCv() < 0.)
    {
        cerr << "\n[CheckForOutOfRange] Negative hCv=" << hCv()
             << "  node=(" << n->x() << ", " << n->y();
        if (dim == 3) cerr << ", "<< n->z();
        cerr<<")"
           << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
    }

    if (with_air && hCa() < 0.)
    {
        cerr << "\n[CheckForOutOfRange] Negative hCa=" << hCa()
             << "  node=(" << n->x() << ", " << n->y();
        if (dim == 3) cerr << ", "<< n->z();
        cerr<<")"
           << "  p=" << p() << " Pa  t=" << t() << " oC" << endl;
    }

    if (xCl() < 0.)
    {
        cerr << "\n[CheckForOutOfRange] Negative xCl=" << xCl()
             << "  node=(" << n->x() << ", " << n->y();
        if (dim == 3) cerr << ", "<< n->z();
        cerr<<")"
           << "  wt=" << wt() << " wt%" << endl;
    }

    if (xCv() < 0.)
    {
        cerr << "\n[CheckForOutOfRange] Negative xCv=" << xCv()
             << "  node=(" << n->x() << ", " << n->y();
        if (dim == 3) cerr << ", "<< n->z();
        cerr<<")"
           << "  wt=" << wt() << " wt%" << endl;
    }

    // ── EOS range and mass sanity — cerr summary + cout full state dump ───────
    if ( /*p() < 100000.0 || */t() < 5.0 )
    {
        cerr << endl << "[CheckForOutOfRange] p or t below EOS minimum";
        cerr << endl << "Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";
        DumpNodeState("CheckForOutOfRange");
        csmp_error.Note( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                         "\nPressure or temperature below minimum values of lookup table!");
    }

    if ( mt() < 0. )
    {
        cerr << endl << "[CheckForOutOfRange] Negative mt: " << mt();
        cerr << endl << " Node=(" << n->x() << ", " << n->y(); if (dim == 3) cerr << ", "<< n->z(); cerr<<")";

        DumpNodeState("CheckForOutOfRange");
        csmp_error.Note( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                         "\nmt() below zero!");
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::CheckBoundaryFlags( Node<dim>* n )
{
    if (t.Flag() == DIRICH && n->AtBoundary() == NOT)
    {
        fixed_temperature = true;
        t_fixed           = t();
    }

    else if (t.Flag() == DIRICH && n->AtBoundary() != NOT)
    {
        if (n->AtBoundary() != INTERNAL && open_boundaries) // if model boundaries are open to flow, t can change!
        {
            //cerr<<" Node is not at an internal boundary and external bounaries ar open, Temp is NOT fixed.";
            fixed_temperature = false;
        }
        else
        {
            fixed_temperature = true;
            t_fixed           = t();
        }
    }

    //Benoit add for OCEAN region which is not a boundary and where i want T ans salinity dirich
    else if (t.Flag() == DIRICH && wt.Flag() == DIRICH)
    {
        fixed_temperature = true;
        t_fixed           = t();
    }

    else
    {
        fixed_temperature = false;
    }
}

template<uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::OpenBoundaries(double ref_wt)
{
    open_boundaries    = true;
    ref_sal            = ref_wt * 0.01;

    cerr << "\n[OpenBoundaries, variable p] Boundary reference conditions:\n"
         << "  ref_wt           = " << ref_wt       << " wt%\n"
         << "  ref_sal          = " << ref_sal       << " [-]\n"
         << "  ref_humidity     = " << ref_humidity  << " [-]\n";

    // Reference enthalpies and ref_frac_vapor are computed per-node in
    // InitialPropertiesFromPTX from the local EOS state at each boundary node,
    // stored as nodal ScalarVariables, and read back at the start of each
    // BoundaryIteration call.  This gives the correct enthalpy for nodes at
    // different elevations / pressures without a single global reference T/P.
    // ref_humidity is set by ActivateAir and used per-node in InitialPropertiesFromPTX.
}

// PW May 2016: function to use improved convergence strategies for thermal equilibration
template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetConvergenceSpeedUpTo(bool boost)
{
    equilibrator.SetConvergenceSpeedUpTo(boost);
}

// PW Sept 2018
template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetAvoidPressureOscillationsAtVLHTo(bool avoid_p_VLH)
{
    avoid_pressure_oscillations_at_VLH = avoid_p_VLH;
}

// PW May 2018: function to use thermodynamic density in gravity term
template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetThermodynamicDensityInGravityTermTo(bool thermodynamic_density)
{
    thermodymamic_density_in_gravity_term = thermodynamic_density;
}

// PW May 2018: function to use thermodynamic density in gravity term
template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetThermodynamicDensityInCompressibilityTermTo(bool thermodynamic_density)
{
    thermodynamic_density_in_compressibility_term = thermodynamic_density;
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::WithRockLiquidusSolidus(bool with_rock_liquidus_solidus_)
{
    with_rock_liquidus_solidus = with_rock_liquidus_solidus_;
    if(with_rock_liquidus_solidus)
    {
        equilibrator.WithRockLiquidusSolidus(with_rock_liquidus_solidus);

        tl_key = pref.StorageKey("liquidus temperature");
        ts_key = pref.StorageKey("solidus temperature");

        //The hardcoded variables "liquidus temperature" and "solidus temperature" are nodal scalars that should be defined externally
        //Typically those temperatures are composition and pressure dependent.

        if (tl_key.type != SCALAR || tl_key.place != NODE)
            throw Exception(ERROR, "NaClH2OPropertiesVisitorPHX<dim>::WithRockLiquidusSolidus",
                            "liquidus temperature", " must be an nodal scalar property." );
        if (ts_key.type != SCALAR || ts_key.place != NODE)
            throw Exception(ERROR, "NaClH2OPropertiesVisitorPHX<dim>::WithRockLiquidusSolidus",
                            "solidus temperature", " must be an nodal scalar property." );
    }
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::ActivateAir(
        double humidity,
        bool treat_air_and_vapor_as_mixture,
        bool treat_air_and_vapor_as_hydraulic_mixture,
        double rain_mm_per_year)
{
    with_air = true;

    equilibrator.WithAirPhase(with_air);

    rain_mm_per_year_top = rain_mm_per_year;

    ref_humidity = std::max(0., std::min(1., humidity));

    full_gas_mixture = treat_air_and_vapor_as_mixture;

    // Full mixture mode requires with_air = true. Activates:
    //   - partial pressure for air density (p * xa)
    //   - partial pressure for Equilibrator (p * xv)
    //   - mixture viscosity for gas mobilities
    //   - mixture compressibility for CT
    //   - shared Darcy velocity for air and vapor transport
    // Need much more work for the thermodynamic side! Use hydraulic mixture mode for now!

    hydraulic_gas_mixture = treat_air_and_vapor_as_hydraulic_mixture;

    // ── Hydraulic gas mixture mode ─────────────────────────────────────────────
    // Vapor and air share the same Darcy velocity (coupled transport) but
    // thermodynamics remain fully decoupled:
    //   - No partial pressures: p_current_ = p_total; Equilibrator and air density
    //     both use total pore pressure
    //   - No moisture transport or condensation between phases
    //   - ComputeGasMixtureProperties called with p_total for both phases
    //     (all mixture diagnostics computed; ra unchanged at p_total)
    //   - Vapor handled by H2O-NaCl Equilibrator at full pore pressure
    //   - Air handled as ideal gas at full pore pressure
    // Transport quantities are coupled via ComputeGasMixtureProperties (called in
    // UpdateCSMPVariables) with p_partial_air = p_partial_vapor = p_total:
    //   rho_gas_mix = rv + ra         (both at p_total)
    //   mu_gas_mix  = volume-fraction weighted blend of muv and mua
    //   xv_mix, xa_mix, h_gas_mix, beta_gas_mix also populated for diagnostics
    // In VolumeFactorComputations, rv_transport and ra_transport are individually
    // rescaled so each phase gets the correct density-fraction mass flux allocation.
    // This is appropriate when vapor and air physically co-occupy the gas phase
    // and travel together, but thermodynamic interaction is negligible
    // (e.g. dry magmatic vapor meeting dry atmospheric air).

    if(full_gas_mixture && hydraulic_gas_mixture)
    {
        cerr<<endl<<"Both full_gas_mixture and hydraulic_gas_mixture are set to true, full_gas_mixture supersedes, disabling hydraulic_gas_mixture mode!";
        hydraulic_gas_mixture = false;
    }
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetBrooksCoreyKr(bool b, double lambda)
{
    use_brooks_corey_kr = b;
    if (b)
    {
        lambda_bc = lambda;
        cerr << "\n[SetBrooksCoreyKr] Brooks-Corey kr pair enabled"
             << "  lambda=" << lambda_bc
             << "  res_sl=" << res_sl
             << "  res_sv=" << res_sv
             << "  res_sa=" << res_sa
             << "\n  kr_l   = Se^" << ((2.0 + 3.0*lambda_bc)/lambda_bc)
             << "\n  kr_gas = (1-Se)^2 * (1 - Se^" << ((2.0 + lambda_bc)/lambda_bc) << ")"
             << endl;
    }
    else
    {
        cerr << "\n[SetBrooksCoreyKr] Brooks-Corey kr disabled"
             << "  (reverting to linear kr_l = Se, kr_gas = 1 - Se)"
             << endl;
    }
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetPrecipitationRecharge(Model<dim>& model)
{
    constexpr double rho_rain = 1000.;   // kg/m³
    const double rain_mass_flux = rain_mm_per_year_top * 1.e-3
            / (365.25 * 24. * 3600.)
            * rho_rain;   // [kg/m²/s]

    //Initialize
    model.InputPropertyValue("rain recharge top", ScalarVariable(ANY, 0.));

    Boundary<dim>& top = model.Boundary("TOP");

    for (auto eit = top.CellsBegin(); eit != top.CellsEnd(); eit++)
    {
        auto unrml = (*eit)->UnitNormal();
        double horiz_proj  = 0.;

        if constexpr(dim == 3U) horiz_proj = std::max(0., dotProduct(unrml, Point<3U>(0., 1., 0.)));
        if constexpr(dim == 2U) horiz_proj = std::max(0., dotProduct(unrml, Point<2U>(0., 1.)));


        //cerr<<endl<<"Cell: "<<(*eit)->Idx()<<", area: "<<(*eit)->Area()<<", horiz_proj: "<<horiz_proj;

        const double node_area = (*eit)->Area() / (*eit)->Nodes();
        for (auto i = 0U; i < (*eit)->Nodes(); i++)
        {
            (*eit)->N(i)->Read(cpv_key, cpv);
            (*eit)->N(i)->Read(rain_recharge_top_key, rain_recharge_top);

            rain_recharge_top() += (cpv() > 0.) ? rain_mass_flux * node_area * horiz_proj / cpv() : 0.;

            (*eit)->N(i)->Store(rain_recharge_top_key, rain_recharge_top);
        }
    }
    cerr << "\n[SetPrecipitationRecharge]"
         << "  rain=" << rain_mm_per_year_top << " mm/yr"
         << "  flux=" << rain_mass_flux << " kg/m2/s";
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::AddFluidContributionToHeatCapacity(bool add_fluid_contribution_to_heat_capacity_)
{
    add_fluid_contribution_to_heat_capacity = add_fluid_contribution_to_heat_capacity_;
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::AttemptToSurviveFluidPropertiesError(bool attempt_to_survive_fluid_properties_error_)
{
    equilibrator.AttemptToSurviveFluidPropertiesError(attempt_to_survive_fluid_properties_error_);
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetRockHeatCapacity( double mini_cp )
{
    rock.        SetRockHeatCapacity( mini_cp );
    equilibrator.SetRockHeatCapacity( mini_cp );
}

template< uint32_t dim>
void NaClH2OPropertiesVisitorPHX<dim>::SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                                                   double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_)
{
    rock.        SetRockCrystallizationCurve(nu_coefficient, sigma1_coefficient, latent_heat_of_fusion, b_coefficient, crystallization_curve_);
    equilibrator.SetRockCrystallizationCurve(nu_coefficient, sigma1_coefficient, latent_heat_of_fusion, b_coefficient, crystallization_curve_);
}

template class NaClH2OPropertiesVisitorPHX<1>;
template class NaClH2OPropertiesVisitorPHX<2>;
template class NaClH2OPropertiesVisitorPHX<3>;

}