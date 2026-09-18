// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ThreePhaseTransportPHX.h"
#include "Model.h"
#include "Region.h"
#include "ExplicitFiniteVolumeTransportPHX.h"
#include "UpwindControlVisitor.h"
#include "SplitRegionLeakage.h"
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

namespace csmp {

// ── Transport sub-cycling toggles (all in one place) ─────────────────────────
// AdvectMassConserved sub-cycles internally so a per-cell over-drain sizes an
// internal sub-step instead of cutting the outer timestep. Toggle/tune here.
//   tp_full_subcycling: false -> LEGACY: CheckDry over-drain returns an outer cut
//                     (control_dt*factor), no internal sub-stepping.
//   tp_verbose      : full per-substep trace vs one summary line per call.
//   tp_max_substeps : N. Every sub-step advances by at least control_dt/N (the
//                     SUBSTEP FLOOR), so the stage completes in <= N sub-steps by
//                     construction. N is an accuracy dial, not just a cap: a cell
//                     whose safe drain step is shorter than control_dt/N is
//                     clamped instead of resolved exactly (see AdvectMassConserved
//                     header), with a front-position error ~1/N per such cell.
//                     Raising N costs nothing on healthy steps — the floor only
//                     binds when a cell demands a smaller sub-step.

static constexpr bool   tp_verbose      = false;

// DRY-OUT REDO (tp_dryout_stage_end) — REINSTATED 2026-07 after validation
// runs showed subcycling can walk a stage into states the frozen field cannot
// recover from (cells fully drained -> negative pressures -> tiny-dt spiral).
// STAGE_END semantics (the validated ~13-day variant): after the sub-cycled
// stage completes, any cell that was wet at stage start and is STILL dry at
// stage close (i.e. every refill opportunity failed) rejects the stage:
// reset + retry at 0.5 x control_dt. Transient dries that the relay refills
// are forgiven — only proven-permanent holes cut. "Newly dry" = above
// DRYOUT_END_MAX at stage start, below it now (eternal dust cannot trigger).
// Universal: no node exclusions (Dirichlet vents included).
static constexpr bool   tp_dryout_stage_end = false;

// DRY-OUT DETECTOR (tp_dryout_detect): the SAME wet-at-start / dry-at-close test
// as above, but it only MEASURES — no cut, no redo, no dt change. Sweeps every
// node and phase (the redo path stops at the first hit, so it can never say how
// many are involved) and reports a count plus the worst offender as a
// stage_events_ CSV row.
//
// SCOPE: Region("Model"), which INCLUDES the split-boundary in/mid/out nodes —
// CheckDry iterates the same region and finds them in the leakage drain map, and
// UpdateTransportVariables stores to the very keys read here ("fluid mass
// liquid"/"vapor"). So a mid node drained by ExchangeSubstep counts exactly like
// one drained by advection. The snapshot is taken before UpdateProjection and
// PrepareStageRates, so m0 is true stage-start state.
//
// Costs one nodal sweep per phase per STAGE (not per beat, unlike the ledger).
static constexpr bool   tp_dryout_detect    = true;

static constexpr double DRYOUT_START_MIN = 1.0e-9;  // was wet at stage start
//static constexpr double DRYOUT_START_MIN = 1.0;  // was DEFINITELY wet at stage start
static constexpr double DRYOUT_END_MAX   = 1.0e-9;  // is dry at stage close


static constexpr int    tp_max_substeps = 100;

// ══════════════════════════════════════════════════════════════════════════
//  MASTER SWITCH — tp_full_subcycling
// ══════════════════════════════════════════════════════════════════════════
//  true  (NEW scheme): transport is sub-cycled AND the split-boundary (SB)
//        leakage runs INSIDE the sub-cycle loop, one exchange per beat, on the
//        same frozen-velocity contract as the facet fluxes, sized by the same
//        CheckDry, kept positive by per-beat budget clamps. No dt cuts from
//        transport or leakage; bounded exits only (tp_max_substeps,
//  false (LEGACY): no sub-cycling; transport takes the full step and, on
//        over-drain, requests an outer dt CUT (old behavior). Leakage runs
//        one-shot in the SCHEME after transport and likewise CUTS on
//        over-drain (mandatory) and on the CFL thresholds (optional, see
//        sb_drain_threshold / sb_stuff_threshold). This is the regression
//        reference against the pre-sub-cycling code.
//  MESHES WITHOUT A SPLIT BOUNDARY: nothing to configure. SetInterfaceExchange
//  is never called, leakage_ stays nullptr, every SB clause above is inert,
//  and 'true' simply means sub-cycled transport (the config banner confirms:
//  "leakage attached: no").

static constexpr bool   tp_full_subcycling = true;

// CLAMP LEDGER (tp_clamp_ledger): per-stage measurement of clamp parking.
// Per sub-step, the ledger sums WITHHELD = requested-but-clamped outflow over
// all cells and phases; that sum IS the queue standing at clamped cells at
// that instant. Reported at stage end: carryover (the LAST sub-step's queue —
// the parking that never recovered, the honest per-stage error proxy) and
// max_queue (the worst transient queue). Units: per-pore-volume, kg per m3 of
// pore space — the LHS units — so ratios are meaningful, absolute values are
// relative monitors, not literal kilograms. Costs one CheckDry-sized sweep
// per sub-step. Measurement only; no dt control is wired to it.

static constexpr bool   tp_clamp_ledger  = false;

/** custom constructor */
template<uint32_t dim>
ThreePhaseTransportPHX<dim>::ThreePhaseTransportPHX( Model<dim>& model,
                                                    UpwindControlVisitor<dim>& upwind_visitor,
                                                    ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                                                    ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                                                    ExplicitFiniteVolumeTransportPHX<dim>& fv_air)
    : model_ref( model ),
    UpwindVisitor( upwind_visitor ),
    fv_transport_vapor( fv_vapor ),
    fv_transport_liquid( fv_liquid ),
    fv_transport_air( fv_air ),
    with_air_(false),//default, needs to be activated via ActivateAir(), otherwise can function as TwoPhaseTransportPHX without air
    control_dt(0.0),
    upper_shell_T(450.),
    lower_shell_T(350.)
{
    Sf_key   = model.Database().StorageKey( "salt flux"                          );
    ff_key   = model.Database().StorageKey( "fluid flux"                         );
    ffi_key  = model.Database().StorageKey( "fluid flux in"                      );
    Ef_key   = model.Database().StorageKey( "energy flux"                        );
    mff_key  = model.Database().StorageKey( "magmatic fluid flux"                );
    smff_key = model.Database().StorageKey( "shell magmatic fluid flux"          );

    SfI_key   = model.Database().StorageKey( "salt flux integral"                );
    ffI_key   = model.Database().StorageKey( "fluid flux integral"               );
    EfI_key   = model.Database().StorageKey( "energy flux integral"              );
    mffI_key  = model.Database().StorageKey( "magmatic fluid flux integral"      );
    smffI_key = model.Database().StorageKey( "shell magmatic fluid flux integral");

    pvn_key  = model.Database().StorageKey( "nodal pore velocity"                );
    cffd_key = model.Database().StorageKey( "fluid flux with direction integral" );
    ffd_key  = model.Database().StorageKey( "fluid flux with direction"          );

    T_key    = model.Database().StorageKey( "temperature"                        );
    // (pv_key removed: it was added for an early CheckDry-voting design that
    //  ended up using per-pore-volume rates directly; never read. The unused
    //  pv_key / pv members can also be dropped from the header.)

    cout << endl << "ThreePhaseTransportPHX<" << typeid(double).name() << "," << dim;
    cout << ">: Constructed successfully." << endl;
}

/** default destructor */
template<uint32_t dim>
ThreePhaseTransportPHX<dim>::~ThreePhaseTransportPHX()
{
}

/** set maximum time step for all three phases */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::SetLargestTimeStep( const double& max_time_step )
{
    fv_transport_liquid.SetMaximumTimeStep( max_time_step );
    fv_transport_vapor .SetMaximumTimeStep( max_time_step );
    if(with_air_) fv_transport_air   .SetMaximumTimeStep( max_time_step );
}

/** main function to coordinate three-phase advection */
template<uint32_t dim>
bool ThreePhaseTransportPHX<dim>::HandlesLeakage() const
{
    // true only in the new scheme; legacy leaves leakage to the scheme (cuts)
    return leakage_ != nullptr && tp_full_subcycling;
}

// ============================================================================
//  AdvectMassConserved — explicit FV transport over the full requested dt,
//  made mass-safe by INTERNAL sub-cycling instead of cutting the outer step.
//
//  Per-cell over-drain (CheckDry factor<1) sizes an internal sub-step rather
//  than escalating to an outer reset + pressure re-solve. We advance by the safe
//  sub-step, re-read state, repeat until the sub-steps sum to control_dt. The
//  outer dt is then governed by CFL alone; CheckDry only sizes sub-steps.
//  (Set tp_full_subcycling=false at the top of the file for the legacy outer-cut path.)
//
//  Mass-conservative to machine precision, two combined guarantees:
//    (1) positivity: AdjustFluxOut clamps each cell's outflow to its own holdup
//        (read before inflow), so no cell exports more than it has.
//    (2) balance: a clamped cell's shortfall scales the matching downstream
//        inflow (mass_balance_vectors / CalculateFluxIn), so neighbours receive
//        exactly what was exported. Each sub-step is a full balanced advection.
//
//  SUBSTEP FLOOR (termination + accuracy model). Sizing every sub-step so the
//  worst cell drains EXACTLY to empty does not terminate: a small cell refilled
//  by its neighbours re-binds with a geometrically decaying factor (a Zeno
//  sequence). So sub-steps are floored at control_dt/tp_max_substeps. A cell
//  whose safe factor falls below the floor no longer binds: its REQUESTED
//  outflow exceeds its holdup within the sub-step, AdjustFluxOut clamps the
//  export to the holdup (the cell simply drains empty), and the mass_balance
//  factor scales the matching downstream inflow. Consequences:
//    * Mass: exact at any floor (positivity + balance are untouched).
//    * Completion: <= tp_max_substeps sub-steps, guaranteed — each advances by at
//      least control_dt/tp_max_substeps or finishes the remainder.
//    * Distribution: a clamped cell exports one sub-step LATE — inflow arriving
//      during sub-step j parks and is exported in j+1 (it IS the holdup then).
//      Parking does not accumulate; per clamped cell the residual is the final
//      sub-step's inflow, plus one sub-step of transit delay per clamped cell a
//      front crosses. Concretely: the floor is control_dt/tp_max_substeps, so one
//      sub-step's worth of flux is ~1/tp_max_substeps of that cell's throughput
//      over the stage — with tp_max_substeps=100, each clamped cell's residual is
//      bounded by ~1% of what it transported this stage.
//  Cells with drain time > control_dt/tp_max_substeps are still resolved exactly.
//
//  MASS_FLOOR vs SUBSTEP FLOOR — complementary, both required. MASS_FLOOR (see
//  CheckDry) excludes dust cells (holdup <= 1e-3) from sizing so they cannot
//  drag the factor DOWN to the substep floor and force all tp_max_substeps
//  sub-steps every stage; the substep floor stops real-but-small cells from
//  grinding BELOW it and guarantees termination. MASS_FLOOR must stay far below
//  any real holdup (see CheckDry for the through-route argument).
//
//  FROZEN RATES + CONTENT REBUILD (no toggle; the only path). The MASS flux
//  rates are assembled once per stage (DetermineFacetFluxRatesOnce) — their
//  ingredients (pressure-projected velocities and rhs mobilities) are written
//  once per outer step, so re-assembling per sub-step would produce identical
//  mass fluxes anyway. Content fluxes (enthalpy, salt, ...) are NOT kept
//  frozen: RescaleFacetFluxToTimestep rebuilds them every sub-step from the
//  mass facet fluxes times the donor's CURRENT content/mass ratio, so every
//  content stays an exact scaled copy of the mass field as holdups evolve —
//  this is what keeps mass and energy/salt from decoupling at clamped cells.
//
//  Termination: the substep floor guarantees completion in <= tp_max_substeps
//  where the seconds floor exceeds the substep floor); both return a clean
//  outer cut, never a tiny dt.
// ============================================================================
template<uint32_t dim>
double ThreePhaseTransportPHX<dim>::AdvectMassConserved( const double& time_increment )
{
    // Sub-cycling toggles are file-scope (top of file) and used directly:

    // ── one-time config banner + sanity checks (first call only) ──────────
    // ════════════════════════════════════════════════════════════════════════
    //  §1  ONE-TIME CONFIG BANNER + SANITY GUARDS
    // ════════════════════════════════════════════════════════════════════════
    static bool config_announced = false;
    if ( !config_announced )
    {
        config_announced = true;
        cerr << "\n========================= TRANSPORT / SB CONFIG ========================="
             << "\n  tp_full_subcycling : " << ( tp_full_subcycling ? "TRUE  (subcycled transport + SB leakage IN LOOP, clamps, no cuts)"
                                                                   : "FALSE (LEGACY: no subcycling, one-shot SB leakage, cuts on overdrain/CFL)" )
             << "\n  leakage attached   : " << ( leakage_ != nullptr ? "yes" : "no" )
             << "  -> transport " << ( ( leakage_ && tp_full_subcycling ) ? "RUNS leakage per beat"
                                                                       : "leaves leakage to the scheme" )
             << "\n  tp_clamp_ledger    : " << ( tp_clamp_ledger ? "on" : "off" )
             << "\n  tp_verbose         : " << ( tp_verbose ? "on" : "off" );

        // ── nonsensical-configuration guards ─────────────────────────────
        // With the toggle collapse these should be unreachable, but assert the
        // invariants so a future edit that reintroduces a bad combo fails loud.
        if ( tp_max_substeps < 1 )
            throw csmp::Exception( ERROR, "ThreePhaseTransportPHX::AdvectMassConserved",
                                  "tp_max_substeps < 1 — the substep floor is undefined." );
        cerr << "\n========================================================================\n";
    }

    // ════════════════════════════════════════════════════════════════════════
    //  §2  STAGE SETUP — reset state, freeze rates, capture interface
    //      stage rates
    // ════════════════════════════════════════════════════════════════════════
    control_dt = time_increment;
    transport_cut_reason_.clear();        // reason channel for transport-initiated
    stage_events_.clear();                // event channel (rebuilt in §5)

    // dry-out support: stage-start per-phase holdup snapshots (model reads —
    // transport property vectors are stale before the first Rescale).
    const size_t n_model_nodes = model_ref.Region("Model").Nodes();
    const Index dry_key_l = model_ref.Database().StorageKey( "fluid mass liquid" );
    const Index dry_key_v = model_ref.Database().StorageKey( "fluid mass vapor"  );
    const Index dry_key_a = with_air_ ? model_ref.Database().StorageKey( "fluid mass air" )
                                      : dry_key_l;   // unused when !with_air_
    std::vector<double> dry0_l, dry0_v, dry0_a;
    if ( tp_dryout_detect || tp_dryout_stage_end )
    {
        dry0_l.resize( n_model_nodes ); dry0_v.resize( n_model_nodes );
        if ( with_air_ ) dry0_a.resize( n_model_nodes );
        for ( size_t k = 0; k < n_model_nodes; ++k )
        {
            dry0_l[k] = model_ref.Region("Model").N(k)->Read( dry_key_l );
            dry0_v[k] = model_ref.Region("Model").N(k)->Read( dry_key_v );
            if ( with_air_ ) dry0_a[k] = model_ref.Region("Model").N(k)->Read( dry_key_a );
        }
    }
    // dry-out detector results (filled after the loop, read by the §5 report)
    size_t      dryout_n           = 0;     // node-phase pairs wet at start, dry at close
    size_t      dryout_worst_idx   = 0;
    const char* dryout_worst_phase = "";
    double      dryout_worst_m0    = 0.;    // how wet the worst one was at stage start

    // cuts (set by the §4 bounded exits)
    constexpr double NEGLIGIBLE_REMAINING_TIME = 1.0e-6; // absolute time remaining [s] below which the stage is done


    UpdateProjection();

    DetermineFacetFluxRatesOnce();   // frozen mass rates for the whole stage


    // ── leakage-in-the-loop stage setup ──────────────────────────────────
    // New scheme only: capture frozen interface rates once, then exchange per
    // beat inside the loop. Legacy (tp_full_subcycling=false) leaves leakage to
    // the scheme, which runs the one-shot cutting pass after transport.
    const bool run_leakage_in_loop = ( leakage_ != nullptr ) && tp_full_subcycling;
    if ( run_leakage_in_loop )
    {
        leakage_->ResetInterfaceLedger();
        leakage_->PrepareStageRates();   // frozen interface rates + drain-rate map
        cerr << "\n[LeakLoop] stage rates captured: triples="
             << leakage_->StageRateCount()
             << "  voting nodes=" << leakage_->NodeDrainRatesPerPV().size();
    }
    leak_vote_dt = 0.;                   // set per sub-step when leakage votes in CheckDry

    double t_done       = 0.0;
    int    substeps     = 0;
    double first_factor = 1.0;

    // ── clamp ledger (tp_clamp_ledger) — per-substep standing queue ────────
    // WITHHELD per sub-step (requested - applied outflow, all cells+phases) =
    // mass parked at clamped cells. Units per-PV (kg/m3 pore space); read
    // ratios, not absolutes. Stage-end report:
    //   carryover = LAST sub-step's queue -> parking that never recovered
    //               (the error proxy, an upper bound);
    //   max_queue = worst transient (max >> carryover: drained, GOOD;
    //               carryover ~ max: persisted, BAD);
    //   clamped-cell-substeps = grind gauge (floor-cadence volume).
    double ledger_queue_last = 0., ledger_queue_max = 0., ledger_advected = 0.;
    bool   stage_clamped = false;   // did ANY sub-step clamp/floor this stage?
    // gates the flush split (no queue -> no flush)
    double ledger_worst_cell = 0.; size_t ledger_worst_idx = 0; const char* ledger_worst_phase = "";
    size_t ledger_clamped_cell_substeps = 0;

    auto queue_measure = [&]( double factor )
    {
        double substep_withheld = 0.;
        auto sweep_phase = [&]( auto& tp, const char* ph )
        {
            for ( typename vector<Node<dim>*>::const_iterator
                     fvit = model_ref.Region("Model").NodesBegin();
                 fvit != model_ref.Region("Model").NodesEnd(); fvit++ )
            {
                const size_t idx = (*fvit)->Idx();
                const double req = tp.GetFluxOut( idx ) * factor;
                if ( req <= 0. ) continue;
                const double avail   = std::max( 0., tp.GetMainPropertyLHS( idx ) );
                const double applied = std::min( req, avail );
                ledger_advected += applied;
                const double w = req - applied;
                if ( w > 0. )
                {
                    substep_withheld += w;
                    ++ledger_clamped_cell_substeps;
                    if ( w > ledger_worst_cell )
                    { ledger_worst_cell = w; ledger_worst_idx = idx; ledger_worst_phase = ph; }
                }
            }
        };

        sweep_phase( fv_transport_liquid, "liquid" );
        sweep_phase( fv_transport_vapor,  "vapor"  );
        if ( with_air_ ) sweep_phase( fv_transport_air, "air" );

        ledger_queue_last = substep_withheld;   // overwritten each sub-step; the
        // stage-end survivor IS the carryover
        if ( substep_withheld > ledger_queue_max ) ledger_queue_max = substep_withheld;
    };

    // Loop while there is still a physically meaningful remainder
    // ════════════════════════════════════════════════════════════════════════
    //  §3  SUB-CYCLE LOOP — per beat: rescale -> size (CheckDry) -> floors
    //      -> flush split -> advect -> interface beat -> ledger
    // ════════════════════════════════════════════════════════════════════════
    while ( (control_dt - t_done) > NEGLIGIBLE_REMAINING_TIME && substeps < tp_max_substeps )
    {
        const double dt_remaining = control_dt - t_done;

        RescaleFacetFluxes( dt_remaining );   // mass: frozen rates scaled;
        // contents: rebuilt at current ratios

        // leakage voting: CheckDry adds the frozen interface drain (per-PV,
        // scaled to dt_remaining) to the transport outflow of interface nodes
        leak_vote_dt = run_leakage_in_loop ? dt_remaining : 0.;

        // RESOLVABLE-ONLY SIZING: a cell whose own safe factor is already
        // below the sub-step floor is clamped NO MATTER what sub_dt we pick —
        // letting it vote in the sizing minimum only drags every sub-step to
        // floor cadence (the observed 20-80-substep grind on steady fault
        // cells) at zero accuracy benefit: it is clamped either way, and the
        // ledger accounts its parking. So CheckDry ignores sub-floor cells
        // when SIZING; they still transport, clamp, and appear in the ledger.
        // Legacy passes 0 (all cells vote) so its cut stays faithful.
        // To A/B against the old sizing behavior, set this to 0.
        sizing_floor_factor_ = tp_full_subcycling
                                   ? std::min( 1.0, ( control_dt / tp_max_substeps ) / dt_remaining )
                                   : 0.;
        CheckDryFiniteVolumesAndAdjustTimestep();

        if ( substeps == 0 )
        {
            first_factor = time_step_factor;
            if ( tp_verbose && time_step_factor < 1.0 )
                cerr << "\n[SubCycle] CheckDry would have cut: factor=" << time_step_factor
                     << "  phase=" << binding_phase_ << "  node=" << binding_idx_
                     << "  LHS=" << binding_lhs_
                     << "  (old code -> outer reset to dt=" << control_dt * time_step_factor << ")";
        }

        // LEGACY (tp_full_subcycling=false): on the first over-drain, hand the cut
        // straight to the outer scheme (control_dt*factor) instead of sub-stepping.
        // Returns exactly what the pre-sub-cycling code returned. No-op when factor==1.
        if ( !tp_full_subcycling && time_step_factor < 1.0 )
        {
            if ( tp_verbose )
                cerr << "\n[SubCycle] legacy (no subcycling) -> outer cut dt="
                     << control_dt * time_step_factor << "\n";
            return control_dt * time_step_factor;
        }

        // SUBSTEP FLOOR — never advance by less than control_dt/tp_max_substeps.
        // Raising time_step_factor (rather than sub_dt alone) is what keeps this
        // safe: PerformFacetFluxes scales the fluxes by the SAME factor, so the
        // floored cell's requested outflow matches the time actually advanced,
        // and AdjustFluxOut's holdup clamp + mass_balance absorb the over-drain.
        // first_factor (the outer-cut fallback) is captured above, pre-floor.
        const double min_sub_dt      = control_dt / tp_max_substeps;
        const double unfloored_factor = time_step_factor;
        if ( unfloored_factor < 1.0 ) stage_clamped = true;   // this stage built a queue
        if ( time_step_factor * dt_remaining < min_sub_dt )
            time_step_factor = std::min( 1.0, min_sub_dt / dt_remaining );

        double sub_dt = dt_remaining * time_step_factor;

        // ── FLUSH SPLIT — never end the stage on a LARGE sub-step ──────────
        // Carryover (the error that survives the stage) is whatever the FINAL
        // sub-step withholds, which scales with that sub-step's size: a big
        // snap-to-finish beat parks a lot with no later beat to relay it. So
        // when this beat would finish the stage and the remainder is bigger
        // than two floors, take (remainder - one floor) now and leave exactly
        // one floor-sized finishing beat: the queue gets one extra fine-
        // resolution relay pass, and the stage always closes on a small beat.
        // Cost: at most one extra sub-step per stage.
        // GATE: only when this stage actually CLAMPED at some point — i.e. it
        // built a standing queue that a fine finisher can relay. A stage that
        // sails through unclamped (factor=1 throughout) has zero carryover
        // regardless, so splitting its single finishing beat into 99%+1% would
        // add a sub-step for nothing (the observed quiescent-stage case).
        if ( stage_clamped
            && sub_dt >= dt_remaining * (1. - 1.e-12)
            && dt_remaining > 2. * min_sub_dt )
        {
            time_step_factor = ( dt_remaining - min_sub_dt ) / dt_remaining;
            sub_dt = dt_remaining * time_step_factor;
            if ( tp_verbose )
                cerr << "\n[SubCycle] flush split: finishing beat capped, floor-sized finisher follows";
        }


        // Safe to advance — measure this sub-step's standing queue
        if ( tp_clamp_ledger ) queue_measure( time_step_factor );
        PerformFacetFluxes();
        TrackFluxes();

        // one interface relay beat: frozen rates x sub_dt against the state
        // transport just wrote; per-beat budgets own positivity
        if ( run_leakage_in_loop )
            leakage_->ExchangeSubstep( sub_dt );

        t_done += sub_dt;

        ++substeps;

        if ( tp_verbose && unfloored_factor < 1.0 && binding_transport_ )
        {
            const double m_old = binding_lhs_;
            const double out   =  binding_transport_->GetFluxOut( binding_idx_ );
            const double in    = -binding_transport_->GetFluxIn(  binding_idx_ );
            const double m_new =  binding_transport_->GetMainPropertyLHS( binding_idx_ );

            cerr << "\n[SubCycle]   step " << substeps
                 << ": sub_dt=" << sub_dt << "  factor=" << time_step_factor;
            if ( unfloored_factor < time_step_factor )
                cerr << " (FLOORED from " << unfloored_factor << ")";
            cerr << "  node=" << binding_idx_ << " (" << binding_phase_ << ")"
                 << "  m_old=" << m_old << "  out=" << out << "  in=" << in << "  m_new=" << m_new;
            // m_new above is post-TRANSPORT (property vectors); the leakage
            // beat has since written the MODEL. Read it so this line chains
            // with the next CheckDry LHS (which sees the post-beat state).
            if ( run_leakage_in_loop )
            {
                const char* bname = ( binding_transport_ == &fv_transport_liquid ) ? "fluid mass liquid"
                                    : ( binding_transport_ == &fv_transport_vapor  ) ? "fluid mass vapor"
                                                                                  :                                                  "fluid mass air";
                const Index bkey = model_ref.Database().StorageKey( bname );
                cerr << "  after_beat=" << model_ref.Region("Model").N( binding_idx_ )->Read( bkey );
            }
            cerr << "  (" << 100.0 * t_done / control_dt << "% done)\n";
        }
        else if ( tp_verbose && unfloored_factor >= 1.0 )   // factor = 1 → no cut needed
        {
            cerr << "\n[SubCycle]   step " << substeps
                 << ": sub_dt=" << sub_dt << "  factor=1.0  "
                 << ( last_floored_cells_ > 0
                         ? "(no RESOLVABLE over-drain, finishing; drained cells excluded from sizing)  "
                         : "(no over‑drain, finishing)  " )
                 << 100.0 * t_done / control_dt << "% done\n";
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    //  §4  BOUNDED EXIT — substep cap hit: hand a clean, known-completable
    //      cut to the scheme (reset-retry). With the substep floor each beat
    //      advances >= control_dt/tp_max_substeps, so this is reachable only
    //      via floating-point dust or the flush split consuming the last
    //      slot — rare and benign (one extra pressure solve). No transport-
    //      side abort: if dt collapses, the scheme's own minimum-dt abort is
    //      the single authority (seconds floor, early-stop flag, and both
    //      transport FATALs retired as redundant with it).
    // ════════════════════════════════════════════════════════════════════════
    if ( control_dt - t_done > NEGLIGIBLE_REMAINING_TIME )
    {
        // Retry at t_done (a dt we KNOW completed safely), 0.99 margin for the
        // pressure field shifting on the reset; first-beat failure falls back
        // to the stage's first CheckDry factor (the old code's cut).
        const double ret_dt = ( t_done > 0.0 ) ? t_done * 0.99
                                             : control_dt * first_factor;
        {
            std::ostringstream r;
            r << "transport bounded exit (cap_hit): " << substeps << " substeps, "
              << 100.0 * t_done / control_dt << "% done; retry at " << ret_dt;
            transport_cut_reason_ = r.str();
        }
        cerr << "\n[SubCycle][CUT] " << transport_cut_reason_ << "\n";
        return ret_dt;
    }

    // ── DRY-OUT: cells wet at stage start, still dry at stage close ────────
    // ONE sweep serves both consumers. Unlike the old redo-only version it does
    // NOT stop at the first hit, so dryout_n is a real count. Detection is
    // unchanged: m0 >= DRYOUT_START_MIN and now < DRYOUT_END_MAX. Recovery
    // WITHIN the stage is forgiven — only the stage-close value is read.
    // Includes split-boundary in/mid/out nodes (see tp_dryout_detect above).
    if ( tp_dryout_detect || tp_dryout_stage_end )
    {
        auto sweep_dry = [&]( const Index key, std::vector<double>& m0, const char* ph )
        {
            for ( size_t k = 0; k < n_model_nodes; ++k )
            {
                if ( m0[k] < DRYOUT_START_MIN ) continue;               // never wet
                if ( model_ref.Region("Model").N(k)->Read( key ) >= DRYOUT_END_MAX )
                    continue;                                           // recovered
                ++dryout_n;
                if ( m0[k] > dryout_worst_m0 )                          // wettest one lost
                { dryout_worst_m0 = m0[k]; dryout_worst_idx = k; dryout_worst_phase = ph; }
            }
        };
        sweep_dry( dry_key_l, dry0_l, "liquid" );
        sweep_dry( dry_key_v, dry0_v, "vapor"  );
        if ( with_air_ ) sweep_dry( dry_key_a, dry0_a, "air" );

        // The REDO is a separate decision from the measurement and stays off by
        // default: it returns a blind 0.95 factor that the scheme scales again by
        // 0.95 (~0.9025 effective) — a ratchet, not a sized cut. It cannot size
        // one, because it only knows the stage-close state, not when the node dried.
        if ( tp_dryout_stage_end && dryout_n > 0 )
        {
            std::ostringstream r;
            r << "stage DRY-OUT redo: node " << dryout_worst_idx << " ("
              << dryout_worst_phase << ") m " << dryout_worst_m0
              << " -> dry at stage close, never refilled; " << dryout_n
              << " node-phase pair(s) affected; retry dt=" << 0.95 * control_dt;
            transport_cut_reason_ = r.str();
            cerr << "\n[SubCycle][dryout] " << transport_cut_reason_ << "\n";
            return 0.95 * control_dt;    // scheme resets + retries (existing path)
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    //  §5  STAGE REPORTS — parking ledger (transport + interface)
    // ════════════════════════════════════════════════════════════════════════
    // ── clamp ledger report ────────────────────────────────────────────────
    // See the queue_measure block above for what carryover / max_queue mean
    // and how to read them (max >> carryover = recovered; ~equal = persisted).
    if ( tp_clamp_ledger && ledger_queue_max > 0. )
        cerr << "[SubCycle] parked [transport]: carryover=" << ledger_queue_last
             << "  (carryover/advected=" << ( ledger_advected > 0. ? ledger_queue_last / ledger_advected : 0. ) << ")"
             << "  max_queue=" << ledger_queue_max
             << "  clamped-cell-substeps=" << ledger_clamped_cell_substeps
             << "  worst single withhold: node=" << ledger_worst_idx << " (" << ledger_worst_phase
             << ") w=" << ledger_worst_cell << "\n";

    // interface (leakage): same two numbers for the fault's budget clamp.
    // NOTE: this line being ABSENT from a log is itself information — it means
    // the per-beat budgets never withheld anything, which is the EXPECTED
    // steady state under the velocity freeze: a leg's request scales with the
    // donor's fresh holdup (vel * m_donor * sub_dt), so it exceeds the
    // (1 - sb_drain_margin) budget only when vel*sub_dt > (1 - margin) — i.e.
    // only for genuinely fast events (a flash filling the fault), not for
    // ordinary draining. Silence = the fault never parked this stage.
    if ( tp_clamp_ledger && run_leakage_in_loop )
    {
        const auto& il = leakage_->GetInterfaceLedger();
        if ( il.parked_max > 0. )
            cerr << "[SubCycle] parked [interface]: carryover=" << il.parked_now
                 << "  (carryover/applied=" << ( il.applied > 0. ? il.parked_now / il.applied : 0. ) << ")"
                 << "  max_queue=" << il.parked_max
                 << "  worst single withhold: mid=" << il.worst_idx << " w=" << il.worst << "\n";
    }

    // ── stage events -> CSV (via the scheme, same channel style as cuts) ────
    // Noteworthy-but-not-cut conditions worth a spreadsheet row. Thresholds
    // deliberately loose: rows should be RARE (events, not telemetry).
    //
    // NOTE: this block used to sit INSIDE the (tp_clamp_ledger && run_leakage_in_loop)
    // guard above, so with the ledger off — as now — or on a mesh with no split
    // boundary, stage_events_ was never even CLEARED, and a stale string from an
    // earlier stage could ride into a later CSV row. Only the interface clause
    // belongs behind those guards.
    stage_events_.clear();
    {
        std::ostringstream ev;

        if ( tp_clamp_ledger )
        {
            const double co_ratio = ( ledger_advected > 0. )
            ? ledger_queue_last / ledger_advected : 0.;
            if ( co_ratio > 1.e-3 )
                ev << "carryover/advected=" << co_ratio << "; ";

            if ( run_leakage_in_loop )
            {
                const auto& il = leakage_->GetInterfaceLedger();
                if ( il.parked_now > 0. )
                    ev << "interface carryover=" << il.parked_now << "; ";
            }
        }

        // Cells wet at stage start and still dry at stage close — the "drained
        // and never came back" population, matrix and fault nodes alike.
        if ( dryout_n > 0 )
            ev << "dried-out=" << dryout_n << " node-phase (worst node "
               << dryout_worst_idx << " " << dryout_worst_phase
               << ", m0=" << dryout_worst_m0 << "); ";

        if ( substeps >= tp_max_substeps / 2 )
            ev << "substeps=" << substeps << " (>= half cap); ";

        if ( !ev.str().empty() )
            stage_events_ = "stage events: " + ev.str();
    }

    return control_dt;
}

// (DetermineFacetFluxes(dt) removed: the per-sub-step full re-assembly path is
//  gone — its inputs (velocities, rhs mobilities) are stage-frozen model
//  variables, so it produced byte-identical MASS fluxes to the rates-once +
//  rescale path while LACKING the per-sub-step content rebuild that keeps
//  energy/salt consistent with mass. Drop the declaration from the header too.)

template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::CheckDryFiniteVolumesAndAdjustTimestep()
{
    // Local logging switch — keep in sync with the one in AdvectMassConserved.
    // verbose = print the per-call binding diagnostic below; non-verbose = silent
    // (AdvectMassConserved prints only the one-line per-call summary).
    constexpr bool verbose = false;

    time_step_factor = 1.0;

    // Min factor over ALL cells with outflow — MASS_FLOOR-excluded dust
    // included (their factor is max(0,LHS)/Outflow). This does NOT size the
    // sub-step; it answers one question for the growth cap in
    // AdvectMassConserved: "would any cell clamp at the candidate sub-step?"
    // (time_step_factor_all < candidate <=> yes).
    time_step_factor_all = 1.0;

    // A cell whose holdup is below MASS_FLOOR is treated as drained: it has
    // essentially nothing left to export (AdjustFluxOut clamps its outflow to its
    // holdup, i.e. to ~zero), so it CANNOT go negative and must not be allowed to
    // constrain the sub-step. Without this floor a cell drained to floating-point
    // dust (~1e-17) divided by a finite *requested* outflow yields a spurious
    // near-zero factor, which drives sub_dt -> 0 and wedges the sub-cycle.
    //
    // IMPORTANT — the floor only affects SIZING, not transport. A floored cell
    // still drains normally in PerformFacetFluxes; it is merely removed from the
    // sub_dt minimum. Mass is therefore conserved exactly whether or not a cell
    // is floored (positivity from AdjustFluxOut's clamp, balance from the
    // mass_balance factor that scales the matching inflow).
    //
    // WHY THE FLOOR MUST STAY SMALL. The floor is safe ONLY for cells whose
    // holdup is essentially all the mass they would move this step (dust cells,
    // no through-flux). A LARGE floor would also skip THROUGH-ROUTE cells — cells
    // that receive inflow and re-drain repeatedly within one control_dt, moving
    // several times their instantaneous holdup. Sub-cycling resolves those
    // (drain, refill, drain, ...); a single clamped step would cap them at one
    // holdup, and the mass_balance factor would throttle the matching inflow.
    // That stays mass-conserving but advances the front too slowly — a
    // distribution error. So MASS_FLOOR sits in the wide gap between numerical
    // dust (~1e-16) and the smallest physically meaningful vapor holdup
    // (~1e-3..1e-1 here): far above noise, far below any real, transportable mass.
    //
    // EDGE CASE — a near-empty cell ABOUT TO RECEIVE A LARGE INFLOW. Suppose a
    // cell holds, say, 5e-4 (below MASS_FLOOR, so it does NOT bind the sub-step),
    // but this step a large inflow is heading into it. Is it safe to have skipped
    // it in sizing? Yes — and here is exactly why:
    //
    //   * The floor only removes the cell from the sub_dt MINIMUM. It does NOT
    //     stop the cell from receiving inflow or from draining. PerformFacetFluxes
    //     runs over every cell regardless of floor status.
    //
    //   * Within the sub-step (which may be large, since this cell didn't bind),
    //     AdjustFluxOut clamps the cell's OUTFLOW to its PRE-inflow holdup (5e-4).
    //     So it can export at most 5e-4 this step — it CANNOT over-drain, no
    //     matter how big sub_dt is or that it wasn't sized for.
    //
    //   * ComposeAdvection then ADDS the large inflow and SUBTRACTS the clamped
    //     (<=5e-4) outflow. The cell ends the step at roughly inflow - 5e-4, i.e.
    //     it FILLS UP. The inflow is not lost; it accumulates in the holdup.
    //
    //   * CheckDry re-runs every sub-step. On the NEXT sub-step the cell now holds
    //     well above MASS_FLOOR, is no longer floored, and binds normally if its
    //     outflow demands it — so it is sized correctly the moment its mass is
    //     real.
    //
    //   Net: the floored-then-flooded cell over-drains nothing (clamp), loses
    //   nothing (inflow accumulates, mass_balance keeps neighbours consistent),
    //   and is throttled by at most ONE sub-step of delay before it is sized
    //   properly. Mass-conservative to machine precision; at worst a one-sub-step
    //   deferral of throughput, never an error.
    constexpr double MASS_FLOOR = 1.e-3;

    // Track the most restrictive node for diagnostics
    double worst_factor = 1.0;
    double worst_LHS = 0., worst_Outflow = 0., worst_leak = 0.;
    const char* worst_phase = nullptr;
    double worst_x = 0., worst_y = 0., worst_z = 0.;

    double LHS, Outflow, temp_factor;
    size_t floored_cells = 0;   // cells skipped as drained (diagnostic only)
    size_t subfloor_cells = 0;  // cells below the sizing floor (clamped regardless; no sizing vote)

    binding_transport_ = nullptr;
    binding_phase_     = nullptr;
    binding_idx_       = 0;
    binding_lhs_       = 0.;

    //TEST ONLY!!!
    double force_control_dt_factor = 1.;

    for ( typename vector<Node<dim>*>::const_iterator
             fvit = model_ref.Region("Model").NodesBegin();
         fvit != model_ref.Region("Model").NodesEnd(); fvit++ )
    {
        const size_t idx = (*fvit)->Idx();

        // ── leakage voting (per_substep mode) ───────────────────────────────
        // Interface nodes' sizing must see the FULL demand on their holdup:
        // transport outflow PLUS the frozen leakage drain (per-PV, scaled to
        // this sub-step's dt_remaining). Everything downstream is inherited:
        // MASS_FLOOR excludes interface dust from sizing (budgets clamp it),
        // the substep floor bounds interface tails, the ledger counts them.
        // leak_vote_dt == 0 outside per_substep mode -> zero overhead.
        double leak_l = 0., leak_v = 0., leak_a = 0.;
        if ( leak_vote_dt > 0. && leakage_ != nullptr )
        {
            const auto& drains = leakage_->NodeDrainRatesPerPV();
            const auto  it     = drains.find( *fvit );
            if ( it != drains.end() )
            {
                leak_l = it->second.l * leak_vote_dt;
                leak_v = it->second.v * leak_vote_dt;
                leak_a = it->second.a * leak_vote_dt;
            }
        }

        // ── liquid ─────────────────────────────────────────────────────────
        LHS     = fv_transport_liquid.GetMainPropertyLHS(idx);
        Outflow = fv_transport_liquid.GetFluxOut(idx) + leak_l;
        if ( Outflow > 0. && LHS <= MASS_FLOOR ) ++floored_cells;
        if ( Outflow > 0. )
            time_step_factor_all = std::min( time_step_factor_all, std::max(0., LHS) / Outflow );
        temp_factor = ( Outflow > 0. && LHS > MASS_FLOOR ) ? LHS / Outflow * force_control_dt_factor : 1.;
        if ( temp_factor < sizing_floor_factor_ ) ++subfloor_cells;   // clamped regardless — no sizing vote
        else if ( temp_factor < time_step_factor )
        {
            time_step_factor = temp_factor;
            worst_factor = temp_factor;
            worst_LHS = LHS; worst_Outflow = Outflow; worst_leak = leak_l;
            worst_phase = "liquid";
            worst_x = (*fvit)->x(); worst_y = (*fvit)->y(); worst_z = (*fvit)->z();
            binding_transport_ = &fv_transport_liquid;
            binding_idx_       = idx;
            binding_lhs_       = LHS;
            binding_phase_     = "liquid";
        }

        // ── vapor ──────────────────────────────────────────────────────────
        LHS     = fv_transport_vapor.GetMainPropertyLHS(idx);
        Outflow = fv_transport_vapor.GetFluxOut(idx) + leak_v;
        if ( Outflow > 0. && LHS <= MASS_FLOOR ) ++floored_cells;
        if ( Outflow > 0. )
            time_step_factor_all = std::min( time_step_factor_all, std::max(0., LHS) / Outflow );
        temp_factor = ( Outflow > 0. && LHS > MASS_FLOOR ) ? LHS / Outflow * force_control_dt_factor : 1.;
        if ( temp_factor < sizing_floor_factor_ ) ++subfloor_cells;   // clamped regardless — no sizing vote
        else if ( temp_factor < time_step_factor )
        {
            time_step_factor = temp_factor;
            worst_factor = temp_factor;
            worst_LHS = LHS; worst_Outflow = Outflow; worst_leak = leak_v;
            worst_phase = "vapor";
            worst_x = (*fvit)->x(); worst_y = (*fvit)->y(); worst_z = (*fvit)->z();
            binding_transport_ = &fv_transport_vapor;
            binding_idx_       = idx;
            binding_lhs_       = LHS;
            binding_phase_     = "vapor";
        }

        if(with_air_)
        {
            // ── air ────────────────────────────────────────────────────────────
            LHS     = fv_transport_air.GetMainPropertyLHS(idx);
            Outflow = fv_transport_air.GetFluxOut(idx) + leak_a;
            if ( Outflow > 0. && LHS <= MASS_FLOOR ) ++floored_cells;
            if ( Outflow > 0. )
                time_step_factor_all = std::min( time_step_factor_all, std::max(0., LHS) / Outflow );
            temp_factor = ( Outflow > 0. && LHS > MASS_FLOOR ) ? LHS / Outflow : 1.;
            if ( temp_factor < sizing_floor_factor_ ) ++subfloor_cells;   // clamped regardless — no sizing vote
            else if ( temp_factor < time_step_factor )
            {
                time_step_factor = temp_factor;
                worst_factor = temp_factor;
                worst_LHS = LHS; worst_Outflow = Outflow; worst_leak = leak_a;
                worst_phase = "air";
                worst_x = (*fvit)->x(); worst_y = (*fvit)->y(); worst_z = (*fvit)->z();
                binding_transport_ = &fv_transport_air;
                binding_idx_       = idx;
                binding_lhs_       = LHS;
                binding_phase_     = "air";
            }
        }
    }

    last_floored_cells_ = floored_cells;   // for the step line's honesty

    // Report if a significant timestep cut was triggered.
    // Outflow is split into its two demands: the transport facet outflow and
    // the frozen leakage vote (fault drain on interface nodes) — the leak part
    // is applied by the BEAT (after the step line), which is why a node's next
    // LHS can be far below the step line's m_new. Node = id + coordinates.
    if ( tp_verbose && time_step_factor < 0.99999 )
        cerr << "\n[CheckDryFiniteVolumes] Timestep cut: factor=" << time_step_factor
             << "  phase=" << worst_phase
             << "  LHS=" << worst_LHS
             << "  Outflow=" << worst_Outflow
             << " (transport " << ( worst_Outflow - worst_leak )
             << " + leak " << worst_leak << ")"
             << "  node=" << binding_idx_
             << " (" << worst_x << ", " << worst_y << ", " << worst_z << ")"
             << "  dt_in=" << control_dt
             << "  dt_out=" << control_dt * time_step_factor
             << "  (floored " << floored_cells << " drained, " << subfloor_cells << " sub-floor cells)" << endl;



}

/** apply time-step factor and execute advection for all three phases */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::PerformFacetFluxes()
{
    fv_transport_liquid.AdjustAndPerformFacetFlux( time_step_factor );
    fv_transport_vapor .AdjustAndPerformFacetFlux( time_step_factor );
    if(with_air_) fv_transport_air   .AdjustAndPerformFacetFlux( time_step_factor );
}

/** activate gravity for all three phases */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::WithGravityComponentAllPhases()
{
    fv_transport_liquid.WithGravityComponent();
    fv_transport_vapor .WithGravityComponent();
    if(with_air_) fv_transport_air   .WithGravityComponent();
}

/** project velocities onto facet normals for all phases */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::UpdateProjection()
{
    fv_transport_liquid.UpdateProjection();
    fv_transport_vapor .UpdateProjection();
    if(with_air_) fv_transport_air   .UpdateProjection();
}

/** accumulate per-timestep fluxes into integral trackers */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::TrackFluxes()
{
    // Hard-coded property indices inside fv_transport_liquid / vapor / air:
    //   [0] = mass
    //   [1] = energy (enthalpy)
    //   [2] = salt                   (liquid + vapor only — air has no salt)
    //   [3] = magmatic fluid mass    (liquid + vapor only — air has no magmatic mass)
    //   [4] = magmatic salt          (liquid + vapor only)
    constexpr unsigned int IDX_ENERGY         = 1;
    constexpr unsigned int IDX_SALT           = 2;
    constexpr unsigned int IDX_MAGMATIC_MASS  = 3;

    for ( typename vector<Node<dim>*>::const_iterator
             fvit = model_ref.Region("Model").NodesBegin();
         fvit != model_ref.Region("Model").NodesEnd(); fvit++ )
    {
        const size_t idx = (*fvit)->Idx();

        (*fvit)->Read( T_key,     temperature );
        (*fvit)->Read( SfI_key,   SfI         );
        (*fvit)->Read( EfI_key,   EfI         );
        (*fvit)->Read( ffI_key,   ffI         );
        (*fvit)->Read( mffI_key,  mffI        );
        (*fvit)->Read( smffI_key, smffI       );
        (*fvit)->Read( pvn_key,   pvn         );
        (*fvit)->Read( cffd_key,  cffd        );

        Sf()   = 0.;
        ff()   = 0.;  ffi() = 0.;
        Ef()   = 0.;
        mff()  = 0.;
        smff() = 0.;

        // ── Salt flux: liquid + vapor only (air carries no salt) ──────────
        Sf += fv_transport_liquid.GetFluxOut(idx, IDX_SALT);
        Sf += fv_transport_vapor .GetFluxOut(idx, IDX_SALT);

        // ── Total fluid mass flux: all three phases ───────────────────────
        ff += fv_transport_liquid.GetFluxOut(idx);
        ff += fv_transport_vapor .GetFluxOut(idx);
        if(with_air_) ff += fv_transport_air   .GetFluxOut(idx);   // air mass flux

        ffi += fv_transport_liquid.GetFluxIn(idx);
        ffi += fv_transport_vapor .GetFluxIn(idx);
        if(with_air_) ffi += fv_transport_air   .GetFluxIn(idx);   // air mass flux in

        // ── Energy flux: all three phases ─────────────────────────────────
        Ef += fv_transport_liquid.GetFluxOut(idx, IDX_ENERGY);
        Ef += fv_transport_vapor .GetFluxOut(idx, IDX_ENERGY);
        if(with_air_) Ef += fv_transport_air   .GetFluxOut(idx, IDX_ENERGY);   // air enthalpy flux

        // ── Magmatic fluid flux: liquid + vapor only ──────────────────────
        mff += fv_transport_liquid.GetFluxOut(idx, IDX_MAGMATIC_MASS);
        mff += fv_transport_vapor .GetFluxOut(idx, IDX_MAGMATIC_MASS);

        if ( temperature() < upper_shell_T && temperature() > lower_shell_T )
        {
            smff += fv_transport_liquid.GetFluxOut(idx, IDX_MAGMATIC_MASS);
            smff += fv_transport_vapor .GetFluxOut(idx, IDX_MAGMATIC_MASS);
        }

        // ── Directional fluid flux ────────────────────────────────────────
        double norm_pv = sqrt( pvn(0)*pvn(0) + pvn(1)*pvn(1) + pvn(2)*pvn(2) );
        ffd(0) = ffd(1) = ffd(2) = 0.;

        if ( norm_pv > 0. )
        {
            ffd(0) = pvn(0) / norm_pv * ff();
            ffd(1) = pvn(1) / norm_pv * ff();
            ffd(2) = pvn(2) / norm_pv * ff();

            cffd(0) += ffd(0);
            cffd(1) += ffd(1);
            cffd(2) += ffd(2);
        }

        // ── Update integral accumulators ──────────────────────────────────
        SfI   += Sf();
        ffI   += ff();
        EfI   += Ef();
        mffI  += mff();
        smffI += smff();

        (*fvit)->Store( Sf_key,   Sf   );
        (*fvit)->Store( Ef_key,   Ef   );
        (*fvit)->Store( ff_key,   ff   );
        (*fvit)->Store( ffi_key,  ffi  );
        (*fvit)->Store( mff_key,  mff  );
        (*fvit)->Store( smff_key, smff );
        (*fvit)->Store( SfI_key,   SfI   );
        (*fvit)->Store( EfI_key,   EfI   );
        (*fvit)->Store( ffI_key,   ffI   );
        (*fvit)->Store( mffI_key,  mffI  );
        (*fvit)->Store( smffI_key, smffI );
        (*fvit)->Store( cffd_key,  cffd  );
        (*fvit)->Store( ffd_key,   ffd   );
    }
}

//Benoit add
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::DetermineFacetFluxRatesOnce()
{
    fv_transport_liquid.DetermineFacetFluxRatesOnce( UpwindVisitor.UpwindMatrices( fv_transport_liquid.GetDensityKey()) );
    fv_transport_vapor .DetermineFacetFluxRatesOnce( UpwindVisitor.UpwindMatrices( fv_transport_vapor .GetDensityKey()) );
    if(with_air_) fv_transport_air.DetermineFacetFluxRatesOnce( UpwindVisitor.UpwindMatrices( fv_transport_air.GetDensityKey()) );
}

template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::RescaleFacetFluxes(double dt)
{
    fv_transport_liquid.RescaleFacetFluxToTimestep(dt);
    fv_transport_vapor .RescaleFacetFluxToTimestep(dt);
    if(with_air_) fv_transport_air.RescaleFacetFluxToTimestep(dt);
}

/** add advection variable carried by H2O-NaCl phases (liquid + vapor) only.
    Air does not carry salt, magmatic mass, lithium, or tracer species. */
template<uint32_t dim>
void ThreePhaseTransportPHX<dim>::AddAdvectionVariable( const char* new_lhs_liquid, const char* new_rhs_liquid,
                                                       const char* new_lhs_vapor,  const char* new_rhs_vapor )
{
    fv_transport_liquid.AddAdvectionVariable( new_lhs_liquid, new_rhs_liquid );
    fv_transport_vapor .AddAdvectionVariable( new_lhs_vapor,  new_rhs_vapor  );
    // fv_transport_air: intentionally not extended — air carries no H2O-NaCl
    // dissolved species (salt, magmatic mass, lithium, tracer).
}


template class ThreePhaseTransportPHX<1U>;
template class ThreePhaseTransportPHX<2U>;
template class ThreePhaseTransportPHX<3U>;

} // end namespace csmp