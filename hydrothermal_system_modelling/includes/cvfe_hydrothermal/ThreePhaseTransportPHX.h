// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef THREE_PHASE_TRANSPORT_PHX_H
#define THREE_PHASE_TRANSPORT_PHX_H

/*   Changelog
     February 2015, Philipp Weis:
     - initial two-phase port to CSMP++.
     2025, Benoit LC:
     - extended to three phases by adding an independent air phase.
     - renamed TwoPhaseTransportPHX → ThreePhaseTransportPHX.
*/

#include "ExplicitFiniteVolumeTransportPHX.h"
#include "UpwindControlVisitor.h"
#include <string>

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class UpwindControlVisitor;
template<uint32_t> class ExplicitFiniteVolumeTransportPHX;
template<uint32_t> class SplitRegionLeakage;

/**
     @class ThreePhaseTransportPHX ThreePhaseTransportPHX.h

     @author Philipp Weis (original), Benoit LC (air phase)
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @section motivation Motivation
      Orchestrator for three-phase (liquid / vapor / air) FV advection within
      the CVFEM scheme (Weis et al., Geofluids, 2014). AdvectMassConserved
      sub-cycles transport (and SB leakage, per beat) inside one outer pressure
      step on stage-frozen velocities: stability comes from per-beat drain
      clamps, accuracy from the scheme-side pore-velocity CFL.

     @section usage Usage
      Replaces TwoPhaseTransportPHX.  Accepts a third
      ExplicitFiniteVolumeTransportPHX instance for the air phase.  Air
      contributes to fluid-mass flux (ff) and energy flux (Ef) but NOT to salt
      flux (Sf) or magmatic-fluid flux (mff/smff) — air carries neither salt
      nor magmatic mass.

     @section dependencies Dependencies
      Needs ExplicitFiniteVolumeTransportPHX and UpwindControlVisitor.
  */

template<uint32_t dim>
class ThreePhaseTransportPHX {
public:

    ThreePhaseTransportPHX( Model<dim>& model,
                           UpwindControlVisitor<dim>& upwind_visitor,
                           ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                           ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                           ExplicitFiniteVolumeTransportPHX<dim>& fv_air/*,
                           bool with_air = false*/);

    ~ThreePhaseTransportPHX();

    double AdvectMassConserved( const double& time_increment );

    // ── Split-boundary leakage handoff ───────────────────────────────────
    // Give transport a pointer to the SB leakage. Whether transport RUNS it
    // (new scheme, per beat) or leaves it to the scheme (legacy, one-shot with
    // cuts) is decided by tp_full_subcycling (file-scope, .cpp) and reported by
    // HandlesLeakage(): true  -> transport runs leakage per beat, scheme skips;
    //                    false -> scheme runs the one-shot cutting leakage.
    void SetInterfaceExchange( SplitRegionLeakage<dim>* leakage ) { leakage_ = leakage; }
    bool HandlesLeakage() const;   // true when the scheme must NOT call leakage itself
    void   SetLargestTimeStep( const double& max_time_step );
    void   WithGravityComponentAllPhases();
    void   UpdateProjection();
    void   AddAdvectionVariable( const char* new_lhs_liquid, const char* new_rhs_liquid,
                              const char* new_lhs_vapor,  const char* new_rhs_vapor );
    void   ActivateAir() { with_air_ = true; }

    /** reason string for a transport-initiated reset-cut (the bounded exits:
        substep cap hit / seconds-floor stop); empty when the returned dt was
        not a transport rejection. The scheme forwards it to the dt-cut CSV. */
    const std::string& TransportCutReason() const { return transport_cut_reason_; }
    /** noteworthy-but-not-cut stage conditions (large carryover, interface
        carryover, near-cap substep count); empty on unremarkable stages.
        The scheme logs it as a non-cut row in the dt-cut CSV. */
    const std::string& StageEvents() const { return stage_events_; }

private:

    Model<dim>&                            model_ref;
    UpwindControlVisitor<dim>&             UpwindVisitor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_vapor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_liquid;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_air;

    bool with_air_;

    void CheckDryFiniteVolumesAndAdjustTimestep();
    void PerformFacetFluxes();
    void TrackFluxes();

    //Benoit add
    void DetermineFacetFluxRatesOnce();
    void RescaleFacetFluxes(double dt);

    double control_dt, time_step_factor, time_step_factor_all;
    double leak_vote_dt = 0.;   // dt_remaining for CheckDry's leakage voting (0 = voting off)
    SplitRegionLeakage<dim>* leakage_ = nullptr;
    double upper_shell_T, lower_shell_T;

    // flux tracking variables
    csmp::Index  Sf_key,   ff_key,   ffi_key,  Ef_key,
        SfI_key,  ffI_key,  EfI_key,
        mff_key,  mffI_key,
        smff_key, smffI_key,
        T_key,
        pvn_key,  cffd_key, ffd_key;

    ScalarVariable       Sf, ff, ffi, Ef, SfI, ffI, EfI, mff, mffI, smff, smffI;
    ScalarVariable       temperature;
    VectorVariable<dim>  pvn, cffd, ffd;

    // CheckDry's binding cell — the node/phase whose drain budget produced the
    // smallest safe factor this pass; consumed by the per-substep log lines.
    size_t                                  binding_idx_;
    double                                  binding_lhs_;
    ExplicitFiniteVolumeTransportPHX<dim>*  binding_transport_;
    const char*                             binding_phase_;
    // Resolvable-only sizing: cells whose own safe factor is already below
    // this floor cannot be helped by a smaller sub-step (they clamp anyway),
    // so they do not vote in CheckDry's sizing minimum. Set per call by the
    // sub-cycle loop (control_dt-relative); 0 = legacy-faithful (all vote).
    double sizing_floor_factor_ = 0.;
    size_t last_floored_cells_ = 0;   // cells excluded from sizing as drained (log honesty)

    std::string transport_cut_reason_;
    std::string stage_events_;
};

} // end namespace csmp

#endif