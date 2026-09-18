// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONSERVATION_CHECK_H
#define CONSERVATION_CHECK_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"

#include <fstream>
#include <string>

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**
    @class ConservationCheck

    Global mass / salt / air / energy balance over the model, including the
    boundary flow that the equilibrator applies in BoundaryIteration.

    ── WHAT IT MEASURES ─────────────────────────────────────────────────────
    For each conserved quantity Q it forms

        residual(Q) = [ storage(Q, now) - storage(Q, previous pass) ]
                    + [ boundary out(Q) this pass ]

    and reports it both absolutely and relative to the current storage. With no
    other sources active, residual should be ~0 to round-off. A residual that
    GROWS with timestep size points at an operator-splitting error; one that is
    a fixed fraction of throughput points at a flux/clamp bug.

    ── UNITS (verified against NaClH2OPropertiesVisitorPHX) ─────────────────
      mt, ms, ma, bfm, bfs, bfa   per PORE volume   -> x pore_volume
      Htp, bfe                    per BULK volume   -> x bulk_volume
    Salt and energy are read from "previous mass salt" / "previous total
    enthalpy" because the equilibrator closes the step with msp() = ms() and
    Htp() = H_current_ — so once Apply() has returned, those ARE the current
    values, and no separate "mass salt" property exists.
    The equilibrator applies `mt() -= bfm()` and `H_current_ -= bfe()` in one
    block, so the per-step boundary flows line up with the storage change.
    Positive bf* means mass/energy LEAVING the model.

    ── WHAT IS *NOT* ACCOUNTED FOR ──────────────────────────────────────────
    These are real sources/sinks that this visitor does not see. If any is
    active, expect a residual and subtract it yourself via AddExternalSource():
      - well exchange (enters through mtp in ComputeWellSourceTerms)
      - PointInjection
      - rain / "fluid source wt" recharge
      - magmatic fluid production
      - rock <-> fluid energy is INTERNAL to Htp, so it cancels (by design)
      - split-boundary leakage is INTERNAL (mid receives minus what the matrix
        sides receive), so it must cancel — a residual localised at the fault
        would itself be a finding

    ── USAGE ────────────────────────────────────────────────────────────────
        ConservationCheck<DIM> cons(model, "run_label");
        // once, after the initial equilibration, to set the baseline:
        cons.BeginPass();  model.Accept(cons);  cons.EndPass(0., 0.);

        // then in the time loop, AFTER CVFEM_PHX->Apply() returns:
        cons.BeginPass();  model.Accept(cons);  cons.EndPass(model_time, dt);

    BeginPass/EndPass are explicit rather than relying on Visit(Model*) firing
    in a particular order.
 */
template<size_t dim>
class ConservationCheck : public Visitor<dim> {

public:
    ConservationCheck(Model<dim>& model,
                      const std::string& label = "conservation",
                      bool with_air = false);

    virtual ~ConservationCheck();

    virtual void Visit(Node<dim>*   node);
    virtual void Visit(Model<dim>*  m);
    virtual void Visit(Region<dim>* r);

    /// zero the per-pass accumulators; call immediately before model.Accept()
    void BeginPass();

    /// close the pass: form residuals, write a CSV row, optionally print
    void EndPass(double model_time, double dt);

    /// register a known external source (+) or sink (-) for this pass, in SI
    /// (kg for mass/salt/air, J for energy) so it can be subtracted out
    void AddExternalSource(double d_mass, double d_energy,
                           double d_salt = 0., double d_air = 0.);

    void SetVerbose(bool v) { verbose_ = v; }

    // ── accessors (current pass) ─────────────────────────────────────────
    double MassStorage()      const { return mass_now_; }
    double EnergyStorage()    const { return energy_now_; }
    double MassResidual()     const { return mass_res_; }
    double EnergyResidual()   const { return energy_res_; }
    double MassResidualRel()  const { return mass_res_rel_; }
    double EnergyResidualRel()const { return energy_res_rel_; }

    /// worst single node by |bfm| this pass — where the boundary is working
    size_t WorstBoundaryNode() const { return worst_bf_node_; }

private:
    void OpenLog();

    Model<dim>&                  model_ref_;
    const PropertyDatabase<dim>& prop_ref_;

    csmp::Index
        mt_key_, ms_key_, ma_key_, Htp_key_,
        bfm_key_, bfe_key_, bfs_key_, bfa_key_,
        pore_volume_key_, bulk_volume_key_;

    ScalarVariable
        mt_, ms_, ma_, Htp_,
        bfm_, bfe_, bfs_, bfa_,
        pore_volume_, bulk_volume_;

    bool with_air_, verbose_, have_baseline_;

    // storage, this pass
    double mass_now_, salt_now_, air_now_, energy_now_;
    // storage, previous pass
    double mass_prev_, salt_prev_, air_prev_, energy_prev_;
    // boundary out, this pass
    double bf_mass_, bf_energy_, bf_salt_, bf_air_;
    // boundary out, cumulative
    double bf_mass_cum_, bf_energy_cum_, bf_salt_cum_, bf_air_cum_;
    // externally registered sources, this pass
    double ext_mass_, ext_energy_, ext_salt_, ext_air_;
    // residuals
    double mass_res_, energy_res_, salt_res_, air_res_;
    double mass_res_rel_, energy_res_rel_;
    // diagnostics
    double worst_bf_abs_;
    size_t worst_bf_node_;
    size_t n_nodes_, pass_;

    std::ofstream log_;
    std::string   label_;
};

} // namespace csmp

#endif
