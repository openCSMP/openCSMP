// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ConservationCheck.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

using namespace std;

namespace csmp {

template<size_t dim>
ConservationCheck<dim>::ConservationCheck(Model<dim>& model,
                                          const std::string& label,
                                          bool with_air)
    : Visitor<dim>(MODEL, NODE),
      model_ref_(model),
      prop_ref_(model.Database()),

      mt_key_          (model.Database().StorageKey("fluid density")),
      // NOTE: there is no stored "mass salt". The equilibrator ends the step
      // with msp() = ms() (exactly as it does Htp() = H_current_), so after
      // Apply() returns the "previous" store holds THIS step's value.
      ms_key_          (model.Database().StorageKey("previous mass salt")),
      ma_key_          (model.Database().StorageKey("fluid mass air")),
      Htp_key_         (model.Database().StorageKey("previous total enthalpy")),

      bfm_key_         (model.Database().StorageKey("boundary flow mass")),
      bfe_key_         (model.Database().StorageKey("boundary flow enthalpy")),
      bfs_key_         (model.Database().StorageKey("boundary flow salt")),
      bfa_key_         (model.Database().StorageKey("boundary flow air")),

      pore_volume_key_ (model.Database().StorageKey("pore volume")),
      bulk_volume_key_ (model.Database().StorageKey("bulk volume")),

      with_air_(with_air), verbose_(true), have_baseline_(false),

      mass_now_(0.), salt_now_(0.), air_now_(0.), energy_now_(0.),
      mass_prev_(0.), salt_prev_(0.), air_prev_(0.), energy_prev_(0.),
      bf_mass_(0.), bf_energy_(0.), bf_salt_(0.), bf_air_(0.),
      bf_mass_cum_(0.), bf_energy_cum_(0.), bf_salt_cum_(0.), bf_air_cum_(0.),
      ext_mass_(0.), ext_energy_(0.), ext_salt_(0.), ext_air_(0.),
      mass_res_(0.), energy_res_(0.), salt_res_(0.), air_res_(0.),
      mass_res_rel_(0.), energy_res_rel_(0.),
      worst_bf_abs_(0.), worst_bf_node_(0), n_nodes_(0), pass_(0),

      label_(label)
{
    OpenLog();
}

template<size_t dim>
ConservationCheck<dim>::~ConservationCheck()
{
    if (log_.is_open()) log_.close();
}

template<size_t dim>
void ConservationCheck<dim>::OpenLog()
{
    time_t     now = time(nullptr);
    struct tm* lt  = localtime(&now);
    char stamp[32];
    strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", lt);

    const std::string fname = "conservation_" + label_ + "_" + stamp + ".csv";
    log_.open(fname.c_str());
    if (log_.is_open())
    {
        log_ << "pass,model_time_s,model_time_yr,dt_s,"
             << "mass_storage_kg,mass_d_storage_kg,mass_bf_out_kg,mass_ext_kg,"
             << "mass_residual_kg,mass_residual_rel,"
             << "energy_storage_J,energy_d_storage_J,energy_bf_out_J,energy_ext_J,"
             << "energy_residual_J,energy_residual_rel,"
             << "salt_storage_kg,salt_residual_kg,"
             << "air_storage_kg,air_residual_kg,"
             << "bf_mass_cum_kg,bf_energy_cum_J,worst_bf_node\n";
        log_.flush();
        cerr << "\n[ConservationCheck] logging to " << fname;
    }
    else
        cerr << "\n[ConservationCheck][WARNING] could not open " << fname;
}

template<size_t dim>
void ConservationCheck<dim>::Visit(Model<dim>*)  {}

template<size_t dim>
void ConservationCheck<dim>::Visit(Region<dim>*) {}

template<size_t dim>
void ConservationCheck<dim>::BeginPass()
{
    mass_now_ = salt_now_ = air_now_ = energy_now_ = 0.;
    bf_mass_  = bf_energy_ = bf_salt_ = bf_air_    = 0.;
    worst_bf_abs_  = 0.;
    worst_bf_node_ = 0;
    n_nodes_       = 0;
}

template<size_t dim>
void ConservationCheck<dim>::AddExternalSource(double d_mass, double d_energy,
                                               double d_salt, double d_air)
{
    ext_mass_   += d_mass;
    ext_energy_ += d_energy;
    ext_salt_   += d_salt;
    ext_air_    += d_air;
}

/** Accumulate one node's storage and its boundary flow.

    Positive bf* = leaving the model (the equilibrator does mt -= bfm,
    H_current_ -= bfe), so the balance is

        d(storage) + boundary_out - external = 0
*/
template<size_t dim>
void ConservationCheck<dim>::Visit(Node<dim>* node)
{
    node->Read(pore_volume_key_, pore_volume_);
    node->Read(bulk_volume_key_, bulk_volume_);

    const double pv = pore_volume_();
    const double bv = bulk_volume_();

    // ── storage ──────────────────────────────────────────────────────────
    node->Read(mt_key_,  mt_);
    node->Read(ms_key_,  ms_);
    node->Read(Htp_key_, Htp_);

    mass_now_   += mt_()  * pv;     // per pore volume  -> kg
    salt_now_   += ms_()  * pv;     // per pore volume  -> kg
    energy_now_ += Htp_() * bv;     // per bulk volume  -> J

    if (with_air_)
    {
        node->Read(ma_key_, ma_);
        air_now_ += ma_() * pv;
    }

    // ── boundary flow applied THIS step by BoundaryIteration ─────────────
    // Reset to 0 by the equilibrator on every visit, so these are per-step
    // and are zero on every non-boundary node.
    node->Read(bfm_key_, bfm_);
    node->Read(bfe_key_, bfe_);
    node->Read(bfs_key_, bfs_);

    bf_mass_   += bfm_() * pv;
    bf_energy_ += bfe_() * bv;
    bf_salt_   += bfs_() * pv;

    if (with_air_)
    {
        node->Read(bfa_key_, bfa_);
        bf_air_ += bfa_() * pv;
    }

    const double abs_bf = fabs(bfm_() * pv);
    if (abs_bf > worst_bf_abs_)
    {
        worst_bf_abs_  = abs_bf;
        worst_bf_node_ = node->Idx();
    }

    ++n_nodes_;
}

template<size_t dim>
void ConservationCheck<dim>::EndPass(double model_time, double dt)
{
    const double d_mass   = mass_now_   - mass_prev_;
    const double d_salt   = salt_now_   - salt_prev_;
    const double d_air    = air_now_    - air_prev_;
    const double d_energy = energy_now_ - energy_prev_;

    if (have_baseline_)
    {
        // storage gained + what left through the boundary - what was injected
        mass_res_   = d_mass   + bf_mass_   - ext_mass_;
        salt_res_   = d_salt   + bf_salt_   - ext_salt_;
        air_res_    = d_air    + bf_air_    - ext_air_;
        energy_res_ = d_energy + bf_energy_ - ext_energy_;

        // scale-free forms: relative to storage, which is the only scale that
        // means the same thing on a 0.05 m3 matrix node and a fault node
        mass_res_rel_   = (mass_now_   != 0.) ? mass_res_   / fabs(mass_now_)   : 0.;
        energy_res_rel_ = (energy_now_ != 0.) ? energy_res_ / fabs(energy_now_) : 0.;
    }
    else
    {
        mass_res_ = salt_res_ = air_res_ = energy_res_ = 0.;
        mass_res_rel_ = energy_res_rel_ = 0.;
        have_baseline_ = true;
    }

    bf_mass_cum_   += bf_mass_;
    bf_energy_cum_ += bf_energy_;
    bf_salt_cum_   += bf_salt_;
    bf_air_cum_    += bf_air_;

    if (log_.is_open())
    {
        log_ << pass_ << ","
             << setprecision(12) << model_time << ","
             << model_time / 3.1557600e7 << ","
             << dt << ","
             << setprecision(10)
             << mass_now_   << "," << d_mass   << "," << bf_mass_   << "," << ext_mass_   << ","
             << mass_res_   << "," << mass_res_rel_   << ","
             << energy_now_ << "," << d_energy << "," << bf_energy_ << "," << ext_energy_ << ","
             << energy_res_ << "," << energy_res_rel_ << ","
             << salt_now_   << "," << salt_res_ << ","
             << air_now_    << "," << air_res_  << ","
             << bf_mass_cum_ << "," << bf_energy_cum_ << ","
             << worst_bf_node_ << "\n";
        log_.flush();
    }

    if (verbose_)
    {
        cerr << "\n[Conservation] pass " << pass_
             << "  nodes=" << n_nodes_
             << "\n    mass   : storage " << mass_now_   << " kg"
             << "  d " << d_mass   << "  bf_out " << bf_mass_
             << "  RESID " << mass_res_   << "  (rel " << mass_res_rel_   << ")"
             << "\n    energy : storage " << energy_now_ << " J"
             << "  d " << d_energy << "  bf_out " << bf_energy_
             << "  RESID " << energy_res_ << "  (rel " << energy_res_rel_ << ")";
        if (worst_bf_abs_ > 0.)
            cerr << "\n    worst boundary node " << worst_bf_node_
                 << "  |bfm| " << worst_bf_abs_ << " kg";
    }

    mass_prev_   = mass_now_;
    salt_prev_   = salt_now_;
    air_prev_    = air_now_;
    energy_prev_ = energy_now_;

    ext_mass_ = ext_energy_ = ext_salt_ = ext_air_ = 0.;
    ++pass_;
}

template class ConservationCheck<1U>;
template class ConservationCheck<2U>;
template class ConservationCheck<3U>;

} // namespace csmp
