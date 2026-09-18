// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef UPWIND_CONTROL_VISITOR_H
#define UPWIND_CONTROL_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "ExplicitFiniteVolumeTransportPHX.h"

/*   Changelog
     February 2015, Philipp Weis:
     - initial port to CSMP++ and strong simplification as compared to the csmp5-version.
     2025, Benoit LC:
     - added air as a third independent phase (fv_transport_air reference,
       three-way velocity dispatch, "upwind control air" key).
     September 2026, JK / Benoit LC:
     - velocity output of lower-dimensional elements: gravity projected onto the
       element plane (visualisation only).
*/

namespace csmp {

/**
     @class UpwindControlVisitor UpwindControlVisitor.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes

     @section motivation Motivation
      Specialised visitor to define upwind nodes for the three-phase CVFEM scheme.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
     Calculates upwind nodes and stores information in a per-phase matrix.
     Also calculates velocity and pore-velocity vectors for visualisation.
     @endcode

     @section dependencies Dependencies
     The functions are tailored for use in ExplicitFiniteVolumeTransportPHX.
     This visitor is needed in ThreePhaseTransportPHX and upwind PDE operators.

     @section issues Known issues

     @section testing Testing
     Testing was done in the period before publication in 2014; air phase
     added 2026.

  */

template<uint32_t dim>
class UpwindControlVisitor : public Visitor<dim> {
public:
    // Full constructor — with air (fv_air provided explicitly)
    UpwindControlVisitor(Model<dim>& model,
                         ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                         ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                         ExplicitFiniteVolumeTransportPHX<dim>* fv_air,
                         const char* permeability,
                         const char* porosity,
                         std::vector<std::string>& densities,
                         std::vector<std::string>& relperm_vis,
                         std::vector<std::string>& saturations,
                         std::vector<std::string>& cfl_variables,
                         std::vector<std::string>& velocities,
                         std::vector<std::string>& pore_velocities);

    // No-air constructor — fv_air omitted entirely
    UpwindControlVisitor(Model<dim>& model,
                         ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                         ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                         const char* permeability,
                         const char* porosity,
                         std::vector<std::string>& densities,
                         std::vector<std::string>& relperm_vis,
                         std::vector<std::string>& saturations,
                         std::vector<std::string>& cfl_variables,
                         std::vector<std::string>& velocities,
                         std::vector<std::string>& pore_velocities);

    ~UpwindControlVisitor();
    void Verbose();

    virtual void Visit(Model<dim>* m);
    virtual void Visit(Region<dim>* r);
    virtual void Visit(Element<dim>* e);

    /** Which donor map a pass writes. Selected implicitly: Reset() -> Pressure,
        WithVelocity() -> Transport. See WithVelocity() for the full rationale. */
    enum DonorMap { Pressure, Transport };

    void Reset();
    /** Flip CHECK: when on, a facet whose donor REVERSED relative to D0 (the
        direction pass's pre-solve choice) is refused rather than re-decided —
        UPWIND_FLIPPED is written instead of a donor, and that facet's velocity
        is zeroed, so it drops out of the nodal velocity sum and contributes no
        CFL candidate. This is the legacy behaviour, off by default.

        Note it does NOT enable "flipping" in any positive sense: nothing flips.
        It enables a CHECK that neutralises reversing facets. With it off, the
        reversed donor is simply written and the facet keeps its full velocity. */
    void ZeroVelocityOnFlip( bool on );
    void SetVerbose(bool verbose);
    void Gravity( bool with_gravity );
    /** Enable the velocity/CFL computation AND select the TRANSPORT donor map.

        There are two donor maps, identically shaped, because the pressure and
        transport equations want different things from the donor choice:

          Upwinder          (D0)  written by the DIRECTION pass (after Reset()),
                                  which runs BEFORE the solve. Read by the upwind
                                  PDE operators via UpwindMatrix(). The scheme may
                                  skip that pass on dt-cut retries, in which case
                                  D0 stays fixed for the step and every retry
                                  assembles the SAME pressure matrix.

          UpwinderTransport (D1)  written by the VELOCITY pass (after
                                  WithVelocity()), which runs AFTER
                                  UpdateProjection on EVERY solve. Read by
                                  transport via UpwindMatrices(). NEVER frozen —
                                  it is rebuilt from the velocities each solve
                                  produced, retries included, so the donor map
                                  handed to AdvectMassConserved always matches
                                  the field its flux rates come from.

        The selection is implicit and not settable from outside: Reset() targets
        D0, WithVelocity() targets D1. There is deliberately no way to freeze D1.

        The flip check (ZeroVelocityOnFlip) applies to whichever map is written,
        and always compares against D0 — so "reversed" keeps its meaning: THIS
        solve turned the facet around. */
    void WithVelocity();

    void SetLargestTimeStep( double timestep );
    void Adjust_CFL_Criterion( double scale_factor, bool take_pore_velocity);

    void ActivateAir(Model<dim>& model,
                     ExplicitFiniteVolumeTransportPHX<dim>& fv_air,
                     const std::string& density_air,
                     const std::string& relperm_visc_air,
                     const std::string& saturation_air,
                     const std::string& cfl_air,
                     const std::string& velocity_air,
                     const std::string& pore_velocity_air);

    /// D0 — pressure operators (CVFEM_Upwind_NumIntegral_dNT_op_dN_dV).
    DenseMatrix<DM_MIN>& UpwindMatrix(csmp::Index rho_index, size_t eidx);
    /// D1 — transport (ThreePhaseTransportPHX::DetermineFacetFluxRatesOnce).
    std::vector<DenseMatrix<DM_MIN> >& UpwindMatrices(csmp::Index rho_index);

    double   WorstCfl()        const { return worst_cfl_; }
    size_t   WorstCflElement()  const { return worst_cfl_element_; }
    uint32_t WorstCflPhase()   const { return worst_cfl_phase_; }
    double   WorstCflSat()     const { return worst_cfl_sat_; }
    double   WorstCflVel()     const { return worst_cfl_vel_; }
    double   WorstCflPoreVel() const { return worst_cfl_porevel_; }
    double WorstCflRelperm() const { return worst_cfl_relperm_; }

    void ResetCflDiagnostics()
    {
        worst_cfl_         = largest_time_step;
        worst_cfl_element_ = 0;
        worst_cfl_phase_   = 0;
        worst_cfl_sat_ = worst_cfl_vel_ = worst_cfl_porevel_ = 0.;
        worst_cfl_relperm_ = 0.;
    }

private:
    void DetermineUpwindNodes( Element<dim>& e );

    // ── Phase FV transport objects ────────────────────────────────────────────
    // Phase index convention matches densities[] / relperm_vis[] / saturations[]:
    //   p == 0  →  liquid    (fv_transport_liquid)
    //   p == 1  →  vapor     (fv_transport_vapor)
    //   p == 2  →  air       (fv_transport_air)
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_vapor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_liquid;

    ExplicitFiniteVolumeTransportPHX<dim>* fv_transport_air = nullptr;

    std::vector<csmp::Index> rho_key, relperm_visc_key, S_key;
    std::vector<csmp::Index> uc_key, cfl_key, pore_vel_key, vel_key;
    csmp::Index k_key, phi_key, KgradP_key, tot_pore_vel_key, tot_vel_key, p_key, sh_key;
    csmp::Index thickness_key;

    ScalarVariable                              k, uc_scal, phi;
    ScalarVariable                              thickness;
    std::vector<ScalarVariable>                 sh;
    VectorVariable<dim>                         KgradP, total_pore_velocity, total_velocity;
    std::vector<VectorVariable<dim> >           pore_phase_velocity, phase_velocity;

    std::vector<std::vector<ScalarVariable> >   rho, relperm_visc, S;

    // D0: written by the direction pass, read by the pressure operators.
    std::vector<std::vector<DenseMatrix<DM_MIN> > >  Upwinder;
    // D1: written by the velocity pass, read by transport. Same shape as Upwinder.
    std::vector<std::vector<DenseMatrix<DM_MIN> > >  UpwinderTransport;
    DenseMatrix<DM_MIN> facet_normals, projection, grad, inversed_facet_normals;

    std::vector<ScalarVariable> cfl;

    std::vector<std::vector<double> >  facet_pore_velocity;
    std::vector<double>  facet_normal;

    uint32_t VERTICAL_AXIS;
    double gravity, pot_crit;
    double velocity_inside, velocity_outside;
    double velocity, abs_velocity, norm, mobility, density, sat, relperm_visc_up;
    double normal_component, g;
    double distance, largest_time_step;
    double cfl_scaling;

    uint32_t xyz;
    uint32_t inside_node_, outside_node_;

    bool zero_velocity_on_flip_;   // legacy flip check — see ZeroVelocityOnFlip()
    DonorMap donor_map_target_;    // Pressure -> Upwinder (D0); Transport -> UpwinderTransport (D1)
    bool verbose, grav, with_velocity, facet_cfl;
    bool cfl_with_pore_velocity;

    uint32_t phases, facets;

    //Benoit add:
    double   worst_cfl_;          // smallest cfl_dt seen this visit
    size_t   worst_cfl_element_;
    uint32_t worst_cfl_phase_;
    double   worst_cfl_sat_;      // saturation at the binding facet
    double   worst_cfl_vel_;      // raw (Darcy) velocity
    double   worst_cfl_porevel_;  // velocity/(sat*phi) actually used
    double worst_cfl_relperm_;
};

} // csmp

#endif