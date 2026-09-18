// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "UpwindControlVisitor.h"
#include "Region.h"
#include <iostream>
#include <chrono>

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim>
UpwindControlVisitor<dim>::UpwindControlVisitor( Model<dim>& model,
                                                ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                                                ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                                                ExplicitFiniteVolumeTransportPHX<dim>* fv_air,      // ← pointer, optional
                                                const char* permeability,
                                                const char* porosity,
                                                std::vector<std::string>& densities,
                                                std::vector<std::string>& relperm_vis,
                                                std::vector<std::string>& saturations,
                                                std::vector<std::string>& cfl_variables,
                                                std::vector<std::string>& velocities,
                                                std::vector<std::string>& pore_velocities)
    : fv_transport_vapor( fv_vapor ),
    fv_transport_liquid( fv_liquid ),
    fv_transport_air( fv_air ),
    phases( static_cast<uint32_t>(densities.size())),   // 2 (no air) or 3 (with air)
    gravity(-9.80665),
    pot_crit(1.e-20),
    xyz(dim),
    zero_velocity_on_flip_(false),
    donor_map_target_(Pressure),
    verbose(false),
    grav(true),
    with_velocity(false),                   // off until WithVelocity() is called
    largest_time_step(60.*60.*24.*365.),    // CFL ceiling: 1 year
    VERTICAL_AXIS( (dim==1u) ? 0u : 1u ),   // gravity axis: x in 1D, y in 2D/3D
    cfl_scaling(0.1),                       // default safety fraction; overridden via Adjust_CFL_Criterion
    cfl_with_pore_velocity(true)           // default = true; see DetermineUpwindNodes
{
    if (dim==3U) xyz = 2;
    // ── Vertical-axis index for the gravity term ─────────────────────────────
    // Benoit: xyz-1 and VERTICAL_AXIS agree in every
    // dimension — they both give 0,1,1 (x,y,y) for dim 1,2,3:
    //   dim=1: xyz=1 -> xyz-1 = 0 ; VERTICAL_AXIS = 0
    //   dim=2: xyz=2 -> xyz-1 = 1 ; VERTICAL_AXIS = 1
    //   dim=3: xyz=2 -> xyz-1 = 1 ; VERTICAL_AXIS = 1
    // xyz-1 indexes the facet-normal component (GetFacetNormalComponent);
    // VERTICAL_AXIS indexes the velocity-vector component the gravity term adds to.
    // Both land on the SAME axis, so there is no inconsistency — just two members
    // feeding two different APIs. Convention worth knowing: the vertical is y(=1),
    // NOT z(=2), even in 3D, so a 3D mesh must align gravity with the y axis. The
    // pair is redundant (could be one member) but correct; left as-is to avoid
    // touching the two call sites.

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);

    k_key            = model.Database().StorageKey(permeability);
    phi_key          = model.Database().StorageKey(porosity);
    KgradP_key       = model.Database().StorageKey("KgradP");
    tot_pore_vel_key = model.Database().StorageKey("pore velocity");
    tot_vel_key      = model.Database().StorageKey("velocity");
    p_key            = model.Database().StorageKey("fluid pressure");
    sh_key           = model.Database().StorageKey("saturation halite");

    facet_pore_velocity.resize(phases);
    pore_vel_key.resize(phases);
    pore_phase_velocity.resize(phases);
    vel_key.resize(phases);
    phase_velocity.resize(phases);
    facet_normal.resize(dim);

    if ( k_key.type != SCALAR )
        throw csmp::Exception( ERROR, "UpwindControlVisitor::(constructor)",
                              permeability, " must be a scalar property." );

    if ( phi_key.type != SCALAR )
        throw csmp::Exception( ERROR, "UpwindControlVisitor::(constructor)",
                              porosity, " must be a scalar property." );

    rho_key.resize(phases);
    rho.resize(phases);
    cfl_key.resize(phases);
    cfl.resize(phases);
    relperm_visc_key.resize(phases);
    relperm_visc.resize(phases);
    S_key.resize(phases);
    S.resize(phases);

    for ( uint32_t i = 0; i < phases; i++)
    {
        rho_key[i]         = model.Database().StorageKey(densities[i].c_str());
        relperm_visc_key[i] = model.Database().StorageKey(relperm_vis[i].c_str());
        S_key[i]           = model.Database().StorageKey(saturations[i].c_str());
        cfl_key[i]         = model.Database().StorageKey(cfl_variables[i].c_str());
        vel_key[i]         = model.Database().StorageKey(velocities[i].c_str());
        pore_vel_key[i]    = model.Database().StorageKey(pore_velocities[i].c_str());

        if ( rho_key[i].place != NODE || rho_key[i].type != SCALAR )
            throw csmp::Exception( ERROR, "UpwindControlVisitor::(constructor)",
                                  densities[i].c_str(), " must be a scalar nodal property." );

        if ( relperm_visc_key[i].place != NODE || relperm_visc_key[i].type != SCALAR )
            throw csmp::Exception( ERROR, "UpwindControlVisitor::(constructor)",
                                  relperm_vis[i].c_str(), " must be a scalar nodal property." );

        if ( S_key[i].place != NODE || S_key[i].type != SCALAR )
            throw csmp::Exception( ERROR, "UpwindControlVisitor::(constructor)",
                                  saturations[i].c_str(), " must be a scalar nodal property." );
    }

    typename std::vector<Element<dim>* >::const_iterator it;

    // Two donor maps, identically shaped (see the header for who reads which):
    //   Upwinder          — written by the DIRECTION pass  -> D0 -> pressure operators
    //   UpwinderTransport — written by the VELOCITY  pass  -> D1 -> transport
    Upwinder.resize(phases);
    UpwinderTransport.resize(phases);
    for ( uint32_t i = 0; i < phases; i++)
    {
        Upwinder[i].resize( model.Region("Model").Cells() );
        UpwinderTransport[i].resize( model.Region("Model").Cells() );
        for ( it = model.Region("Model").CellsBegin(); it < model.Region("Model").CellsEnd(); it++)
        {
            // Benoit 10/09/2026: DenseMatrix has no operator=(double); use Zero()
            // AFTER Resize() so the zeroed extent matches the active rows/cols.
            Upwinder[i][(*(*it)).Idx()].Resize((*(*it)).Nodes(), (*(*it)).Nodes());
            Upwinder[i][(*(*it)).Idx()].Zero();
            UpwinderTransport[i][(*(*it)).Idx()].Resize((*(*it)).Nodes(), (*(*it)).Nodes());
            UpwinderTransport[i][(*(*it)).Idx()].Zero();
        }
    }

    // ── Upwind control keys — one per phase, in phase-index order ─────────────
    // Phase index convention: 0 = liquid, 1 = vapor, 2 = air.
    // This MUST stay in sync with the densities[] order passed in.
    //
    // In production the no-air delegating constructor is used (the scheme
    // builds liquid+vapor only — Coupled_..._Scheme constructs upwind_control with
    // no fv_air), so here phases==2 and only [0],[1] are set. The air slot [2] is
    // created later by ActivateAir() (the live path when air is enabled). So the
    // air key is genuinely registered in ONE place per run: here if a caller ever
    // constructs 3-phase directly, or in ActivateAir() for the no-air-then-activate
    // production path. The two must register the SAME air key name ("upwind control
    // air") — they currently do; keep them in sync if that name changes. This is
    // documented rather than refactored to avoid touching the working build path.

    uc_key.resize(phases);
    uc_key[0] = model.Database().StorageKey("upwind control liquid");
    uc_key[1] = model.Database().StorageKey("upwind control vapor");

    if (phases == 3)
        uc_key[2] = model.Database().StorageKey("upwind control air");
}

// No-air delegating constructor
template<uint32_t dim>
UpwindControlVisitor<dim>::UpwindControlVisitor(
    Model<dim>& model,
    ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
    ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
    const char* permeability,
    const char* porosity,
    std::vector<std::string>& densities,
    std::vector<std::string>& relperm_vis,
    std::vector<std::string>& saturations,
    std::vector<std::string>& cfl_variables,
    std::vector<std::string>& velocities,
    std::vector<std::string>& pore_velocities)
    : UpwindControlVisitor(model, fv_liquid, fv_vapor,
                           nullptr,                     // ← fv_air = nullptr
                           permeability, porosity,
                           densities, relperm_vis, saturations,
                           cfl_variables, velocities, pore_velocities)
{}

/** default destructor */
template<uint32_t dim>
UpwindControlVisitor<dim>::~UpwindControlVisitor()
{}

template<uint32_t dim>
void UpwindControlVisitor<dim>::Visit(Model<dim>* m)
{}

/** visit function for region — unused hook (per-element work is in Visit(Element)). */
template<uint32_t dim>
void UpwindControlVisitor<dim>::Visit(Region<dim>* r)
{}

/** visit function for elements.
    Per element: load nodal properties, (re)initialise per-element CFL, compute the
    upwind directions + per-phase velocities/CFL in DetermineUpwindNodes, then store
    the per-node upwind-control scalar, velocities and CFL back to the mesh. */
template<uint32_t dim>
void UpwindControlVisitor<dim>::Visit(Element<dim>* e)
{
    facets = e->FV()->Facets();

    e->Read( k_key, k );                              // element permeability (scalar)
    e->NodePropertyVector( sh_key, sh );              // nodal halite saturation (for plug check)

    // Nodal phase properties needed for upwinding / velocity (all phases).
    for ( uint32_t i = 0; i < phases; i++ )
    {
        e->NodePropertyVector( rho_key[i],         rho[i]         );
        e->NodePropertyVector( relperm_visc_key[i], relperm_visc[i] );
        e->NodePropertyVector( S_key[i],           S[i]           );
    }

    if ( with_velocity )
    {
        e->Read( phi_key,   phi   );
        e->Read( KgradP_key, KgradP );

        // Per-element CFL reset to the ceiling
        for ( uint32_t i = 0; i < phases; i++ )
        {
            cfl[i] = largest_time_step;
            facet_pore_velocity[i].resize(facets);
        }
    }

    DetermineUpwindNodes( *e );

    for ( uint32_t p = 0; p < phases; ++p )
    {
        uc_scal() = 0;
        for ( uint32_t i{0}; i < e->FV()->Facets(); i++ )
        {
            e->FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
            uc_scal() += Upwinder[p][e->Idx()](inside_node_, outside_node_);
        }
        e->Store( uc_key[p], uc_scal );
    }

    if ( with_velocity )
    {
        // Write per-phase velocity, pore velocity and CFL back to the mesh, plus totals.
        for ( uint32_t p = 0; p < phases; ++p )
        {
            e->Store( vel_key[p],      phase_velocity[p]      );
            e->Store( pore_vel_key[p], pore_phase_velocity[p] );
            e->Store( cfl_key[p],      cfl[p]                 );
        }

        e->Store( tot_pore_vel_key, total_pore_velocity );
        e->Store( tot_vel_key,      total_velocity      );
    }
}


/** define upwind nodes for CVFEM scheme — and calculate velocities if requested.
    For each phase and each facet of the element:
      1. build the phase velocity across the facet (K·gradP + gravity, then ×kr/μ),
      2. decide the upwind direction from the signs of the two one-sided velocities,
      3. record it in the per-phase Upwinder matrix,
      4. (if with_velocity) derive the facet's CFL contribution.
    Then per phase, assemble the nodal phase velocity and pore velocity.

    Upwinder matrix encoding (per element, per phase, [inside][outside]):
      0  = UPWIND_NONE     no net flow across this facet
      1  = UPWIND_INSIDE   upwind node is the inside node
      2  = UPWIND_OUTSIDE  upwind node is the outside node
      10 = UPWIND_FLIPPED  donor reversed vs D0 while the flip check was on
                           (legacy path only; never written while frozen) */
template<uint32_t dim>
void UpwindControlVisitor<dim>::DetermineUpwindNodes( Element<dim>& e )
{
    constexpr uint32_t UPWIND_NONE    = 0;
    constexpr uint32_t UPWIND_INSIDE  = 1;
    constexpr uint32_t UPWIND_OUTSIDE = 2;
    constexpr uint32_t UPWIND_FLIPPED = 10;

    // gravity component used to visualise velocity vectors
    std::vector<double> gravity_component;
    if (dim == 1U)
        gravity_component = {1.0};
    else if (dim == 2U)
        gravity_component = {0.0, 1.0};
    else if (dim == 3U)
        gravity_component = {0.0, 1.0, 0.0};

    // initialise total velocity to zero
    total_velocity = total_pore_velocity = 0.;

    for ( uint32_t p = 0; p < phases; ++p )
    {
        // Flux-weighted accumulators for this phase's nodal velocity (reset per phase).
        norm = mobility = density = sat = relperm_visc_up = 0.;

        for ( uint32_t i = 0; i < facets; ++i )
        {
            if (verbose && p==0) cerr << endl << "Facet: " << i;

            facet_cfl = true;                                  // cleared if facet is a DIRICH boundary or has no net flow
            e.FV()->FacetEdgeNodes(i, inside_node_, outside_node_);

            double rho_in  = rho[p][inside_node_]();
            double rho_out = rho[p][outside_node_]();

            double relperm_visc_in  = relperm_visc[p][inside_node_]();
            double relperm_visc_out = relperm_visc[p][outside_node_]();

            // Which donor map does THIS pass write?
            //   direction pass (after Reset())        -> Upwinder          (D0)
            //   velocity  pass (after WithVelocity()) -> UpwinderTransport (D1)
            // Everything below is unchanged; only the target array differs.
            // WRITE target: the direction pass writes D0, the velocity pass D1.
            // Set implicitly — Reset() selects Pressure, WithVelocity() selects
            // Transport — so the transport map can never be frozen by mistake.
            DenseMatrix<DM_MIN>& write_map =
                ( donor_map_target_ == Transport ) ? UpwinderTransport[p][e.Idx()]
                                                 : Upwinder[p][e.Idx()];
            auto& upwinder_inside_outside = write_map(inside_node_, outside_node_);
            auto& upwinder_outside_inside = write_map(outside_node_, inside_node_);

            // READ reference for the flip CHECK: ALWAYS D0, i.e. the
            // donor this step's DIRECTION pass chose from the pre-solve velocities.
            // So "flipped" keeps its original meaning — this solve reversed the
            // facet — regardless of which map is being written.
            const double prior_donor = Upwinder[p][e.Idx()](inside_node_, outside_node_);

            // ── Phase velocity dispatch ───────────────────────────────────────
            // Each phase has its own ExplicitFiniteVolumeTransportPHX whose base
            // velocity field is KgradP; gravity and relperm_visc are applied here.
            // Phase index convention: 0 = liquid, 1 = vapor, 2 = air.
            if      ( p == 0 ) velocity = fv_transport_liquid.GetFacetNormalVelocity(e.Idx(), i);
            else if ( p == 1 ) velocity = fv_transport_vapor .GetFacetNormalVelocity(e.Idx(), i);
            else               velocity = fv_transport_air   ->GetFacetNormalVelocity(e.Idx(), i);

            velocity_inside = velocity_outside = velocity;

            // ── Gravity / buoyancy term, added per side with that side's density ──
            // (each side uses its own nodal density, so the two directions can differ)
            if ( grav )
            {
                normal_component  = fv_transport_vapor.GetFacetNormalComponent( e.Idx(), i, xyz-1 );
                g                 = gravity * normal_component;
                velocity_inside  += k() * rho_in  * g;
                velocity_outside += k() * rho_out * g;
            }

            if (verbose && p==0) cerr << " {Init} v_in: " << velocity_inside << ", v_out: " << velocity_outside;
            if (verbose && p==0)
            {
                if(upwinder_inside_outside == UPWIND_NONE)    cerr << endl << "upwinder_inside_outside: UPWIND_NONE";
                if(upwinder_inside_outside == UPWIND_INSIDE)  cerr << endl << "\033[92m upwinder_inside_outside: UPWIND_INSIDE \033[0m";
                if(upwinder_inside_outside == UPWIND_OUTSIDE) cerr << endl << "\033[93m upwinder_inside_outside: UPWIND_OUTSIDE \033[0m";
                if(upwinder_inside_outside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_inside_outside: UPWIND_FLIPPED \033[0m";
                if(upwinder_outside_inside == UPWIND_NONE)    cerr << endl << "upwinder_outside_inside: UPWIND_NONE";
                if(upwinder_outside_inside == UPWIND_INSIDE)  cerr << endl << "\033[93m upwinder_outside_inside: UPWIND_INSIDE \033[0m";
                if(upwinder_outside_inside == UPWIND_OUTSIDE) cerr << endl << "\033[92m upwinder_outside_inside: UPWIND_OUTSIDE \033[0m";
                if(upwinder_outside_inside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_outside_inside: UPWIND_FLIPPED \033[0m";
            }

            // ── Directional gating ───────────────────────────────────────────
            // Keep only the physically admissible one-sided component:
            //   inside→outside flow must be >0, outside→inside flow must be <0.
            if ( velocity_inside  < 0 ) velocity_inside  = 0.0;
            if ( velocity_outside > 0 ) velocity_outside = 0.0;

            // No fluid on a side (non-positive density) ⇒ that side cannot supply flow.
            if ( rho_in  <= 0 ) velocity_inside  = 0.0;
            if ( rho_out <= 0 ) velocity_outside = 0.0;

            // ── Mobility weighting ───────────────────────────────────────────
            // Multiply by kr/μ so the velocities become PHASE fluxes, not raw Darcy.
            // From here `velocity` (once selected below) carries kr/μ. A low-relperm
            // phase therefore moves slowly even under a large pressure gradient.
            velocity_inside  *= relperm_visc_in;
            velocity_outside *= relperm_visc_out;

            if (verbose && p==0) cerr << endl << " {Step 1} v_in: " << velocity_inside << ", v_out: " << velocity_outside;

            if ( velocity_inside > 0 && velocity_outside < 0 )
            {
                if ( abs(velocity_inside) <= abs(velocity_outside) ) velocity_inside  = 0.0;
                if ( abs(velocity_inside) >= abs(velocity_outside) ) velocity_outside = 0.0;
            }

            if ( velocity_inside > 0 && velocity_outside >= 0 )
            {
                velocity = velocity_inside;

                if ( !zero_velocity_on_flip_ )
                {
                    upwinder_inside_outside = UPWIND_INSIDE;
                    upwinder_outside_inside = UPWIND_OUTSIDE;
                }
                else
                {
                    if ( prior_donor != UPWIND_FLIPPED )
                    {
                        if ( prior_donor != UPWIND_INSIDE )
                        {
                            upwinder_inside_outside = UPWIND_FLIPPED;
                            upwinder_outside_inside = UPWIND_FLIPPED;
                            velocity = 0.;
                        }
                        else
                        {
                            upwinder_inside_outside = UPWIND_INSIDE;
                            upwinder_outside_inside = UPWIND_OUTSIDE;
                        }
                    }
                }

                if (verbose && p==0) cerr << endl << " {A} velocity_inside > 0 && velocity_outside >= 0, velocity: " << velocity;
                if (verbose && p==0)
                {
                    if(upwinder_inside_outside == UPWIND_NONE)    cerr << endl << "upwinder_inside_outside: UPWIND_NONE";
                    if(upwinder_inside_outside == UPWIND_INSIDE)  cerr << endl << "\033[92m upwinder_inside_outside: UPWIND_INSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_OUTSIDE) cerr << endl << "\033[93m upwinder_inside_outside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_inside_outside: UPWIND_FLIPPED \033[0m";
                    if(upwinder_outside_inside == UPWIND_NONE)    cerr << endl << "upwinder_outside_inside: UPWIND_NONE";
                    if(upwinder_outside_inside == UPWIND_INSIDE)  cerr << endl << "\033[93m upwinder_outside_inside: UPWIND_INSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_OUTSIDE) cerr << endl << "\033[92m upwinder_outside_inside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_outside_inside: UPWIND_FLIPPED \033[0m";
                }

                abs_velocity = abs(velocity);
                density  += rho_in * abs_velocity;
                mobility += relperm_visc_in * abs_velocity;
                norm     += abs_velocity;
                sat       = S[p][inside_node_]();
                relperm_visc_up = relperm_visc_in;

                if ( e.N(outside_node_)->Status(p_key) == DIRICH ) facet_cfl = false;
            }

            else if ( velocity_inside <= 0 && velocity_outside < 0 )
            {
                velocity = velocity_outside;

                if ( !zero_velocity_on_flip_ )
                {
                    upwinder_inside_outside = UPWIND_OUTSIDE;
                    upwinder_outside_inside = UPWIND_INSIDE;
                }
                else
                {
                    if ( prior_donor != UPWIND_FLIPPED )
                    {
                        if ( prior_donor != UPWIND_OUTSIDE )
                        {
                            upwinder_inside_outside = UPWIND_FLIPPED;
                            upwinder_outside_inside = UPWIND_FLIPPED;
                            velocity = 0.;
                        }
                        else
                        {
                            upwinder_inside_outside = UPWIND_OUTSIDE;
                            upwinder_outside_inside = UPWIND_INSIDE;
                        }
                    }
                }

                if (verbose && p==0) cerr << endl << " {B} velocity_inside <= 0 && velocity_outside < 0, velocity: " << velocity;
                if (verbose && p==0)
                {
                    if(upwinder_inside_outside == UPWIND_NONE)    cerr << endl << "upwinder_inside_outside: UPWIND_NONE";
                    if(upwinder_inside_outside == UPWIND_INSIDE)  cerr << endl << "\033[92m upwinder_inside_outside: UPWIND_INSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_OUTSIDE) cerr << endl << "\033[93m upwinder_inside_outside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_inside_outside: UPWIND_FLIPPED \033[0m";
                    if(upwinder_outside_inside == UPWIND_NONE)    cerr << endl << "upwinder_outside_inside: UPWIND_NONE";
                    if(upwinder_outside_inside == UPWIND_INSIDE)  cerr << endl << "\033[93m upwinder_outside_inside: UPWIND_INSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_OUTSIDE) cerr << endl << "\033[92m upwinder_outside_inside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_outside_inside: UPWIND_FLIPPED \033[0m";
                }

                abs_velocity = abs(velocity);
                density  += rho_out * abs_velocity;
                mobility += relperm_visc_out * abs_velocity;
                norm     += abs_velocity;
                sat       = S[p][outside_node_]();
                relperm_visc_up = relperm_visc_out;

                if ( e.N(inside_node_)->Status(p_key) == DIRICH ) facet_cfl = false;
            }

            else
            {
                velocity = 0.;

                if ( !zero_velocity_on_flip_ )
                {
                    upwinder_inside_outside = UPWIND_NONE;
                    upwinder_outside_inside = UPWIND_NONE;
                }
                else
                {
                    if ( prior_donor != UPWIND_FLIPPED )
                    {
                        if ( prior_donor != UPWIND_NONE )
                        {
                            upwinder_inside_outside = UPWIND_FLIPPED;
                            upwinder_outside_inside = UPWIND_FLIPPED;
                        }
                        else
                        {
                            upwinder_inside_outside = UPWIND_NONE;
                            upwinder_outside_inside = UPWIND_NONE;
                        }
                    }
                }

                if (verbose && p==0) cerr << endl << " {C} v_in and v_out have different signs or are null, velocity: " << velocity;
                if (verbose && p==0)
                {
                    if(upwinder_inside_outside == UPWIND_NONE)    cerr << endl << "upwinder_inside_outside: UPWIND_NONE";
                    if(upwinder_inside_outside == UPWIND_INSIDE)  cerr << endl << "\033[92m upwinder_inside_outside: UPWIND_INSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_OUTSIDE) cerr << endl << "\033[93m upwinder_inside_outside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_inside_outside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_inside_outside: UPWIND_FLIPPED \033[0m";
                    if(upwinder_outside_inside == UPWIND_NONE)    cerr << endl << "upwinder_outside_inside: UPWIND_NONE";
                    if(upwinder_outside_inside == UPWIND_INSIDE)  cerr << endl << "\033[93m upwinder_outside_inside: UPWIND_INSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_OUTSIDE) cerr << endl << "\033[92m upwinder_outside_inside: UPWIND_OUTSIDE \033[0m";
                    if(upwinder_outside_inside == UPWIND_FLIPPED) cerr << endl << "\033[94m upwinder_outside_inside: UPWIND_FLIPPED \033[0m";
                }

                facet_cfl = false;
            }

            // Halite blocks all phases. Applies to whichever map this pass writes.
            if ( sh[inside_node_]()  >= 1.0 ) { upwinder_inside_outside = UPWIND_OUTSIDE; upwinder_outside_inside = UPWIND_INSIDE;  }
            if ( sh[outside_node_]() >= 1.0 ) { upwinder_outside_inside = UPWIND_OUTSIDE; upwinder_inside_outside = UPWIND_INSIDE;  }

            if ( with_velocity )
            {
                // ════════════════════════════════════════════════════════════
                //  CFL CRITERION: why the Darcy (bulk-flux) form is correct here,
                //  not the pore-velocity (saturation-front) form.
                // ════════════════════════════════════════════════════════════
                //  `velocity` here is the MOBILITY-WEIGHTED phase flux (kr/mu and
                //  gravity already applied in the upwind branches above) — NOT raw
                //  Darcy. Two candidate criteria:
                //
                //   (A) pore-velocity / saturation-front:  v / (sat·phi)
                //       The textbook IMPES criterion — correct WHEN saturation is
                //       the advected variable moving along a fractional-flow curve.
                //       THIS SCHEME DOES NOT DO THAT: it advects MASS, and
                //       saturation is a flash OUTPUT, not a marched variable. There
                //       is no saturation wave for v/(sat·phi) to bound, so it limits
                //       a wave the scheme never evolves. Worse, with (near-)linear
                //       vapor relperm kr_v = S_v the saturation cancels
                //       (kr_v/(sat·phi) = 1/(mu_v·phi)), so the "pore velocity"
                //       collapses to the intrinsic vapor Darcy speed INDEPENDENT of
                //       how much vapor is present: a 0.2%-vapor cell then pins dt at
                //       the same speed as a 27%-vapor cell. Verified directly —
                //       relperm_visc/sat is ~constant across that whole saturation
                //       range. That is a division artifact, not a physical limit.
                //
                //   (B) Darcy / bulk-flux:  v
                //       Bounds how far MASS moves per step — the quantity this
                //       scheme actually advects — so it is the correct stability
                //       limit. The other thing a phase CFL might guard, per-phase
                //       over-draining, is already handled exactly and separately by
                //       CheckDry / transport sub-cycling, so dropping the sat
                //       division loses no protection.
                //
                //  => Use (B). cfl_with_pore_velocity should be FALSE for this
                //     scheme. (A) is retained behind the flag for reference / for a
                //     future scheme that genuinely advects saturation.
                //
                //  CAVEAT: this rests on kr being ~linear in saturation (true for
                //  the near-critical steam-water front here). For Corey-type
                //  kr ∝ S^n (n>1) — e.g. a cold air-water front — kr/sat → 0 as
                //  S → 0 and the trade-off changes; revisit if such a relperm is
                //  introduced. (facet_pore_velocity is still computed below because
                //  it is also stored as an output field, independent of the CFL.)

                if ( sat * phi() != 0. )
                    facet_pore_velocity[p][i] = velocity / (sat * phi());
                else
                    facet_pore_velocity[p][i] = 0.;

                distance = e.N(inside_node_)->Coordinate().DistanceTo(e.N(outside_node_)->Coordinate());

                if ( facet_pore_velocity[p][i] != 0. && facet_cfl )
                {
                    double cfl_candidate;
                    if ( cfl_with_pore_velocity )
                    {
                        // (A) pore-velocity / saturation-front criterion (tighter).
                        // Carries the 1/sat artifact above — NOT recommended here.
                        cfl_candidate = distance * cfl_scaling / fabs(facet_pore_velocity[p][i]);
                    }
                    else
                    {
                        // (B) Darcy / bulk-flux criterion — the correct one for this
                        // mass-advecting scheme. No sat/phi division.
                        // FIX (history): the previous form was
                        //   distance*cfl_scaling / (facet_pore_velocity * phi)
                        // which only cancelled one phi and left an orphan
                        // velocity/sat — neither Darcy nor pore velocity (an
                        // incomplete edit from when the sat factor was introduced).
                        cfl_candidate = distance * cfl_scaling / fabs(velocity);
                    }

                    cfl[p]() = std::min( cfl[p](), cfl_candidate );

                    // Track the single facet that sets the global CFL minimum, for
                    // diagnostics (reset once per pass). worst_cfl_vel_ is the
                    // MOBILITY-WEIGHTED phase flux (kr/mu applied), not raw Darcy;
                    // recover raw Darcy as worst_cfl_vel_ / worst_cfl_relperm_.
                    if ( cfl_candidate < worst_cfl_ )
                    {
                        worst_cfl_         = cfl_candidate;
                        worst_cfl_element_ = e.Idx();
                        worst_cfl_phase_   = p;
                        worst_cfl_sat_     = sat;
                        worst_cfl_vel_     = velocity;                    // mobility-weighted phase flux
                        worst_cfl_porevel_ = facet_pore_velocity[p][i];  // sat-divided front speed
                        worst_cfl_relperm_ = relperm_visc_up;            // kr/mu at the upwind node
                    }
                }
            }



        } // end facet loop

        if ( with_velocity )
        {
            // ── Nodal phase velocity assembly (flux-weighted over this element's facets) ──
            // phase_velocity = mobility-weighted KgradP (+gravity), averaged by the
            // total flux magnitude `norm`; pore velocity divides by porosity.
            if ( norm != 0. )
            {
                phase_velocity[p] = KgradP;


                // JK: updated velocity vectors: needed to visualise in-plane flow of LDE regions
                // NOTE: this is used for visualisation only
                uint32_t n_nodes = e.Nodes();
                bool is_lde = (n_nodes == dim);

                if ( is_lde )
                {
                    // Project model gravity onto LDE : g_proj = g - (g·n)n
                    auto unrml = e.UnitNormal();

                    if( dim == 2U )
                    {
                        double gdotn = gravity_component[0]*unrml[0] + gravity_component[1]*unrml[1];
                        gravity_component[0] -= gdotn * unrml[0];
                        gravity_component[1] -= gdotn * unrml[1];
                    }
                    else if( dim == 3U )
                    {
                        double gdotn = gravity_component[0]*unrml[0] + gravity_component[1]*unrml[1] + gravity_component[2]*unrml[2];
                        gravity_component[0] -= gdotn * unrml[0];
                        gravity_component[1] -= gdotn * unrml[1];
                        gravity_component[2] -= gdotn * unrml[2];
                    }
                }

                // gravity contribution uses flux-weighted mean density (density/norm);
                // gravity_component already reduces to the VERTICAL_AXIS-only case for
                // non-LDE elements; -> this is the only gravity addition needed.
                if ( grav && norm != 0. ) {
                    for (size_t i = 0; i < gravity_component.size(); i++)
                        phase_velocity[p](i) += gravity_component[i] * gravity * density * k() / norm;
                }
                phase_velocity[p] *= mobility / norm;          // flux-weighted mobility average

                if ( phi() != 0. ) pore_phase_velocity[p] = phase_velocity[p] / phi();
                else               pore_phase_velocity[p] = 0.;
            }
            else
            {
                // no flux contributions on any facet ⇒ zero velocity for this phase
                phase_velocity[p]      = 0.;
                pore_phase_velocity[p] = 0.;
            }

            total_velocity      += phase_velocity[p];
            total_pore_velocity += pore_phase_velocity[p];
        }

    } // end phase loop

} // end DetermineUpwindNodes


/** reset booleans (call before reconfiguring; does NOT clear CFL diagnostics) */
template<uint32_t dim>
void UpwindControlVisitor<dim>::Reset()
{
    zero_velocity_on_flip_ = false;
    with_velocity      = false;
    donor_map_target_ = Pressure;   // the direction pass always writes D0
}

/** set boolean flip — enables the UPWIND_FLIPPED hysteresis (set by WithVelocity) */
template<uint32_t dim>
void UpwindControlVisitor<dim>::ZeroVelocityOnFlip(bool on)
{
    zero_velocity_on_flip_ = on;
}

/** set boolean verbose (phase-0 diagnostic prints in DetermineUpwindNodes) */
template<uint32_t dim>
void UpwindControlVisitor<dim>::SetVerbose(bool verb)
{
    verbose = verb;
}

/** access to upwind matrix of specified phase (looked up by density key) and element */
template<uint32_t dim>
DenseMatrix<DM_MIN>& UpwindControlVisitor<dim>::UpwindMatrix(csmp::Index rho_index, size_t eidx)
{
    // D0 — the pressure operators' map, written only by the direction pass.
    for (uint32_t i = 0; i < rho_key.size(); i++)
        if (rho_index == rho_key[i])
            return Upwinder[i][eidx];

    throw csmp::Exception( ERROR, "UpwindControlVisitor<dim>::UpwindMatrix:", "density index never identified." );
}

/** access to vector of upwind matrices for a phase (looked up by density key).
    Used by ThreePhaseTransportPHX::DetermineFacetFluxes to get the per-phase
    upwind operator. */
template<uint32_t dim>
std::vector<DenseMatrix<DM_MIN> >& UpwindControlVisitor<dim>::UpwindMatrices(csmp::Index rho_index)
{
    // D1 — transport's map, refreshed by the velocity pass after every solve.
    for (uint32_t i = 0; i < rho_key.size(); i++)
        if (rho_index == rho_key[i])
            return UpwinderTransport[i];

    throw csmp::Exception( ERROR, "UpwindControlVisitor<dim>::UpwindMatrices:", "density index never identified." );
}

/** activate or deactivate gravity component */
template<uint32_t dim>
void UpwindControlVisitor<dim>::Gravity(bool with_gravity)
{
    grav = with_gravity;
}

/** activate velocity + CFL calculations, and target the transport map (D1).
    The flip check is NOT enabled here — call ZeroVelocityOnFlip() for that. */
template<uint32_t dim>
void UpwindControlVisitor<dim>::WithVelocity()
{
    with_velocity     = true;
    donor_map_target_ = Transport;   // the velocity pass ALWAYS writes D1
}

/** set the CFL ceiling (the "no constraint" dt returned when nothing binds) */
template<uint32_t dim>
void UpwindControlVisitor<dim>::SetLargestTimeStep(double timestep)
{
    largest_time_step = timestep;
}

/** modify CFL criterion.
    scale_factor      → cfl_scaling (safety fraction of the transit time)
    take_pore_velocity→ true: pore-velocity (saturation-front) CFL, the 1/(sat·phi)
                        form that is tighter and blows up at trace saturation;
                        false: Darcy (bulk-flux) CFL, looser. See DetermineUpwindNodes
                        for the full caveat on which is appropriate for this scheme. */
template<uint32_t dim>
void UpwindControlVisitor<dim>::Adjust_CFL_Criterion(double scale_factor, bool take_pore_velocity)
{
    cfl_scaling            = scale_factor;
    cfl_with_pore_velocity = take_pore_velocity;
}

/** Add the air phase after construction (no-air → with-air upgrade).
    Performs the first-time air allocation for the production path: the scheme
    builds this visitor 2-phase (no-air constructor, phases==2) and calls
    ActivateAir() when air is enabled (Coupled_..._Scheme::ActivateAir). So this is
    LIVE in any air run, not dead code. It extends every per-phase vector to 3,
    registers the air keys, allocates the air Upwinder matrices, and adds the air
    upwind-control key — i.e. it reproduces, for the new phase-2 slot, the per-phase
    setup the constructor does for phases 0 and 1.
    MAINTENANCE: this and the constructor's per-phase setup are NOT factored into a
    shared helper, so any change to how a phase is registered (a new key, a resize,
    Upwinder init) must be made in BOTH places or an air run will be under-set while
    a no-air run stays fine (both compile). Caller contract: the air key-name
    strings must already exist — the scheme calls names.ActivateAir() first. */
template<uint32_t dim>
void UpwindControlVisitor<dim>::ActivateAir(
    Model<dim>& model,
    ExplicitFiniteVolumeTransportPHX<dim>& fv_air,
    const std::string& density_air,
    const std::string& relperm_visc_air,
    const std::string& saturation_air,
    const std::string& cfl_air,
    const std::string& velocity_air,
    const std::string& pore_velocity_air)
{
    // Guard
    if (phases == 3) return;  // already active

    fv_transport_air = &fv_air;

    // 1. Extend all per-phase key vectors by one slot
    phases = 3;
    rho_key        .resize(3);    rho        .resize(3);
    relperm_visc_key.resize(3);   relperm_visc.resize(3);
    S_key          .resize(3);    S          .resize(3);
    cfl_key        .resize(3);    cfl        .resize(3);
    vel_key        .resize(3);
    pore_vel_key   .resize(3);
    facet_pore_velocity.resize(3);
    pore_phase_velocity.resize(3);
    phase_velocity .resize(3);

    // 2. Register the air storage keys from the model database
    rho_key[2]          = model.Database().StorageKey(density_air.c_str());
    relperm_visc_key[2] = model.Database().StorageKey(relperm_visc_air.c_str());
    S_key[2]            = model.Database().StorageKey(saturation_air.c_str());
    cfl_key[2]          = model.Database().StorageKey(cfl_air.c_str());
    vel_key[2]          = model.Database().StorageKey(velocity_air.c_str());
    pore_vel_key[2]     = model.Database().StorageKey(pore_velocity_air.c_str());

    // 3. Allocate the Upwinder matrices for air — same loop as constructor
    // Benoit 10/09/2026: UpwinderTransport (D1) was NOT extended here while
    // Upwinder (D0) was, so the velocity pass wrote UpwinderTransport[2] out of
    // range as soon as air was active. Both maps are extended together now —
    // they must stay identically shaped (see header).
    Upwinder.resize(3);
    UpwinderTransport.resize(3);
    Upwinder[2].resize(model.Region("Model").Cells());
    UpwinderTransport[2].resize(model.Region("Model").Cells());
    for (auto it = model.Region("Model").CellsBegin();
         it < model.Region("Model").CellsEnd(); ++it)
    {
        Upwinder[2][(*(*it)).Idx()].Resize((*(*it)).Nodes(), (*(*it)).Nodes());
        Upwinder[2][(*(*it)).Idx()].Zero();
        UpwinderTransport[2][(*(*it)).Idx()].Resize((*(*it)).Nodes(), (*(*it)).Nodes());
        UpwinderTransport[2][(*(*it)).Idx()].Zero();
    }

    // 4. Register the upwind control key for air
    uc_key.resize(3);
    uc_key[2] = model.Database().StorageKey("upwind control air");
}

template class UpwindControlVisitor<1U>;
template class UpwindControlVisitor<2U>;
template class UpwindControlVisitor<3U>;

} // csmp
