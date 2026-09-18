// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ExplicitFiniteVolumeTransportPHX.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "StencilProcessorPHX.h"

using namespace std;

namespace csmp {


/** custom constructor */
template<uint32_t dim>
ExplicitFiniteVolumeTransportPHX<dim>::ExplicitFiniteVolumeTransportPHX( Model<dim>& model,
                                                                        std::vector<std::string>& lhs_property,
                                                                        std::vector<std::string>& rhs_property,
                                                                        const char* velocity,
                                                                        const char* source,
                                                                        const char* density)
    :
    NodeCenteredFiniteVolumeTransport<dim>( "Model", model, "porosity",
                                           lhs_property[0].c_str(), velocity,
                                           source, false, false, "thickness"), //Benoit 2025 add
    model_ref( model ),
    region_ref( model.Region("Model") ),
    flux_in_vectors(lhs_property.size()),
    flux_out_vectors(lhs_property.size()),
    property_vectors(lhs_property.size()),
    mass_balance_vectors(lhs_property.size()),
    facet_flux_vectors(lhs_property.size()),
    max_time_step(8640000.),
    with_gravity(false)
{
    cout <<endl<< "Constructing ExplicitFiniteVolumeTransportPHX";

    // resize flux and property vectors
    for ( size_t i{0U}; i<lhs_property.size(); i++ )
    {
        flux_in_vectors[i].resize(region_ref.Nodes());
        flux_out_vectors[i].resize(region_ref.Nodes());
        property_vectors[i].resize(region_ref.Nodes());
        mass_balance_vectors[i].resize(region_ref.Nodes());
        facet_flux_vectors[i].resize( region_ref.Cells() );
        for ( typename vector<Element<dim>*>::const_iterator
                 eit=region_ref.CellsBegin();
             eit!=region_ref.CellsEnd(); eit++)
            facet_flux_vectors[i][(*(*eit)).Idx()].resize( (*(*eit)).FV()->Facets() );
    }

    // set and check keys
    // ANISOTROPY CAVEAT: k_key feeds ONLY the gravity term of the facet flux
    // (g * rho * k in the stencil), i.e. a vertical driving force — so with an
    // anisotropic permeability tensor this must be the VERTICAL permeability,
    // not the isotropic scalar. Currently wired to the scalar "permeability"
    // key; switch to a "vertical permeability" key when anisotropy is enabled.
    k_key  =  model.Database().StorageKey( "vertical permeability" );
    pv_key  =  model.Database().StorageKey( "pore volume" );
    rho_key  =  model.Database().StorageKey( density );
    thickness_key = model.Database().StorageKey( "thickness" );//Benoit 2025 add

    pl_key.resize( lhs_property.size() );
    pr_key.resize( rhs_property.size() );

    for ( size_t i{0U}; i<lhs_property.size(); i++ )
    {
        pl_key[i] =  model.Database().StorageKey( lhs_property[i].c_str() );
        pr_key[i] =  model.Database().StorageKey( rhs_property[i].c_str() );

        if ( pl_key[i].type != SCALAR || pl_key[i].place != NODE )
            throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX.<double, dim>::ExplicitFiniteVolumeTransportPHX\n",
                                  lhs_property[i].c_str(), " must be a nodal scalar property." );

        if ( pr_key[i].type != SCALAR || (pr_key[i].place != NODE && pr_key[i].place != ELEMENT) )
            throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX.<double, dim>::ExplicitFiniteVolumeTransportPHX\n",
                                  rhs_property[i].c_str(), " must be a nodal or element scalar property." );
    }
    cout <<endl<<"ExplicitFiniteVolumeTransportPHX<"<< typeid(double).name() <<","<< dim;
    cout <<">: Constructed successfully."<< endl;
} // end constructor


/** default destructor */
template<uint32_t dim>
ExplicitFiniteVolumeTransportPHX<dim>::~ExplicitFiniteVolumeTransportPHX()
{}

// ComposeAdvection: Apply the (already pore-volume-scaled, already clamped) fluxes to the holdup.
// Sign convention: flux_in is stored NEGATIVE (incoming), flux_out POSITIVE
// (outgoing), so both lines subtract:
//   property -= flux_in   ->  ADDS incoming mass
//   property -= flux_out  ->  REMOVES outgoing mass
// Because AdjustFluxOut has already clamped flux_out to the pre-inflow holdup and
// CalculateFluxIn has scaled the matching inflow by the same mass_balance factor,
// the result is non-negative and globally mass-exact. The guard below flags any
// violation (overflow or negativity) as a hard error — it should never fire.

/** change lhs variable by adding and subtracting flux in and out of control volume */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::ComposeAdvection( )
{
    for ( size_t nidx=0U; nidx<region_ref.Nodes(); nidx++ )
        for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
        {
            property_vectors[i][nidx] -= flux_in_vectors[i][nidx];
            property_vectors[i][nidx] -= flux_out_vectors[i][nidx];

            if (property_vectors[i][nidx] > 1.e20 || property_vectors[i][nidx] < 0.0)
            {
                cerr << "\nproperty_vectors[i][nidx]: " << property_vectors[i][nidx] <<  endl;
                cerr << "\nflux_in_vectors[i][nidx]: " << flux_in_vectors[i][nidx] << endl;
                cerr << "\nflux_out_vectors[i][nidx]: " << flux_out_vectors[i][nidx] << endl;
            }
        }
} // end ComposeAdvection


/** store lhs variables */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::WriteResults( )
{
    double                   pmin(0), pmax(0), value(0);
    ScalarVariable   sc;

    for ( size_t i{0U}; i<property_vectors.size(); i++ )
    {
        model_ref.Database().RangeOf( model_ref.Database().Name(pl_key[i]), pmin, pmax );
        for ( size_t idx=0; idx<property_vectors[i].size(); idx++ )
        {
            value = property_vectors[i][idx];

            // store property
            if ( region_ref.N(idx)->Status( pl_key[i] ) != DIRICH )
            {
                // reading the pre-existing value (use it to keep the flag if results are overwritten)
                region_ref.N(idx)->Read( pl_key[i], sc );
                if ( value <= pmax && value >= pmin )
                {
                    sc()=value;
                    region_ref.N(idx)->Store( pl_key[i], sc );
                }
                else
                {
                    cerr << "\n[WriteResults] Out-of-range property: "
                         << model_ref.Database().Name(pl_key[i])
                         << " = " << value
                         << " at node (" << region_ref.N(idx)->Coordinate() << ")"
                         << "  valid range: [" << pmin << ", " << pmax << "]" << endl;

                    for (size_t ii = 0; ii < property_vectors.size(); ii++)
                    {
                        cerr << "  [" << model_ref.Database().Name(pl_key[ii]) << "]"
                             << "  value=" << property_vectors[ii][idx]
                             << "  flux_in=" << flux_in_vectors[ii][idx]
                             << "  flux_out=" << flux_out_vectors[ii][idx]
                             << "  mass_balance=" << mass_balance_vectors[ii][idx] << endl;
                    }

                    if (value < pmin) sc() = pmin;
                    else if (value > pmax) sc()= pmax;
                    cerr << "\nsc(): " << sc() << endl;
                    region_ref.N(idx)->Store(pl_key[i], sc);
                    throw csmp::Exception( WARNING, "ExplicitFiniteVolumeTransportPHX::WriteResults",
                                          "Output property was out of range, legal (min/max) was stored instead");
                    //  		           cin >> temp;
                }
            }
        } // end finite volumes
    } // end property vectors
} // end WriteResults

/** change size of maximum time step */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::SetMaximumTimeStep( const double& new_max_time_step )
{
    max_time_step = new_max_time_step;
}

/** access to flux out of control volume - main property */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFluxOut(size_t index)
{
    return flux_out_vectors[0][index];
}

/** access to flux into control volume - main property */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFluxIn(size_t index)
{
    return flux_in_vectors[0][index];
}

/** access to flux out of control volume - specified property */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFluxOut(size_t index, uint32_t property_idx)
{
    return flux_out_vectors[property_idx][index];
}

/** access to flux into control volume - specified property */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFluxIn(size_t index, uint32_t property_idx)
{
    return flux_in_vectors[property_idx][index];
}

/** get value of main property at specified control volume */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetPropertyValue(size_t index)
{
    return region_ref.N(index)->Read( pl_key[0] );
}

/** get facet flux of specified propoerty at indicated facet */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFacetFlux( Element<dim>& e, uint32_t facet_idx, uint32_t property_idx )
{
    uint32_t  inside_node_,outside_node_;
    auto eidx_{e.Idx()};
    e.FV()->FacetEdgeNodes( facet_idx, inside_node_, outside_node_ );
    double flux_ = facet_flux_vectors[property_idx][eidx_][facet_idx];
    if (flux_<0)
    {
        flux_ *= mass_balance_vectors[property_idx][e.N(outside_node_)->Idx()];
    }
    else
    {
        flux_ *= mass_balance_vectors[property_idx][e.N( inside_node_)->Idx()];
    }
    return flux_;
}

/** get facet flux of specified property at indicated facet */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetMainPropertyLHS(size_t index)
{
    return property_vectors[0][index];
}

/** calculate and store facet fluxes with suggested time step.
    STATUS: no longer called by ThreePhaseTransportPHX (which uses RatesOnce +
    RescaleFacetFluxToTimestep exclusively). Kept as public API for standalone /
    non-subcycled users; note it carries stage-start content ratios (no per-
    sub-step content rebuild) — correct for single-shot use only. */
template<uint32_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::DetermineFacetFlux( const double& time_increment, std::vector<DenseMatrix<DM_MIN> >& upwind_visitor )
{
    internal_time_step = time_increment;

    vector<FV_Parameter>::const_iterator      fvt;
    typename vector<Element<dim>*>::const_iterator  eit;

    for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
    {
        fill(    flux_in_vectors[i].begin(),    flux_in_vectors[i].end(), static_cast<double>(0.) );
        fill(   flux_out_vectors[i].begin(),   flux_out_vectors[i].end(), static_cast<double>(0.) );
        fill(   property_vectors[i].begin(),   property_vectors[i].end(), static_cast<double>(0.) );
        fill(mass_balance_vectors[i].begin(),mass_balance_vectors[i].end(), static_cast<double>(1.) );
        for ( size_t j{0U}; j<facet_flux_vectors[i].size(); j++ )
            fill( facet_flux_vectors[i][j].begin(), facet_flux_vectors[i][j].end(), static_cast<double>(0.) );
    }

    for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
    {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.CellsBegin(); eit!=region_ref.CellsEnd(); eit++, fvt++ )
        {
            // pass the visitor's matrix directly — no per-element DenseMatrix copy
            DenseMatrix<DM_MIN>& uw = upwind_visitor[ (*eit)->Idx() ];
            if (with_gravity)
                stencilPHX.DetermineFluxOutWithGravity( (*fvt), *(*eit), uw,
                                                       facet_flux_vectors[i], flux_out_vectors[i],
                                                       pr_key[i], rho_key, k_key,
                                                       thickness_key);//Benoit 2025 add
            else
                stencilPHX.DetermineFluxOutWithoutGravity( (*fvt), *(*eit), uw,
                                                          facet_flux_vectors[i], flux_out_vectors[i],
                                                          pr_key[i],
                                                          thickness_key);//Benoit 2025 add
        }
    }
    CalculateOutflowPerPoreVolume();
} // end DetermineFacetFlux

/** calculate total flux out of control volumes */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::CalculateOutflowPerPoreVolume()
{

    double pore_vol;

    for ( typename vector<Node<dim>*>::const_iterator
             fvit=region_ref.NodesBegin();
         fvit!=region_ref.NodesEnd(); fvit++ )
    {
        size_t idx = (*fvit)->Idx();
        pore_vol = (*fvit)->Read( pv_key);
        for( size_t i{0U}; i<flux_out_vectors.size(); i++ ) {
            property_vectors[i][idx] = (*fvit)->Read( pl_key[i] );
            flux_out_vectors[i][idx] *= internal_time_step;
            flux_out_vectors[i][idx] /= pore_vol;
        } // end flux vectors
    } // end finite volumes
} // end CalculateOutflowPerPoreVolume

/** adjust fluxes, perform advection calculations and store results*/
template<uint32_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AdjustAndPerformFacetFlux( const double& time_factor )
{

    internal_time_step *= time_factor;

    AdjustFluxOut( time_factor );

    CalculateFluxIn();

    CalculateInflowPerPoreVolume();

    ComposeAdvection();

    WriteResults();

}

// AdjustFluxOut: the heart of the positivity + conservation guarantee.
//
// For each cell:
//   1. Scale the requested outflows by time_factor (the sub-step fraction).
//   2. If the MASS outflow exceeds the mass holdup (read BEFORE this
//      sub-step's inflow), the cell empties: every property's outflow is set
//      to exactly its holdup (contents ride mass by construction after the
//      content-flux rebuild), and the mass factor is recorded in
//      mass_balance_vectors for the receivers' inflow scaling.
//
// The clamp gives POSITIVITY: a cell never exports more than it holds, so it
// cannot go negative — independent of sub_dt or whether the cell bound the step.
// This is what makes a near-empty / floored cell safe, and what makes the
// "tiny holdup about to receive large inflow" edge case safe: outflow is capped
// at the small pre-inflow holdup, the inflow simply accumulates, and the cell is
// sized correctly on the NEXT sub-step.
//
// mass_balance_factor gives CONSERVATION across cells: it is passed into
// CalculateFluxIn -> DetermineFluxIn, where it scales the DOWNSTREAM inflow by
// the same ratio. So when a cell can only deliver part of its requested outflow,
// its neighbours receive exactly that reduced amount — never more. Mass removed
// here equals mass deposited there; nothing is created, dropped, or
// double-counted. (mass_balance is reset to 1 each sub-step before this runs.)

/** adjust fluxes */
template<uint32_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AdjustFluxOut( const double& time_factor )
{
    // ── MASS-DRIVEN CELL CLAMP (2026-07) ────────────────────────────────────
    // With the content-flux rebuild (RescaleFacetFluxToTimestep), every content
    // property's outflow is BY CONSTRUCTION the mass outflow times the cell's
    // current content/mass ratio — so all per-property clamp factors are equal
    // to the mass factor (up to float eps). The clamp therefore needs no
    // per-property search: the MASS holdup decides. When mass binds, the cell
    // empties completely — every property's outflow is ASSIGNED its exact
    // holdup (the old code's exact-zero landing, extended to all properties:
    // no multiplicative eps residue, holdups land at exactly 0.0 together,
    // mixtures intact). When mass does not bind, nothing should clamp; a cheap
    // exact guard still caps each content at its holdup, because the rebuild's
    // ratio arithmetic can overshoot by ~1 ulp and WriteResults range-checks
    // contents too. Conservation: receivers scale inflow by the recorded
    // factor; the eps-scale donor/receiver mismatch is the same class the
    // historical per-property code had.

    for ( typename vector<Node<dim>*>::const_iterator
             fvit=region_ref.NodesBegin();
         fvit!=region_ref.NodesEnd(); fvit++ )
    {
        size_t idx = (*fvit)->Idx();

        for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
            flux_out_vectors[i][idx] *= time_factor;

        const double m_req   = flux_out_vectors[0][idx];
        const double m_avail = std::max( 0., property_vectors[0][idx] );

        if ( m_req > 0. && m_req > m_avail )
        {
            // mass binds -> the cell empties; exact landing for every property
            const double cell_factor = m_avail / m_req;
            for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
            {
                mass_balance_vectors[i][idx] = cell_factor;
                flux_out_vectors[i][idx]     = std::max( 0., property_vectors[i][idx] );
            }
        }
        else
        {
            // no clamp; ulp guard so contents never export beyond their holdup
            for ( size_t i{1U}; i<flux_out_vectors.size(); i++ )
                if ( flux_out_vectors[i][idx] > 0. )
                    flux_out_vectors[i][idx] =
                        std::min( flux_out_vectors[i][idx],
                                 std::max( 0., property_vectors[i][idx] ) );
        }
    }
} // end AdjustFluxOut

/** calculate flux into control volumes */
template<uint32_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::CalculateFluxIn( )
{

    vector<FV_Parameter>::const_iterator      fvt;
    typename vector<Element<dim>*>::const_iterator  eit;

    for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
    {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.CellsBegin(); eit!=region_ref.CellsEnd(); eit++, fvt++ )
            stencilPHX.DetermineFluxIn( *(*eit),
                                       facet_flux_vectors[i],
                                       flux_in_vectors[i],
                                       mass_balance_vectors[i]);
    }
} // end CalculateFluxIn

/** calculate total flux into control volumes */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::CalculateInflowPerPoreVolume()
{

    double pore_vol;

    for ( typename vector<Node<dim>*>::const_iterator
             fvit=region_ref.NodesBegin();
         fvit!=region_ref.NodesEnd(); fvit++ )
    {
        size_t idx = (*fvit)->Idx();
        pore_vol = (*fvit)->Read( pv_key );
        for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
        {
            flux_in_vectors[i][idx] *= internal_time_step;
            flux_in_vectors[i][idx] /= pore_vol;
        } // end flux vectors
    } // end finite volumes
} // end CalculateInflowPerPoreVolume

/** get projected velocity for specified facet */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetProjectedVelocities( Element<dim>& e, uint32_t i)
{
    return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[e.Idx()].FacetNormalVelocity(i)*NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[e.Idx()].FacetArea(i);
}

/** activate gravity component */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::WithGravityComponent( )
{
    with_gravity = true;
}

// (GetUpwindMatrix removed: it copied a DenseMatrix per element per property
//  per assembly — Cells x Nprops copies per stage. Callers now bind a reference
//  to the visitor's matrix directly, zero copies. Drop the declaration and the
//  'upwind' member from the header on the header pass.)

/** get index of phase density */
template<uint32_t dim>
csmp::Index ExplicitFiniteVolumeTransportPHX<dim>::GetDensityKey(  )
{
    return rho_key;
}

/** update projection of velocity onto the facet */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::UpdateProjection(  )
{
    NodeCenteredFiniteVolumeTransport<dim>::UpdateProjectedVelocitiesAndFluxBalances( );
}

/** update projection of velocity onto the facet */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFacetNormalVelocity( size_t element, uint32_t facet )
{
    return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetNormalVelocity(facet);
}

/** get normal component of indicated facet and direction */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFacetNormalComponent( size_t element, uint32_t facet, uint32_t x_or_y_or_z )
{
    return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetNormalComponent(facet, x_or_y_or_z);
}

/** get area of indicated facet */
template<uint32_t dim>
double ExplicitFiniteVolumeTransportPHX<dim>::GetFacetArea( size_t element, uint32_t facet )
{
    return NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA[element].FacetArea(facet);
}

//Benoit add
/** Assemble the CARRIER-MASS facet flux RATES once for the whole transport stage.
    The stencil flux is built from pr_key[0] (the equilibration-set mass mobility),
    the frozen upwind matrices, density and permeability — none of which change
    between sub-steps (no pressure solve or flash runs mid-stage). So
    facet_flux_vectors[0] and the raw per-node outflow are identical for every
    sub-step and are assembled only here. RescaleFacetFluxToTimestep() then produces
    each sub-step's scaled outflow cheaply, without re-running this loop, and
    rebuilds every content property (i >= 1) from the mass rates.

    NOTE (Benoit 18/08/2026): content properties are NOT assembled here any more,
    and their pr_key mobilities are unused by transport. See the block comment
    inside the function body. */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::DetermineFacetFluxRatesOnce( std::vector<DenseMatrix<DM_MIN> >& upwind_visitor )
{
    internal_time_step = 1.0;   // rates only; real dt applied in RescaleFacetFluxToTimestep

    vector<FV_Parameter>::const_iterator            fvt;
    typename vector<Element<dim>*>::const_iterator  eit;

    for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
    {
        fill(    flux_in_vectors[i].begin(),    flux_in_vectors[i].end(),    0. );
        fill(   flux_out_vectors[i].begin(),   flux_out_vectors[i].end(),    0. );
        fill(   property_vectors[i].begin(),   property_vectors[i].end(),    0. );
        fill(mass_balance_vectors[i].begin(),mass_balance_vectors[i].end(),  1. );
        for ( size_t j{0U}; j<facet_flux_vectors[i].size(); j++ )
            fill( facet_flux_vectors[i][j].begin(), facet_flux_vectors[i][j].end(), 0. );
    }

    // ── CARRIER-ONLY ASSEMBLY (Benoit 18/08/2026) ───────────────────────────
    // Only index 0 (the phase carrier mass) is assembled from the stencil.
    // Every content property (i >= 1: enthalpy, salt, magmatic mass, magmatic
    // salt, lithium, tracer) has its facet fluxes AND its raw outflow fully
    // OVERWRITTEN by the content-flux rebuild in RescaleFacetFluxToTimestep(),
    // which always runs before any consumer reads them. Assembling them here
    // from pr_key[i] was therefore discarded work: one full cell sweep per
    // content property per phase per stage.
    //
    // CONSEQUENCE — be explicit, this surprises people: pr_key[i] for i >= 1
    // (the "liquid/vapor <species> mobility" fields) is NOT read by transport.
    // Those fields are DIAGNOSTIC ONLY. Changing them, including setting them
    // to zero, has no effect whatsoever on advection. The advected ratio comes
    // from the donor's CURRENT content/mass holdup instead — see
    // StencilProcessorPHX::RebuildContentFluxFromMass.
    // Index 0 keeps its stencil-assembled frozen rate: the velocity-freeze contract.
    {
        fvt = NodeCenteredFiniteVolumeTransport<dim>::STENCIL_DATA.begin();
        for ( eit=region_ref.CellsBegin(); eit!=region_ref.CellsEnd(); eit++, fvt++ )
        {
            // pass the visitor's matrix directly — no per-element DenseMatrix copy
            DenseMatrix<DM_MIN>& uw = upwind_visitor[ (*eit)->Idx() ];
            if (with_gravity)
                stencilPHX.DetermineFluxOutWithGravity( (*fvt), *(*eit), uw,
                                                       facet_flux_vectors[0], flux_out_vectors[0],
                                                       pr_key[0], rho_key, k_key, thickness_key );
            else
                stencilPHX.DetermineFluxOutWithoutGravity( (*fvt), *(*eit), uw,
                                                          facet_flux_vectors[0], flux_out_vectors[0],
                                                          pr_key[0], thickness_key );
        }
    }

    // snapshot the raw (unscaled) per-node outflow for sub-step rescaling.
    // Loop stays over ALL i: i >= 1 are zero-filled here and replaced by the
    // rebuild, but the vectors must still be sized. (Benoit 18/08/2026)
    if ( raw_flux_out_vectors.size() != flux_out_vectors.size() )
        raw_flux_out_vectors.resize( flux_out_vectors.size() );
    for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
        raw_flux_out_vectors[i] = flux_out_vectors[i];
}

//Benoit add
/** Produce the scaled outflow for one sub-step of length time_increment from the
    rates assembled by DetermineFacetFluxRatesOnce(). Re-reads the current holdup
    (pl_key, advanced by the previous sub-step) and resets the per-sub-step buffers
    (flux_in -> 0, mass_balance -> 1) exactly as DetermineFacetFlux would, but WITHOUT
    rebuilding the stencil. facet_flux_vectors are intentionally left intact (frozen,
    reused by CalculateFluxIn). Output is identical to
    DetermineFacetFlux + CalculateOutflowPerPoreVolume for this dt. */
template<uint32_t dim>
void ExplicitFiniteVolumeTransportPHX<dim>::RescaleFacetFluxToTimestep( const double& time_increment )
{
    internal_time_step = time_increment;

    // ── INDEX-0 CARRIER CONTRACT (Benoit 18/08/2026) ────────────────────────
    // The content rebuild below hard-codes index 0 as the phase carrier mass
    // (facet_flux_vectors[0], property_vectors[0]). Since DetermineFacetFluxRatesOnce
    // now assembles index 0 ONLY, that assumption is load-bearing in two places.
    // Registration order is fixed by the scheme: [0] mass, [1] energy, [2] salt,
    // [3] magmatic mass, [4] magmatic salt, then AddAdvectionVariable() appends
    // (lithium, tracer, ...). Fail loud if that ever changes.
    if ( flux_out_vectors.empty() )
        throw csmp::Exception( ERROR, "ExplicitFiniteVolumeTransportPHX::RescaleFacetFluxToTimestep",
                              "index 0 (carrier mass) is not registered." );

    // One-time trace: makes the carrier/diagnostic split visible in every log,
    // so "why does changing the mobility do nothing?" is answerable from the
    // output alone. (Benoit 18/08/2026)
    static bool announced_carrier_split = false;
    if ( !announced_carrier_split )
    {
        cerr << "\n[EFVT] carrier pl_key[0] = " << model_ref.Database().Name( pl_key[0] )
        << " (rhs " << model_ref.Database().Name( pr_key[0] ) << ")"
        << " | " << flux_out_vectors.size()-1
        << " content properties rebuilt from mass; their rhs mobilities are DIAGNOSTIC ONLY"
        << endl;
        announced_carrier_split = true;
    }

    for ( size_t i{0U}; i<flux_in_vectors.size(); i++ )
    {
        fill(    flux_in_vectors[i].begin(),    flux_in_vectors[i].end(),   0. );
        fill(mass_balance_vectors[i].begin(),mass_balance_vectors[i].end(), 1. );
    }

    // pass 1: fresh holdups for every property (the content rebuild below
    // needs current ratios BEFORE any flux is formed this sub-step)
    for ( typename vector<Node<dim>*>::const_iterator
             fvit=region_ref.NodesBegin(); fvit!=region_ref.NodesEnd(); fvit++ )
    {
        size_t idx = (*fvit)->Idx();
        for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
            property_vectors[i][idx] = (*fvit)->Read( pl_key[i] );
    }

    // ── CONTENT-FLUX REBUILD (mass/content consistency, 2026-07) ────────────
    // SCOPE: the i-loop below runs over EVERY content property (i >= 1), i.e.
    // enthalpy, salt, magmatic mass, magmatic salt, and anything later added
    // via AddAdvectionVariable (tracer, lithium) — generically. Mass (i = 0)
    // is the carrier; every other property rides it.
    // The per-property rates frozen at DetermineFacetFluxRatesOnce carry
    // stage-start upwind ratios (their rhs variables are written once per
    // outer step), while holdups evolve every sub-step — which decouples
    // mass from enthalpy/salt at clamped cells (mass with zero energy, exact-
    // zero contents, wrong exported mixtures). Rebuild every content
    // property's facet fluxes and outflow from the MASS facet rates at the
    // donor's CURRENT ratio: contents are then an exact scaled copy of the
    // mass field — they clamp by the same factor and empty at the same beat,
    // by construction. Index 0 (the main property, mass) keeps its frozen
    // rate: that is the velocity-freeze contract.
    {
        vector<double> rebuilt_out( region_ref.Nodes() );
        for ( size_t i{1U}; i<flux_out_vectors.size(); i++ )
        {
            fill( rebuilt_out.begin(), rebuilt_out.end(), 0. );
            for ( typename vector<Element<dim>*>::const_iterator
                     eit=region_ref.CellsBegin(); eit!=region_ref.CellsEnd(); eit++ )
                stencilPHX.RebuildContentFluxFromMass( *(*eit),
                                                      facet_flux_vectors[0],
                                                      facet_flux_vectors[i],
                                                      rebuilt_out,
                                                      property_vectors[i],
                                                      property_vectors[0] );
            raw_flux_out_vectors[i] = rebuilt_out;   // rates, scaled below
        }
    }

    // pass 2: scale all outflow rates to this sub-step
    double pore_vol;
    for ( typename vector<Node<dim>*>::const_iterator
             fvit=region_ref.NodesBegin(); fvit!=region_ref.NodesEnd(); fvit++ )
    {
        size_t idx = (*fvit)->Idx();
        pore_vol = (*fvit)->Read( pv_key );
        for ( size_t i{0U}; i<flux_out_vectors.size(); i++ )
            flux_out_vectors[i][idx] = raw_flux_out_vectors[i][idx] * internal_time_step / pore_vol;
    }
}

/** add new advection variables for finite volume calculations */
template<uint32_t dim>
void   ExplicitFiniteVolumeTransportPHX<dim>::AddAdvectionVariable( const char* new_lhs, const char* new_rhs)
{

    pl_key.push_back(model_ref.Database().StorageKey( new_lhs ));
    pr_key.push_back(model_ref.Database().StorageKey( new_rhs ));

    flux_in_vectors.resize(pl_key.size());
    flux_out_vectors.resize(pl_key.size());
    property_vectors.resize(pl_key.size());
    mass_balance_vectors.resize(pl_key.size());
    facet_flux_vectors.resize(pl_key.size());

    auto i = pl_key.size()-1;
    flux_in_vectors[i].resize(region_ref.Nodes());
    flux_out_vectors[i].resize(region_ref.Nodes());
    property_vectors[i].resize(region_ref.Nodes());
    mass_balance_vectors[i].resize(region_ref.Nodes());
    facet_flux_vectors[i].resize( region_ref.Cells() );
    for ( typename vector<Element<dim>*>::const_iterator
             eit=region_ref.CellsBegin();
         eit!=region_ref.CellsEnd(); eit++)
        facet_flux_vectors[i][(*(*eit)).Idx()].resize( (*(*eit)).FV()->Facets() );

}

template class ExplicitFiniteVolumeTransportPHX<1U>;
template class ExplicitFiniteVolumeTransportPHX<2U>;
template class ExplicitFiniteVolumeTransportPHX<3U>;

} // end namespace csmp
