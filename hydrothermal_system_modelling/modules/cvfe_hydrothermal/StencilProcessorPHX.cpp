// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "StencilProcessorPHX.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "Element.h"
#include <algorithm>

using namespace std;

namespace csmp {


/** custom destructor */
template<uint32_t dim>
StencilProcessorPHX<dim>::StencilProcessorPHX( )
{
    std::vector<double> g;
    g.resize(dim,0.0);
    if (dim==1) g[0] = -9.80665;
    else g[1] = -9.80665;
    VectorVariable<dim> gvv( g );
    gravity_ = gvv;
} // end StencilProcessorPHX


/** default destructor */
template<uint32_t dim>
StencilProcessorPHX<dim>::~StencilProcessorPHX()
{
} // end destructor


/** Calculating fluxes out of the control volumes, including gravity */
template<uint32_t dim>
void StencilProcessorPHX<dim>::DetermineFluxOutWithGravity( const FV_Parameter& param,
                                                           const Element<dim>& e,
                                                           DenseMatrix<DM_MIN>& upwind,
                                                           std::vector<std::vector<double> >& facet_flux,
                                                           std::vector<double>& flux_out,
                                                           csmp::Index rhs_key,
                                                           csmp::Index rho_key,
                                                           csmp::Index k_key,
                                                           csmp::Index thickness_key)// Benoit 2025 add
{
    eidx_ = e.Idx();
    e.Read( k_key, perm_ );
    e.Read( thickness_key, thickness_ );//Benoit 2025 add

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
    {
        grav_ = param.FacetNormalProjection(i,gravity_);
        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        if ( upwind(inside_node_,outside_node_) == 1 )
        {
            e.N(inside_node_)->Read( rho_key, density_ );
            g_component_ = grav_ * density_() * perm_();
            vel_ =  param.FacetNormalVelocity(i) + g_component_;

            if ( vel_ < 0.0 ) vel_ = 0.0;

            e.N(inside_node_)->Read( rhs_key, rhs_property_ );
        }
        else if ( upwind(inside_node_,outside_node_) == 2 )
        {
            e.N(outside_node_)->Read( rho_key, density_ );
            g_component_ = grav_ * density_() * perm_();
            vel_ =  param.FacetNormalVelocity(i) + g_component_;

            if ( vel_ > 0.0 ) vel_ = 0.0;

            e.N(outside_node_)->Read( rhs_key, rhs_property_ );
        }
        else vel_ =0.0;

        if ( vel_ != 0.0 )
        {
            facet_flux[eidx_][i] = vel_ * param.FacetArea(i) * rhs_property_();
            facet_flux[eidx_][i] *= thickness_();//Benoit 2025 add

            if ( facet_flux[eidx_][i] < 0 )
                flux_out[ e.N(outside_node_)->Idx() ] -= facet_flux[eidx_][i];
            else
                flux_out[ e.N( inside_node_)->Idx() ] += facet_flux[eidx_][i];
        }
    }

} // end DetermineFluxOutWithGravity


/** Calculating fluxes out of the control volumes without a gravity component */
template<uint32_t dim>
void StencilProcessorPHX<dim>::DetermineFluxOutWithoutGravity( const FV_Parameter& param,
                                                              const Element<dim>& e,
                                                              DenseMatrix<DM_MIN>& upwind,
                                                              std::vector<std::vector<double> >& facet_flux,
                                                              std::vector<double>& flux_out,
                                                              csmp::Index rhs_key,
                                                              csmp::Index thickness_key)// Benoit 2025 add
{

    eidx_ = e.Idx();
    e.Read( thickness_key, thickness_ );//Benoit 2025 add

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
    {

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
        vel_ =  param.FacetNormalVelocity(i);

        if (upwind(inside_node_,outside_node_) == 1)
        {
            if (vel_ < 0.0) vel_ = 0.0;
            e.N( inside_node_ )->Read( rhs_key, rhs_property_ );
        }
        else if (upwind(inside_node_,outside_node_) == 2)
        {
            if ( vel_ > 0.0 ) vel_ = 0.0;
            e.N( outside_node_ )->Read( rhs_key, rhs_property_ );
        }
        else vel_ = 0.0;

        if ( vel_ != 0.0 )
        {
            facet_flux[eidx_][i] = vel_ * param.FacetArea(i) * rhs_property_();
            facet_flux[eidx_][i] *= thickness_();//Benoit 2025 add

            if (facet_flux[eidx_][i]<0)
                flux_out[ e.N( outside_node_ )->Idx() ] -= facet_flux[eidx_][i];
            else
                flux_out[ e.N( inside_node_ )->Idx() ] += facet_flux[eidx_][i];
        }
    }

} // end DetermineFluxOutWithoutGravity

/** Calculating fluxes into the control volumes */
template<uint32_t dim>
void StencilProcessorPHX<dim>::DetermineFluxIn( const Element<dim>& e,
                                               std::vector<std::vector<double> >& facet_flux,
                                               std::vector<double>& flux_in,
                                               std::vector<double>& mass_balance)
{

    eidx_ = e.Idx();

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
    {

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        if (facet_flux[eidx_][i]<0)
            flux_in [ e.N( inside_node_ )->Idx() ] += facet_flux[eidx_][i]*mass_balance[ e.N( outside_node_ )->Idx() ];
        else
            flux_in [ e.N( outside_node_ )->Idx() ] -= facet_flux[eidx_][i]*mass_balance[ e.N( inside_node_ )->Idx() ];

    }
} // end DetermineFluxIn


/** Rebuild a content property's facet fluxes + node outflow from the mass
    facet fluxes at the donor's CURRENT content/mass ratio (see header).

    (Benoit 18/08/2026) This REPLACES the mobility-driven route for every content
    property. The two are algebraically identical when the mobility is defined
    consistently — for a donor d,
        mf * C[d]/ml[d] = vel*A*t * (krl*rl/mul) * C/(rl*sl) = vel*A*t * krl*C_V/mul
    i.e. exactly vel*A*t*mobility — but the rebuild evaluates the ratio at the
    CURRENT sub-step holdup instead of the stage-start mobility, and gives every
    content the same donor and the same clamp beat as the mass it rides on
    (content_flux_out = ratio * mass_flux_out exactly, so it clamps by the same
    factor). This is what keeps mass and enthalpy/salt coupled at clamped cells.

    ASSUMPTION / TRIPWIRE: the species is carried at the donor's BULK ratio. Any
    species whose effective velocity differs from the carrier phase cannot be
    represented here — sorption/retardation, a species-specific dispersion, an
    immobile fraction stored inside the same content variable, or a
    concentration-dependent mobility would all be silently advected at the
    carrier's speed. Nothing currently does this (lithium's precipitated phase
    lives in a separate "lithium solid" variable, so the transported content is
    fully mobile). If that changes, this function is the place to revisit. */
template<uint32_t dim>
void StencilProcessorPHX<dim>::RebuildContentFluxFromMass( const Element<dim>& e,
                                                          const std::vector<std::vector<double> >& mass_facet_flux,
                                                          std::vector<std::vector<double> >& content_facet_flux,
                                                          std::vector<double>& content_flux_out,
                                                          const std::vector<double>& content_holdup,
                                                          const std::vector<double>& mass_holdup )
{
    eidx_ = e.Idx();

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
    {
        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        const double mf = mass_facet_flux[eidx_][i];
        if ( mf == 0.0 ) { content_facet_flux[eidx_][i] = 0.0; continue; }

        // donor convention as in DetermineFluxIn: mf<0 -> outside donates
        const size_t donor = ( mf < 0.0 ) ? e.N( outside_node_ )->Idx()
                                        : e.N( inside_node_ )->Idx();
        const double m = mass_holdup[donor];
        const double ratio = ( m > 0.0 ) ? std::max( 0.0, content_holdup[donor] ) / m : 0.0;

        const double cf = mf * ratio;
        content_facet_flux[eidx_][i] = cf;

        if ( cf < 0.0 )
            content_flux_out[ e.N( outside_node_ )->Idx() ] -= cf;
        else
            content_flux_out[ e.N( inside_node_ )->Idx() ] += cf;
    }
} // end RebuildContentFluxFromMass

// /** Calculating fluxes out of the control volumes without pre-defined upwind nodes */
// template<uint32_t dim>
// void StencilProcessorPHX<dim>::DetermineFluxOut(const FV_Parameter& param,
//                                                 const Element<dim>& e,
//                                                 std::vector<std::vector<double> >& facet_flux,
//                                                 std::vector<double>& flux_out,
//                                                 csmp::Index rhs_key)
// {

//     const double zero{0.};
//     ScalarVariable rhs_property;
//     eidx_ = e.Idx();

//     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
//     {

//         e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

//         if ( param.FacetNormalVelocity(i) != zero ) {

//             const double n_x_A(param.FacetNormalVelocity(i) * param.FacetArea(i));
//             // identifying the upstream node and initialize variables for it
//             const uint32_t  nidx( (param.FacetNormalVelocity(i) < zero) ? outside_node_ : inside_node_ );

//             if (rhs_key.place == NODE) e.N(nidx)->Read( rhs_key, rhs_property );
//             else if (rhs_key.place == ELEMENT) e.Read( rhs_key, rhs_property );
//             facet_flux[eidx_][i] = n_x_A * rhs_property();

//             if (facet_flux[eidx_][i]<0)
//                 flux_out[ e.N(outside_node_)->Idx() ] -= facet_flux[eidx_][i];
//             else
//                 flux_out[ e.N( inside_node_)->Idx() ] += facet_flux[eidx_][i];

//         }
//     }

// } // end DetermineFluxOut

template class StencilProcessorPHX<1U>;
template class StencilProcessorPHX<2U>;
template class StencilProcessorPHX<3U>;


} // namespace csmp