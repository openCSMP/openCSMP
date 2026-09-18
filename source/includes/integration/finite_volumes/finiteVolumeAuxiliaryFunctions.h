// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_FINITE_VOLUME_AUXILIARY_FUNCTIONS_H
#define CSMP_FINITE_VOLUME_AUXILIARY_FUNCTIONS_H

#include "CSMP_mathUtilities.h"

namespace csmp {

/**
@file finiteVolumeAuxiliaryFunctions.h
*/

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t dim, template<uint32_t> class CELL> class ModelSubDomain;
template<uint32_t> class Region;
template<uint32_t> class Model;


/**
@addtogroup CSMPglobalFunctions
*/

    /// identifies "halo" elements/faces/interfaces, i.e. which contribute to domain FVs, but are outside of domain, returns number
//template<uint32_t dim, template<uint32_t> class CELL> TODO: implement efficient version
//size_t collectHaloStencils( const ModelSubDomain<dim,CELL>&, std::vector<CELL<dim>*>& halo_stencils );

/**
    taking into account element thickness and total velocity, initialises:
      finite volume, sector volume, sector pore volume, FV pore volume, facet area, facet normal,
      facet flux, flux balance.
*/
template<uint32_t dim> void initializeFiniteVolumeProperties( Model<dim>&, Region<dim>&, bool initialize_flux );

/// taking into account thickness, initialising FV sector volume, FV pore volume, facet area, facet normal
template<uint32_t dim> void initializeBasicFiniteVolumeProperties( Model<dim>&, Region<dim>& );

/// sums precomputed fluxes over the facets that surround the FE - FV sector; thickness of dim-1 elements is taken into account
template<uint32_t dim> double sectorFlux( const Element<dim>* const eptr, uint32_t sector, const csmp::Index& flux_key );

template<uint32_t dim>
double fluxThroughFiniteVolume( Node<dim> const& node, Index const& velocityKey );

template<uint32_t dim, typename ForwardIt>
double fluxThroughFiniteVolumes( ForwardIt nodesBegin, ForwardIt nodesEnd, Index const& velocityKey );


double  diffusionVelocity( const Region<1>& model_domain,
                             const Node<1U>* const nd,
                             const std::set<uint32_t>& ngraph_entry,
                             const csmp::Index& advected_var_key );

double  diffusionVelocity( const Region<2>&,
                             const Node<2U>* const,
                             const std::set<uint32_t>&,
                             const csmp::Index& );

double  diffusionVelocity( const Region<3>&,
                             const Node<3U>* const,
                             const std::set<uint32_t>&,
                             const csmp::Index& );

double delta_X_FromFV_Volume( double FV_volume, uint32_t dim );

double omega( double dt, double pore_vol, double src );

double NVD_Function( double xi, double U_f, double U_c );

/// minmod limiter with coefficient xi that controls maximum permissible gradient
double limitProperty( double psi_hat_c, double psi_hat_d, double psi_facet,
                        const std::pair<double,double>& SMINMAX_at_upstream_node,
                        double limiter_value_xi=2. );

double limitProperty( double psi_hat_c, double psi_hat_u,
                        double min, double max );

/// limiting with node-based property gradients, see Geiger et al. (2004, Geofluids)
template <uint32_t dim>
void limitProperty_LSMGRAD( const Element<dim>& e,
                            const csmp::Index& mass_center_key,
                            const csmp::Index& grad_sn_key,
                            const csmp::Index& grad_sn_limiter_key,
                            uint32_t inside_node, uint32_t outside_node, uint32_t iFacet,
                            const double sn_inside_node, const double sn_outside_node,
                            double& limited_sn_inside_node, double& limited_sn_outside_node );

/// limiting with node-based property gradients, see Geiger et al. (2004, Geofluids)
template <uint32_t dim>
double limitProperty_LSMGRAD( const Element<dim>& e,
                            const csmp::Index& mass_center_key,
                            const csmp::Index& grad_sn_key,
                            const csmp::Index& grad_sn_limiter_key,
                            uint32_t ustream_node, uint32_t iFacet,
                            const double sn_upstream_node);

/// Matthai et al. 2009 (TIPM) explains theta limiting concept
double thetaLimiter( double hf_t0, double hf_t1,
                       double vol_upstr, double vol_dwstr,
                       double psi_hat_c_t1_minus_psi_hat_c_t0,
                       double psi_hat_d_t1_minus_psi_hat_d_t0,
                       double dt ); 

/// gradient limiter for the solution of diffusion equations with the FVM method
double diffusionLimiter( double grad_psi_O1, double grad_psi_O2 );

/// computes the diameter  and vertical of a surface finite volume using midpoints to connected nodes on the same surface
template<uint32_t dim>
std::pair<double,double>  diameterAndVerticalExtentOfLowerDimensional_FV( const Node<dim>* const );

/**
@}
*/




} // end namespace csmp

#endif
