#ifndef FINITE_VOLUME_AUXILIARY_FUNCTIONS_H
#define FINITE_VOLUME_AUXILIARY_FUNCTIONS_H

#include "CSMP_mathUtilities.h"

namespace csmp {

/**
@file finiteVolumeAuxiliaryFunctions.h
*/

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;

/**
@addtogroup CSMPglobalFunctions
*/

double64  diffusionVelocity( const Region<1>& super_group,
                             const Node<1U>* const nd,
                             const std::set<size_t>& ngraph_entry,
                             const csmp::Index& advected_var_key );

double64  diffusionVelocity( const Region<2>&,
                             const Node<2U>* const,
                             const std::set<size_t>&,
                             const csmp::Index& );

double64  diffusionVelocity( const Region<3>&,
                             const Node<3U>* const,
                             const std::set<size_t>&,
                             const csmp::Index& );

template<size_t dim>
double64 fluxThroughFiniteVolume( Node<dim> const& node, Index const& velocityKey );

template<size_t dim, typename ForwardIt>
double64 fluxThroughFiniteVolumes( ForwardIt nodesBegin, ForwardIt nodesEnd, Index const& velocityKey );

double64 delta_X_FromFV_Volume( double64 FV_volume, size_t dim );

double64 omega( double64 dt, double64 pore_vol, double64 src );

double64 NVD_Function( double64 xi, double64 U_f, double64 U_c );

/// minmod limiter with coefficient xi that controls maximum permissible gradient
double64 limitProperty( double64 psi_hat_c, double64 psi_hat_d, double64 psi_facet,
                        const std::pair<double64,double64>& SMINMAX_at_upstream_node,
                        double64 limiter_value_xi=2. );

double64 limitProperty( double64 psi_hat_c, double64 psi_hat_u,
                        double64 min, double64 max );

/// limiting with node-based property gradients, see Geiger et al. (2004, Geofluids)
template <size_t dim>
void limitProperty_LSMGRAD( const Element<dim>& e,
                            const csmp::Index& mass_center_key,
                            const csmp::Index& grad_sn_key,
                            const csmp::Index& grad_sn_limiter_key,
                            size_t inside_node, size_t outside_node, size_t iFacet,
                            const double64 sn_inside_node, const double64 sn_outside_node,
                            double64& limited_sn_inside_node, double64& limited_sn_outside_node );

/// limiting with node-based property gradients, see Geiger et al. (2004, Geofluids)
template <size_t dim>
double64 limitProperty_LSMGRAD( const Element<dim>& e,
                            const csmp::Index& mass_center_key,
                            const csmp::Index& grad_sn_key,
                            const csmp::Index& grad_sn_limiter_key,
                            size_t ustream_node, size_t iFacet,
                            const double64 sn_upstream_node);

/// Matthai et al. 2009 (TIPM) explains theta limiting concept
double64 thetaLimiter( double64 hf_t0, double64 hf_t1,
                       double64 vol_upstr, double64 vol_dwstr,
                       double64 psi_hat_c_t1_minus_psi_hat_c_t0,
                       double64 psi_hat_d_t1_minus_psi_hat_d_t0,
                       double64 dt ); 

/// gradient limiter for the solution of diffusion equations with the FVM method
double64 diffusionLimiter( double64 grad_psi_O1, double64 grad_psi_O2 );

/**
@}
*/




} // end namespace csmp

#endif
