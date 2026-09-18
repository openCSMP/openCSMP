// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STENCIL_PROCESSOR_PHX_H
#define STENCIL_PROCESSOR_PHX_H

/*   Changelog
     February 2015, Philipp Weis:
     - initial port to CSMP++.
*/

#include "Model.h"
#include "FV_Parameter.h"

namespace csmp {

template<uint32_t> class Element;

/**
     @class StencilProcessorPHX StencilProcessorPHX.h

     Stencil processor for finitive volume calculations within CVFEM scheme (Weis et al., Geofluids, 2014).

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes

     @section motivation Motivation
      Specialized stencil processor for PHX scheme.

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
     Performs finite volume fluxes across facets.

     @endcode

     @section dependencies Dependencies
     The functions are tailored for use in ExplicitFiniteVolumeTransportPHX

     @section issues Known issues

     @section testing Testing
     testing was done in the period before publication in 2014.

  */

template<uint32_t dim>
class StencilProcessorPHX {
public:
    StencilProcessorPHX( );

    void DetermineFluxOutWithGravity( const FV_Parameter& param,
                                     const Element<dim>& e,
                                     DenseMatrix<DM_MIN>& upwind, //matrix indicating pre-defined upwind nodes
                                     std::vector<std::vector<double> >& facet_flux,  //vector storing flux across each facet
                                     std::vector<double>& flux_out, // vector storing the flux out of the control volumes
                                     csmp::Index rhs_key, // key for the right-hand side variable of finite volume calculations
                                     csmp::Index rho_key, // density of the transported phase
                                     csmp::Index k_key, // element permeability, Benoit: this has to be vertical permeability if we have anisotropic permeability
                                     csmp::Index thickness_key);  //Benoit 2025 add

    void DetermineFluxOutWithoutGravity( const FV_Parameter& param,
                                        const Element<dim>& e,
                                        DenseMatrix<DM_MIN>& upwind, //matrix indicating pre-defined upwind nodes
                                        std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
                                        std::vector<double>& flux_out, // vector storing the flux out of the control volumes
                                        csmp::Index rhs_key,  // density of the transported phase
                                        csmp::Index thickness_key);  //Benoit 2025 add

    void DetermineFluxIn( const Element<dim>& e,
                         std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
                         std::vector<double>& flux_in, // vector storing the flux into of the control volumes
                         std::vector<double>& mass_balance); // vector for mass balance corrections

    // ── content-flux rebuild (2026-07, mass/energy consistency) ───────────
    // Rebuild a CONTENT property's facet fluxes and node outflow from the MASS
    // facet fluxes and the donor's CURRENT content/mass holdup ratio:
    //   content_facet = mass_facet * content_holdup(donor)/mass_holdup(donor)
    // Donor convention identical to DetermineFluxIn: facet_flux < 0 -> donor is
    // the outside node, > 0 -> inside node. Rationale: the rhs variables that
    // feed DetermineFluxOut* are written once per outer step, so per-property
    // flux ratios are pinned at stage-start while holdups evolve over the
    // sub-cycle — mass and enthalpy/salt then clamp and drain inconsistently //Benoit: are we sure it can drain inconsistently??
    // (mass with zero energy, etc.). Rebuilding contents from the mass flux at
    // the CURRENT donor ratio keeps every property an exact scaled copy of the
    // mass field: contents empty exactly when mass empties, and the unified
    // AdjustFluxOut factor is equal across properties by construction.
    // Zero-mass donor -> zero content flux (no phantom content export).
    void RebuildContentFluxFromMass( const Element<dim>& e,
                                    const std::vector<std::vector<double> >& mass_facet_flux,
                                    std::vector<std::vector<double> >& content_facet_flux,
                                    std::vector<double>& content_flux_out,      // accumulated (rate units)
                                    const std::vector<double>& content_holdup,  // property_vectors[i]
                                    const std::vector<double>& mass_holdup );   // property_vectors[0]

    // void DetermineFluxOut( const FV_Parameter& param,
    //                         const Element<dim>& e,
    //                         std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
    //                         std::vector<double>& flux_out, // vector storing the flux out of the control volumes
    //                         csmp::Index rhs_key); // density of the transported phase

    ~StencilProcessorPHX();

private:

    VectorVariable<dim> gravity_;

    size_t            eidx_;
    mutable uint32_t  inside_node_, outside_node_;

    ScalarVariable    rhs_property_, density_, perm_;
    double            grav_, g_component_, vel_;

    ScalarVariable    thickness_;

};

} // namespace csmp

#endif