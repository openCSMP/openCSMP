// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_FLUX_EVALUATOR_H
#define CSMP_FLUX_EVALUATOR_H

#include "VectorVariable.h"
#include "DenseMatrix.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t dim,template<uint32_t> class CELL> class ModelSubDomain;

/**
@class FluxEvaluator FluxEvaluator  "reservoir_simulator/FluxEvaluator.h"

\brief     Policy class that computes fluxes across finite-volume facets.
\details   Part of trial implementation of Colleoli transport scheme.
\author    Stephan K. Matthai
\version   1a
\date      21/2/2013
\pre       high-level class depending on CSMP++ API
\copyright Stephan K. Matthai

The transported variable on current and future time level is denoted as C0 and C1, respectively.

All interim results are stored on the model.

@todo   double  FluxBalanceFromGradP_K_mu_C0( Node<dim>* const ) const;
@todo   double  FluxBalanceFromGradP_KK_mu_C0( Node<dim>* const ) const;
@todo   include dispersion and options for diffusive fluxes

*/
template<uint32_t dim, template<uint32_t> class USER>
class FluxEvaluator {
  public:
    FluxEvaluator() { /* assumes external initialisation of variables in base class */ };
  
    // 1st-ORDER FLUX PROJECTIONS AND INTEGRATIONS

    /// computes F = (A_i vD . n_i) and its products with the upstream C0 on all facets, also computes FV flux balances
    void   Advective_O1_FluxesAndBalances( Element<dim>* const ) const;
  
    /// computes F = (A_i vD . n_i)  and its products with transported variables, also computes FV flux balances
    void  Advective_O1_FluxesAndBalances( Node<dim>* const ) const;

// TODO: tensor permeability fluxes and higher-order in space approximations of fluxes

    /// computes  volumetric flows and (chemical) fluxes across facets using FacetFlux (facet flux) and stores them in target region
    template<template<uint32_t> class CELL>
    void VolumetricFlowAndTransportVariableFluxBalances( ModelSubDomain<dim,CELL>&, 
                                                         std::vector<CELL<dim>*>& halo_stencils );

    // COMPUTATIONS DEPENDENT ON PRECOMPUTED FACET FLUXES

    /// computes flux balances on boundary FVs, using precomputed facet fluxes; volumetric balances are set to zero
    void  FluxBalancesAtBoundary( Node<dim>* const ) const;

    /// computes first-order facet fluxes, F = ff * C0;  ff=(pre-computed) volumetric flow, C0=transport variable
    void  TransportVariableFluxes( Element<dim>* const ) const;

    /// like TransportVariableFluxes(), but computes FV balances of F as well, storing them in 'accumulation'.
    void  TransportVariableFluxesAndBalances( Element<dim>* const ) const;

    /// computes balance of transport variable fluxes from precomputed facet fluxes and transport variable
    void  FluxBalancesFromFacetFluxes( Node<dim>* const ) const;

     /// stores FV flux balances computed from current facet fluxes and C0 concentrations; returns outflow from cell
    double  FluxBalanceAndOutFlow( Node<dim>* const ) const;
  
    /// computes the volumetric flow into the current cell; @return inflow which is always positive
    double  InFlow( const Node<dim>* const ) const;
  
    /// computes the volumetric flow outside of the current cell; @return outflow which is always positive
    double  OutFlow( const Node<dim>* const ) const;

    /// computes the volumetric flow outside of the current cell; @return outflow which is always positive
    double  VolumetricFlowBalance( const Node<dim>* const ) const;

  private:
    /// shorthand for accessing the class that FluxEvaluator is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
  private:
    mutable ScalarVariable       sc_;
    mutable VectorVariable<dim>  vD_, nrml_, gradC_; ///< for efficient re-use
    mutable Point<dim>           pt_;
    mutable DenseMatrix<DM_MIN>  DN_;
};

} // end csmp

#endif /* CSMP_FLUX_EVALUATOR_H */
