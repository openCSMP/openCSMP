#ifndef FACET_FLUX_TRACER_TRANSFER_EXPLICIT_H
#define FACET_FLUX_TRACER_TRANSFER_EXPLICIT_H

#include "VectorVariable.h"
#include "DenseMatrix.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> struct VariableSet_TracerTransferExplicit;



/**

@brief  Policy class that computes fluxes across finite-volume facets;
for application to individual FVs using terms computed on linear FEs.

\details   Part of trial implementation of Colleoli transport scheme.
\author    Stephan K. Matthai
\version   1a
\date      21/2/2013
\pre       high-level class depending on CSMP++ API
\bug       NONE
\warning   NONE
\copyright Stephan K. Matthai

*/
template<size_t dim, template<size_t> class USER>
class FacetFlux_TracerTransferExplicit {
  public:
    FacetFlux_TracerTransferExplicit() { /* assumes external initialisation of variables in base class */ };
  
    /// initialisation
    void FacetNormalPermeabilities();
    void UpwindDirection();
  
    /// computes (A_i vD . n_i) * upstream C on all facets in element stencil and stores them there
    void      Advective_O1_FluxesInterior( bool reuse_previous_velocity, Element<dim>* ) const;
    double64  Advective_O1_FluxesAtBoundary( Node<dim>* ) const;
    void      Advective_O2_FluxesInterior( bool reuse_previous_velocity, Element<dim>* ) const;
    double64  Advective_O2_FluxesAtBoundary( Node<dim>* ) const;
  
    /// stores and returns FV flux balance computed from current facet fluxes
    double64  FluxBalance( Node<dim>* const ) const;

    /// returns in- and outflow from FV computed from the current facet fluxes
    std::pair<double64,double64>  InflowOutflow( const Node<dim>* const ) const;
  
//TODO: include dispersion and options for diffusive fluxes

  private:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
  private:
    mutable VectorVariable<dim>  nrml_, gradC_; ///< for efficient re-use
    mutable Point<dim>           pt_;
    mutable DenseMatrix<DM_MIN>  DN_;
};



} // end csmp

#endif /* FACET_FLUX_TRACER_TRANSFER_EXPLICIT_H */
