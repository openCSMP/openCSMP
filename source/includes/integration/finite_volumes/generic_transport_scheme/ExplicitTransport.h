#ifndef CSMP_EXPLICIT_TRANSPORT_H
#define CSMP_EXPLICIT_TRANSPORT_H

#include "VariableSet_TracerTransfer.h"
#include "FluxEvaluator.h"
#include "TimeStepEvaluator.h"

namespace csmp {

template<size_t> class Region;
template<size_t> class Model;

// TODO: gradient calculation: compare different implementations (ExtrapolateElementPropertyToNode...)
// TODO: make the transported variable a template as well: Scalar, Array, FlaggedArray...
/**
    Single-phase flow tracer, heat or other transport, taking into account the thickness attribute
    of lower-dimensional elements in the flux calculations.
 
    @attention terminology is super important:
    - flux refers to the total quantity that passes through the surface, i.e.  (v . n_i) * A_i,
      where subscript i refers to the finite volume facet
    - volume flux  - refers to volumetric flow
    - flux - refers to products of the advected quantity and the volumetric flow
 
    flux balance - sum over all finite-volume facet fluxes
 
    @author SKM
    @date 25/6/2015
 
    @copyright Stephan K. Matthai
 
    TODO: account for the effects of potential element based 'fluid volume source' terms.
    TODO: return rates of solute influx/outflux into/from the model
    double InFlux() const;
    double OutFlux() const;
    TODO: add choice of transport scheme: 1st versus 2nd order in space
    TODO: Generalise scheme so that it can handle Face and InterFace objects
*/
template<size_t dim>
class ExplicitTransport : public variables::VariableSet_TracerTransfer,
                          public FluxEvaluator<dim,ExplicitTransport>,
                          public TimeStepEvaluator<dim,ExplicitTransport> {
  public:
    friend class ExplicitTransport_Test;
    
    /// shorthand to get to variable names
    const VariableSet_TracerTransfer&  Notation;

    /// constructor for target region; by default all driving forces are considered; @attention if there is a velocity field, it is used to initialise facet fluxes
    ExplicitTransport( Model<dim>&, const std::string& target_region );
    
    /// initialises transport class for current pressure/velocity/transport variable field
    void UpdateFluxesAndFluxBalances();
  
    /// computes the time constraint
    double TimeIncrement() const;
    
    /// executes incremental time-stepping (a suitable time increment is computed by scheme)
    void AdvectVariable( double time_interval );
    
    double IncomingVolumetricFlow() const;
    double OutgoingVolumetricFlow() const;

  private:
    /// identifies "halo elements", i.e. which contribute to domain FVs, but are outside of domain, returns number
    size_t CollectHaloStencils();
    
    /// true if the computational domain has parts that are removed from the model boundary
    bool HasHaloStencils() const { return !halo_elmts_.empty(); }
    
    /// 1.a computations of  volumetric flows and (chemical) fluxes across facets using FacetFlux (facet flux) policy
    void VolumetricFlowAndTransportVariableFluxBalances();
    
    /// 1.b recomputed transport-variable flux balance, stored in 'accumulation'; assumes that transport velocity did not change
    void TransportVariableFluxBalances();
    
    /// 2. calculates optimal time increment, flux balance, and in- and out flows for each FV
    double TimeIncrement_CFL_Outflow( double max_time_increment ) const;
    
    /// 4. composes 'new concentration': C^t+1 = C^t - dt/(phi Vi) * sum_j^faces Aj n . [vi]
    void Assemble1stOrderSolution( double delta_t, bool enforce_divergence_free_vt_field );
    
    /// 5. transfer results updating concentration, zeroing out 'new concentration' values, and performing range checks; returns error
    double VerifyAndAssignResults( bool show_range, bool do_range_check ) const;
    
  private:
    // DEPRACATE ?
    /// if we know beforehand that velocity field will be divergence free, this method compensates for small abberations from this
    void AdjustResultsAssumingDivergenceFreeVelocityField( double time_interval );
    
  private:
    Region<dim>&                subdomain_;   ///< region to which transport algorithm is applied
    std::vector<Element<dim>*>  halo_elmts_;  ///< elements outside of subdomain, contributing stencils to subdomain FVs
    double  upper_limit_, lower_limit_;     ///< range in which the result is allowed to vary
};


/**
@class ExplicitTransport ExplicitTransport "reservoir_simulator/TransportEquationSolver.h"

\brief     For implicit transport calculations
\details   Part of the Colleoli transport scheme.
\author    Stephan K. Matthai
\version   1a
\date      21/2/2013
\pre       high-level class depending on CSMP++ API
\bug
\warning
\copyright Stephan K. Matthai

*/

} // end csmp


#endif /* CSMP_EXPLICIT_TRANSPORT_H */
