#ifndef CSMP_EXPLICIT_TRANSPORT_H
#define CSMP_EXPLICIT_TRANSPORT_H

#include "Variables_TracerTransfer.h"
#include "FacetFlux_TracerTransferExplicit.h"
#include "TimeStepEvaluator.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;
template<size_t> class TwoPhaseModel;
template<size_t> class ExplicitTransport;

// TODO: gradient calculation: compare different implementations (ExtrapolateElementPropertyToNode computing and averaging element gradients vs. node by node approach)
// TODO: make the transported variable a template as well: Scalar, Array, FlaggedArray...
template<size_t dim>
  class ExplicitTransport : public variables::Variables_TracerTransfer,
                          public FacetFlux_TracerTransferExplicit<dim,ExplicitTransport>,
                          public TimeStepEvaluator<dim,ExplicitTransport> {
  public:
    // TODO: add choice of transport scheme: 1st versus 2nd order in space
    /// constructor for target region; by default all driving forces are considered
    ExplicitTransport( Model<dim>&, const char* target_region, bool second_order_in_space );
  
    /// computes the time constraint
    double64 TimeIncrement();
    
    /// executes incremental time-stepping (a suitable time increment is computed by scheme)
    void AdvectVariable( double64 time_interval );

    /// gets the model
    const Model<dim>& GetModel() const;
  
  private:    
    /// 1.b computations of piecewise constatn element velocities and facet fluxes using FacetFlux (facet flux) policy
    void PostProcessVelocityAndUpdateFacetFluxes();
    
    /// 2. calculates optimal time increment, flux balance, and in- and out flows for each FV
    double64 TimeIncrementAndFluxBalance( double64 max_time_increment );
    
    /// 3. inflow and outflow flux compensation
    void AccumulateBoundaryConditions();
    
    /// 4. composes 'new concentration': C^t+1 = C^t - dt/(phi Vi) * sum_j^faces Aj n . [vi]
    void AssembleSolution( double64 delta_t, bool enforce_divergence_free_vt_field );

    /// if we know beforehand that velocity field will be divergence free, this method compensates for small abberations from this
    void AdjustResultsAssumingDivergenceFreeVelocityField( double64 time_interval );
    
    /// 5. transfer results updating concentration, zeroing out 'new concentration' values, and performing range checks; returns error
    double64 VerifyAndAssignResults( bool show_range, bool do_range_check ) const;
  
  private:
    Model<dim>& model_;
    Region<dim>& gref_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    bool second_order_in_space_;
};


/**
\brief     For implicit transport calculations
\details   Part of the Colleoli transport scheme.
\author    Stephan K. Matthai
\version   1a
\date      21/2/2013
\pre       high-level class depending on CSMP++ API
\bug
\warning
\copyright Stephan K. Matthai

@section motivation Motivation


@section design Design Intent


@section applicability Applicability


@section collaborations Collaborations


@section implementation Implementation


@section examples Application Examples

@code

@endcode

*/

} // end csmp


#endif /* CSMP_EXPLICIT_TRANSPORT_H */
