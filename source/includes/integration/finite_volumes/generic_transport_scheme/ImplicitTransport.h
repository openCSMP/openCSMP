#ifndef CSMP_IMPLICIT_TRANSPORT_H
#define CSMP_IMPLICIT_TRANSPORT_H

#include "VariableSet_TracerTransferImplicit.h"
#include "FacetFlux_TracerTransferExplicit.h"
#include "TimeStepEvaluator.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "LinearSolver.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Region;
template<size_t> class Model;
template<size_t> class TwoPhaseModel;

template<size_t dim>
class ImplicitTransport : public VariableSet_TracerTransferImplicit,
                          public FacetFlux_TracerTransferExplicit<dim,ImplicitTransport>,
                          public TimeStepEvaluator<dim,ImplicitTransport> {
  public:
    // TODO: add choice of transport scheme: 1st versus 2nd order in space
    /// constructor for target region; by default all driving forces are considered
    ImplicitTransport( Model<dim>&, const char* target_region, bool second_order );
  
    /// computes the time constraint
    double64 TimeIncrement();
    
    /// executes incremental time-stepping (a suitable time increment is computed by scheme)
    void AdvectVariable( double64 time_interval );

    /// gets the model
    const Model<dim>& GetModel() const;
  
  private:
    /// 1.a computations of facet flux using FacetFlux (facet flux) policy
    void UpdateFacetFluxes(bool reuse_previous_velocity);
    void UpdateFacetFluxes_O1(bool reuse_previous_velocity);
    void UpdateFacetFluxes_O2(bool reuse_previous_velocity);

    /// 1.b computations of piecewise constatn element velocities and facet fluxes using FacetFlux (facet flux) policy
    void PostProcessVelocityAndUpdateFacetFluxes();

    /// 2. calculates optimal time increment, flux balance, and in- and out flows for each FV
    double64 TimeIncrementAndFluxBalance( double64 max_time_increment );

    /// 3. process 
    void AccumulateSystem();

    /// 4. inflow and outflow flux compensation
    void AccumulateBoundaryConditions();

    /// 5. composes 'new concentration': C^t+1 = C^t - dt/(phi Vi) * sum_j^faces Aj n . [vi]
    void Solve();

    /// 6. transfer results updating concentration, zeroing out 'new concentration' values, and performing range checks; returns error
    double64 VerifyAndAssignResults( bool show_range, bool do_range_check ) const;

  private:
    Model<dim>& model_;
    Region<dim>& gref_;
    double64 upper_limit_, lower_limit_; ///< range in which the result is allowed to vary
    bool second_order_;
};


/**
@class ImplicitTransport ImplicitTransport "integration/finite_volumes/generic_transport_scheme/ImplicitTransport.h"

\brief     Implicit transport calculations
\details   Part of the Colleoli transport scheme.
\author    Andrew J. Bromage
\version   0a
\date      28/11/2017
\pre       high-level class depending on CSMP++ API
\bug
\warning
\copyright The University of Melbourne

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


#endif /* CSMP_IMPLICIT_TRANSPORT_H */
