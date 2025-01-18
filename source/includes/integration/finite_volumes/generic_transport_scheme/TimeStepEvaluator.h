#ifndef TIME_STEP_EVALUATOR_H
#define TIME_STEP_EVALUATOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;

/// for application to individual FVs
template<uint32_t dim, template<uint32_t> class USER>
class TimeStepEvaluator {
  public:
    explicit TimeStepEvaluator( double step_size_reduction_factor=0.5, double max_time_increment=86400. * 365. );

    /// returns robust delta_t criterion in the presence of fluid sources and sinks; also computes and stores flux balance
    double OutFlowLessThanContentIncrement( Node<dim>* const, double outflow ) const;
  
    /// returns robust delta_t criterion in the presence of fluid sources and sinks at model boundary (no variables are touched)
    double OutFlowLessThanContentIncrementBoundary( const Node<dim>* const, double outflow ) const;

    /// anisotropic CFL for the strictly hyperbolic case; also computes and stores flux balance
    double StreamlineCFL( const Element<dim>* const ) const;
  
    /// adjust setting for the solve (default=0.5)
    void StepSizeReductionFactor( double=0.5 );
    double StepSizeReductionFactor() const;
  
    /// return maximum time-increment that is applied if there is ~no flow
    double MaxTimeIncrement() const;
    
  private:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

  private:
    const double max_time_increment_;    ///< set through constructor argument
    double step_size_reduction_factor_;  ///< CFL multiplier set through constructor argument
};

/**
@class TimeStepEvaluator TimeStepEvaluator  "reservoir_simulator/TimeStepEvaluator.h"

\brief     Face-by-face calculation of CFL like timestepping criteria.
\details   Part of Colleoli transport scheme.
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


#endif /* TIME_STEP_EVALUATOR_H */
