#ifndef TIME_STEP_EVALUATOR_H
#define TIME_STEP_EVALUATOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Node;

/**

@brief Local time-increment calculator for application to individual FVs.

Face-by-face calculation of CFL like timestepping criteria.

\details   Part of Colleoli transport scheme.
\author    Stephan K. Matthai
\version   1a
\date      21/2/2013
\pre       high-level class depending on CSMP++ API
\copyright Stephan K. Matthai

*/
template<size_t dim, template<size_t> class USER>
class TimeStepEvaluator {
  public:
    explicit TimeStepEvaluator( double64 step_size_reduction_factor=0.5, double64 max_time_increment=86400. * 365. );

    /// returns robust delta_t criterion in the presence of fluid sources and sinks; also computes and stores flux balance
    double64 OutFlowLessThanContentIncrement( Node<dim>* const ) const;
  
    /// returns robust delta_t criterion in the presence of fluid sources and sinks at model boundary (no variables are touched)
    double64 OutFlowLessThanContentIncrementBoundary( const Node<dim>* const ) const;

    /// anisotropic CFL for the strictly hyperbolic case; also computes and stores flux balance
    double64 StreamlineCFL( Node<dim>* const ) const;
  
    /// adjust setting for the solve (default=0.5)
    void StepSizeReductionFactor( double64=0.5 );
    double64 StepSizeReductionFactor() const;
  
    /// return maximum time-increment that is applied if there is ~no flow
    double64 MaxTimeIncrement() const;
    
  private:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

  private:
    const double64 max_time_increment_;
    double64       step_size_reduction_factor_;
};

} // end csmp


#endif /* TIME_STEP_EVALUATOR_H */
