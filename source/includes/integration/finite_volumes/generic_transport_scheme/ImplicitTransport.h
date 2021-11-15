#ifndef CSMP_IMPLICIT_TRANSPORT_H
#define CSMP_IMPLICIT_TRANSPORT_H

#include "IntegralEquation.h"
#include "LinearAlgebraicSystem.h"
#include "Accumulator.h"
#include "Integrator.h"
#include "PostProcessor.h"
#include "FluxEvaluator.h"
#include "TimeStepEvaluator.h"

#include "VariableSet_TracerTransfer.h"

namespace csmp {

template<size_t> class Region;
template<size_t> class Model;

// alias template (since C++11)
template<size_t dim> 
using GoverningEquation = IntegralEquation<dim,variables::VariableSet_TracerTransfer>;

/**
TODO: include thickness attributes to support computations with lower-dimensional elements
TODO: make this algorithm general so that it can also be applied to Boundary and SplitBoundary objects
TODO: Implement bijective grid-to-grid mapping to loose grid orientation effects completely!
*/
template<size_t dim>
class ImplicitTransport : public GoverningEquation<dim>,
                          public Accumulator<dim,ImplicitTransport>,        // needs base class to get to the INDEX keys
                          public Integrator<dim,ImplicitTransport>,
                          public PostProcessor<dim,ImplicitTransport>,      // anything we may need
                          public FluxEvaluator<dim,ImplicitTransport>,
                          public TimeStepEvaluator<dim,ImplicitTransport> {
  public:
    /// constructor for target region; by default all driving forces are considered; default equation is used; model non-const because FVs may have to be created
    ImplicitTransport( Model<dim>&, const std::string& target_region, 
                       const GoverningEquation<dim>&, 
                       bool second_order_in_space );
 
    // TODO: review whether this functio should be replaced by: VolumetricFlowAndTransportVariableFluxBalances()
    double TimeIncrementAndFluxBalance( double max_time_increment ) const; 
    
    /// computes the time constraint
    double TimeIncrement() const;
    
    /// initialises transport class for current pressure/velocity/transport variable field
    void UpdateFluxesAndFluxBalances();

    /// executes incremental time-stepping (a suitable time increment is computed by scheme) and writes the results back to model
    void AdvectVariable( double time_interval );

    double IncomingVolumetricFlow() const;
    double OutgoingVolumetricFlow() const;


    // MEMBER ACCESS
    
    Region<dim>&                  ComputationDomain()       { return subdomain_; }
    const Region<dim>&            ComputationDomain() const { return subdomain_; }
    LinearAlgebraicSystem&        LinearSystem()            { return linalg_sys_; }
    const LinearAlgebraicSystem&  LinearSystem() const      { return linalg_sys_; }
    
    /// true if the computational domain has parts that are removed from the model boundary
    bool HasHaloStencils() const { return !halo_cells_.empty(); }
    
    typename std::vector<Element<dim>*>::const_iterator HaloCellsBegin() const { return halo_cells_.begin(); }
    typename std::vector<Element<dim>*>::const_iterator HaloCellsEnd() const { return halo_cells_.end(); }

    void Verbose( bool value ) { verbose_ = value; }
    bool Verbose() const { return verbose_; }
    
  private:
    // note that the sequence of these classes is critical for correct construction of ImplicitTransport
    Region<dim>&                subdomain_;   ///< computational region
    LinearAlgebraicSystem       linalg_sys_;  ///<  A x = b
    std::vector<Element<dim>*>  halo_cells_;  ///< elements outside of subdomain, contributing stencils to subdomain FVs
    bool second_order_ = false, verbose_ = true;
    
    // Solver is in the Integrator base class
};

} // end csmp


#endif /* CSMP_IMPLICIT_TRANSPORT_H */
