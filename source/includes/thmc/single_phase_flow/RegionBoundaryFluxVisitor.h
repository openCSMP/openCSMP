//
//  RegionBoundaryFluxVisitor.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/30/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_API_library2014__RegionBoundaryFluxVisitor_h
#define CSMP_API_library2014__RegionBoundaryFluxVisitor_h

#include "Visitor.h"
#include "VectorVariable.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Element;
template<size_t,template<size_t> class> class ModelSubDomain;

/// using the FVM this visitor computes the cumulative single-phase influx into the target region
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class RegionBoundaryFluxVisitor : public Visitor<dim> {
  public:
    /// choice between Region, Boundary and SplitBoundary is desired
    RegionBoundaryFluxVisitor( Model<dim>& model,
                               const char*  target_region,
                               const char*  Darcy_velocity );
  
    virtual ~RegionBoundaryFluxVisitor() {};

    /// application target = node-centered finite volume
    virtual void Visit( Node<dim>* );
  
    /// returns the flux through the region boundary 
    double64 InFlux() const;
    void ResetFlux();
    void TimeIncrement( double64 );
 
  protected:
    ModelSubDomain<dim,Element>*  domain_ptr_; // either of:
    csmp::Index                   flux_key_;   // for the user-defined flux variable
    double64                      FVinflux_,   // flux summation variable
                                  delta_t_;    // time increment used in the flux integration
    VectorVariable<dim>           vt_;
};

} // end csmp

#endif /* defined(__CSMP_API_library2014__RegionBoundaryFluxVisitor__) */
