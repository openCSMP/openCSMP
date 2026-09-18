// Copyright (c) 2014 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_API_library2014__RegionBoundaryFluxVisitor_h
#define CSMP_API_library2014__RegionBoundaryFluxVisitor_h

#include "Visitor.h"
#include "VectorVariable.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Element;
template<uint32_t,template<uint32_t> class> class ModelSubDomain;

/// using the FVM this visitor computes the cumulative single-phase influx into the target region
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
class RegionBoundaryFluxVisitor final : public Visitor<dim> {
  public:
    /// choice between Region, Boundary and SplitBoundary is desired
    RegionBoundaryFluxVisitor( Model<dim>&,
                               const char*  target_region,
                               const char*  Darcy_velocity );

    /// application target = node-centered finite volume
    void Visit( Node<dim>* ) override final;
  
    /// returns the flux through the region boundary 
    double InFlux() const;
    void ResetFlux();
    void TimeIncrement( double );
 
  protected:
    ModelSubDomain<dim,Element>*  domain_ptr_; ///< either of:
    csmp::Index                   flux_key_;   ///< for the user-defined flux variable
    double                        FVinflux_,   ///< flux summation variable
                                  delta_t_;    ///< time increment used in the flux integration
    VectorVariable<dim>           vt_;
};

} // end csmp

#endif /* defined(__CSMP_API_library2014__RegionBoundaryFluxVisitor__) */
