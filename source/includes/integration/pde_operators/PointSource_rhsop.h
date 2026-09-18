// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POINT_SOURCE_RHSOP_H
#define POINT_SOURCE_RHSOP_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/** @brief Assignment of scalar nodal (Neumann) source or sink terms.

    This operator applies (absolute) nodal source terms to the RHS vector as is. 
    These values must be precomputed elsewhere. One possibility to do this
    is via the method

    @code
      NodeCenteredFiniteVolumeTranport::TransformScalarBoundaryValuesIntoNeumannConditions()
    @endcode

    It computes gradients normal to the outer model boundary from influxes (sources)
    assigned and flagged as NEUMANN elsewhere. To do this it needs material parameters
    from the corresponding flow law, for instance the hydraulic conductivity
    in the case of Darcy's law.
 
    @attention the nodal term gets applied as many times as the node is a member of an element, i.e.
    n.Parents() times. Therefore the value should be divided by the number of parent elements of the node.
    This is already accounted for by the method GetOperands().
 
    @attention when applying this term, remember that model must have 
    at least one Dirichlet constraint assigned to make solution unique.

@date 2000 
@author Stephan K. Matthaei

*/      
template<uint32_t dim, template<uint32_t> class CELL=Element>
class PointSource_rhsop : public MathOperatorRHS<dim,CELL> {
  public:
    PointSource_rhsop( const PropertyDatabase<dim>&, const char* nodal_src, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    PointSource_rhsop<dim,CELL>* clone() const override { return new PointSource_rhsop<dim,CELL> (*this); }

private:
    std::vector<ScalarVariable>  SRC_;
};

} // csmp

#endif
















