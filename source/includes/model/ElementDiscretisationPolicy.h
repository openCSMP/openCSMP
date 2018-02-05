#ifndef CSMP_ELEMENT_DISCRETISATION_POLICY_H
#define CSMP_ELEMENT_DISCRETISATION_POLICY_H

#include "VariablePlacement.h"

namespace csmp {

template<size_t> class Element;

/**

@brief Policy class to support calculations involving variables with different placements,
and iteration over integration points, etc.

@author A.J. Bromage
@date 2018

@section motivation Motivation

@section design Design Intent

@section Applicability

@section structure Structure

@section collaborations Collaborations
 
@section consequences Consequences

@section implementation Implementation

@section application Application Examples

*/

template<size_t dim>
class ElementDiscretisationPolicy
{
public:
    ElementPlacement<dim,ELEMENT> AtBarycenter();

    ElementPlacementCollection<dim,NODE> AllNodes();
    ElementPlacementCollection<dim,ELEMENT_INTEGRATION_POINT> AllElementIntegrationPoints();
    ElementPlacementCollection<dim,FACET_INTEGRATION_POINT> AllFacetIntegrationPoints();
    ElementPlacementCollection<dim,SECTOR_INTEGRATION_POINT> AllSectorIntegrationPoints();
};

} // csmp

#endif



