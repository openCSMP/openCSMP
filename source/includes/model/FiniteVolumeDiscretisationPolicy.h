#ifndef CSMP_FINITE_VOLUME_DISCRETISATION_POLICY_H
#define CSMP_FINITE_VOLUME_DISCRETISATION_POLICY_H

#include "VariablePlacement.h"

namespace csmp {

template<size_t> class Node;

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
class FiniteVolumeDiscretisationPolicy
{
public:
    FiniteVolumePlacement<dim,NODE> AtNode();

    NeighbourNodeCollection<dim> AllNeighbourNodes();

    FiniteVolumePlacementCollection<dim,ELEMENT> AllElements();
    FiniteVolumePlacementCollection<dim,FACET_INTEGRATION_POINT> AllFacetIntegrationPoints();
    FiniteVolumePlacementCollection<dim,SECTOR_INTEGRATION_POINT> AllSectorIntegrationPoints();
};

} // csmp

#endif



