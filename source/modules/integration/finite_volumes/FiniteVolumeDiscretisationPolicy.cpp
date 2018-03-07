#include "Node.h"
#include "VariablePlacement.h"


namespace csmp {

template<size_t dim>
FiniteVolumePlacement<dim,NODE>
FiniteVolumeDiscretisationPolicy<dim>::AtNode()
{
    Node<dim>& n = *static_cast<Node<dim>*>(this);
    return FiniteVolumePlacement<dim,NODE>(n, 0, 0, 0);
}

template<size_t dim>
FiniteVolumePlacementCollection<dim,ELEMENT>
FiniteVolumeDiscretisationPolicy<dim>::AllElements()
{
    Node<dim>& n = *static_cast<Node<dim>*>(this);
    return FiniteVolumePlacementCollection<dim,ELEMENT>(n);
}

template<size_t dim>
FiniteVolumePlacementCollection<dim,FACET_INTEGRATION_POINT>
FiniteVolumeDiscretisationPolicy<dim>::AllFacetIntegrationPoints()
{
    Node<dim>& n = *static_cast<Node<dim>*>(this);
    return FiniteVolumePlacementCollection<dim,FACET_INTEGRATION_POINT>(n);
}

template<size_t dim>
FiniteVolumePlacementCollection<dim,SECTOR_INTEGRATION_POINT>
FiniteVolumeDiscretisationPolicy<dim>::AllSectorIntegrationPoints()
{
    Node<dim>& n = *static_cast<Node<dim>*>(this);
    return FiniteVolumePlacementCollection<dim,SECTOR_INTEGRATION_POINT>(n);
}


template class FiniteVolumeDiscretisationPolicy<1u>;
template class FiniteVolumeDiscretisationPolicy<2u>;
template class FiniteVolumeDiscretisationPolicy<3u>;

}
