#include "Element.h"
#include "VariablePlacement.h"

namespace csmp {

template<size_t dim>
FiniteElementPlacement<dim,ELEMENT>
FiniteElementDiscretisationPolicy<dim>::AtBarycenter()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return FiniteElementPlacement<dim,ELEMENT>(e, 0, 0);
}


template<size_t dim>
FiniteElementPlacementCollection<dim,NODE>
FiniteElementDiscretisationPolicy<dim>::AllNodes()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return FiniteElementPlacementCollection<dim,NODE>(e);
}


template<size_t dim>
FiniteElementPlacementCollection<dim,ELEMENT_INTEGRATION_POINT>
FiniteElementDiscretisationPolicy<dim>::AllElementIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return FiniteElementPlacementCollection<dim,ELEMENT_INTEGRATION_POINT>(e);
}


template<size_t dim>
FiniteElementPlacementCollection<dim,FACET_INTEGRATION_POINT>
FiniteElementDiscretisationPolicy<dim>::AllFacetIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return FiniteElementPlacementCollection<dim,FACET_INTEGRATION_POINT>(e);
}


template<size_t dim>
FiniteElementPlacementCollection<dim,SECTOR_INTEGRATION_POINT>
FiniteElementDiscretisationPolicy<dim>::AllSectorIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return FiniteElementPlacementCollection<dim,SECTOR_INTEGRATION_POINT>(e);
}


template class FiniteElementDiscretisationPolicy<1u>;
template class FiniteElementDiscretisationPolicy<2u>;
template class FiniteElementDiscretisationPolicy<3u>;

}
