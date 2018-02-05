#include "Element.h"
#include "VariablePlacement.h"

namespace csmp {

template<size_t dim>
ElementPlacement<dim,ELEMENT>
ElementDiscretisationPolicy<dim>::AtBarycenter()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return ElementPlacement<dim,ELEMENT>(e, 0, 0);
}


template<size_t dim>
ElementPlacementCollection<dim,NODE>
ElementDiscretisationPolicy<dim>::AllNodes()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return ElementPlacementCollection<dim,NODE>(e);
}


template<size_t dim>
ElementPlacementCollection<dim,ELEMENT_INTEGRATION_POINT>
ElementDiscretisationPolicy<dim>::AllElementIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return ElementPlacementCollection<dim,ELEMENT_INTEGRATION_POINT>(e);
}


template<size_t dim>
ElementPlacementCollection<dim,FACET_INTEGRATION_POINT>
ElementDiscretisationPolicy<dim>::AllFacetIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return ElementPlacementCollection<dim,FACET_INTEGRATION_POINT>(e);
}


template<size_t dim>
ElementPlacementCollection<dim,SECTOR_INTEGRATION_POINT>
ElementDiscretisationPolicy<dim>::AllSectorIntegrationPoints()
{
    Element<dim>& e = *static_cast<Element<dim>*>(this);
    return ElementPlacementCollection<dim,SECTOR_INTEGRATION_POINT>(e);
}


template class ElementDiscretisationPolicy<1u>;
template class ElementDiscretisationPolicy<2u>;
template class ElementDiscretisationPolicy<3u>;

}
