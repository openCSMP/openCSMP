#include "FiniteVolumePlacement.h"
#include "VariablePlacement.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"


namespace csmp {

    template<size_t dim>
    Element<dim>&
    FiniteVolumePlacementOperations<dim,FACET_INTEGRATION_POINT>::TheElement()
    {
      auto user = User();
      return *user->n_.Parent(user->idx1_);
    }

    template Element<1u>& FiniteVolumePlacementOperations<1u,FACET_INTEGRATION_POINT>::TheElement();
    template Element<2u>& FiniteVolumePlacementOperations<2u,FACET_INTEGRATION_POINT>::TheElement();
    template Element<3u>& FiniteVolumePlacementOperations<3u,FACET_INTEGRATION_POINT>::TheElement();
    

#define INSTANTIATE_FV_PROPERTY_INTERPOLATOR(interp) \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<1,SCALAR>(Index const&, Node<1> const*, size_t, size_t, size_t, VariableTypeTraits<1,SCALAR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<2,SCALAR>(Index const&, Node<2> const*, size_t, size_t, size_t, VariableTypeTraits<2,SCALAR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<3,SCALAR>(Index const&, Node<3> const*, size_t, size_t, size_t, VariableTypeTraits<3,SCALAR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<1,VECTOR>(Index const&, Node<1> const*, size_t, size_t, size_t, VariableTypeTraits<1,VECTOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<2,VECTOR>(Index const&, Node<2> const*, size_t, size_t, size_t, VariableTypeTraits<2,VECTOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<3,VECTOR>(Index const&, Node<3> const*, size_t, size_t, size_t, VariableTypeTraits<3,VECTOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<1,TENSOR>(Index const&, Node<1> const*, size_t, size_t, size_t, VariableTypeTraits<1,TENSOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<2,TENSOR>(Index const&, Node<2> const*, size_t, size_t, size_t, VariableTypeTraits<2,TENSOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<3,TENSOR>(Index const&, Node<3> const*, size_t, size_t, size_t, VariableTypeTraits<3,TENSOR>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<1,ARRAY>(Index const&, Node<1> const*, size_t, size_t, size_t, VariableTypeTraits<1,ARRAY>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<2,ARRAY>(Index const&, Node<2> const*, size_t, size_t, size_t, VariableTypeTraits<2,ARRAY>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<3,ARRAY>(Index const&, Node<3> const*, size_t, size_t, size_t, VariableTypeTraits<3,ARRAY>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<1,FLAGGEDARRAY>(Index const&, Node<1> const*, size_t, size_t, size_t, VariableTypeTraits<1,FLAGGEDARRAY>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<2,FLAGGEDARRAY>(Index const&, Node<2> const*, size_t, size_t, size_t, VariableTypeTraits<2,FLAGGEDARRAY>::VariableType&); \
template void FiniteVolumePropertyInterpolator<interp>::Interpolate<3,FLAGGEDARRAY>(Index const&, Node<3> const*, size_t, size_t, size_t, VariableTypeTraits<3,FLAGGEDARRAY>::VariableType&);



  template<>
  FiniteVolumePropertyInterpolator<FV_READ_NODE>::FiniteVolumePropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  FiniteVolumePropertyInterpolator<FV_READ_NODE>::Interpolate( const Index& prop, const Node<dim>* n, size_t idx1, size_t idx2, size_t idx3, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    n->Read(prop, var);
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_READ_NODE)

  template<>
  FiniteVolumePropertyInterpolator<FV_READ_ELMT>::FiniteVolumePropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  FiniteVolumePropertyInterpolator<FV_READ_ELMT>::Interpolate( const Index& prop, const Node<dim>* n, size_t idx1, size_t idx2, size_t idx3, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    n->Parent(idx1)->Read(prop, var);
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_READ_ELMT)
  
  template<>
  FiniteVolumePropertyInterpolator<FV_READ_FIP>::FiniteVolumePropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  FiniteVolumePropertyInterpolator<FV_READ_FIP>::Interpolate( const Index& prop, const Node<dim>* n, size_t idx1, size_t idx2, size_t idx3, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    auto eptr = n->Parent(idx1);
    size_t pnid = n->ParentNodeNumber(idx1);
    auto fv = eptr->FV();
    size_t facet = fv->FacetSurroundingSector(pnid, idx2);
    
    eptr->Read(facet, 0U, prop, var);
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_READ_FIP)  

  template<>
  FiniteVolumePropertyInterpolator<FV_NODE_TO_FIP>::FiniteVolumePropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  FiniteVolumePropertyInterpolator<FV_NODE_TO_FIP>::Interpolate( const Index& prop, const Node<dim>* n, size_t idx1, size_t idx2, size_t idx3, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    auto eptr = n->Parent(idx1);
    size_t pnid = n->ParentNodeNumber(idx1);
    auto fv = eptr->FV();
    size_t facet = fv->FacetSurroundingSector(pnid, idx2);

    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    calculateN(*eptr, fv->FacetIntegrationPoint( facet, idx3 ), &coeff[0]);
    
    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_NODE_TO_FIP)

  template<>
  FiniteVolumePropertyInterpolator<FV_NODE_TO_SIP>::FiniteVolumePropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  FiniteVolumePropertyInterpolator<FV_NODE_TO_SIP>::Interpolate( const Index& prop, const Node<dim>* n, size_t idx1, size_t idx2, size_t idx3, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    auto eptr = n->Parent(idx1);
    size_t pnid = n->ParentNodeNumber(idx1);
    auto fv = eptr->FV();

    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    calculateN(*eptr, fv->SectorIntegrationPoint( pnid, idx3 ), &coeff[0]);
    
    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_NODE_TO_SIP)


  template<size_t dim>
  Point<dim>
  FiniteVolumePlacementOperations<dim,FACET_INTEGRATION_POINT>::DirectedArea() const
  {
    auto user = User();
    auto& e = *user->n_.Parent(user->idx1_);
    size_t pnid = user->n_.ParentNodeNumber(user->idx1_);
    auto fv = e.FV();
    size_t iFacet = fv->FacetSurroundingSector(pnid, user->idx2_);
    return directedAreaOfFacet(e, iFacet);
  }

  template Point<1u> FiniteVolumePlacementOperations<1u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<2u> FiniteVolumePlacementOperations<2u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<3u> FiniteVolumePlacementOperations<3u,FACET_INTEGRATION_POINT>::DirectedArea() const;

  template<size_t dim>
  Point<dim>
  FiniteVolumePlacementOperations<dim,FACET_INTEGRATION_POINT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop) const
  {
    auto user = User();
    auto& e = *user->n_.Parent(user->idx1_);
    
    auto fv = e.FV();
    auto fe = e.FE();

    const size_t element_dim = (size_t)fv->Geometry();
    const size_t num_nodes = e.Nodes();

    e.CoordinateMatrix();

    std::vector<double64> DN[dim];
    for (size_t i = 0; i < element_dim; ++i) {
      DN[i].resize(num_nodes);
    }
    calculateDN(e, fv->Barycenter(), DN);

    Point<dim> grad(0.);
    for ( size_t i=0U; i<num_nodes; ++i ) {
      const double64 value_at_node(e.N(i)->Read(prop));
      for (size_t j = 0; j < element_dim; ++j) {
        grad[j] += DN[j][i] * value_at_node;
      }
    }
    fe->JacobianInverse();

    return Point<dim>(fe->JINV * grad.Coordinates());
  }
  
  template Point<1u> FiniteVolumePlacementOperations<1u,FACET_INTEGRATION_POINT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop) const;
  template Point<2u> FiniteVolumePlacementOperations<2u,FACET_INTEGRATION_POINT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop) const;
  template Point<3u> FiniteVolumePlacementOperations<3u,FACET_INTEGRATION_POINT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop) const;   
  
}

