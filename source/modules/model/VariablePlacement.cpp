#include "VariablePlacement.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"

namespace {
  using namespace csmp;
  
  template<size_t dim>
  size_t elementDimension(const Element<dim>& e)
  {
    if (e.IsLineElement()) {
      return 1;
    }
    else if (e.IsSurfaceElement()) {
      return 2;
    }
    else if (e.IsVolumeElement()) {
      return 3;
    }
    throw Exception(ERROR, "elementDimension", "element type unknown");
  }

  template<size_t dim>
  void CalculateN(const Element<dim>* eptr, size_t element_dim, const Point<dim>& p, double64* coeff);

  template<>
  void CalculateN<1u>(const Element<1u>* eptr, size_t element_dim, const Point<1u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 1:
        fe->Nr( p[0], coeff );
    }
  }

  template<>
  void CalculateN<2u>(const Element<2u>* eptr, size_t element_dim, const Point<2u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 1:
        fe->Nr( p[0], coeff );
        break;

      case 2:
        fe->Nrs( p[0], p[1], coeff );
        break;
    }
  }

  template<>
  void CalculateN<3u>(const Element<3u>* eptr, size_t element_dim, const Point<3u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 3:
        fe->Nrst( p[0], p[1], p[2], coeff );
        break;
      case 2:
        fe->Nrs( p[0], p[1], coeff );
        break;
      case 1:
        fe->Nr( p[0], coeff );
        break;
    }
  }

  template<size_t dim>
  void
  CalculateDN(const Element<dim>& e, const Point<dim>& p, std::vector<double64>* DN);

  template<>
  void
  CalculateDN(const Element<1u>& e, const Point<1u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (elementDimension(e)) {
      case 1:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
    }
  }


  template<>
  void
  CalculateDN(const Element<2u>& e, const Point<2u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (elementDimension(e)) {
      case 1:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
        break;

      case 2:
        fe->dNr(p[0], p[1], DN[0]);
        fe->dNs(p[0], p[1], DN[1]);
        fe->Jacobian( DN[0], DN[1] );
        break;
    }
  }


  template<>
  void
  CalculateDN(const Element<3u>& e, const Point<3u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (elementDimension(e)) {
      case 3:
        fe->dNr(p[0], p[1], p[2], DN[0]);
        fe->dNs(p[0], p[1], p[2], DN[1]);
        fe->dNt(p[0], p[1], p[2], DN[2]);
        fe->Jacobian( DN[0], DN[1], DN[2] );
        break;
      case 2:
        fe->dNs(p[0], p[1], DN[1]);
        fe->Jacobian( DN[0], DN[1] );
        break;
      case 1:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
        break;
    }
  }
}



namespace csmp {

#define INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(interp) \
template void ElementPropertyInterpolator<interp>::Interpolate<1,SCALAR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,SCALAR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<2,SCALAR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,SCALAR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<3,SCALAR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,SCALAR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<1,VECTOR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,VECTOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<2,VECTOR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,VECTOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<3,VECTOR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,VECTOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<1,TENSOR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,TENSOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<2,TENSOR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,TENSOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<3,TENSOR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,TENSOR>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<1,ARRAY>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,ARRAY>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<2,ARRAY>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,ARRAY>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<3,ARRAY>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,ARRAY>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<1,FLAGGEDARRAY>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,FLAGGEDARRAY>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<2,FLAGGEDARRAY>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,FLAGGEDARRAY>::VariableType&); \
template void ElementPropertyInterpolator<interp>::Interpolate<3,FLAGGEDARRAY>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,FLAGGEDARRAY>::VariableType&);


#if 0

  template<>
  ElementPropertyInterpolator<READ_MODEL>::ElementPropertyInterpolator()
  {
  }


  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_MODEL>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    throw csmp::Exception(ERROR, "ElementPropertyInterpolator<READ_MODEL>::Interpolate", "READ_MODEL NYI");
  }


  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_MODEL)
#endif


  template<>
  ElementPropertyInterpolator<READ_ELMT>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_ELMT>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_ELMT)

  template<>
  ElementPropertyInterpolator<READ_NODE>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_NODE>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->N(idx1)->Read(prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_NODE)

  template<>
  ElementPropertyInterpolator<READ_EIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_EIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_EIP)

  template<>
  ElementPropertyInterpolator<READ_FIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_FIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, idx2, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_FIP)

  template<>
  ElementPropertyInterpolator<READ_SIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<READ_SIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, idx2, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(READ_SIP)


  template<>
  ElementPropertyInterpolator<NODE_TO_BCTR>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<NODE_TO_BCTR>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    CalculateN(eptr, elementDimension(*eptr), eptr->FV()->Barycenter(), &coeff[0]);

    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(NODE_TO_BCTR)

  template<>
  ElementPropertyInterpolator<NODE_TO_EIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<NODE_TO_EIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    CalculateN(eptr, elementDimension(*eptr), eptr->IntegrationPoint(idx1), &coeff[0]);

    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(NODE_TO_EIP)

  template<>
  ElementPropertyInterpolator<NODE_TO_FIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<NODE_TO_FIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    CalculateN(eptr, elementDimension(*eptr), eptr->FV()->FacetIntegrationPoint(idx1, idx2), &coeff[0]);

    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(NODE_TO_FIP)

  template<>
  ElementPropertyInterpolator<NODE_TO_SIP>::ElementPropertyInterpolator()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyInterpolator<NODE_TO_SIP>::Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    CalculateN(eptr, elementDimension(*eptr), eptr->FV()->SectorIntegrationPoint(idx1, idx2), &coeff[0]);

    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_INTERPOLATOR(NODE_TO_SIP)





  template<size_t dim>
  Point<dim>
  ElementPlacementOperations<dim,ELEMENT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop)
  {
    auto user = User();
    auto& e = user->e_;
    auto fv = e.FV();
    auto fe = e.FE();

    const size_t element_dim = elementDimension(e);
    const size_t num_nodes = e.Nodes();

    e.CoordinateMatrix();

    std::vector<double64> DN[dim];
    for (size_t i = 0; i < element_dim; ++i) {
      DN[i].resize(num_nodes);
    }
    CalculateDN(e, fv->Barycenter(), DN);

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

  template Point<1ul> csmp::ElementPlacementOperations<1ul, (csmp::PLACEMENT)5>::Gradient(csmp::INDEX<(csmp::VARIABLE_TYPE)1, (csmp::PLACEMENT)17> const&);
  template Point<2ul> csmp::ElementPlacementOperations<2ul, (csmp::PLACEMENT)5>::Gradient(csmp::INDEX<(csmp::VARIABLE_TYPE)1, (csmp::PLACEMENT)17> const&);
  template Point<3ul> csmp::ElementPlacementOperations<3ul, (csmp::PLACEMENT)5>::Gradient(csmp::INDEX<(csmp::VARIABLE_TYPE)1, (csmp::PLACEMENT)17> const&);

  template<size_t dim>
  Point<dim>
  ElementPlacementOperations<dim,FACET_INTEGRATION_POINT>::DirectedArea() const
  {
    auto user = User();
    auto& e = user->e_;
    auto fv = e.FV();
    size_t iFacet(user->idx1_);

    switch (elementDimension(e)) {
      case 1:
      {
        Point<dim> normal;
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
          const Point<dim> n(e.N(iNode)->Coordinate());
          auto weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          normal += weights.first * n;
        }
        return normal;
      }

      case 2:
      {
        Point<dim> tangent;
        Point<dim> bitangent;
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
          const Point<dim> n(e.N(iNode)->Coordinate());
          auto weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          tangent += weights.first * n;
          bitangent += weights.second * n;
        }
        double64 length = exteriorProductLength(tangent, bitangent);
        tangent.NormalizeLengthTo(1.0);
        Point<dim> normal = bitangent - dotProduct(tangent,bitangent) * tangent;
        normal.NormalizeLengthTo(length);
        return normal;
      }

      case 3:
      {
        Point<dim> v0(0.0);
        Point<dim> v1(0.0);
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0; iNode < iNrNodes; ++iNode) {
          auto xform_weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          const Point<dim> n(e.N(iNode)->Coordinate());
          v0 += xform_weights.first * n;
          v1 += xform_weights.second * n;
        }
        return crossProduct(v1, v0);
      }
    }
    throw csmp::Exception(ERROR, "FiniteElementHelper::NormalOfFacet", "Element dimension must be 1, 2, or 3");
  }

  template Point<1u> ElementPlacementOperations<1u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<2u> ElementPlacementOperations<2u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<3u> ElementPlacementOperations<3u,FACET_INTEGRATION_POINT>::DirectedArea() const;





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
    CalculateN(eptr, elementDimension(*eptr), fv->FacetIntegrationPoint( facet, idx3 ), &coeff[0]);
    
    var = 0.;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      var += eptr->N(i)->Read( prop ) * coeff[i];
    }
  }

  INSTANTIATE_FV_PROPERTY_INTERPOLATOR(FV_NODE_TO_FIP)

#if 0
  template<>
  struct ElementPropertyInterpolator<NODE_TO_ELMT>
  {
    const ElementInterpolatorType type_ = NODE_TO_ELMT;

    static const bool needs_recalculation_ = true;

    std::vector<double64> coeff_;

    PropertyInterpolator()
    {
      coeff_.resize(DM_MAX);
    }

    ElementInterpolatorType Type() const { return type_; }

  };

  template<>
  struct ElementPropertyInterpolator<NODE_TO_FIP>
  {
    const ElementInterpolatorType type_ = NODE_TO_FIP;

    static const bool needs_recalculation_ = true;

    size_t facets_, ips_per_facet_, num_nodes_;
    std::vector<double64> coeff_;

    ElementPropertyInterpolator()
    {
      coeff_.resize(DM_MAX);
    }

    ElementInterpolatorType Type() const { return type_; }

    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim )
    {
      auto fv = eptr->FV();
      num_nodes_ = eptr->Nodes();
      facets_ = fv->Facets();
      ips_per_facet_ = fv->IntegrationPointsPerFacet();

      const size_t coeff_size = facets_ * ips_per_facet_ * num_nodes_;

      if (coeff_.size() < coeff_size) {
        coeff_.resize(coeff_size);
      }
      size_t offset = 0;
      for (size_t iFacet = 0; iFacet < facets_; ++iFacet) {
        for (size_t iFip = 0; iFip < ips_per_facet_; ++iFip) {
          CalculateN(eptr, element_dim, fv->FacetIntegrationPoint(iFacet, iFip), &coeff_[offset]);
          offset += num_nodes_;
        }
      }
    }

    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t facet, size_t fip, typename VariableTypeTraits<dim, ty>::VariableType& var )
    {
      auto fv = eptr->FV();
      const size_t num_nodes = eptr->Nodes();
      const size_t offset = (facet * fv->IntegrationPointsPerFacet() + fip) * num_nodes;
      var = 0.;
      for ( size_t i=0; i<num_nodes; i++ ) {
        var += eptr->N(i)->Read( prop ) * coeff_[offset+i];
      }
    }
  };

  struct PropertyInterpolators
  {
    std::bitset<ELMT_INTERPOLATOR_COUNT> interp_valid_;

    ElementPropertyInterpolator<READ_ELMT> read_elmt_;
    ElementPropertyInterpolator<READ_NODE> read_node_;
    // ElementPropertyInterpolator<NODE_TO_ELMT> node_elmt_;
    // ElementPropertyInterpolator<NODE_TO_EIP> node_eip_;
    ElementPropertyInterpolator<NODE_TO_FIP> node_fip_;
    // ElementPropertyInterpolator<NODE_TO_SIP> node_sip_;
  };
#endif

  /*

   template<size_t dim>
   struct FiniteElementHelper<dim>::Impl : public PropertyInterpolators
   {
   Element<dim>* eptr_;
   size_t element_dim_;
   size_t num_nodes_;
   std::vector<double64> DN_bctr_[dim];

   void FiniteElement( Element<dim>* eptr )
   {
   eptr_ = eptr;
   if (eptr_->IsLineElement()) {
   element_dim_ = 1;
   }
   else if (eptr_->IsSurfaceElement()) {
   element_dim_ = 2;
   }
   else if (eptr_->IsVolumeElement()) {
   element_dim_ = 3;
   }

   auto fe = eptr_->FE();
   num_nodes_ = fe->Nodes();

   eptr->CoordinateMatrix();
   interp_valid_.reset();
   }
   };





   template<size_t dim>
   FiniteElementHelper<dim>::FiniteElementHelper()
   : pimpl_(new FiniteElementHelper::Impl())
   {
   }


   template<size_t dim>
   void FiniteElementHelper<dim>::FiniteElement( Element<dim>* eptr )
   {
   pimpl_->FiniteElement(eptr);
   }


   template<size_t dim>
   FiniteElementHelper<dim>::~FiniteElementHelper()
   {
   }


   template<size_t dim>
   Element<dim>* FiniteElementHelper<dim>::FiniteElement( )
   {
   return pimpl_->eptr_;
   }



   template<size_t dim> template<VARIABLE_TYPE ty,PLACEMENT pl>
   void FiniteElementHelper<dim>::ReadAtBarycenter( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var )
   {
   typedef typename InterpolatorDispatch<pl,ELEMENT>::Interpolator Interpolator;
   Interpolator& interpolator = InterpolatorDispatch<pl,ELEMENT>::GetInterpolator(*pimpl_);
   if (Interpolator::needs_recalculation_ && !pimpl_->interp_valid_[interpolator.type_]) {
   pimpl_->interp_valid_[interpolator.type_] = true;
   interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
   }

   interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, 0, 0, var );
   }

   template<size_t dim>
   template<VARIABLE_TYPE ty,PLACEMENT pl>
   void FiniteElementHelper<dim>::ReadAtNode( const csmp::INDEX<ty,pl>& prop, size_t n, typename VariableTypeTraits<dim,ty>::VariableType& var )
   {
   typedef typename InterpolatorDispatch<pl,ELEMENT>::Interpolator Interpolator;
   Interpolator& interpolator = InterpolatorDispatch<pl,ELEMENT>::GetInterpolator(*pimpl_);
   if (Interpolator::needs_recalculation_ && !pimpl_->interp_valid_[interpolator.type_]) {
   pimpl_->interp_valid_[interpolator.type_] = true;
   interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
   }

   interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, n, 0, var );
   }


   template<size_t dim>
   template<VARIABLE_TYPE ty,PLACEMENT pl>
   void FiniteElementHelper<dim>::ReadAtFacetIntegrationPoint( const csmp::INDEX<ty,pl>& prop, size_t facet, size_t fip, typename VariableTypeTraits<dim,ty>::VariableType& var )
   {
   typedef typename InterpolatorDispatch<pl,ELEMENT>::Interpolator Interpolator;
   Interpolator& interpolator = InterpolatorDispatch<pl,ELEMENT>::GetInterpolator(*pimpl_);
   if (Interpolator::needs_recalculation_ && !pimpl_->interp_valid_[interpolator.type_]) {
   pimpl_->interp_valid_[interpolator.type_] = true;
   interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
   }

   interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, facet, fip, var );
   }




   template<size_t dim>
   Point<dim>
   FiniteElementHelper<dim>::ReadGradientAtBarycenter( const csmp::INDEX<SCALAR,NODE>& prop )
   {
   const size_t num_nodes = pimpl_->num_nodes_;
   const size_t element_dim = pimpl_->element_dim_;
   auto& DN_bctr = pimpl_->DN_bctr_;
   auto eptr = pimpl_->eptr_;
   for (unsigned i = 0; i < element_dim; ++i) {
   if (DN_bctr[i].size() < num_nodes) {
   DN_bctr[i].resize(num_nodes);
   }
   }
   CalculateDN(eptr->FV()->Barycenter(), DN_bctr);

   Point<dim> grad(0.);
   for ( size_t i=0U; i<num_nodes; ++i ) {
   const double64 value_at_node(eptr->N(i)->Read(prop));
   for (size_t j = 0; j < element_dim; ++j) {
   grad[j] += DN_bctr[j][i] * value_at_node;
   }
   }
   eptr->FE()->JacobianInverse();

   return Point<dim>(eptr->FE()->JINV * grad.Coordinates());
   }



   template<>
   void
   FiniteElementHelper<1u>::CalculateDN(const Point<1u>& p, std::vector<double64>* DN)
   {
   auto fe = pimpl_->eptr_->FE();
   switch (pimpl_->element_dim_) {
   case 1:
   fe->dNr(p[0], DN[0]);
   fe->Jacobian( DN[0] );

   }
   }


   template<>
   void
   FiniteElementHelper<2u>::CalculateDN(const Point<2u>& p, std::vector<double64>* DN)
   {
   auto fe = pimpl_->eptr_->FE();
   switch (pimpl_->element_dim_) {
   case 1:
   fe->dNr(p[0], DN[0]);
   fe->Jacobian( DN[0] );
   break;

   case 2:
   fe->dNr(p[0], p[1], DN[0]);
   fe->dNs(p[0], p[1], DN[1]);
   fe->Jacobian( DN[0], DN[1] );
   break;
   }
   }


   template<>
   void
   FiniteElementHelper<3u>::CalculateDN(const Point<3u>& p, std::vector<double64>* DN)
   {
   auto fe = pimpl_->eptr_->FE();
   switch (pimpl_->element_dim_) {
   case 3:
   fe->dNr(p[0], p[1], p[2], DN[0]);
   fe->dNs(p[0], p[1], p[2], DN[1]);
   fe->dNt(p[0], p[1], p[2], DN[2]);
   fe->Jacobian( DN[0], DN[1], DN[2] );
   break;
   case 2:
   fe->dNs(p[0], p[1], DN[1]);
   fe->Jacobian( DN[0], DN[1] );
   break;
   case 1:
   fe->dNr(p[0], DN[0]);
   fe->Jacobian( DN[0] );
   break;
   }
   }





   */

}
