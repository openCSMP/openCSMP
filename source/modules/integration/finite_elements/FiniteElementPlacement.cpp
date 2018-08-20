#include "VariablePlacement.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"


namespace csmp {


#define INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(interp) \
template void ElementPropertyObtainer<interp>::Obtain<1,SCALAR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,SCALAR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<2,SCALAR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,SCALAR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<3,SCALAR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,SCALAR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<1,VECTOR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,VECTOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<2,VECTOR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,VECTOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<3,VECTOR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,VECTOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<1,TENSOR>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,TENSOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<2,TENSOR>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,TENSOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<3,TENSOR>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,TENSOR>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<1,ARRAY>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,ARRAY>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<2,ARRAY>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,ARRAY>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<3,ARRAY>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,ARRAY>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<1,FLAGGEDARRAY>(Index const&, Element<1> const*, size_t, size_t, VariableTypeTraits<1,FLAGGEDARRAY>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<2,FLAGGEDARRAY>(Index const&, Element<2> const*, size_t, size_t, VariableTypeTraits<2,FLAGGEDARRAY>::VariableType&); \
template void ElementPropertyObtainer<interp>::Obtain<3,FLAGGEDARRAY>(Index const&, Element<3> const*, size_t, size_t, VariableTypeTraits<3,FLAGGEDARRAY>::VariableType&);



  template<>
  ElementPropertyObtainer<FE_READ_ELMT>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_READ_ELMT>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_READ_ELMT)

  template<>
  ElementPropertyObtainer<FE_READ_NODE>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_READ_NODE>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->N(idx1)->Read(prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_READ_NODE)

  template<>
  ElementPropertyObtainer<FE_READ_EIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_READ_EIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_READ_EIP)

  template<>
  ElementPropertyObtainer<FE_READ_FIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_READ_FIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, idx2, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_READ_FIP)

  template<>
  ElementPropertyObtainer<FE_READ_SIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_READ_SIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    eptr->Read(idx1, idx2, prop, var);
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_READ_SIP)


  template<>
  ElementPropertyObtainer<FE_NODE_TO_BCTR>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_NODE_TO_BCTR>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    calculateN(*eptr, eptr->FV()->Barycenter(), &coeff[0]);

    var = 0.;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      eptr->N(i)->Read( prop, v );
      var += v * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_NODE_TO_BCTR)

  template<>
  ElementPropertyObtainer<FE_NODE_TO_EIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_NODE_TO_EIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    //calculateN(*eptr, eptr->IntegrationPoint(idx1), &coeff[0]);
    eptr->FE()->N_AtIntegrationPoint( idx1, coeff );

    var = 0.;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      eptr->N(i)->Read( prop, v );
      var += v * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_NODE_TO_EIP)

  template<>
  ElementPropertyObtainer<FE_NODE_TO_FIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_NODE_TO_FIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    calculateN(*eptr, eptr->FV()->FacetIntegrationPoint(idx1, idx2), &coeff[0]);

    var = 0.;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      eptr->N(i)->Read( prop, v );
      var += v * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_NODE_TO_FIP)

  template<>
  ElementPropertyObtainer<FE_NODE_TO_SIP>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_NODE_TO_SIP>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrNodes(eptr->Nodes());
    std::vector<double64> coeff(iNrNodes);
    calculateN(*eptr, eptr->FV()->SectorIntegrationPoint(idx1, idx2), &coeff[0]);

    var = 0.;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t i=0; i<iNrNodes; i++ ) {
      eptr->N(i)->Read( prop, v );
      var += v * coeff[i];
    }
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_NODE_TO_SIP)

  template<>
  ElementPropertyObtainer<FE_FIP_TO_ELMT>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_FIP_TO_ELMT>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrFacets(eptr->Facets());

    var = 0.;
    size_t count = 0;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t iFacet=0; iFacet<iNrFacets; iFacet++ ) {
        const size_t iNrIps = eptr->IntegrationPointsPerFacet();
        for ( size_t iIp=0; iIp<iNrIps; iIp++ ) {
            eptr->Read(iFacet, iIp, prop, v);
            var += v;
            ++count;
        }
    }
    var *= 1.0 / (double64)count;
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_FIP_TO_ELMT)

  template<>
  ElementPropertyObtainer<FE_SIP_TO_ELMT>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_SIP_TO_ELMT>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrSectors(eptr->Nodes());

    var = 0.;
    size_t count = 0;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t iSector=0; iSector<iNrSectors; iSector++ ) {
        const size_t iNrIps = eptr->IntegrationPointsPerSector();
        for ( size_t iIp=0; iIp<iNrIps; iIp++ ) {
            eptr->Read(iSector, iIp, prop, v);
            var += v;
            ++count;
        }
    }
    var *= 1.0 / (double64)count;
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_SIP_TO_ELMT)

  template<>
  ElementPropertyObtainer<FE_EIP_TO_ELMT>::ElementPropertyObtainer()
  {
  }

  template<>
  template<size_t dim,VARIABLE_TYPE ty>
  void
  ElementPropertyObtainer<FE_EIP_TO_ELMT>::Obtain( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    const size_t iNrIps(eptr->IntegrationPoints());

    var = 0.;
    size_t count = 0;
    typename VariableTypeTraits<dim,ty>::VariableType v;
    for ( size_t iIp=0; iIp<iNrIps; iIp++ ) {
        eptr->Read(iIp, prop, v);
        var += v;
        ++count;
    }
    var *= 1.0 / (double64)count;
  }

  INSTANTIATE_ELEMENT_PROPERTY_OBTAINER(FE_EIP_TO_ELMT)


  template<size_t dim>
  Point<dim>
  FiniteElementPlacementOperations<dim,ELEMENT>::Gradient(csmp::INDEX<SCALAR, NODE> const& prop)
  {
  
    /*  
    auto user = User();
    auto& e = user->e_;
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
    */
  
    auto user = User();
    auto& e = user->e_;
    
    const size_t num_nodes = e.Nodes();
    
    DenseMatrix<DM_MIN> DN;
    e.dN_AtBaryCenter( DN );
    
    Point<dim> grad (0.);
    for ( size_t i=0U; i<num_nodes; ++i ) {
      const double64 value_at_node(e.N(i)->Read(prop));
      for (size_t j = 0; j < dim; ++j) {
        grad[j] += DN(j,i) * value_at_node;
      }
    }   
    
    return grad;
  }
  
  template Point<1ul> csmp::FiniteElementPlacementOperations<1ul, ELEMENT>::Gradient(csmp::INDEX<SCALAR, NODE> const&);
  template Point<2ul> csmp::FiniteElementPlacementOperations<2ul, ELEMENT>::Gradient(csmp::INDEX<SCALAR, NODE> const&);
  template Point<3ul> csmp::FiniteElementPlacementOperations<3ul, ELEMENT>::Gradient(csmp::INDEX<SCALAR, NODE> const&);  

 

  template<size_t dim>
  Point<dim>
  FiniteElementPlacementOperations<dim,FACET_INTEGRATION_POINT>::DirectedArea() const
  {
    auto user = User();
    auto& e = user->e_;
    size_t iFacet(user->idx1_);
    return directedAreaOfFacet(e, iFacet);
  }

  template Point<1u> FiniteElementPlacementOperations<1u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<2u> FiniteElementPlacementOperations<2u,FACET_INTEGRATION_POINT>::DirectedArea() const;
  template Point<3u> FiniteElementPlacementOperations<3u,FACET_INTEGRATION_POINT>::DirectedArea() const;

}
