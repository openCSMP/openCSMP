#ifndef CSMP_FINITE_ELEMENT_PLACEMENT_H
#define CSMP_FINITE_ELEMENT_PLACEMENT_H

#include "VariableTypeTraits.h"

namespace csmp {

  template<size_t dim> class Element;
  template<size_t dim> class Node;
  template<size_t dim, PLACEMENT pl> class FiniteElementPlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteElementPlacementCollection;
  template<size_t dim, PLACEMENT pl> class FiniteVolumePlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteVolumePlacementCollection;

  enum FiniteElementInterpolatorType {
    FE_READ_MODEL,
    FE_READ_REGION,
    FE_READ_ELMT,
    FE_READ_NODE,
    FE_READ_EIP,
    FE_READ_FIP,
    FE_READ_SIP,
    FE_NODE_TO_BCTR,
    FE_NODE_TO_EIP,
    FE_NODE_TO_FIP,
    FE_NODE_TO_SIP,
    FE_FIP_TO_ELMT,
    FE_SIP_TO_ELMT,
    FE_EIP_TO_ELMT,
    FE_EIP_TO_FIP,
    FE_EIP_TO_SIP,
    FE_SIP_TO_FIP,
    FE_SIP_TO_EIP,
    FE_FIP_TO_EIP,
    FE_FIP_TO_SIP,
    FE_INTERPOLATOR_COUNT
  };

  template<FiniteElementInterpolatorType interp>
  struct ElementPropertyInterpolator
  {
    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2,
                  typename VariableTypeTraits<dim,ty>::VariableType& var );

    ElementPropertyInterpolator();
  };

  template<PLACEMENT from,PLACEMENT to>
  struct ElementInterpolatorDispatch
  {
  };

#define CSMP_FE_INTERPOLATOR_DISPATCH(from,to,interp) \
template<> struct ElementInterpolatorDispatch<from,to> { static constexpr FiniteElementInterpolatorType TYPE = interp; };

  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,MODEL,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,REGION,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,ELEMENT,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,NODE,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,ELEMENT_INTEGRATION_POINT,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,FACET_INTEGRATION_POINT,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(MODEL,SECTOR_INTEGRATION_POINT,FE_READ_MODEL)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,REGION,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,ELEMENT,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,NODE,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,ELEMENT_INTEGRATION_POINT,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,FACET_INTEGRATION_POINT,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(REGION,SECTOR_INTEGRATION_POINT,FE_READ_REGION)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT,ELEMENT,FE_READ_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT,NODE,FE_READ_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT,ELEMENT_INTEGRATION_POINT,FE_READ_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT,FACET_INTEGRATION_POINT,FE_READ_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT,SECTOR_INTEGRATION_POINT,FE_READ_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(NODE,NODE,FE_READ_NODE)
  CSMP_FE_INTERPOLATOR_DISPATCH(NODE,ELEMENT,FE_NODE_TO_BCTR)
  CSMP_FE_INTERPOLATOR_DISPATCH(NODE,FACET_INTEGRATION_POINT,FE_NODE_TO_FIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(NODE,ELEMENT_INTEGRATION_POINT,FE_NODE_TO_EIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(NODE,SECTOR_INTEGRATION_POINT,FE_NODE_TO_SIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT_INTEGRATION_POINT,ELEMENT_INTEGRATION_POINT,FE_READ_EIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(FACET_INTEGRATION_POINT,FACET_INTEGRATION_POINT,FE_READ_FIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(SECTOR_INTEGRATION_POINT,SECTOR_INTEGRATION_POINT,FE_READ_SIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(FACET_INTEGRATION_POINT,ELEMENT,FE_FIP_TO_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(SECTOR_INTEGRATION_POINT,ELEMENT,FE_SIP_TO_ELMT)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT_INTEGRATION_POINT,ELEMENT,FE_EIP_TO_ELMT)
# if 0
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT_INTEGRATION_POINT,FACET_INTEGRATION_POINT,FE_EIP_TO_FIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(ELEMENT_INTEGRATION_POINT,SECTOR_INTEGRATION_POINT,FE_EIP_TO_SIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(SECTOR_INTEGRATION_POINT,FACET_INTEGRATION_POINT,FE_SIP_TO_FIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(SECTOR_INTEGRATION_POINT,ELEMENT_INTEGRATION_POINT,FE_SIP_TO_EIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(FACET_INTEGRATION_POINT,ELEMENT_INTEGRATION_POINT,FE_FIP_TO_EIP)
  CSMP_FE_INTERPOLATOR_DISPATCH(FACET_INTEGRATION_POINT,SECTOR_INTEGRATION_POINT,FE_FIP_TO_SIP)
#endif

  template<size_t dim, PLACEMENT pl>
  struct FiniteElementPlacementOperations
  {
  };


  template<size_t dim>
  struct FiniteElementPlacementOperations<dim,NODE>
  {
  private:
    typedef FiniteElementPlacement<dim,NODE> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:
    size_t NodeIdx() const
    {
      auto user = User();
      return user->e_.N(user->idx1_)->Idx();
    }

    double64 SectorVolume() const
    {
      auto user = User();
      return user->e_.SectorVolume(user->idx1_);
    }

    NeighbourNodeCollection<dim> AllNeighbourNodes() const
    {
      auto user = User();
      return NeighbourNodeCollection<dim>(*user->e_.N(user->idx1_));
    }
  };

  template<size_t dim>
  struct FiniteElementPlacementOperations<dim,ELEMENT_INTEGRATION_POINT>
  {
  private:
    typedef FiniteElementPlacement<dim,ELEMENT_INTEGRATION_POINT> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:
    double64 IntegrationWeight() const
    {
      return User()->e_.IntegrationWeight();
    }
  };

  template<size_t dim>
  struct FiniteElementPlacementOperations<dim,FACET_INTEGRATION_POINT>
  {
  private:
    typedef FiniteElementPlacement<dim,FACET_INTEGRATION_POINT> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:
    double64 IntegrationWeight() const
    {
      auto user = User();
      auto fv = user->e_.FV();
      return fv->FacetIntegrationWeight(user->idx1_, user->idx2_);
    }

    Point<dim> DirectedArea() const;

    double64 ProjectOntoDirectedArea(const VectorVariable<dim>& v) const
    {
      return dotProduct(v.P(), DirectedArea());
    }

    double64 ProjectOntoDirectedArea(const Point<dim>& v) const
    {
      return dotProduct(v, DirectedArea());
    }

    template<PLACEMENT from>
    double64 ProjectOntoDirectedArea(const csmp::INDEX<VECTOR,from>& prop) const
    {
      VectorVariable<dim> v;
      User()->Interpolate(prop, v);
      return ProjectOntoDirectedArea(v);
    }

    Point<dim> FacetNormal() const
    {
        Point<dim> norm = DirectedArea();
        norm.NormalizeLengthTo(1.0);
        return norm;
    }

    double64 ProjectOntoFacetNormal(const VectorVariable<dim>& v) const
    {
      return dotProduct(v.P(), FacetNormal());
    }

    double64 ProjectOntoFacetNormal(const Point<dim>& v) const
    {
      return dotProduct(v, FacetNormal());
    }

    template<PLACEMENT from>
    double64 ProjectOntoFacetNormal(const csmp::INDEX<VECTOR,from>& prop) const
    {
      VectorVariable<dim> v;
      User()->Interpolate(prop, v);
      return ProjectOntoFacetNormal(v);
    }

    FiniteElementPlacement<dim,NODE> InsideNode() const
    {
      auto user = User();
      auto nid = user->e_.FV()->InsideNode(user->idx1_);
      return FiniteElementPlacement<dim,NODE>(user->e_, nid, 0);
    }

    FiniteElementPlacement<dim,NODE> OutsideNode() const
    {
      auto user = User();
      auto nid = user->e_.FV()->OutsideNode(user->idx1_);
      return FiniteElementPlacement<dim,NODE>(user->e_, nid, 0);
    }

    FiniteElementPlacement<dim,NODE> UpstreamNode(double64 ff) const
    {
      return ff < 0 ? OutsideNode() : InsideNode();
    }

    template<PLACEMENT from>
    FiniteElementPlacement<dim,NODE> UpstreamNode(const csmp::INDEX<VECTOR,from>& vel) const
    {
      VectorVariable<dim> v;
      User()->Interpolate(vel, v);
      return UpstreamNode(ProjectOntoDirectedArea(v));
    }

    FiniteElementPlacement<dim,NODE> UpstreamNode(const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>& ff) const
    {
      return UpstreamNode(this->Interpolate(ff));
    }

    FiniteElementPlacement<dim,NODE> DownstreamNode(double64 ff) const
    {
      return ff > 0 ? OutsideNode() : InsideNode();
    }

    template<PLACEMENT from>
    FiniteElementPlacement<dim,NODE> DownstreamNode(const csmp::INDEX<VECTOR,from>& vel) const
    {
      return DownstreamNode(ProjectOntoDirectedArea(vel));
    }

    FiniteElementPlacement<dim,NODE> DownstreamNode(const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>& ff) const
    {
      return DownstreamNode(this->Interpolate(ff));
    }
  };


  template<size_t dim, PLACEMENT pl>
  struct FiniteElementPlacementWorker
  {
    size_t FeaturesPerElement(Element<dim>& e);

    size_t PlacementsPerFeature(Element<dim>& e);

    VARIABLE_FLAG Status(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop);

    template<class VarType>
    void Store(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, const VarType& var);

    template<class VarType>
    void Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, VarType& var);

    double64 Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop);
  };


  template<size_t dim>
  struct FiniteElementPlacementWorker<dim,ELEMENT>
  {
    size_t FeaturesPerElement(Element<dim>& e)
    {
      return 1;
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return 1;
    }

    VARIABLE_FLAG Status(Element<dim>& e, size_t, size_t, const csmp::Index& prop)
    {
      return e.Status(prop);
    }

    template<class VarType>
    void Store(Element<dim>& e, size_t, size_t, const csmp::Index& prop, const VarType& var)
    {
      e.Store(prop, var);
    }

    template<class VarType>
    void Read(Element<dim>& e, size_t, size_t, const csmp::Index& prop, VarType& var)
    {
      e.Read(prop, var);
    }

    double64 Read(Element<dim>& e, size_t, size_t, const csmp::Index& prop)
    {
      return e.Read(prop);
    }
  };


  template<size_t dim>
  struct FiniteElementPlacementWorker<dim,NODE>
  {
    size_t FeaturesPerElement(Element<dim>& e)
    {
      return e.Nodes();
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return 1;
    }

    VARIABLE_FLAG Status(Element<dim>& e, size_t idx1, size_t, const csmp::Index& prop)
    {
      return e.N(idx1)->Status(prop);
    }

    template<class VarType>
    void Store(Element<dim>& e, size_t idx1, size_t, const csmp::Index& prop, const VarType& var)
    {
      e.N(idx1)->Store(prop, var);
    }

    template<class VarType>
    void Read(Element<dim>& e, size_t idx1, size_t, const csmp::Index& prop, VarType& var)
    {
      e.N(idx1)->Read(prop, var);
    }

    double64 Read(Element<dim>& e, size_t idx1, size_t, const csmp::Index& prop)
    {
      return e.N(idx1)->Read(prop);
    }
  };

  template<size_t dim>
  struct FiniteElementPlacementWorker<dim,ELEMENT_INTEGRATION_POINT>
  {
    size_t FeaturesPerElement(Element<dim>& e)
    {
      return 1;
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return e.IntegrationPoints();
    }

    VARIABLE_FLAG Status(Element<dim>& e, size_t, size_t idx2, const csmp::Index& prop)
    {
      return e.Status(idx2, prop);
    }

    template<class VarType>
    void Store(Element<dim>& e, size_t, size_t idx2, const csmp::Index& prop, const VarType& var)
    {
      e.Store(idx2, prop, var);
    }

    template<class VarType>
    void Read(Element<dim>& e, size_t, size_t idx2, const csmp::Index& prop, VarType& var)
    {
      e.Read(idx2, prop, var);
    }

    double64 Read(Element<dim>& e, size_t, size_t idx2, const csmp::Index& prop)
    {
      return e.Read(idx2, prop);
    }
  };

  template<size_t dim>
  struct FiniteElementPlacementWorker<dim,FACET_INTEGRATION_POINT>
  {
    size_t FeaturesPerElement(Element<dim>& e)
    {
      return e.Facets();
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return e.IntegrationPointsPerFacet();
    }

    VARIABLE_FLAG Status(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop)
    {
      return e.Status(idx1, idx2, prop);
    }

    template<class VarType>
    void Store(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, const VarType& var)
    {
      e.Store(idx1, idx2, prop, var);
    }

    template<class VarType>
    void Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, VarType& var)
    {
      e.Read(idx1, idx2, prop, var);
    }

    double64 Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop)
    {
      return e.Read(idx1, idx2, prop);
    }
  };

  template<size_t dim>
  struct FiniteElementPlacementWorker<dim,SECTOR_INTEGRATION_POINT>
  {
    size_t FeaturesPerElement(Element<dim>& e)
    {
      return e.Nodes();
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return e.IntegrationPointsPerSector();
    }

    VARIABLE_FLAG Status(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop)
    {
      return e.Status(idx1, idx2, prop);
    }

    template<class VarType>
    void Store(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, const VarType& var)
    {
      e.Store(idx1, idx2, prop, var);
    }

    template<class VarType>
    void Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop, VarType& var)
    {
      e.Read(idx1, idx2, prop, var);
    }

    double64 Read(Element<dim>& e, size_t idx1, size_t idx2, const csmp::Index& prop)
    {
      return e.Read(idx1, idx2, prop);
    }
  };


  template<size_t dim>
  struct FiniteElementPlacementOperations<dim,ELEMENT>
  {
  private:
    typedef FiniteElementPlacement<dim,ELEMENT> user_type;
    
    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }
    
  public:
      Point<dim> Gradient(const csmp::INDEX<SCALAR,NODE>& prop);
  };

  template<size_t dim>
  struct FiniteElementPlacementOperations<dim,SECTOR_INTEGRATION_POINT>
  {
  private:
    typedef FiniteElementPlacement<dim,SECTOR_INTEGRATION_POINT> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:
    double64 IntegrationWeight() const
    {
      auto user = User();
      auto fv = user->e_.FV();
      return fv->SectorIntegrationWeight(user->idx1_, user->idx2_);
    }

    size_t NodeIdx() const
    {
      auto user = User();
      return user->e_.N(user->idx1_)->Idx();
    }

    double64 SectorVolume() const
    {
      auto user = User();
      return user->e_.SectorVolume(user->idx1_);
    }
  };


  template<size_t dim, PLACEMENT pl>
  class FiniteElementPlacement : public FiniteElementPlacementOperations<dim,pl>
  {
  private:
    friend struct FiniteElementPlacementOperations<dim,pl>;

    Element<dim>& e_;
    size_t idx1_, idx2_;

  public:
    FiniteElementPlacement(Element<dim>& e, size_t idx1, size_t idx2)
    : e_(e), idx1_(idx1), idx2_(idx2)
    {
    }

    template<VARIABLE_TYPE ty>
    void Read( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var ) const
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      worker.Read(e_, idx1_, idx2_, prop, var);
    }

    double Read( const csmp::INDEX<SCALAR,pl>& prop ) const
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      return worker.Read(e_, idx1_, idx2_, prop);
    }

    Point<dim> Read( const csmp::INDEX<VECTOR,pl>& prop ) const
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      VectorVariable<dim> var;
      worker.Read(e_, idx1_, idx2_, prop, var);
      return var.P();
    }

    template<VARIABLE_TYPE ty>
    VARIABLE_FLAG Status( const csmp::INDEX<ty,pl>& prop ) const
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      return worker.Status(e_, idx1_, idx2_, prop);
    }

    // When values are stored back, by default flags should not be overwritten.
    template<VARIABLE_TYPE ty>
    void Store( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType const & var )
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      worker.Store(e_, idx1_, idx2_, prop, var);
    }

    template<VARIABLE_TYPE ty, PLACEMENT from>
    void Interpolate( const csmp::INDEX<ty,from>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var ) const
    {
        ElementPropertyInterpolator<ElementInterpolatorDispatch<from,pl>::TYPE> interp;
        interp.template Interpolate<dim,ty>( prop, &e_, idx1_, idx2_, var );
    }

    template<PLACEMENT from>
    double64 Interpolate( const csmp::INDEX<SCALAR,from>& prop ) const
    {
        ElementPropertyInterpolator<ElementInterpolatorDispatch<from,pl>::TYPE> interp;
        ScalarVariable var;
        interp.template Interpolate<dim,SCALAR>( prop, &e_, idx1_, idx2_, var );
        return var();
    }

    template<PLACEMENT from>
    Point<dim> Interpolate( const csmp::INDEX<VECTOR,from>& prop ) const
    {
        ElementPropertyInterpolator<ElementInterpolatorDispatch<from,pl>::TYPE> interp;
        VectorVariable<dim> var;
        interp.template Interpolate<dim,VECTOR>( prop, &e_, idx1_, idx2_, var );
        return var.P();
    }

  };



  template<size_t dim, PLACEMENT pl>
  struct FiniteElementPlacementIterator
  {
    typedef FiniteElementPlacement<dim,pl> value_type;

    Element<dim>& e_;
    size_t idx1_, idx2_, max_idx2_;

    bool operator==(const FiniteElementPlacementIterator<dim,pl>& rhs) const
    {
      return &e_ == &rhs.e_ && idx1_ == rhs.idx1_ && idx2_ == rhs.idx2_;
    }

    bool operator!=(const FiniteElementPlacementIterator<dim,pl>& rhs) const
    {
      return &e_ != &rhs.e_ || idx1_ != rhs.idx1_ || idx2_ != rhs.idx2_;
    }

    value_type operator*() {
      return value_type(e_, idx1_, idx2_);
    }

    const FiniteElementPlacementIterator& operator++()
    {
      if (++idx2_ >= max_idx2_) {
        ++idx1_;
        idx2_ = 0;
      }
      return *this;
    }

    FiniteElementPlacementIterator operator++(int)
    {
      auto result = *this;
      ++(*this);
      return result;
    }

    FiniteElementPlacementIterator(Element<dim>& e, size_t idx1)
    : e_(e), idx1_(idx1), idx2_(0)
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      max_idx2_ = worker.PlacementsPerFeature(e);
    }
  };

  template<size_t dim, PLACEMENT pl>
  struct FiniteElementPlacementCollection
  {
    typedef FiniteElementPlacement<dim,pl> value_type;
    typedef FiniteElementPlacementIterator<dim,pl> iterator;
    typedef FiniteElementPlacementIterator<dim,pl> const_iterator;

    Element<dim>& e_;
    size_t idx_max_;

    FiniteElementPlacementCollection(Element<dim>& e)
    : e_(e)
    {
      FiniteElementPlacementWorker<dim,pl> worker;
      idx_max_ = worker.FeaturesPerElement(e);
    }

    iterator begin() { return iterator(e_, 0); }
    const_iterator begin() const { return const_iterator(e_, 0); }
    const_iterator cbegin() const { return const_iterator(e_, 0); }

    iterator end() { return iterator(e_, idx_max_); }
    const_iterator end() const { return const_iterator(e_, idx_max_); }
    const_iterator cend() const { return const_iterator(e_, idx_max_); }
  };


}

#endif
