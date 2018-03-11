#ifndef CSMP_FINITE_VOLUME_PLACEMENT_H
#define CSMP_FINITE_VOLUME_PLACEMENT_H

#include "VariableTypeTraits.h"

namespace csmp {

  template<size_t dim> class Element;
  template<size_t dim> class Node;
  template<size_t dim, PLACEMENT pl> class FiniteElementPlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteElementPlacementCollection;
  template<size_t dim, PLACEMENT pl> class FiniteVolumePlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteVolumePlacementCollection;
  template<size_t dim> struct NeighbourNodeCollection;

  enum FiniteVolumeInterpolatorType {
    FV_READ_MODEL,
    FV_READ_REGION,
    FV_READ_ELMT,
    FV_READ_NODE,
    FV_READ_EIP,
    FV_READ_FIP,
    FV_READ_SIP,
    FV_NODE_TO_FIP,
    FV_NODE_TO_SIP,

    FV_INTERPOLATOR_COUNT
  };

  template<FiniteVolumeInterpolatorType interp>
  struct FiniteVolumePropertyInterpolator
  {
    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Node<dim>* eptr, size_t idx1, size_t idx2, size_t idx3,
                  typename VariableTypeTraits<dim,ty>::VariableType& var );

    FiniteVolumePropertyInterpolator();
  };

  template<PLACEMENT from,PLACEMENT to>
  struct FiniteVolumeInterpolatorDispatch
  {
  };

#define CSMP_FV_INTERPOLATOR_DISPATCH(from,to,interp) \
template<> struct FiniteVolumeInterpolatorDispatch<from,to> { static constexpr FiniteVolumeInterpolatorType TYPE = interp; };

  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,MODEL,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,REGION,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,ELEMENT,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,NODE,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,ELEMENT_INTEGRATION_POINT,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,FACET_INTEGRATION_POINT,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(MODEL,ELEMENT,FV_READ_MODEL)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,REGION,FV_READ_REGION)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,ELEMENT,FV_READ_REGION)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,NODE,FV_READ_REGION)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,ELEMENT_INTEGRATION_POINT,FV_READ_REGION)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,FACET_INTEGRATION_POINT,FV_READ_REGION)
  // CSMP_FV_INTERPOLATOR_DISPATCH(REGION,SECTOR_INTEGRATION_POINT,FV_READ_REGION)
  CSMP_FV_INTERPOLATOR_DISPATCH(ELEMENT,ELEMENT,FV_READ_ELMT)
  CSMP_FV_INTERPOLATOR_DISPATCH(ELEMENT,NODE,FV_READ_ELMT)
  CSMP_FV_INTERPOLATOR_DISPATCH(ELEMENT,FACET_INTEGRATION_POINT,FV_READ_ELMT)
  CSMP_FV_INTERPOLATOR_DISPATCH(ELEMENT,SECTOR_INTEGRATION_POINT,FV_READ_ELMT)
  CSMP_FV_INTERPOLATOR_DISPATCH(ELEMENT,ELEMENT_INTEGRATION_POINT,FV_READ_ELMT)
  CSMP_FV_INTERPOLATOR_DISPATCH(NODE,NODE,FV_READ_NODE)
  CSMP_FV_INTERPOLATOR_DISPATCH(FACET_INTEGRATION_POINT,FACET_INTEGRATION_POINT,FV_READ_FIP)
  CSMP_FV_INTERPOLATOR_DISPATCH(NODE,FACET_INTEGRATION_POINT,FV_NODE_TO_FIP)
  CSMP_FV_INTERPOLATOR_DISPATCH(NODE,SECTOR_INTEGRATION_POINT,FV_NODE_TO_SIP)

  template<size_t dim, PLACEMENT pl>
  struct FiniteVolumePlacementWorker
  {
    size_t FeaturesPerSector(Element<dim>& e, size_t);

    size_t PlacementsPerFeature(Element<dim>& e);

    VARIABLE_FLAG Status(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop);

    template<class VarType>
    void Store(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, const VarType& var);

    template<class VarType>
    void Read(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, VarType& var);

    double64 Read(Node<dim>& e, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop);
  };

  template<size_t dim>
  struct FiniteVolumePlacementWorker<dim,NODE>
  {
    std::false_type USES_PARENTS;

    size_t FeaturesPerSector(Element<dim>& e, size_t)
    {
      return 1;
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return 1;
    }

    VARIABLE_FLAG Status(Node<dim>& n, size_t, size_t, size_t, const csmp::Index& prop)
    {
      return n.Status(prop);
    }

    template<class VarType>
    void Store(Node<dim>& n, size_t, size_t, size_t, const csmp::Index& prop, const VarType& var)
    {
      n.Store(prop, var);
    }

    template<class VarType>
    void Read(Node<dim>& n, size_t, size_t, size_t, const csmp::Index& prop, VarType& var)
    {
      n.Read(prop, var);
    }

    double64 Read(Node<dim>& n, size_t, size_t, size_t, const csmp::Index& prop)
    {
      return n.Read(prop);
    }
  };

  template<size_t dim>
  struct FiniteVolumePlacementWorker<dim,ELEMENT>
  {
    std::true_type USES_PARENTS;
    
    size_t FeaturesPerSector(Element<dim>& e, size_t)
    {
      return 1;
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return 1;
    }

    VARIABLE_FLAG Status(Node<dim>& n, size_t idx1, size_t, size_t, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Status(prop);
    }

    template<class VarType>
    void Store(Node<dim>& n, size_t idx1, size_t, size_t, const csmp::Index& prop, const VarType& var)
    {
      n.Parent(idx1)->Store(prop, var);
    }

    template<class VarType>
    void Read(Node<dim>& n, size_t idx1, size_t, size_t, const csmp::Index& prop, VarType& var)
    {
      n.Parent(idx1)->Read(prop, var);
    }

    double64 Read(Node<dim>& n, size_t idx1, size_t, size_t, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Read(prop);
    }
  };

  template<size_t dim>
  struct FiniteVolumePlacementWorker<dim,SECTOR_INTEGRATION_POINT>
  {
    std::true_type USES_PARENTS;
    
    size_t FeaturesPerSector(Element<dim>& e, size_t)
    {
      return 1;
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return e.IntegrationPointsPerSector();
    }

    VARIABLE_FLAG Status(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Status(n.ParentNodeNumber(idx1), idx2, idx3, prop);
    }

    template<class VarType>
    void Store(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, const VarType& var)
    {
      n.Parent(idx1)->Store(n.ParentNodeNumber(idx1), idx2, idx3, prop, var);
    }

    template<class VarType>
    void Read(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, VarType& var)
    {
      n.Parent(idx1)->Read(n.ParentNodeNumber(idx1), idx2, idx3, prop, var);
    }

    double64 Read(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Read(n.ParentNodeNumber(idx1), idx2, idx3, prop);
    }

  };

  template<size_t dim>
  struct FiniteVolumePlacementWorker<dim,FACET_INTEGRATION_POINT>
  {
    std::true_type USES_PARENTS;
    
    size_t FeaturesPerSector(Element<dim>& e, size_t sector)
    {
      return e.FV()->FacetsPerSector( sector );
    }

    size_t PlacementsPerFeature(Element<dim>& e)
    {
      return e.IntegrationPointsPerFacet();
    }

    VARIABLE_FLAG Status(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Status(idx2, idx3, prop);
    }

    template<class VarType>
    void Store(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, const VarType& var)
    {
      n.Parent(idx1)->Store(idx2, idx3, prop, var);
    }

    template<class VarType>
    void Read(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop, VarType& var)
    {
      n.Parent(idx1)->Read(idx2, idx3, prop, var);
    }

    double64 Read(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3, const csmp::Index& prop)
    {
      return n.Parent(idx1)->Read(idx2, idx3, prop);
    }

  };

  template<size_t dim, PLACEMENT pl>
  struct FiniteVolumePlacementIterator
  {
    typedef FiniteVolumePlacement<dim,pl> value_type;

    Node<dim>& n_;
    size_t parent_, num_parents_, idx_, max_idx_, ppf_;

    bool operator==(const FiniteVolumePlacementIterator<dim,pl>& rhs) const
    {
      return &n_ == &rhs.n_ && parent_ == rhs.parent_ && idx_ == rhs.idx_;
    }

    bool operator!=(const FiniteVolumePlacementIterator<dim,pl>& rhs) const
    {
      return &n_ != &rhs.n_ || parent_ != rhs.parent_ || idx_ != rhs.idx_;
    }

    value_type operator*() {
      return value_type(n_, parent_, idx_ / ppf_, idx_ % ppf_);
    }

    void ResetIdx()
    {
      idx_ = 0;
      FiniteVolumePlacementWorker<dim,pl> worker;

      if (worker.USES_PARENTS && parent_ < num_parents_) {
        auto& e = *n_.Parent(parent_);
        size_t fpfv = worker.FeaturesPerSector(e, n_.ParentNodeNumber(parent_));
        size_t ppf = worker.PlacementsPerFeature(e);
        ppf_ = ppf;
        max_idx_ = fpfv * ppf;
      }
    }

    const FiniteVolumePlacementIterator& operator++()
    {
      if (++idx_ >= max_idx_) {
        ++parent_;
        ResetIdx();
      }
      return *this;
    }

    FiniteVolumePlacementIterator operator++(int)
    {
      auto result = *this;
      ++(*this);
      return result;
    }

    FiniteVolumePlacementIterator(Node<dim>& n, size_t parent)
    : n_(n), parent_(parent), num_parents_(n.Parents()), ppf_(1)
    {
      ResetIdx();
    }
  };

  template<size_t dim, PLACEMENT pl>
  struct FiniteVolumePlacementOperations
  {
  };


  template<size_t dim>
  struct FiniteVolumePlacementOperations<dim,FACET_INTEGRATION_POINT>
  {
  private:
    typedef FiniteVolumePlacement<dim,FACET_INTEGRATION_POINT> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:

    Element<dim>& TheElement();
    
    Point<dim> Gradient( const csmp::INDEX<SCALAR,NODE>& prop ) const;

    size_t FacetId() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      return fv->FacetSurroundingSector(pnid, user->idx2_);
    }

    size_t FacetIp() const
    {
      auto user = User();
      return user->idx3_;
    }

    double64 IntegrationWeight() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      size_t facet = fv->FacetSurroundingSector(pnid, user->idx2_);
      return fv->FacetIntegrationWeight(facet, user->idx3_);
    }

    bool FromInside() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      size_t facet = fv->FacetSurroundingSector(pnid, user->idx2_);
      return fv->InsideNode(facet) == pnid;
    }

    Point<dim> DirectedArea() const;
    
    double64 FacetArea() const
    {
        auto user = User();
        auto e = user->n_.Parent(user->idx1_);
        auto fv = e->FV();
        auto pnid = user->n_.ParentNodeNumber(user->idx1_);
        size_t facet = fv->FacetSurroundingSector(pnid, user->idx2_);
        return e->FacetArea(facet);
    }    

    Point<dim> FacetNormal() const
    {
        Point<dim> norm = DirectedArea();
        norm.NormalizeLengthTo(1.0);
        return norm;
    }

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

    FiniteVolumePlacement<dim,NODE> InsideNode() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      size_t facet = fv->FacetSurroundingSector(pnid, user->idx2_);
      auto nid = fv->InsideNode(facet);
      return FiniteVolumePlacement<dim,NODE>(*e->N(nid), 0, 0, 0);
    }

    FiniteVolumePlacement<dim,NODE> OutsideNode() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      size_t facet = fv->FacetSurroundingSector(pnid, user->idx2_);
      auto nid = fv->OutsideNode(facet);
      return FiniteVolumePlacement<dim,NODE>(*e->N(nid), 0, 0, 0);
    }

    FiniteVolumePlacement<dim,NODE> UpstreamNode(double64 ff) const
    {
      return ff < 0 ? OutsideNode() : InsideNode();
    }

    template<PLACEMENT from>
    FiniteVolumePlacement<dim,NODE> UpstreamNode(const csmp::INDEX<VECTOR,from>& vel) const
    {
      VectorVariable<dim> v;
      User()->Interpolate(vel, v);
      return UpstreamNode(ProjectOntoDirectedArea(v));
    }

    FiniteVolumePlacement<dim,NODE> UpstreamNode(const csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT>& ff) const
    {
      return UpstreamNode(this->Interpolate(ff));
    }

    FiniteVolumePlacement<dim,NODE> DownstreamNode(double64 ff) const
    {
      return ff > 0 ? OutsideNode() : InsideNode();
    }

    template<PLACEMENT from>
    FiniteVolumePlacement<dim,NODE> DownstreamNode(const csmp::INDEX<VECTOR,from>& vel) const
    {
      return DownstreamNode(ProjectOntoDirectedArea(vel));
    }

    FiniteVolumePlacement<dim,NODE> DownstreamNode(const csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT>& ff) const
    {
      return DownstreamNode(this->Interpolate(ff));
    }
  };


  template<size_t dim>
  struct FiniteVolumePlacementOperations<dim,NODE>
  {
  private:
    typedef FiniteVolumePlacement<dim,NODE> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }
      
  public:
    NeighbourNodeCollection<dim> AllNeighbourNodes() const
    {
       return NeighbourNodeCollection<dim>(User()->n_);
    }

    FiniteVolumePlacementCollection<dim,FACET_INTEGRATION_POINT> AllFacetIntegrationPoints() const
    {
       return FiniteVolumePlacementCollection<dim,FACET_INTEGRATION_POINT>(User()->n_);
    }

    FiniteVolumePlacementCollection<dim,SECTOR_INTEGRATION_POINT> AllSectorIntegrationPoints() const
    {
       return FiniteVolumePlacementCollection<dim,SECTOR_INTEGRATION_POINT>(User()->n_);
    }

  };

  template<size_t dim>
  struct FiniteVolumePlacementOperations<dim,SECTOR_INTEGRATION_POINT>
  {
  private:
    typedef FiniteVolumePlacement<dim,SECTOR_INTEGRATION_POINT> user_type;

    const user_type* User() const
    {
      return static_cast<const user_type*>(this);
    }

  public:
    double64 IntegrationWeight() const
    {
      auto user = User();
      auto e = user->n_.Parent(user->idx1_);
      auto fv = e->FV();
      auto pnid = user->n_.ParentNodeNumber(user->idx1_);
      return fv->SectorIntegrationWeight(pnid, user->idx2_);
    }

    size_t NodeIdx() const
    {
        auto user = User();
        user->n_.Idx();
    }

    double64 SectorVolume() const
    {
        auto user = User();
        auto e = user->n_.Parent(user->idx1_);
        auto pnid = user->n_.ParentNodeNumber(user->idx1_);
        return e->SectorVolume(pnid);
    }
  };

  template<size_t dim, PLACEMENT pl>
  class FiniteVolumePlacement : public FiniteVolumePlacementOperations<dim,pl>
  {
  private:
    friend struct FiniteVolumePlacementOperations<dim,pl>;
    Node<dim>& n_;
    size_t idx1_, idx2_, idx3_;

  public:
    FiniteVolumePlacement(Node<dim>& n, size_t idx1, size_t idx2, size_t idx3)
    : n_(n), idx1_(idx1), idx2_(idx2), idx3_(idx3)
    {
    }

    size_t NodeIdx() const
    {
        return n_.Idx();
    }

    //Point<dim> Gradient( const csmp::INDEX<SCALAR,NODE>& prop ) const;

    template<VARIABLE_TYPE ty>
    void Read( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var ) const
    {
      FiniteVolumePlacementWorker<dim,pl> worker;
      worker.Read(n_, idx1_, idx2_, idx3_, prop, var);
    }

    double64 Read( const csmp::INDEX<SCALAR,pl>& prop ) const
    {
      ScalarVariable var;
      Read(prop, var);
      return var();
    }

    template<VARIABLE_TYPE ty>
    VARIABLE_FLAG Status( const csmp::INDEX<ty,pl>& prop ) const
    {
      FiniteVolumePlacementWorker<dim,pl> worker;
      return worker.Status(n_, idx1_, idx2_, idx3_, prop);
    }

    // When values are stored back, by default flags should not be overwritten.
    template<VARIABLE_TYPE ty>
    void Store( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType const & var )
    {
      FiniteVolumePlacementWorker<dim,pl> worker;
      worker.Store(n_, idx1_, idx2_, idx3_, prop, var);
    }


    template<VARIABLE_TYPE ty, PLACEMENT from>
    void Interpolate( const csmp::INDEX<ty,from>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var ) const
    {
        FiniteVolumePropertyInterpolator<FiniteVolumeInterpolatorDispatch<from,pl>::TYPE> interp;
        interp.template Interpolate<dim,ty>( prop, &n_, idx1_, idx2_, idx3_, var );
    }

    template<PLACEMENT from>
    double64 Interpolate( const csmp::INDEX<SCALAR,from>& prop ) const
    {
        FiniteVolumePropertyInterpolator<FiniteVolumeInterpolatorDispatch<from,pl>::TYPE> interp;
        ScalarVariable var;
        interp.template Interpolate<dim,SCALAR>( prop, &n_, idx1_, idx2_, idx3_, var );
        return var();
    }

    template<PLACEMENT from>
    Point<dim> Interpolate( const csmp::INDEX<VECTOR,from>& prop ) const
    {
        FiniteVolumePropertyInterpolator<FiniteVolumeInterpolatorDispatch<from,pl>::TYPE> interp;
        VectorVariable<dim> var;
        interp.template Interpolate<dim,VECTOR>( prop, &n_, idx1_, idx2_, idx3_, var );
        return var.P();
    }
  };



  template<size_t dim, PLACEMENT pl>
  struct FiniteVolumePlacementCollection
  {
    typedef FiniteVolumePlacement<dim,pl> value_type;
    typedef FiniteVolumePlacementIterator<dim,pl> iterator;
    typedef FiniteVolumePlacementIterator<dim,pl> const_iterator;

    Node<dim>& n_;
    size_t idx_max_;

    FiniteVolumePlacementCollection(Node<dim>& n)
    : n_(n)
    {
      FiniteVolumePlacementWorker<dim,pl> worker;
      idx_max_ = worker.USES_PARENTS ? n_.Parents() : 1;
    }

    iterator begin() { return iterator(n_, 0); }
    const_iterator begin() const { return const_iterator(n_, 0); }
    const_iterator cbegin() const { return const_iterator(n_, 0); }

    iterator end() { return iterator(n_, idx_max_); }
    const_iterator end() const { return const_iterator(n_, idx_max_); }
    const_iterator cend() const { return const_iterator(n_, idx_max_); }
  };

  template<size_t dim>
  struct NeighbourNodeIterator
  {
    typedef FiniteVolumePlacement<dim,NODE> value_type;

    Node<dim>& n_;
    size_t idx_, max_idx_;

    bool operator==(const NeighbourNodeIterator<dim>& rhs) const
    {
      return &n_ == &rhs.n_ && idx_ == rhs.idx_;
    }

    bool operator!=(const NeighbourNodeIterator<dim>& rhs) const
    {
      return &n_ != &rhs.n_ || idx_ != rhs.idx_;
    }

    value_type operator*() {
      return value_type(n_, idx_, 0, 0);
    }

    const NeighbourNodeIterator& operator++()
    {
      ++idx_;
      return *this;
    }

    NeighbourNodeIterator operator++(int)
    {
      auto result = *this;
      ++(*this);
      return result;
    }

    NeighbourNodeIterator(Node<dim>& n, size_t idx)
      : n_(n), idx_(idx)
    {
      max_idx_ = n.Neighbors();
    }
  };

  template<size_t dim>
  struct NeighbourNodeCollection
  {
    typedef FiniteVolumePlacement<dim,NODE> value_type;
    typedef NeighbourNodeIterator<dim> iterator;
    typedef NeighbourNodeIterator<dim> const_iterator;

    Node<dim>& n_;
    size_t idx_max_;

    NeighbourNodeCollection(Node<dim>& n)
      : n_(n)
    {
      idx_max_ = n.Neighbors();
    }

    iterator begin() { return iterator(n_, 0); }
    const_iterator begin() const { return const_iterator(n_, 0); }
    const_iterator cbegin() const { return const_iterator(n_, 0); }

    iterator end() { return iterator(n_, idx_max_); }
    const_iterator end() const { return const_iterator(n_, idx_max_); }
    const_iterator cend() const { return const_iterator(n_, idx_max_); }
  };

}

#endif
