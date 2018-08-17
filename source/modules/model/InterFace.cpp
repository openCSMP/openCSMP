#include "InterFace.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "Visitor.h"
#include "variableOperations.h"

using namespace std;

namespace csmp {


template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f )
 : FiniteElementPolicy<dim,::csmp::InterFace>(f),
   idx_(UINT_MAX),
   interface_connector_(f->Neighbors(),nullptr),
   baseElement_( nullptr ),
   current_side_(INSIDE),
   innerParent_( nullptr ),
   outerParent_( nullptr )
 {
 }



template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f,
                           csmp::FiniteVolumeStencil<dim>* fvs )
 : FiniteElementPolicy<dim,::csmp::InterFace>(f),
   FiniteVolumePolicy<dim,::csmp::InterFace>(fvs),
   idx_(UINT_MAX),
   interface_connector_(f->Neighbors(),nullptr),
   baseElement_( nullptr ),
   current_side_(INSIDE),
   innerParent_( nullptr ),
   outerParent_( nullptr )
 {
 }




/// custom constructor which also builds variable storage; used in most cases
template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f,
                           csmp::FiniteVolumeStencil<dim>* fvs,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,::csmp::InterFace>(f),
    FiniteVolumePolicy<dim,::csmp::InterFace>(fvs),
    idx_(UINT_MAX),
    interface_connector_(f->Neighbors(),nullptr),
    baseElement_( nullptr ),
    current_side_(INSIDE),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



template<size_t dim>
InterFace<dim>::InterFace( size_t index,
                           csmp::FiniteElement* f,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,::csmp::InterFace>(f),
    idx_(index),
    interface_connector_(f->Neighbors(),nullptr),
    baseElement_( nullptr ),
    current_side_(INSIDE),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



template<size_t dim>
InterFace<dim>::~InterFace()
 {
 }


/// copy constructor
template<size_t dim>
InterFace<dim>::InterFace( const InterFace<dim>& ifc )
  : FiniteElementPolicy<dim,::csmp::InterFace>(ifc.FE()),
    FiniteVolumePolicy<dim,csmp::InterFace>(ifc.FV()),
    idx_                    ( ifc.idx_),
    interface_connector_    ( ifc.interface_connector_),
    baseElement_            ( ifc.baseElement_ ),
    current_side_           ( ifc.current_side_ ),
    innerParent_            ( ifc.innerParent_ ),
    outerParent_            ( ifc.outerParent_ ),
    inner_parent_face_id_   ( ifc.inner_parent_face_id_ ),
    outer_parent_face_id_   ( ifc.outer_parent_face_id_ )
  {
    assert( !interface_connector_.empty() /* detected unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( ifc.LVS() );
  }


/// move constructor
template<size_t dim>
InterFace<dim>::InterFace( InterFace<dim>&& ifc )
  : FiniteElementPolicy<dim,::csmp::InterFace>(move(ifc.FE())),
    FiniteVolumePolicy<dim,csmp::InterFace>(move(ifc.FV())),
    idx_(move(ifc.idx_)),
    interface_connector_(move(ifc.interface_connector_)),
    baseElement_(move(ifc.baseElement_)),
    current_side_(move(ifc.current_side_)),
    innerParent_(move(ifc.innerParent_)),
    outerParent_(move(ifc.outerParent_)),
    inner_parent_face_id_(move(ifc.inner_parent_face_id_)),
    outer_parent_face_id_(move(ifc.outer_parent_face_id_))
  {
    assert( !interface_connector_.empty() ); // detected unitialized element
    // variable storage: call of initialization function
    this->LVS( move(ifc.LVS()) );
    
    ifc.innerParent_ = nullptr;
    ifc.outerParent_ = nullptr;

    ifc.AssignFiniteElementNullPtr();
    ifc.AssignFiniteVolumeNullPtr();
  }




/// assignment
template<size_t dim>
InterFace<dim>&  InterFace<dim>::operator=( const InterFace<dim>& ifc )
 {
    if ( &ifc != this )
      {
        FiniteElementPolicy<dim,csmp::InterFace>::Assign(ifc.FE());
        FiniteVolumePolicy<dim,csmp::InterFace>::AssignFiniteVolume(ifc.FV());

        idx_                    = ifc.idx_;
        interface_connector_    = ifc.interface_connector_;
        baseElement_             = ifc.baseElement_;
        current_side_            = ifc.current_side_;
        innerParent_             = ifc.innerParent_;
        outerParent_             = ifc.outerParent_;
        inner_parent_face_id_    = ifc.inner_parent_face_id_;
        outer_parent_face_id_    = ifc.outer_parent_face_id_;

        this->LVS( ifc.LVS() );
      }
    return *this;
 }




template<size_t dim>
InterFace<dim>&  InterFace<dim>::operator=( InterFace<dim>&& ifc )
 {
    assert( &ifc != this );
 
    FiniteElementPolicy<dim,csmp::InterFace>::Assign(move(ifc.FE()));
    FiniteVolumePolicy<dim,csmp::InterFace>::AssignFiniteVolume(move(ifc.FV()));

    idx_                    = ifc.idx_;
    interface_connector_    = ifc.interface_connector_;
    baseElement_            = ifc.baseElement_;
    current_side_           = ifc.current_side_;
    innerParent_            = ifc.innerParent_;
    outerParent_            = ifc.outerParent_;
    inner_parent_face_id_   = ifc.inner_parent_face_id_;
    outer_parent_face_id_   = ifc.outer_parent_face_id_;

    this->LVS( move(ifc.LVS()) );

    ifc.AssignFiniteElementNullPtr();
    ifc.AssignFiniteVolumeNullPtr();

    ifc.innerParent_ = nullptr;
    ifc.outerParent_ = nullptr;

    return *this;
 }







/// Roman, 2014
/// WARNING: this operator is used specifically in the process of creation of particular SplitBoundary.
/// Therefore only important infromation for that process is taken into account in order to distinguish two InterFace's.
/// That must be reference to inner and outer parent Elements, inner and outer parent face ID's

template<size_t dim>
bool  InterFace<dim>::operator==( const InterFace<dim>& ifc )
 {
    if ( &ifc != this )
        if( innerParent_ != ifc.innerParent_ ||
            outerParent_ != ifc.outerParent_ ||
            baseElement_ != ifc.baseElement_ ||
            inner_parent_face_id_ != ifc.inner_parent_face_id_ ||
            outer_parent_face_id_ != ifc.outer_parent_face_id_  )
            return false;

    return true;
 }

// VISITOR

template<size_t dim>
void InterFace<dim>::Accept( csmp::Visitor<dim>& vis )
  {
    if ( vis.ApplicationTarget() == INTER_FACE )
      {
        vis.Visit(this);
        return;
      }
    throw csmp::Exception( ERROR, "InterFace<dim>::Accept", "Target of visitation unresolved." );
  } // end Accept



// CONSTRUCTION

template<size_t dim>
void InterFace<dim>::Assign( size_t i, InterFace<dim>* const ifc_ptr ) // neighbor interface
  {
    assert( this->FE() != nullptr );
    assert( interface_connector_.size() == this->FE()->Neighbors() );
    assert( i < this->Neighbors() );
    interface_connector_[i] = ifc_ptr;
  }



template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const base_elmt )
  {
     assert( base_elmt != nullptr );
     baseElement_ = base_elmt;
  }



/**
    Connects face to the higher-dimensional elements which it is sandwiched between.
    These elements are stored in the local pair  of pointers called parents_.
    
    @attention The node numbering of the Face determines the direction of its normal.
    Here the convention is assumed that the first parent element is that on the inside
    of the Face with regard to the outward pointing normal and the second element is on
    the outside.

    Same method as in Face, but - in addition - connects InterFace to the multiplicated nodes
    on either side.
*/
template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const inner_elmt, Element<dim>* const outer_elmt, bool assign_nodes )
  {
     assert( inner_elmt != nullptr );
     assert( outer_elmt != nullptr );
     innerParent_ = inner_elmt;
     outerParent_ = outer_elmt;
     if ( assign_nodes )
       // connect the nodes of the higher dimensional neighbors to the interface
       InitializeNodeVector();
  }


/**
    As above, but with the extra knowledge of the indices of the element faces that are juxtaposed
*/
template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const inner_elmt, size_t inner_local_face_id,
                             Element<dim>* const outer_elmt, size_t outer_local_face_id,
                             bool assign_nodes )
  {
     assert( inner_elmt != nullptr );
     assert( outer_elmt != nullptr );
     innerParent_ = inner_elmt;
     outerParent_ = outer_elmt;
     assert( inner_local_face_id < innerParent_->Faces() );
     assert( outer_local_face_id < outerParent_->Faces() );
     inner_parent_face_id_ = inner_local_face_id;
     outer_parent_face_id_ = outer_local_face_id;
    
     if ( assign_nodes )
       // connect the nodes of the higher dimensional neighbors to the interface
       InitializeNodeVector( inner_parent_face_id_, outer_parent_face_id_ );
  }




template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const parentElement, size_t faceId, INTERFACE_SIDE side )
  {
    assert( parentElement != nullptr );
    assert( side != MIDDLE );

    if( side == INSIDE ) {
        innerParent_ = parentElement;
        inner_parent_face_id_ = faceId;
    }
    else if( side == OUTSIDE ) {
        outerParent_ = parentElement;
        outer_parent_face_id_ = faceId;
    }
    else if( side == MIDDLE ) {
        baseElement_ = parentElement;
        throw csmp::Exception( ERROR, "csmp::InterFace<dim>::Assign( side, Element, face_id )", "This method is not intended to be used for base element!" );
    }

  }




/**
    Assigns nodes to node_connector_ vector of the InterFace 
    
    @note the nodes on the inner side of the interface are the first in the node connector
    vector, the outer ones follow 
    
    @author SKM 16/6/2016
*/
template<size_t dim>
void InterFace<dim>::Assign( size_t n_local, Node<dim>* nptr, INTERFACE_SIDE side )
 {
    assert( !node_connector_.empty() );
    assert( n_local < node_connector_.size() );
    assert( nptr != nullptr );
    assert( side != MIDDLE );
   
    if ( side == INSIDE ) {
         node_connector_[n_local] = nptr;
         return;
      }
      
    // if the side is OUTSIDE, its nodes follow those on the inside
    node_connector_[ node_connector_.size()/2 + n_local ] = nptr;
   
 } // end Assign node pointers




/** 
    finds the local numbers of the faces of the higher-dimensional element that will be connected by the interface; 
    uses point coordinates that must be matched
*/
template<size_t dim>
std::pair<size_t,size_t>  InterFace<dim>::SharedElementFaces()
 {
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
   
    // 1. building search maps that we will use to find the shared interfaces
    //    key=pointset   face iD
    map<set<Point<dim> >,pair<INTERFACE_SIDE,size_t> >   inner_elmt_faces, outer_elmt_faces;
    vector<size_t>  nids;
    // inner parent element
    for ( size_t face=0U; face<innerParent_->Faces(); ++face ) {
          innerParent_->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j=0U; j<nids.size(); ++j )
            face_key.insert( innerParent_->N(j)->Coordinate() );
          inner_elmt_faces.emplace( make_pair(face_key,make_pair(INSIDE,face) ) );
      }
    // outer parent element
    for ( size_t face=0U; face<outerParent_->Faces(); ++face ) {
          outerParent_->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j=0U; j<nids.size(); ++j )
            face_key.insert( outerParent_->N(j)->Coordinate() );
          outer_elmt_faces.emplace( make_pair(face_key,make_pair(OUTSIDE,face) ) );
      }
   
   // 2. finding the shared faces
   vector<pair<set<Point<dim> >,pair<INTERFACE_SIDE,size_t> > > result(2);
   set_intersection( inner_elmt_faces.begin(), inner_elmt_faces.end(),
                     outer_elmt_faces.begin(), outer_elmt_faces.end(), result.begin() );
 
   // assume that we only got a single result
   auto value(result[0].second);
   size_t inner_face_id = (value.first == INSIDE) ? value.second : result[1].second.second;
   size_t outer_face_id = (value.first == OUTSIDE) ? value.second : result[1].second.second;
   
   return make_pair( inner_face_id, outer_face_id );
 
 } // end SharedElementFaces





/**
    1. From the finite element connected to the Interface we know the topology of the faces of the higher-dimensional
       neighbor elements so that we can assign their nodes. Start with the inside element.
       
      (We rely on the fact that the nodes on either side of the interface are collocated)
     
    2. Find the faces connected to this element and assign them
 
    3. Get its nodes 
 
    4. Assign the nodes in the same order as for that face since the normal points into the Interface
 
    5. Repeat for other side
    
    @attention special provisions are made in 1D where there is no corresponding elemnent type for the Interface.
    Thus the number of nodes is assumed to be 1, duplicated to 2.
 
*/
template<size_t dim>
void InterFace<dim>::InitializeNodeVector()
 {
    // 1. get shared faces and assign them
    std::pair<size_t,size_t> shared_faces = SharedElementFaces();
    inner_parent_face_id_ = shared_faces.first;
    outer_parent_face_id_ = shared_faces.second;
 
    // resizing the nodevector
    if ( node_connector_.empty() ) {
         if ( dim == 1U ) node_connector_.resize(2U);
         else node_connector_.resize(this->FE()->Nodes() * 2U);
      }
    node_connector_.shrink_to_fit();
   
    // 2,3. starting with the inside
    vector<size_t>  nids;
   
    // inside
    innerParent_->FE()->NodesOfFace( shared_faces.first, nids );
    // we retain the order in which the nodes are given to
    if ( dim != 1U ) {
         assert( nids.size() == this->FE()->Nodes() );
         for ( size_t n=0U; n<this->FE()->Nodes(); ++n )
           Assign( n, innerParent_->N( nids[n] ), INSIDE );
      }
    // assuming that the unit normal points from the inside to the outside
    else Assign( 0U, innerParent_->N(1), INSIDE );
   
    // outside
    outerParent_->FE()->NodesOfFace( shared_faces.second, nids );
    // we retain the order in which the nodes are given to
    if ( dim != 1U ) {
         assert( nids.size() == this->FE()->Nodes() );
         for ( size_t n=0U; n<this->FE()->Nodes(); ++n )
           Assign( n, outerParent_->N( nids[n] ), OUTSIDE );
      }
    else Assign( 0U, outerParent_->N(0), OUTSIDE );
   
 } // end




/**
     Same as initialise node vector, but for the case where the shared faces have
     already been established.
 
     @note Nodes() cannot not be used here because the node_connector_ vector is
     just getting initialised
*/
template<size_t dim>
void InterFace<dim>::InitializeNodeVector( size_t inner_elmt_face_id,
                                           size_t outer_elmt_face_id )
 {
    assert( inner_elmt_face_id < innerParent_->Faces() );
    assert( outer_elmt_face_id < outerParent_->Faces() );
   
    // 1. get shared faces and assign them
    inner_parent_face_id_ = inner_elmt_face_id;
    outer_parent_face_id_ = outer_elmt_face_id;
 
    // resizing the nodevector (x2 because there are 2 sides of the interface
    if ( node_connector_.empty() ) {
         if ( dim == 1U ) node_connector_.resize(2U);
         else node_connector_.resize(this->FE()->Nodes() * 2U);
      }
    node_connector_.shrink_to_fit();
   
    // 2,3. starting with the inside
    vector<size_t>  nids;
   
    // inside
    innerParent_->FE()->NodesOfFace( inner_parent_face_id_, nids );
    // we retain the order in which the nodes are given to
    if ( dim != 1U ) {
         assert( nids.size() == this->FE()->Nodes() );
         for ( size_t n=0U; n<this->FE()->Nodes(); ++n )
           Assign( n, innerParent_->N( nids[n] ), INSIDE );
      }
    // assuming that the unit normal points from the inside to the outside
    else Assign( 0U, innerParent_->N(1), INSIDE );
    
    // outside
    outerParent_->FE()->NodesOfFace( outer_parent_face_id_, nids );
    // we retain the order in which the nodes are given to
    if ( dim != 1U ) {
         assert( nids.size() == this->FE()->Nodes() );
         for ( size_t n=0U; n<this->FE()->Nodes(); ++n )
           Assign( n, outerParent_->N( nids[n] ), OUTSIDE );
      }
    else Assign( 0U, outerParent_->N(0), OUTSIDE );
   
 } // end






template<size_t dim>
void  InterFace<dim>::CurrentSide( INTERFACE_SIDE side_to_assign )
 {
    current_side_ = side_to_assign;
 }


// ACCESSORS


template<size_t dim>
void  InterFace<dim>::Idx( size_t idx_to_assign ) const
 {
    idx_ = idx_to_assign;
 }

template<size_t dim>
size_t  InterFace<dim>::Idx() const
 {
    return idx_;
 }

template<size_t dim>
INTERFACE_SIDE  InterFace<dim>::CurrentSide() const
 {
    return current_side_;
 }

template<size_t dim>
typename std::vector<csmp::InterFace<dim>*>&  InterFace<dim>::NeighborElementVector()
 {
    return interface_connector_;
 }


/**

Returns pointers to the nodes on either side of the Interface.

@section input Input Arguments

An integer from 0...n-1, where n is the number of nodes per face of the Element.
These are numbered counterclockwise looking from the outside into the face
of parent element 1.

@param side  side refers to the first or second parent element.

@section implementation Implementation

A range check is performed.

@return A pointer to the Target node.
*/
template<size_t dim>
csmp::Node<dim>*  InterFace<dim>::N( size_t n, INTERFACE_SIDE side ) const
 {
    // the number of nodes on a single side of the interface
    const size_t if_FE_nodes(this->FE()->Nodes());
   
    assert( n < if_FE_nodes );

    if ( side == INSIDE )
      return node_connector_[n];
    else if( side == OUTSIDE )
      return node_connector_[n+if_FE_nodes];

    assert( side == MIDDLE );
    if ( baseElement_ != nullptr )
      return baseElement_->N( n );

    throw csmp::Exception( ERROR, "InterFace<dim>::N( local_id, side )", "Base Element does not exist!" );
    return nullptr;
 }


/**
    Access to all nodes of the interface.
*/
template<size_t dim>
csmp::Node<dim>*  InterFace<dim>::N( size_t n ) const
 {
    assert( n < node_connector_.size() );
    return node_connector_[n];
 }


/**
    Returns the equal dimensional neighbor of the InterFace which also is an interface element.
*/
template<size_t dim>
csmp::InterFace<dim>*  InterFace<dim>::Neighbor( size_t n ) const
 {
    assert( interface_connector_.size() == this->Neighbors() );
    assert( n < this->Neighbors() );
    return interface_connector_[n];
 }


/**
    @note costly method, use judiciously.
*/
template<size_t dim>
size_t  InterFace<dim>::ParentNodeNumber( size_t n, INTERFACE_SIDE side ) const
 {
    assert( n < node_connector_.size() );

    if ( side == INSIDE ) {
        // finding the parent element that corresponds to inside element
        for ( size_t i=0U; i<node_connector_[n]->Parents(); ++i )
          if ( node_connector_[n]->Parent(i) == innerParent_ )
            return node_connector_[n]->ParentNodeNumber(i);
      }
    else if ( side == OUTSIDE ) {
        for ( size_t i=0U; i<node_connector_[n]->Parents(); ++i )
          if ( node_connector_[n]->Parent(i) == outerParent_ )
            return node_connector_[n]->ParentNodeNumber(i);
      }
    else if ( side == MIDDLE )
      return n;

    throw csmp::Exception( ERROR, "InterFace<dim>::ParentNodeNumber:",
                          "Node does not seem to be connected to parent element." );
 }





// watch out if there is no base element this returns a nullptr pointer
template<size_t dim>
Element<dim>*  InterFace<dim>::Parent( INTERFACE_SIDE side ) const
 {
    if ( side == INSIDE )
        return innerParent_;
    else if( side == OUTSIDE)
        return outerParent_;
    return baseElement_;
 }

template<size_t dim>
Element<dim>*  InterFace<dim>::InnerParent() const
 {
    return innerParent_;
 }
template<size_t dim>
Element<dim>*  InterFace<dim>::OuterParent() const
 {
    return outerParent_;
 }

template<size_t dim>
Element<dim>*  InterFace<dim>::BaseElement() const
 {
    return baseElement_;
 }


/// return face ID of inner parent element
template<size_t dim>
size_t  InterFace<dim>::InnerParentFaceID() const
 {
     assert( innerParent_ != nullptr );
     return inner_parent_face_id_;
}

/// return face ID of outer parent element
template<size_t dim>
size_t  InterFace<dim>::OuterParentFaceID() const
 {
     assert( outerParent_ != nullptr );
     return outer_parent_face_id_;
}

/// return face ID of inner or outer parent element depending on index
template<size_t dim>
size_t  InterFace<dim>::ParentFaceID( INTERFACE_SIDE side ) const
{
    if( side == INSIDE )
    {
        assert( innerParent_ != nullptr );
        return inner_parent_face_id_;
    }
    else if( side == OUTSIDE )
    {
        assert( outerParent_ != nullptr );
        return outer_parent_face_id_;
    }


    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( WARNING, "csmp::InterFace<dim>::AssignFaceID:", "ID of outer face could not be determined." );

    return inner_parent_face_id_;

}

// GEOMETRY


template<size_t dim>
double64 InterFace<dim>::Area() const
 {
    return this->Volume();
 }


/** returns unit normal into argument vector variable depending on corresponding parent element side
*/
template<size_t dim>
void  InterFace<dim>::UnitNormal( VectorVariable<dim>& vc, INTERFACE_SIDE side ) const
 {
    if( side == INSIDE )
    {
        assert( innerParent_ != nullptr );
        innerParent_->UnitNormalToFace( inner_parent_face_id_, vc );
        return;
    }
    else if( side == OUTSIDE )
    {
        assert( innerParent_ != nullptr );
        outerParent_->UnitNormalToFace( outer_parent_face_id_, vc );
        return;
    }

    if( baseElement_ != nullptr )
        return baseElement_->UnitNormal( vc );

    throw csmp::Exception( ERROR, "InterFace<dim>::UnitNormal( side )", "Base Element does not exist!" );
    innerParent_->UnitNormalToFace( inner_parent_face_id_, vc );
 }




template<size_t dim>
void  InterFace<dim>::UnitNormal( VectorVariable<dim>& vc ) const
 {
    UnitNormal( vc, current_side_ );
 }




/**
   Initialises VectorVariable with vector between node pair - returns false if overlap, true if distant
*/
template<size_t dim>
bool InterFace<dim>::NodeSpacing( size_t n, VectorVariable<dim>& innerToOuter ) const
  {
    assert( n < node_connector_.size()-Nodes() );
    Point<dim> dxyz = node_connector_[n]->Coordinate() - node_connector_[n+Nodes()]->Coordinate();
    innerToOuter    = dxyz;
    
    if ( dxyz.Length() < numeric_limits<double64>::epsilon() ) return false;

    return true;
  }


/**

Returns a matrix with 'nodes'-rows and 'coordinate-directions' columns.
This matrix defines the positions of the elements nodes for
the finite-element matrix assembly. Since the number of element nodes
may vary among different elements types, the number of rows in XY may
also vary from element to element.

@param XY A DenseMatrix<DM_MIN> class object (value type fT). This matrix is dynamically
resized if necessary but must have been constructed with a finite size
before passing it to CoordinateMatrix().

The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<size_t dim>
void  InterFace<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY, INTERFACE_SIDE side ) const
  {
    const size_t n_nodes(Nodes());
    XY.Resize( n_nodes, dim );
    for ( size_t i=0U; i<n_nodes; ++i )
        XY.AssignRow( i, N(i,side)->Coordinate() );

  } // end CoordinateMatrix







/**

BaryCenter() calculates the node coordinate average for the element. This
coordinate value is equivalent to the center of gravity of the Element
type.

@section input Input Arguments

The barycentre is returned into a CSMP vector variable a reference to which
is supplied as single argument.

@section application Application

The element barycentre could be used for instance to output element material
properties from the Model as point data. This is done if you use the
Model OutputToTextFile() methods.

@return The VARIABLE_FLAG of the returned vector variable will not be changed by
BaryCentre().

@note could be madfe more

@test O.K. SKM25/8/14 after refactoring loop

*/
template<size_t dim>
Point<dim>  InterFace<dim>::BaryCenter() const
  {
    Point<dim>    pt(N(0U)->Coordinate());
    const size_t  n_nodes(node_connector_.size());
    
    // all the nodes on both sides
    for ( size_t i=1U; i<n_nodes; ++i )
      pt += N(i)->Coordinate();

    return pt / static_cast<double64>(Nodes());
  }



/**

Measures the length of the element (in physical space) in a certain direction.

@param vecDirection direction in which the element is to be measured

@return Length of the element.

@sectin implementation Implementation

The main idea is to measure the height of the oriented bounding box
(oriented by the direction of vecDirection) which surrounds the element.
This is done by projecting each node onto the direction vector and calculating
the difference between the largest and smallest magnitude. The distance between both
projections will be given by the absolute value of this difference.

@section application Application

The length of the element in a certain dimension is used to weigh the error of the element,
when evaluating the quality of a certain mesh.

@todo this method could be optimised if only the corner nodes would be used

*/
template<size_t dim>
double64  InterFace<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
  {
    double64 fMinTemp( static_cast<double64>( DBL_MAX) );
    double64 fMaxTemp( static_cast<double64>(-DBL_MAX) );

    // this normalisation is necessary because the vector variable
    // being any physical quantity may have any magnitude
    const double64 fMagnitudeOfDirection(vecDirection.Length());
    // avoid division by zero
    assert( fMagnitudeOfDirection >= numeric_limits<double64>::epsilon() );

    const size_t n_nodes(node_connector_.size());
    for ( size_t i=0; i<n_nodes; ++i ) {
        // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
        double64 fTemp(vecDirection.DotProduct( N(i)->Coordinate() ));
        fTemp /= fMagnitudeOfDirection;

        // update minimum value
        fMinTemp = std::min(fMinTemp, fTemp);
        // update maximum value
        fMaxTemp = std::max(fMaxTemp, fTemp);
    }

    //substract magnitudes
    return fMaxTemp - fMinTemp;
 }




/** 
    returns property values at the nodes
*/
template<size_t dim>
template< class Var>
void  InterFace<dim>::NodePropertyVector( const csmp::Index& idx, std::vector<Var>& V, INTERFACE_SIDE side ) const
 {
    if ( idx.place != NODE ) {
         std::cerr <<"\nInterFace<"<< dim;
         std::cerr <<">::NodePropertyVector: Requested property ";
         std::cerr <<"is not placed on the nodes; property Index: "<< std::endl;
         idx.Out();
         return;
      }

    // resizing V if necessary
    const size_t  n_nodes(Nodes());
    V.resize(n_nodes);

    for ( size_t i=0U; i<n_nodes; i++ )
      N(i,side)->Read( idx, V[i] );
 }

// scalar
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>&, INTERFACE_SIDE ) const;
// vector
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<1U> >&, INTERFACE_SIDE ) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<2U> >&, INTERFACE_SIDE ) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<3U> >&, INTERFACE_SIDE ) const;
// tensor
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<1U> >&, INTERFACE_SIDE ) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<2U> >&, INTERFACE_SIDE ) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<3U> >&, INTERFACE_SIDE ) const;
// array
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>&, INTERFACE_SIDE ) const;
// flagged array
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>&, INTERFACE_SIDE ) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>&, INTERFACE_SIDE ) const;





// SCREEN OUTPUT

template<size_t dim>
void  InterFace<dim>::Out() const
 {
    cout <<"\nInterFace<"<< dim <<">::Out: number: "<< idx_;
    cout <<"\nInternal data: ";
    cout <<"\n\tconnected nodes with boundary flags:  ";
    string str;
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i,INSIDE)->AtBoundary());
         cout << N(i,INSIDE)->Idx() <<":"<< str <<"  ";
      }
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i,OUTSIDE)->AtBoundary());
         cout << N(i,OUTSIDE)->Idx() <<":"<< str <<"  ";
      }
    cout << endl;

    cout <<"\n\tconnected neighbor InterFace types / boundary flags:\n";
    for ( size_t i=0U; i<this->Neighbors(); i++ )
      if ( Neighbor(i) != nullptr ) {
           cout <<"\t\t"<< Idx() <<":";
           cout << parseFiniteElementType( Neighbor(i)->FE_Type()) <<": ";
           //str = parseBoundary(Neighbor(i)->AtBoundary());
           //cout << str;
           cout << endl;
        }
      else cout <<"none.  ";
    cout << endl;

    /// @todo (2-P) Remove typeid, implement name fct
    cout <<"\tInterFace is connected via bridge pattern to: ";
    cout << typeid(this).name() << endl;

    Point<dim>  pt(this->BaryCenter());

    if ( dim == 1U )
       cout <<"\n\tBarycentre at (xyz): "<< pt[0] << endl;
    else if ( dim == 2U )
       cout <<"\n\tBarycentre at (xyz): "<< pt[0] <<", "<< pt[1] << endl;
    else
       cout <<"\n\tBarycentre at (xyz): "<< pt[0] <<", "<< pt[1] <<", "<< pt[2] << endl;

    const size_t ipoints(this->IntegrationPoints());
    if ( ipoints > 0U ) {
         cout <<"\n\tStorage sites for IntegrationPoint properties: "<< ipoints << endl;
      }

    cout <<"\n Connected Node objects, side 1 of interface: ";
    for ( size_t i=0U; i<this->Nodes(); i++ )
      node_connector_[i]->Out();
    cout << endl;

    cout <<"\n Connected Node objects, side 2 of interface: ";
    for ( size_t i=0U; i<this->Nodes(); i++ )
      node_connector_[i]->Out();
    cout << endl;

    /// @todo (2-D) Rm rtti
    cout <<"\n\nParent (higher-dimensional) Element objects:     "<< endl;
    if ( innerParent_ != 0 ) {
         cout <<"\tinward  facing Element: ";
         cout  <<" ("<< typeid(*(this->Parent(INSIDE)->FE())).name() <<")"<< endl;
         this->innerParent_->Out();
      }
    else cout <<"\tnone.\n";
    if ( this->outerParent_ != 0 ) {
         cout <<"\toutward facing Element: ";
         cout <<" ("<< typeid(*(this->Parent(OUTSIDE)->FE())).name() <<")"<< endl;
         this->outerParent_->Out();
      }
    else cout <<"\tnone.\n";

    cout <<"\tUnit Normal:            ";
    VectorVariable<dim> un( PLAIN, 0. );
    UnitNormal( un );
    for ( size_t i=0U; i<dim; i++ ) cout << un[i] <<", ";
    cout << endl;

 } // end Out



template class InterFace<1U>;
template class InterFace<2U>;
template class InterFace<3U>;

} // end namespace csmp
