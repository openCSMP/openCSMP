#include "InterFace.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "Visitor.h"

using namespace std;

namespace csmp {

template<size_t dim>
InterFace<dim>::InterFace()
 : idx_(UINT_MAX),
   fptr_(0),
   fvptr_(0),
   baseElement_( NULL ),
   current_side_(INSIDE),
   innerParent_( NULL ),
   outerParent_( NULL )
 {
 }



template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f )
 : idx_(UINT_MAX),
   fptr_(f),
   fvptr_(0),
   interface_connector_(f->Neighbors(),NULL),
   baseElement_( NULL ),
   current_side_(INSIDE),
   innerParent_( NULL ),
   outerParent_( NULL )
 {
 }



template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f,
                           csmp::FiniteVolumeStencil<dim>* fvs )
 : idx_(UINT_MAX),
   fptr_(f),
   fvptr_(fvs),
   interface_connector_(f->Neighbors(),NULL),
   baseElement_( NULL ),
   current_side_(INSIDE),
   innerParent_( NULL ),
   outerParent_( NULL )
 {
 }




/// custom constructor which also builds variable storage; used in most cases
template<size_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f,
                           csmp::FiniteVolumeStencil<dim>* fvs,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : idx_(UINT_MAX),
    fptr_(f),
    fvptr_(fvs),
    interface_connector_(f->Neighbors(),nullptr),
    baseElement_( nullptr ),
    current_side_(INSIDE),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



template<size_t dim>
InterFace<dim>::InterFace( size_t index,
                           csmp::FiniteElement* f,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : idx_(index),
    fptr_(f),
    fvptr_(nullptr),
    interface_connector_(f->Neighbors(),nullptr),
    baseElement_( nullptr ),
    current_side_(INSIDE),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
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
  : idx_                    ( ifc.idx_                      ),
    fptr_                   ( ifc.fptr_                     ),
    fvptr_                  ( ifc.fvptr_                    ),
    interface_connector_    ( ifc.interface_connector_      ),
    parent_elements_node_connector_   ( ifc.parent_elements_node_connector_  ),
    baseElement_            ( ifc.baseElement_              ),
    current_side_           ( ifc.current_side_             ),
    innerParent_            ( ifc.innerParent_              ),
    outerParent_            ( ifc.outerParent_              ),
    inner_parent_face_id_   ( ifc.inner_parent_face_id_     ),
    outer_parent_face_id_   ( ifc.outer_parent_face_id_     )
  {
    assert( fptr_ != NULL /* detected unitialized element*/ );
    assert( !interface_connector_.empty() /* detected unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( ifc.LVS() );
  }


/// move constructor
template<size_t dim>
InterFace<dim>::InterFace( InterFace<dim>&& ifc )
  : idx_                    { ifc.idx_ },
    fptr_                   { ifc.fptr_ },
    fvptr_                  { ifc.fvptr_ },
    interface_connector_    { ifc.interface_connector_ },
    parent_elements_node_connector_{ ifc.parent_elements_node_connector_ },
    baseElement_            { ifc.baseElement_  },
    current_side_           { ifc.current_side_ },
    innerParent_            { ifc.innerParent_  },
    outerParent_            { ifc.outerParent_ },
    inner_parent_face_id_   { ifc.inner_parent_face_id_ },
    outer_parent_face_id_   { ifc.outer_parent_face_id_ }
  {
    assert( fptr_ != NULL /* detected unitialized element*/ );
    assert( !interface_connector_.empty() /* detected unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( move(ifc.LVS()) );
    ifc.fptr_  = nullptr;
    ifc.fvptr_ = nullptr;
  }



/// assignment
template<size_t dim>
InterFace<dim>&  InterFace<dim>::operator=( const InterFace<dim>& ifc )
 {
    if ( &ifc != this )
      {
        idx_                    = ifc.idx_;
        fptr_                   = ifc.fptr_;
        fvptr_                  = ifc.fvptr_;
        interface_connector_    = ifc.interface_connector_;
        parent_elements_node_connector_    = ifc.parent_elements_node_connector_;

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
void InterFace<dim>::Assign( const csmp::FiniteVolumeStencil<dim>* const stencil_ptr )
  {
    assert( stencil_ptr != NULL );
    fvptr_ = stencil_ptr;
  }

template<size_t dim>
void InterFace<dim>::Assign( size_t i, InterFace<dim>* const ifc_ptr ) // neighbor interface
  {
    assert( this->FE() != NULL );
    assert( interface_connector_.size() == this->Neighbors() );
    assert( i < this->Neighbors() );
    interface_connector_[i] = ifc_ptr;
  }

template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const base_elmt )
  {
     assert( base_elmt != NULL );
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
     assert( inner_elmt != NULL );
     assert( outer_elmt != NULL );
     innerParent_ = inner_elmt;
     outerParent_ = outer_elmt;
     if( assign_nodes )
         InitializeNodeCorrespondanceVector();
  }


template<size_t dim>
void InterFace<dim>::Assign( Element<dim>* const parentElement, size_t faceId, INTERFACE_SIDE side )
  {
    assert( parentElement != NULL );
    assert( side != MIDDLE );

    if( side == INSIDE )
    {
        innerParent_ = parentElement;
        inner_parent_face_id_ = faceId;
        parent_elements_node_connector_.resize( parentElement->FE()->NodesPerFace( faceId) );
    }
    else if( side == OUTSIDE )
    {
        outerParent_ = parentElement;
        outer_parent_face_id_ = faceId;
        parent_elements_node_connector_.resize( parentElement->FE()->NodesPerFace( faceId) );
    }
    else if( side == MIDDLE )
    {
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




template<size_t dim>
void InterFace<dim>::Assign( size_t n_local, size_t parent_node, INTERFACE_SIDE side )
  {
    assert( n_local < parent_elements_node_connector_.size() );
    assert( side != MIDDLE );

    if( side == INSIDE )
    {
        parent_elements_node_connector_[n_local].first = parent_node;
    }
    else if( side == OUTSIDE )
    {
        parent_elements_node_connector_[n_local].second = parent_node;
    }
    else if( side == MIDDLE )
    {
        throw csmp::Exception( ERROR, "csmp::InterFace<dim>::Assign( local_id, side, parent_node )", "This method is not intended to be used for base element!" );
    }
  }

template<size_t dim>
void InterFace<dim>::Assign( size_t n_local, std::pair<size_t,size_t> parent_nodes )
  {
    assert( n_local < parent_elements_node_connector_.size() );
    parent_elements_node_connector_[n_local] = parent_nodes;
  }

template<size_t dim>
void InterFace<dim>::Assign(FiniteElement * fem_ptr )
 {
    assert( fem_ptr != NULL );
    fptr_ = fem_ptr;
 }

/** InitializeNodeCorrespondanceVector

Initialises vector of Node pointers to the nodes on either
side of the split interface. 

The node pointer vector of the base class is left untouched.
Thus, it continues to refer to the face nodes of the first parent element. 

@attention  The method assumes that the parent elements have already been assigned.
Note that this also implies that Interfaces cannot be located on the outer
boundaries of a model.

Internally, the method also assumes that the numbering of faces corresponds to that
of the neighbors.
*/
template<size_t dim>
void  InterFace<dim>::InitializeNodeCorrespondanceVector()
 {
     assert( innerParent_ != NULL );
     assert( outerParent_ != NULL );
     
     pair<size_t,size_t>  numbers_of_shared_face; 
     
     // finding number of face of parent element 1 which is shared between parent elements.
     bool inner_face_found( false );
     for ( size_t i=0U; i<innerParent_->Neighbors(); ++i )
       if ( innerParent_->Neighbor(i) == outerParent_ ) 
         {
           numbers_of_shared_face.first = i;
           inner_face_found = true;
           break;
         }
         
     // finding number of face of parent element 2 which is shared between parent elements.
     bool outer_face_found( false );
     for ( size_t i=0U; i<outerParent_->Neighbors(); i++ )
       if ( outerParent_->Neighbor(i) == innerParent_ ) 
         {
           numbers_of_shared_face.second = i;
           outer_face_found = true;
           break;
         }

     // shared node ids
     std::vector<size_t>  fnids_inside;
     std::vector<size_t>  fnids_outside;
     if( inner_face_found && outer_face_found )
     {
         inner_parent_face_id_ = numbers_of_shared_face.first;
         outer_parent_face_id_ = numbers_of_shared_face.second;

         // initialising node correspondance vector from the face nodes
         innerParent_->FE()->NodesOfFace( numbers_of_shared_face.first,  fnids_inside  );
         outerParent_->FE()->NodesOfFace( numbers_of_shared_face.second, fnids_outside );
     }
     else
         throw csmp::Exception( ERROR,
                                "csmp::InterFace<dim>::InitializeNodeCorrespondanceVector",
                                "Cannot connect elements! Elements can have different nodes or one of it is of low dimensional type!" );

     // check that the shared nodes are found correctly
     assert( fnids_inside.size() == fnids_outside.size() );

     parent_elements_node_connector_.reserve(fnids_inside.size());
     parent_elements_node_connector_.clear();

     for ( size_t i=0U; i<fnids_inside.size(); i++ )
       {
           const size_t  inside_node = fnids_inside[i];
           // finding the matching outside node through coordinate comparison
           for ( size_t j=0U; j<fnids_outside.size(); j++ )
             if ( outerParent_->N(fnids_outside[j])->Coordinate() == innerParent_->N(inside_node)->Coordinate() )
               {
                  parent_elements_node_connector_.push_back( std::make_pair( inside_node, fnids_outside[j] ) );
                  break;
               }
       }

    // free excessive allocated capacity
    std::vector<std::pair<size_t, size_t> >( parent_elements_node_connector_ ).swap( parent_elements_node_connector_ );

    if( fnids_inside.size() != parent_elements_node_connector_.size() )
        throw csmp::Exception( ERROR, "csmp::InterFace<dim>::InitializeNodeCorrepondanceVector", "Cannot connect elements that has different nodes!" );
     
 } // end InitializeNodeCorrespondanceVector




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
INTERFACE_SIDE  InterFace<dim>::CurrentSide( ) const
 {
    return current_side_;
 }

// methods which use subclasses of FiniteElement class via bridge
template<size_t dim>
FiniteElement*  InterFace<dim>::FE() const
 {
    return fptr_;
 }

template<size_t dim>
const FiniteVolumeStencil<dim>*  InterFace<dim>::FV_Stencil() const
 {
    return fvptr_;
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
    assert( n < parent_elements_node_connector_.size() );

    if ( side == INSIDE )
        return innerParent_->N( parent_elements_node_connector_[n].first );
    else if( side == OUTSIDE )
        return outerParent_->N( parent_elements_node_connector_[n].second );

    if( baseElement_ != NULL )
        return baseElement_->N( n );

    throw csmp::Exception( ERROR, "InterFace<dim>::N( local_id, side )", "Base Element does not exist!" );
    return innerParent_->N( parent_elements_node_connector_[n].first );
 }

template<size_t dim>
csmp::Node<dim>*  InterFace<dim>::N( size_t n ) const
 {
    return N( n, current_side_ );
 }

// watch out if there is no neighbor this returns a NULL pointer
template<size_t dim>
csmp::InterFace<dim>*  InterFace<dim>::Neighbor( size_t n ) const
 {
    assert( interface_connector_.size() == this->Neighbors() );
    assert( n < this->Neighbors() );
    return interface_connector_[n];
 }


template<size_t dim>
size_t  InterFace<dim>::ParentNodeNumber( size_t n, INTERFACE_SIDE side ) const
 {
    assert( n < parent_elements_node_connector_.size() );

    if ( side == INSIDE )
        return parent_elements_node_connector_[n].first;
    else if ( side == OUTSIDE )
        return parent_elements_node_connector_[n].second;
    return n;
 }

template<size_t dim>
size_t  InterFace<dim>::ParentNodeNumberOppositeTo( size_t n_from_parent, INTERFACE_SIDE parent_side ) const
 {
    assert ( parent_side != MIDDLE );

    if ( parent_side == INSIDE )
    {
        assert ( n_from_parent < innerParent_->Nodes() );

        for ( size_t i=0; i < parent_elements_node_connector_.size(); ++i )
            if( parent_elements_node_connector_[i].first == n_from_parent )
            {
                return parent_elements_node_connector_[i].second;
            }
    }
    else if ( parent_side == OUTSIDE )
    {
        assert ( n_from_parent < outerParent_->Nodes() );

        for ( size_t i=0; i < parent_elements_node_connector_.size(); ++i )
            if( parent_elements_node_connector_[i].second == n_from_parent )
            {
                return parent_elements_node_connector_[i].first;
            }
    }

    throw csmp::Exception( ERROR, "InterFace<dim>::ParentNodeNumberOppositeTo( n_from_parent, parent_side )", "Base Element does not exist!" );
    return n_from_parent;
 }


// watch out if there is no base element this returns a NULL pointer
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
Element<dim>*  InterFace<dim>::InnerParent( ) const
 {
    return innerParent_;
 }
template<size_t dim>
Element<dim>*  InterFace<dim>::OuterParent( ) const
 {
    return outerParent_;
 }

template<size_t dim>
Element<dim>*  InterFace<dim>::BaseElement( ) const
 {
    return baseElement_;
 }

template<size_t dim>
bool  InterFace<dim>::hasBase( ) const
 {
    if( baseElement_ != NULL )
        return true;
    return false;
 }

/// return face ID of inner parent element
template<size_t dim>
size_t  InterFace<dim>::InnerParentFaceID() const
 {
     assert( innerParent_ != NULL );
     return inner_parent_face_id_;
}

/// return face ID of outer parent element
template<size_t dim>
size_t  InterFace<dim>::OuterParentFaceID() const
 {
     assert( outerParent_ != NULL );
     return outer_parent_face_id_;
}

/// return face ID of inner or outer parent element depending on index
template<size_t dim>
size_t  InterFace<dim>::ParentFaceID( INTERFACE_SIDE side ) const
{
    if( side == INSIDE )
    {
        assert( innerParent_ != NULL );
        return inner_parent_face_id_;
    }
    else if( side == OUTSIDE )
    {
        assert( outerParent_ != NULL );
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
        assert( innerParent_ != NULL );
        innerParent_->UnitNormalToFace( inner_parent_face_id_, vc );
        return;
    }
    else if( side == OUTSIDE )
    {
        assert( innerParent_ != NULL );
        outerParent_->UnitNormalToFace( outer_parent_face_id_, vc );
        return;
    }

    if( baseElement_ != NULL )
        return baseElement_->UnitNormal( vc );

    throw csmp::Exception( ERROR, "InterFace<dim>::UnitNormal( side )", "Base Element does not exist!" );
    innerParent_->UnitNormalToFace( inner_parent_face_id_, vc );
 }

template<size_t dim>
void  InterFace<dim>::UnitNormal( VectorVariable<dim>& vc ) const
 {
    UnitNormal( vc, current_side_ );
 }

/// Sets VectorVariable to vector between node pair - returns false if overlap, true if distant
template<size_t dim>
bool InterFace<dim>::SpacingNode( size_t n_local, VectorVariable<dim>& innerToOuter ) const
  {
    assert( n_local < parent_elements_node_connector_.size() );
    const Point<dim> inner( this->N(n_local,INSIDE)->Coordinate() ), outer( this->N(n_local,OUTSIDE)->Coordinate() );
    VectorVariable<dim> faceUnitNormal( PLAIN, 0. );
    UnitNormal( faceUnitNormal, INSIDE );
    for( size_t d(0); d < dim; ++d )
        innerToOuter(d) = outer[d]-inner[d];
    if( (innerToOuter&faceUnitNormal) < 0. )
      {
        innerToOuter *= -1.;
        return false;
      }
    return true;
  }


// SCREEN OUTPUT

template<size_t dim>
void  InterFace<dim>::Out(std::ostream& os) const
 {
    os <<"\n\n\nInterFace<"<< dim <<">::Out: number: "<< idx_;

    string str;
    //str = parseBoundary(at_boundary_);
    /// @todo (2-P) Remove typeid, implement name fct
    os <<" ("<< typeid(fptr_).name() <<")";
    //if ( at_boundary_ != NOT ) os <<", Boundary flag: "<< str;
    os << endl;

    os <<"\nInternal data: "<< endl;

    os <<"\n\tconnected nodes with boundary flags:  ";
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i,INSIDE)->AtBoundary());
         os << N(i,INSIDE)->Idx() <<":"<< str <<"  ";
      }
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i,OUTSIDE)->AtBoundary());
         os << N(i,OUTSIDE)->Idx() <<":"<< str <<"  ";
      }
    os << endl;

    os <<"\n\tconnected neighbor InterFace types / boundary flags:\n";
    for ( size_t i=0U; i<this->Neighbors(); i++ )
      if ( Neighbor(i) != NULL ) {
           os <<"\t\t"<< Idx() <<":";
           os << parseFiniteElementType( Neighbor(i)->FE_Type()) <<": ";
           //str = parseBoundary(Neighbor(i)->AtBoundary());
           //os << str;
           os << endl;
        }
      else os <<"none.  ";
    os << endl;

    /// @todo (2-P) Remove typeid, implement name fct
    os <<"\tInterFace is connected via bridge pattern to: ";
    os << typeid(fptr_).name() << endl;

    os <<"\n\tAspect ratio (b-box):   "<< this->AspectRatio() << endl;

    Point<dim>  pt(this->BaryCenter());

    if ( dim == 1U )
       os <<"\n\tBarycentre at (xyz): "<< pt[0] << endl;
    else if ( dim == 2U )
       os <<"\n\tBarycentre at (xyz): "<< pt[0] <<", "<< pt[1] << endl;
    else
       os <<"\n\tBarycentre at (xyz): "<< pt[0] <<", "<< pt[1] <<", "<< pt[2] << endl;

    const size_t ipoints(this->IntegrationPoints());
    if ( ipoints > 0U ) {
         os <<"\n\tStorage sites for IntegrationPoint properties: "<< ipoints << endl;
      }

    os <<"\n Connected Node objects, side 1 of interface: ";
    for ( size_t i=0U; i<this->Nodes(); i++ ) os << parent_elements_node_connector_[i].first <<",  ";
    os << endl;

    os <<"\n Connected Node objects, side 2 of interface: ";
    for ( size_t i=0U; i<this->Nodes(); i++ ) os << parent_elements_node_connector_[i].second <<",  ";
    os << endl;

    /// @todo (2-D) Rm rtti
    os <<"\n\nParent (higher-dimensional) Element objects:     "<< endl;
    if ( innerParent_ != 0 ) {
         os <<"\tinward  facing Element: ";
         os  <<" ("<< typeid(*(this->Parent(INSIDE)->FE())).name() <<")"<< endl;
         this->innerParent_->Out(os);
      }
    else os <<"\tnone.\n";
    if ( this->outerParent_ != 0 ) {
         os <<"\toutward facing Element: ";
         os <<" ("<< typeid(*(this->Parent(OUTSIDE)->FE())).name() <<")"<< endl;
         this->outerParent_->Out(os);
      }
    else os <<"\tnone.\n";

    os <<"\tUnit Normal:            ";
    VectorVariable<dim> un( PLAIN, 0. );
    UnitNormal( un );
    for ( size_t i=0U; i<dim; i++ ) os << un[i] <<", ";
    os << endl;

 } // end Out



template class InterFace<1U>;
template class InterFace<2U>;
template class InterFace<3U>;

} // end namespace csmp
