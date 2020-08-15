#include "Face.h"
#include "Node.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "FiniteElementManager.h"
#include "Exception.h"
#include "Visitor.h"

using namespace std;

namespace csmp {

/// default constructor; @todo make private, decimate overall number of constructors
template<size_t dim>
Face<dim>::Face()
 : idx_(NULL_IDX),
   innerParent_( nullptr ),
   outerParent_( nullptr )
 {
 }


/** 
    Constructs face as exact copy of a lower-dimensional element from which it is constructed.

    custom constructor which:
 
    - builds variable storage
 
    - connects face to its nodes
 
    - connects face to its higher-dimensiona parent elements
    
    @note nodes and node numbering will be exactly that 
    of the original lower-dimensional ELement object.
    It is assumed that this Elemen is numbered correctly in terms of the overall
    mesh when the model imported from a third-party tool.
    
    @note for models imported from professional meshing tools, 
    use available 'surface orientation' checks to ascertain that the 
    nodes of lower-dimensional elements are numbered in a consistent way.
    
    @attention the costly part of this constructor is the determination
    of the indices of the faces of the higher dimensional elements  that match this face
    
    @author SKM 
    @date   1/4/2016

*/
template<size_t dim>
Face<dim>::Face( const Element<dim>& elmt,
                 Element<dim>* const inner_parent,
                 Element<dim>* const outer_parent,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::Face>(elmt.FE()),
    idx_(elmt.Idx()),
    node_connector_(elmt.Nodes(),nullptr),
    face_connector_(elmt.Neighbors(),nullptr),
    innerParent_(inner_parent),
    outerParent_(outer_parent)
 {
    if ( dim == 3 ) assert( elmt.IsSurfaceElement() );
    else if ( dim == 2 ) assert( elmt.IsLineElement() );
    assert( innerParent_ != nullptr );
   
    // 0. verification that the lower-dimensional element and the element that will be transformed
    //    into a face have indeed matching nodes
    const size_t     nodes_to_match(elmt.Nodes());
    set<Node<dim>*>  elmt_nodes;
    for ( size_t j=0U; j<nodes_to_match; ++j )
      elmt_nodes.insert( elmt.N(j) );
    // checking
    size_t matching_nodes(0U);
    for ( size_t j=0U; j<innerParent_->Nodes(); ++j ) {
          assert( innerParent_->N(j) != nullptr );
          if ( elmt_nodes.find( innerParent_->N(j) ) != elmt_nodes.end() )
            matching_nodes++;
      }
    assert( matching_nodes == nodes_to_match );

    // 1. creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // 2. connecting the nodes of the face with those of the lower-dimensional element
    //   from which it was created
    const size_t nodes(elmt.Nodes()); // nodes of Element object that is replicated by Face
    for ( size_t i=0U; i<nodes; ++i ) {
         // assignig the node
         assert( elmt.N(i) != nullptr );
         Assign( i, elmt.N(i) );
    }
   
 } // end constructor





/** 
    constructs model-edge line-element face connected with two volumetric elements at model boundary sharing its nodes
    
    The volumetric elements have to be chosen such that that share an edge with the line element 
    and a boundary face which also shares one of its edges with the line element
    
    @param edge_nodes contains node pointers in the sequence in which they appear on the edge of 
    the inner volumetric element that is adjacent to the Face
    
    @attention if the constructor detects that the higher dimensional parent elements
    have a BOX_BOUNDARY flag that does not indicate a boundary location, they are 
    either flagged according to the edge nodes or as IRREGULAR.
*/
template<size_t dim>
Face<dim>::Face( csmp::FiniteElement* FE_type_of_boundary_face,
                 Element<dim>* const inner_parent,
                 Element<dim>* const outer_parent,
                 const std::vector<Node<dim>*>&  edge_nodes,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
 : FiniteElementPolicy<dim,csmp::Face>(FE_type_of_boundary_face),
   idx_(NULL_IDX),
   innerParent_(inner_parent),
   outerParent_(outer_parent),
   node_connector_(edge_nodes),
   face_connector_(2U,nullptr)
 {
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
   
    // creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
      
 } // end edge constructor




/**
    constructs Face object as an exact match of the target face of the supplied inner parent element
 
    @note boundary element from which face is constructed may have multiple 
    boundary faces, therefore the nth_boundary face variable is required
*/
template<size_t dim>
Face<dim>::Face( Element<dim>& e,
                 csmp::FiniteElement* FE_type_for_face,
                 size_t boundary_face,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::Face>(FE_type_for_face),
    idx_(NULL_IDX),
    node_connector_( FE_type_for_face->Nodes(),nullptr),
    face_connector_(FE_type_for_face->Neighbors(),nullptr),
    innerParent_(&e),
    outerParent_(nullptr)
 {
    assert( boundary_face < e.Faces() );
    if ( e.Neighbor(boundary_face) != nullptr ) outerParent_ = e.Neighbor(boundary_face);
    if ( dim == 2 ) assert( e.IsSurfaceElement() );
    if ( dim == 3 ) assert( e.IsVolumeElement() );

    // 1. creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // 2. assigning nodes to face in the same order as the face nodes
    //    of the inner parent element
    vector<size_t> fnids;
    // nodes of the Element object from wich this Face is constructed
    e.FE()->NodesOfFace( boundary_face, fnids );
    for ( size_t j=0U; j<fnids.size(); ++ j ) {
         assert( e.N( fnids[j] ) != nullptr );
         node_connector_[j] = e.N( fnids[j] );
      }
   
 } // end constructor








/**
    constructor that allocates storage for connections and builds variable storage
    
    @note boundary element from which face is constructed may have multiple 
    boundary faces, therefore the nth_boundary face variable is required
*/
/*
template<size_t dim>
Face<dim>::Face( const FiniteElementManager& finiteElementManager,
                 Element<dim>& e,
                 size_t& nth_boudary_face,
                 const LocalVariables&,
                 const IntegrationPointVariables& )
  : idx_(UINT_MAX),
    fptr_(nullptr), // can only be assigned when face is known
    fvptr_(nullptr),
    innerParent_(&e),
    outerParent_(nullptr)
    // face_connector_(e.Neighbors(),nullptr)
 {
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // finding boundary face
    size_t bface_counter(0U);
    assert( e.Faces() == e.Neighbors() );
    const size_t faces(e.Neighbors());
    for ( size_t i=0U; i<faces; ++i )
       // if the face is located at the model boundary
       if ( e.Neighbor(i) == nullptr )
         {
            if ( bface_counter == nth_boudary_face )
              {
                 // assigning parent FE pointer
                 fptr_ = finiteElementManager.E( e->ElementTypeOfFace(i) );
              
                 // assigning nodes to face in the same order as the face nodes
                 // of the inner parent element
                 vector<size_t> fnids;
                 e.FE()->NodesOfFace( i, fnids );
                 node_connector_.reserve( fnids.size() );
                 for ( size_t j=0U; j<fnids.size(); ++ j ) {
                      assert( e.N( fnids[j] ) != nullptr );
                      node_connector_.push_back( e.N( fnids[j] ) );
                   }
              }
            bface_counter++;
         }

    // reporting how many faces are located at model boundary
    nth_boudary_face = bface_counter;
   
 } // end constructor
*/





/**
    constructor that allocates storage for connections and builds variable storage
    
    @note used in most cases
*/
template<size_t dim>
Face<dim>::Face( csmp::FiniteElement* f,
                 csmp::FiniteVolumeStencil<dim>* fvs,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )

  : FiniteElementPolicy<dim,csmp::Face>(f),
    idx_(UINT_MAX),
    node_connector_(f->Nodes(),nullptr),
    face_connector_(f->Neighbors(),nullptr),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    // LVS must be resized here because the finite element pointer has to be initialised before
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



/**
    construct with storage but without connectivity
    
    @note used to reconstruct model from native binary file
*/
template<size_t dim>
Face<dim>::Face( size_t index,
                 csmp::FiniteElement* f,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )

  : FiniteElementPolicy<dim,csmp::Face>(f),
    idx_(index),
    node_connector_(f->Nodes(),nullptr),
    face_connector_(f->Neighbors(),nullptr),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



/// copy constructor
template<size_t dim>
Face<dim>::Face( const Face<dim>& fc )
  : FiniteElementPolicy<dim,csmp::Face>(fc.FE()),
    FiniteVolumePolicy<dim,csmp::Face>(fc.FV()),
    idx_(fc.idx_),
    face_connector_(fc.face_connector_),
    node_connector_( fc.node_connector_),
    innerParent_(fc.innerParent_),
    outerParent_(fc.outerParent_)
  {
    assert( this->FE() != nullptr /* detected unitialized element*/ );
    assert( !face_connector_.empty() /* detected unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( fc.LVS() );
  }




/// move constructor
template<size_t dim>
Face<dim>::Face( Face<dim>&& fc )
  : FiniteElementPolicy<dim,csmp::Face>(move(fc.FE())),
    FiniteVolumePolicy<dim,csmp::Face>(move(fc.FV())),
    idx_(move(fc.idx_)),
    node_connector_(move(fc.node_connector_)),
    face_connector_(move(fc.face_connector_)),
    innerParent_(move(fc.innerParent_)),
    outerParent_(move(fc.outerParent_))
 {
    this->LVS( move(fc.LVS()) );

    fc.AssignFiniteElementNullPtr();
    fc.AssignFiniteVolumeNullPtr();

    fc.innerParent_ = nullptr;
    fc.outerParent_ = nullptr;
 }


// TODO: implement move constructor and assignment operator


template<size_t dim>
Face<dim>::~Face()
 {
 }




/**
    private as there should never be a need to use it
    and the results are undefined in terms of the pointers to parent and neighbor elements
    
    @note assignment operator is used in container::find(Face) operations; why?
*/
template<size_t dim>
Face<dim>&  Face<dim>::operator=( const Face<dim>& fc )
 {
    if ( &fc != this ) {
		if ( fc.FE() ) FiniteElementPolicy<dim,csmp::Face>::Assign(fc.FE());
		if ( fc.FV() ) FiniteVolumePolicy<dim,csmp::Face>::AssignFiniteVolume(fc.FV());
        idx_                  = fc.idx_;
        face_connector_       = fc.face_connector_;
        node_connector_       = fc.node_connector_;
        innerParent_          = fc.innerParent_; // problematic pointer assignment
        outerParent_          = fc.outerParent_; // problematic
        this->LVS( fc.LVS() );
      }
    return *this;
 }




/**
    @note a temporary value cannot be equivalent to lvalue.
*/
template<size_t dim>
Face<dim>&  Face<dim>::operator=( Face<dim>&& fc )
 {
    assert( &fc != this );

    if ( fc.FE() ) FiniteElementPolicy<dim,csmp::Face>::Assign(move(fc.FE()));
    if ( fc.FV() ) FiniteVolumePolicy<dim,csmp::Face>::AssignFiniteVolume(move(fc.FV()));
   
    idx_             = move(fc.idx_ );
    face_connector_  = move(fc.face_connector_);
    node_connector_  = move(fc.node_connector_);
    innerParent_     = move(fc.innerParent_);
    outerParent_     = move(fc.outerParent_);

    this->LVS( move(fc.LVS()) );

    fc.AssignFiniteElementNullPtr();
    fc.AssignFiniteVolumeNullPtr();

    fc.innerParent_ = nullptr;
    fc.outerParent_ = nullptr;

    return *this;
 }




// CONSTRUCTION PROCESS

/**
    Connect Face object with its nodes.
*/
template<size_t dim>
void Face<dim>::Assign( size_t i, csmp::Node<dim>* const nd_ptr )
 {
    assert( nd_ptr != nullptr );
    assert( i < Nodes()*2U );
    assert( !node_connector_.empty() );

    node_connector_[i] = nd_ptr;
 }



template<size_t dim>
void Face<dim>::Assign( size_t i, csmp::Face<dim>* const f_ptr )
 {
    assert( f_ptr != nullptr );
    assert( i < Neighbors() );
    assert( face_connector_.size() == this->FE()->Neighbors() );

    face_connector_[i] = f_ptr;
 }


/**
    unassigns the neighbor face, setting the pointer in the 'face_connector' vector to null
*/
template<size_t dim>
void Face<dim>::Unassign( const Face<dim>* e_ptr ) 
  {
    for ( size_t i = 0U; i < face_connector_.size(); i++ ) {
        if ( e_ptr == nullptr || face_connector_[i] == nullptr )
          continue;
        if ( e_ptr == face_connector_[i] ) {
            face_connector_[i] = nullptr;
            break;
          }
      }
  }


/*
template<size_t dim>
bool Face<dim>::DisconnectNeighbor(Face<dim>* f_ptr)
{
	assert(f_ptr != nullptr);	
	assert(face_connector_.size() == this->FE()->Neighbors());
	for (size_t i(0); i < face_connector_.size(); ++i)
		if (f_ptr == face_connector_[i])
		{
			face_connector_.erase(face_connector_.begin() + i);
			face_connector_.swap(face_connector_);
			return true;
		}
	return false;
} // end DisconnectNeighbor
*/


template<size_t dim>
size_t  Face<dim>::Nodes() const
{
	return node_connector_.size();
}

template<size_t dim>
size_t  Face<dim>::Neighbors() const
{
	return face_connector_.size();
}


template<size_t dim>
size_t  Face<dim>::ConnectedNeighbors() const
{
	size_t nulls(0);
	for (auto f : face_connector_)
		if (!f) nulls++;
	return (face_connector_.size() - nulls);
}


template<size_t dim>
size_t  Face<dim>::Faces() const
{
	return face_connector_.size();
}


/**
    Connects face to the higher-dimensional parent elements which it may be sandwiched between.
    The node pointers are reassigned if they were already pointing to other nodes.
    
    @note method assumes that the Face object already has valid node pointers
    
    @attention the nodes must be assigned to the Face before this method is called
 
    @attention The node numbering of the Face determines the direction of its normal.
    Here the convention is assumed that the first parent element is that on the inside
    of the Face with regard to the outward pointing normal and the second element is on
    the outside. It follows that the node sequences are the same.
    
    TODO: @todo (1) store the number of face or edge of the higher-dimensional inner
    element for future reference; return it from the function InnerParentFaceNumber().
    
    TODO: @todo (1) check code-coverage and remove all constructors and functions 
    that are not used BEFORE working up the InterFace class.
*/
template<size_t dim>
void Face<dim>::Assign( Element<dim>* const innerElement, Element<dim>* const outerElement )
  {
    // 1. argument checks and assignments
    // ----------------------------------
    assert( innerElement != nullptr ); // inner element must be defined
    if      ( dim == 3U ) assert( innerElement->IsVolumeElement() );
    else if ( dim == 2U ) assert( innerElement->IsSurfaceElement() );
    innerParent_ = innerElement;
    
    if ( outerElement != nullptr ) { // outer element is defined if face is in model interior
         if      ( dim == 3U ) assert( outerElement->IsVolumeElement() );
         else if ( dim == 2U ) assert( outerElement->IsSurfaceElement() );
         outerParent_ = outerElement;
      }
    
    // 2. finding the face of the inner element which corresponds to this Face and checking its node numbering
    // -------------------------------------------------------------------------------------------------------
    // creating a search key for Face
    set<Node<dim>*>  face_nds;
    const size_t face_nodes(Nodes());
    for ( size_t i=0U; i<face_nodes; ++i ) {
         // are the nodes there?
         assert( this->N(i) != nullptr );
         face_nds.insert( this->N(i) );
      }

    // 3. matching the lower-dimensional face to a face of the higher-dimensional element
    // -----------------------------------------------------------------------------------
    if ( ((dim == 3U) && this->IsSurfaceElement() ) || // Face is either a triangle or a quadrilateral in 3D
         ((dim == 2U) && this->IsLineElement()) )      // Face is a line in 2D
      {
         // searching the matching Face of the inner parent element
         bool             matching_face_found(false);
         const size_t     faces(innerParent_->Faces());
         vector<size_t>   nodes_of_face;
         set<Node<dim>*>  parent_nds;
        
         for ( size_t i=0U; i<faces; ++i ) {
              innerParent_->FE()->NodesOfFace( i, nodes_of_face );
              for ( size_t j=0U; j<nodes_of_face.size(); ++j )
                parent_nds.insert( innerParent_->N( nodes_of_face[j] ) );
              // checking
              if ( face_nds == parent_nds ) {
                    matching_face_found = true;
                    break;
                }
              parent_nds.clear();
           }         
         return;
      }
    
    // 4. matching the dimension-2 face with an edge of the higher-dimensional element
    // -------------------------------------------------------------------------------
    // (if the higher-dimensional element has 2 more dimensions that the Face, a comparison
    //  with its edges needs to be performed)
    // searching for the matching Segment (Edge) of the inner parent element
    bool             matching_segment_found(false);
    const size_t     segments(innerParent_->Segments());
    vector<size_t>   nodes_of_segm;
    set<Node<dim>*>  parent_nds;
    
     for ( size_t i=0U; i<segments; ++i ) {
          innerParent_->FE()->NodesOfSegment( i, nodes_of_segm );
          for ( size_t j=0U; j<nodes_of_segm.size(); ++j )
            parent_nds.insert( innerParent_->N( nodes_of_segm[j] ) );
          // checking wether line element matches edge dim+2 element
          if ( face_nds == parent_nds ) {
                matching_segment_found = true;
                break;
            }
          parent_nds.clear();
       }

 } // end assign





/// allows to assign a single parent, providing the node ids of face instead of face id. all assigns are dispatched to here.
/*
template<size_t dim>
void  Face<dim>::Assign( Element<dim>* const parentElement, const vector<Node<dim>*>& faceNodes )
  {
    assert( parentElement != NULL );

    // assigning parent ptr
    innerParent_ = parentElement;

    AssignFaceID( faceNodes );

    // SKM FIX - the normals are fine, this is no longer needed
    // CheckNodeOrderingAccordingToUnitNormalOrientation();
  }

SKM NOT SUPPORTED:  NO USE AT THE MOMENT

template<size_t dim>
void Face<dim>::Assign(FiniteElement * fem_ptr )
 {
    assert( fem_ptr != NULL );
    fptr_ = fem_ptr;
 }
*/



/* SKM NOT SUPPORTED:  NO USE AT THE MOMENT

// checking for proper node numbering (unit normal needs to point outward, to outer parent that is)
// both methods works
template<size_t dim>
void  Face<dim>::CheckNodeOrderingAccordingToUnitNormalOrientation()
  {
    //if( !this->isCorrectUnitNormalOrientation( innerParent_ ) )
    if( !this->isCorrectUnitNormalOrientation( innerParent_, inner_parent_face_id_ ) )
    {
        if( outerParent_ != NULL )
        {
            // potential infinite loop in case if orientaion cannot be established
            //this->Flip();
            this->RevertNodeNumbering();
        }else{
            //throw csmp::Exception( ERROR, "Face<dim>::CheckNodeOrderingAccordingToUnitNormalOrientation()", "Face nodes are not ordered correctly!." );
            //std::cerr<<"Face<dim>::CheckNodeOrderingAccordingToUnitNormalOrientation(): Face nodes are not ordered correctly!"<<std::endl;
            this->RevertNodeNumbering();
        }
    }
    this->FE()->CurrentID( FiniteElement::InitialID() );
  }
*/



// VISITOR

/// SKM temporary patch to get Visitor objects to work on boundaries
template<size_t dim>
void Face<dim>::Accept( csmp::Visitor<dim>& vis )
 {
    if ( vis.ApplicationTarget() == FACE ) {
         vis.Visit(this);
         return;
      }
    if ( vis.ApplicationTarget() == NODE ) {
         size_t n_nodes(this->Nodes());
         for ( size_t i = 0U; i< n_nodes; i++ )
             this->N(i)->Accept( vis );
         return;
      }
    throw logic_error("Face<dim>::Accept: target of visitation unresolved.");

 } // end Accept




  

// ACCESSORS

/**
     Reports the local number (0..faces-1) of the face of the inner higher-dimensional parent element that
     matches the nodes of the Face object (also in terms of the sequence of these nodes).
     
     @return local face number of UINT_MAX if no index could be found.

     @attention assumes that the Face has valid nodes and its inner parent element is connected
*/
template<size_t dim>
size_t Face<dim>::InnerParentFaceNumber() const
 {
    assert( !node_connector_.empty() );
   
     // 1. creating a unique key from the nodes of the Face
    set<size_t>  face_key;
    for ( auto nit=node_connector_.begin(); nit!=node_connector_.end(); ++nit ) {
         assert( (*nit) != nullptr );
         face_key.insert( (*nit)->Idx() );
      }
    
    // 2. creating face keys for the inner parent element and trying to match them
    //    with the one created for the current face
    vector<size_t> fnids;
    set<size_t>    face_key_n;
    const size_t faces(innerParent_->Faces());
    for ( size_t i=0U; i<faces; ++i ) {
         innerParent_->FE()->NodesOfFace( i, fnids );
         for ( size_t j=0U; j<fnids.size(); ++j )
           face_key_n.insert( innerParent_->N( fnids[j] )->Idx() );
         // has a matching face been found?
         if ( face_key == face_key_n ) return i;
         else face_key_n.clear();
      }

    return UINT_MAX;
 }





/// returns the current index of this face assuming that a meaningful value was assigned earlier
template<size_t dim>
void  Face<dim>::Idx( size_t idx_to_assign ) const
  {
    idx_ = idx_to_assign;
  }

template<size_t dim>
size_t   Face<dim>::Idx() const
  {
    return idx_;
  }


template<size_t dim>
typename std::vector<csmp::Face<dim>*>&  Face<dim>::NeighborElementVector()
  {
    return face_connector_;
  }


template<size_t dim>
csmp::Node<dim>*  Face<dim>::N( size_t n ) const
  {
     assert( n < Nodes() );
     return node_connector_[n];
  }


/**
    watch out if there is no neighbor this returns a NULL pointer
*/
template<size_t dim>
csmp::Face<dim>*  Face<dim>::Neighbor( size_t n ) const
 {
    assert( n < Neighbors() );
    return face_connector_[n];
 }


/**
    Returns either the INSIDE or OUTSIDE elements
    connected to this face if they exist.
*/
template<size_t dim>
Element<dim>*  Face<dim>::Parent( INTERFACE_SIDE side ) const
 {
    assert( side != MIDDLE );
    assert( side != 0 );
    assert( innerParent_ != nullptr );
    return ( side == INSIDE ) ? innerParent_ : outerParent_;
 }


/**
    inner parent is always initialised for a valid face
*/
template<size_t dim>
Element<dim>*  Face<dim>::InnerParent() const
 {
    assert( innerParent_ != nullptr );
    return innerParent_;
 }


/**
    @attention outer parent will not be initialised if the Face lies on the outside boundary of the model
*/
template<size_t dim>
Element<dim>*  Face<dim>::OuterParent() const
 {
    return outerParent_;
 }

/*  RECREATE IF NEEDED USING THE NODE TO PARENT ELEMENT CAPABILITY
template<size_t dim>
size_t  Face<dim>::ParentNodeNumber( size_t n ) const
 {
    assert( n < parent_element_node_ids_.size() );
    return parent_element_node_ids_[n];
 }
*/

// GEOMETRY

/**
     Returns a vector (not a unit normal) connecting the barycenter of the Face
     with that of the higher-dimensional parent element on its inside.
*/
template<size_t dim>
void Face<dim>::VectorToInnerElementBaryCenter( VectorVariable<dim>& vectorToInner ) const
{
   assert( Parent(INSIDE) != nullptr );
   vectorToInner = Parent(INSIDE)->BaryCenter() - this->BaryCenter();
}


/// returns area (3D) or length (2D) of face
template<size_t dim>
double64 Face<dim>::Area() const
 {
    return this->Volume();
 }




// FUNCTIONALITY

/**

Returns a matrix with 'nodes'-rows and 'coordinate-directions' columns.
This matrix defines the positions of the elements nodes for
the finite-element matrix assembly. Since the number of element nodes
may vary among different elements types, the number of rows in XY may
also vary from element to element.

@param XY A DenseMatrix<DM_MIN> class object (value type fT). This matrix is dynamically
resized if necessary but must have been constructed with a finite size
before passing it to CoordinateMatrix().

@return The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<size_t dim>
void  Face<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
  {
    const size_t n_nodes( Nodes());
    XY.Resize( n_nodes, dim );
    for ( size_t i=0U; i<n_nodes; ++i )
        XY.AssignRow( i, N(i)->Coordinate() );

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
Point<dim>  Face<dim>::BaryCenter() const
  {
    Point<dim>    pt(N(0U)->Coordinate());
    const size_t  n_nodes(Nodes());
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
double64  Face<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
  {
    double64 fMinTemp( static_cast<double64>( DBL_MAX) );
    double64 fMaxTemp( static_cast<double64>(-DBL_MAX) );

    // this normalisation is necessary because the vector variable
    // being any physical quantity may have any magnitude
    const double64 fMagnitudeOfDirection(vecDirection.Length());
    // avoid division by zero
    assert( fMagnitudeOfDirection >= numeric_limits<double64>::epsilon() );

    const size_t n_nodes(Nodes());
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
void  Face<dim>::NodePropertyVector( const csmp::Index& idx, std::vector<Var>& V ) const
 {
    if ( idx.place != NODE ) {
         std::cerr <<"\nFace<"<< dim;
         std::cerr <<">::NodePropertyVector: Requested property ";
         std::cerr <<"is not placed on the nodes; property Index: "<< std::endl;
         idx.Out();
         return;
      }

    // resizing V if necessary
    const size_t  n_nodes(Nodes());
    V.resize(n_nodes);

    for ( size_t i=0U; i<n_nodes; i++ )
      N(i)->Read( idx, V[i] );
 }

// scalar
template void  Face<1U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
template void  Face<2U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
template void  Face<3U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
// vector
template void  Face<1U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<1U> >& ) const;
template void  Face<2U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<2U> >& ) const;
template void  Face<3U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<3U> >& ) const;
// tensor
template void  Face<1U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<1U> >& ) const;
template void  Face<2U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<2U> >& ) const;
template void  Face<3U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<3U> >& ) const;
// array
template void  Face<1U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
template void  Face<2U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
template void  Face<3U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
// flagged array
template void  Face<1U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;
template void  Face<2U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;
template void  Face<3U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;





// OUTPUT

///  outputs local variables of Face(Element) overriding corresponding method of base class
template<size_t dim>
void  Face<dim>::Out() const
 {
    cout <<"\n\n\nFace<"<< dim <<">::Out: number: "<< idx_;

    string str;
    cout <<" ("<< parseFiniteElementType(this->FE_Type()) <<")";
    if ( outerParent_ == nullptr ) cout <<", face is located at model boundary.";
    cout << endl;

    cout <<"\nInternal data: "<< endl;

    cout <<"\n\tconnected nodes with boundary flags:  ";
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i)->AtBoundary());
         cout << N(i)->Idx() <<":"<< str <<"  ";
      }
    cout << endl;

    cout <<"\n\tconnected neighbor Face types / boundary flags:\n";
    for ( size_t i=0U; i<this->Neighbors(); i++ )
      if ( Neighbor(i) != NULL ) {
           cout <<"\t\t"<< Idx() <<":";
           cout << parseFiniteElementType(Neighbor(i)->FE_Type()) <<": ";
           cout << endl;
        }
      else cout <<"none.  ";
    cout << endl;

    cout <<"\tFace is connected via bridge pattern to: ";
    cout << parseFiniteElementType(this->FE_Type()) << endl;

    cout <<"\n\tAspect ratio (b-box):   "<< this->AspectRatio() << endl;

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

    cout <<"\nParent (higher-dimensional) Element objects:\n";
    if ( innerParent_ != nullptr ) {
         cout <<"\tinward  facing Element: "<< this->innerParent_->Idx();
         cout  <<" ("<< parseFiniteElementType(this->Parent(INSIDE)->FE_Type()) <<")"<< endl;
      }
    else cout <<"\tnone.\n";
    if ( this->outerParent_ != nullptr ) {
         cout <<"\toutward facing Element: "<< this->outerParent_->Idx();
         cout <<" ("<< parseFiniteElementType(this->Parent(OUTSIDE)->FE_Type()) <<")"<< endl;
      }
    else cout <<"\tnone.\n";

 } // end Out


// ==========================================================================================================
// TO DEPRECATE
// ==========================================================================================================


/// Roman, 2014
/// WARNING: this operator is used specifically in the process of creation of particular Boundary.
/// Therefore only important infromation for that process is taken into account in order to distinguish two Face's.
/// That must be reference to inner and outer parent Elements and inner parent face ID
template<size_t dim>
bool  Face<dim>::operator==( const Face<dim>& fc ) const
 {
	 if (&fc != this) {
		 if (node_connector_.size() != fc.node_connector_.size()) return false;
		 if (innerParent_ != fc.innerParent_ || outerParent_ != fc.outerParent_) return false;
		 for (size_t i = 0U; i < node_connector_.size(); i++)
			 if (node_connector_[i]->Idx() != fc.node_connector_[i]->Idx()) return false;
	 }
     return true;
 }





// explicits
template class Face<1U>;
template class Face<2U>;
template class Face<3U>;


} // end namespace csmp
