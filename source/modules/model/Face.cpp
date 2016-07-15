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
   fptr_(nullptr),
   fvptr_(nullptr),
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
  : idx_(elmt.Idx()),
    fptr_(elmt.FE()),
    fvptr_(elmt.FV_Stencil()),
    node_connector_(elmt.Nodes(),nullptr),
    face_connector_(elmt.Neighbors(),nullptr),
    innerParent_(inner_parent),
    outerParent_(outer_parent)
 {
    if ( dim == 3 ) assert( elmt.IsSurfaceElement() );
    else if ( dim == 2 ) assert( elmt.IsLineElement() );
    assert( innerParent_ != nullptr );
   
    // 1. creating local storage for face and face integration point variables
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // 2. connecting the nodes of the face with those of the lower-dimensional element
    //   from which it was created
    const size_t nodes(elmt.Nodes());
    const size_t inner_elmt_nodes(innerParent_->Nodes());
    for ( size_t i=0U; i<nodes; ++i ) {
         // assignig the node
         assert( elmt.N(i) != nullptr );
         Assign( i, elmt.N(i) );
    
#ifndef NDEBUG // the node match with the higher-dimensional parent element is verified
         bool matching_node_exists(false);
         for ( size_t j=0U; j<inner_elmt_nodes; ++j ) {
              assert( innerParent_->N(j) != nullptr );
              // assigning parent element node numbers when the pointers match
              if ( elmt.N(i) == innerParent_->N(j) ) {
                  matching_node_exists = true;
                  break;
                }
           }
         assert( matching_node_exists );
#endif
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
 : idx_(NULL_IDX),
   fptr_(FE_type_of_boundary_face),
   fvptr_(nullptr),
   innerParent_(inner_parent),
   outerParent_(outer_parent),
   node_connector_(edge_nodes),
   face_connector_(2U,nullptr)
 {
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
    assert( fptr_ != nullptr );
   
    // creating local storage for face and face integration point variables
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
   
    // checking and adjusting the BOX_BOUNDARY element flagging of
    // the higher-dimensional neighbor elements if these are not already
    // flagged as boundary
    if ( innerParent_->AtBoundary() == NOT )
      for ( const auto nit : node_connector_ )
        if ( isEdge( nit->AtBoundary() ) ) {
             innerParent_->AtBoundary( nit->AtBoundary() );
             break;
          }
    // if no suitable edge was determined the element is nevertheless marked as a boundary one
    if ( innerParent_->AtBoundary() == NOT )
      innerParent_->AtBoundary( IRREGULAR );
    
    if ( outerParent_->AtBoundary() == NOT )
      for ( const auto nit : node_connector_ )
        if ( isEdge( nit->AtBoundary() ) ) {
             outerParent_->AtBoundary( nit->AtBoundary() );
             break;
          }
    if ( outerParent_->AtBoundary() == NOT )
      outerParent_->AtBoundary( IRREGULAR );
   
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
  : idx_(NULL_IDX),
    fptr_(FE_type_for_face),
    fvptr_(nullptr),
    node_connector_( FE_type_for_face->Nodes(),nullptr),
    face_connector_(FE_type_for_face->Neighbors(),nullptr),
    innerParent_(&e),
    outerParent_(nullptr)
 {
    assert( boundary_face < e.Faces() );
    assert( e.Neighbor(boundary_face) == nullptr );

    // 1. creating local storage for face and face integration point variables
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // 2. assigning nodes to face in the same order as the face nodes
    //    of the inner parent element
    vector<size_t> fnids;
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
                 fptr_ = finiteElementManager.E( e.FE()->ElementTypeOfFace(i) );
              
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

  : idx_(UINT_MAX),
    fptr_(f),
    fvptr_(fvs),
    node_connector_(f->Nodes(),nullptr),
    face_connector_(f->Neighbors(),nullptr),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    // LVS must be resized here because the finite element pointer has to be initialised before
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
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

  : idx_(index),
    fptr_(f),
    fvptr_(nullptr),
    node_connector_(f->Nodes(),nullptr),
    face_connector_(f->Neighbors(),nullptr),
    innerParent_( nullptr ),
    outerParent_( nullptr )
 {
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
 }



/// copy constructor
template<size_t dim>
Face<dim>::Face( const Face<dim>& fc )
  : idx_                  ( fc.idx_                  ),
    fptr_                 ( fc.fptr_                 ),
    fvptr_                ( fc.fvptr_                ),
    face_connector_       ( fc.face_connector_       ),
    node_connector_ ( fc.node_connector_ ),
    innerParent_          ( fc.innerParent_          ),
    outerParent_          ( fc.outerParent_          )
  {
    assert( fptr_ != nullptr /* detected unitialized element*/ );
    assert( !face_connector_.empty() /* detected unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( fc.LVS() );
  }


/// move constructor
template<size_t dim>
Face<dim>::Face( Face<dim>&& fc )
  : idx_                  { fc.idx_},
    fptr_                 { fc.fptr_},
    fvptr_                { fc.fvptr_},
    node_connector_{ fc.node_connector_},
    face_connector_       { fc.face_connector_},
    innerParent_          { fc.innerParent_},
    outerParent_          { fc.outerParent_}
  {
    this->LVS( move(fc.LVS()) );
    fc.fptr_  = nullptr;
    fc.fvptr_ = nullptr;
  }



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
        idx_                  = fc.idx_;
        fptr_                 = fc.fptr_;  // O.K.
        fvptr_                = fc.fvptr_; // O.K.
        face_connector_       = fc.face_connector_;
        node_connector_       = fc.node_connector_;
        innerParent_          = fc.innerParent_; // problematic pointer assignment
        outerParent_          = fc.outerParent_; // problematic
        this->LVS( fc.LVS() );
      }
    return *this;
 }


// CONSTRUCTION PROCESS


/**
    Connect Face object with its nodes.
*/
template<size_t dim>
void Face<dim>::Assign( size_t i, csmp::Node<dim>* const nd_ptr )
 {
    assert( i < this->Nodes() );
    assert( nd_ptr != NULL );
    assert( fptr_ != NULL );
    assert( node_connector_.size() == this->Nodes() );

    node_connector_[i] = nd_ptr;
 }



/// connect face to corresponding finite-volume stecncil
template<size_t dim>
void Face<dim>::Assign( const csmp::FiniteVolumeStencil<dim>* const stencil_ptr )
 {
    assert( stencil_ptr != NULL );
    fvptr_ = stencil_ptr;
 }



/// connect face to its neighbors
template<size_t dim>
void Face<dim>::Assign( size_t i, Face<dim>* const fc_ptr ) // neighbor face
 {
    assert( this->FE() != NULL );
    assert( face_connector_.size() == this->Neighbors() );
    assert( i < this->Neighbors() );
    face_connector_[i] = fc_ptr;
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
    assert( fptr_ != nullptr );
    set<Node<dim>*>  face_nds;
    const size_t face_nodes(fptr_->Nodes());
    for ( size_t i=0U; i<face_nodes; ++i ) {
         // are the nodes there?
         assert( this->N(i) != nullptr );
         face_nds.insert( this->N(i) );
      }

    // searching for the Face in inner parent element
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
    assert( matching_face_found );

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
    assert( fptr_ != nullptr );
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

// methods which use subclasses of FiniteElement class via bridge
template<size_t dim>
FiniteElement*  Face<dim>::FE() const
  {
    return fptr_;
  }

template<size_t dim>
const FiniteVolumeStencil<dim>*  Face<dim>::FV_Stencil() const
  {
    return fvptr_;
  }

template<size_t dim>
typename std::vector<csmp::Face<dim>*>&  Face<dim>::NeighborElementVector()
  {
    return face_connector_;
  }


template<size_t dim>
csmp::Node<dim>*  Face<dim>::N( size_t n ) const
  {
     assert( node_connector_.size() == fptr_->Nodes() );
     assert( n < this->Nodes() );
     return node_connector_[n];
  }


/**
    watch out if there is no neighbor this returns a NULL pointer
*/
template<size_t dim>
csmp::Face<dim>*  Face<dim>::Neighbor( size_t n ) const
 {
    assert( face_connector_.size() == fptr_->Neighbors() );
    assert( n < this->Neighbors() );
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
    cout << parseFiniteElementType(fptr_->ElementType()) << endl;

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
    if ( &fc != this )
      if ( innerParent_ != fc.innerParent_ || outerParent_ != fc.outerParent_ ) return false;
    return true;
 }





// explicits
template class Face<1U>;
template class Face<2U>;
template class Face<3U>;


} // end namespace csmp
