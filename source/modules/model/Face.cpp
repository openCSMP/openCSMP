#include "Face.h"
#include "Node.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "Exception.h"
#include "Visitor.h"

using namespace std;

namespace csmp {

/** 
    Constructs face as exact copy of a lower-dimensional element from which it is constructed.
    The lower-dimensional element gets attached to the middle element pointer.

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
    
    @note used by MeshManager::ReplaceElementByFace

    @author SKM
    @date   1/4/2016
    
*/
template<uint32_t dim>
Face<dim>::Face( const Element<dim>& elmt,
                 Element<dim>* const inner_parent,
                 Element<dim>* const outer_parent,
                 uint32_t inner_parent_face_id,
                 uint32_t outer_parent_face_id,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::Face>(elmt.FE()),
    FiniteVolumePolicy<dim,csmp::Face>(elmt.FV()),
    idx_(elmt.Idx()),
    node_connector_(elmt.Nodes(),nullptr),
    face_connector_(elmt.Neighbors(),nullptr),
    innerParent_(inner_parent),
    outerParent_(outer_parent),
    inner_parent_face_id_(inner_parent_face_id),
    outer_parent_face_id_(outer_parent_face_id)
 {
    if constexpr ( dim == 3 ) assert( elmt.IsSurface() );
    if constexpr ( dim == 2 ) assert( elmt.IsLine() );
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
    assert( innerParent_->Neighbor(inner_parent_face_id_) == outerParent_ );
    assert( outerParent_->Neighbor(outer_parent_face_id_) == innerParent_ );
   
    // 0. verification that the lower-dimensional element and the element that will be transformed
    //    into a face have indeed matching nodes
#ifdef DEBUG
    const uint32_t         nodes_to_match(elmt.Nodes());
    set<const Node<dim>*>  elmt_nodes;
    for ( auto j{0U}; j<nodes_to_match; ++j )
      elmt_nodes.insert( elmt.N(j) );
    // checking inner parent
    uint32_t matching_nodes(0U);
    for ( auto j{0U}; j<innerParent_->Nodes(); ++j ) {
          assert( innerParent_->N(j) != nullptr );
          if ( elmt_nodes.find( innerParent_->N(j) ) != elmt_nodes.end() )
            matching_nodes++;
      }
    assert( matching_nodes == nodes_to_match );
    // checking outer parent
    if ( outerParent_ != nullptr ) {
        matching_nodes = 0U;
        for ( auto j{0U}; j<innerParent_->Nodes(); ++j ) {
              assert( innerParent_->N(j) != nullptr );
              if ( elmt_nodes.find( innerParent_->N(j) ) != elmt_nodes.end() )
                matching_nodes++;
          }
        assert( matching_nodes == nodes_to_match );
     }
#endif

    // 1. creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
      this->ResizePropertyStorage( ep, ip );
    else
      this->ResizePropertyStorage( ep );


    // 2. Assigning the nodes from that of the Inner Parents face.  Attention! Not the lower-dim element
    vector<uint32_t> fnids = inner_parent->FE()->NodesOfFace( inner_parent_face_id_ );
    uint32_t i_node{0};
    for ( auto i : fnids )
      node_connector_[i_node++] = inner_parent->N(i);

    /* Old Implementation - Doesnt guarantee nodes of face = InnerParent->NodesOfFace() are the same (e.g could be rotated by one node).
    // 2. connecting the nodes of the face with those of the lower-dimensional element
    //   from which it was created
    const auto nodes(elmt.Nodes()); // nodes of Element object that is replicated by Face
    for ( auto i{0U}; i<nodes; ++i ) {
         // assignig the node
         assert( elmt.N(i) != nullptr );
         Assign( i, elmt.N(i) );
      }
      */
   
 } // end constructor





/**
      Constructs face on the inside of the model, auto detecting their common face
                                     ^^^^^^^^^^^^^
      and shared nodes. Only a pointer to the finite element matching the face has to be supplied.
      
    @note the node numbering of the inner parent element will be adopted for the new Face.
    
    @author SKM
    @date 22/10.2021
*/
template<uint32_t dim>
Face<dim>::Face( const FiniteElementManager& fem_manager,
                 const FiniteVolumeStencilManager<dim>* fvm_manager,
                 Element<dim>* const inner_parent,
                 Element<dim>* const outer_parent,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : idx_(NULL_IDX),
    innerParent_(inner_parent),
    outerParent_(outer_parent),
    inner_parent_face_id_(NULL_IDX),
    outer_parent_face_id_(NULL_IDX)
 {
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
    if ( innerParent_ == outerParent_ ) {
         cerr <<"\n\nFace<"<< dim <<">(ctor: of face between parents): supplied pointers point to the same element: ";
         cerr << inner_parent->Idx();
         inner_parent->Out();
         innerParent_ = nullptr;
         outerParent_ = nullptr;
         assert( innerParent_ != outerParent_ );
         return;
      }
    
    // finding the face which is shared
    // --------------------------------
    bool  matching_face_found{false};
    const size_t n_faces_inner{innerParent_->Faces()};
    for ( auto i{0U}; i<n_faces_inner; ++i )
      // if the faces match
      if ( inner_parent->Neighbor(i) == outer_parent ) {
           inner_parent_face_id_ = i;
           // assigning the finite element
           const CSMP_FEM_TYPE etype = outer_parent->FE()->ElementTypeOfFace(i);
           FiniteElementPolicy<dim,csmp::Face>::Assign( fem_manager.E(etype) );
           if ( fvm_manager ) FiniteVolumePolicy<dim,csmp::Face>::AssignFiniteVolume( fvm_manager->Stencil(etype) );
           node_connector_.resize(this->FE()->Nodes(),nullptr);
           face_connector_.resize(this->FE()->Faces(),nullptr);
           // assigning the nodes
           uint32_t n_count{0U};
           for ( const auto& k : inner_parent->FE()->NodesOfFace(i) )
             node_connector_[n_count++] = inner_parent->N(k);
           // finding the number of the shared face in the outer element
           const auto n_faces_outer{outerParent_->Faces()};
           for ( auto j{0U}; j<n_faces_outer; ++j )
             if ( outer_parent->Neighbor(j) == inner_parent ) {
                  outer_parent_face_id_ = j;
                  break;
               }
           matching_face_found = true;
           break;
        }
 
    // reporting the failed construction
    if ( !matching_face_found ) {
         cerr <<"\n\nFace<"<< dim <<">(ctor: face between parents): parent elements "<< inner_parent->Idx();
         cerr <<" and "<< outer_parent->Idx() <<" do not seem to share a face:\n";
         inner_parent->Out();
         outer_parent->Out();
         innerParent_ = nullptr;
         outerParent_ = nullptr;
         return;
      }

    // pointers should point to the adjacent element
    assert( innerParent_->Neighbor(inner_parent_face_id_) == outerParent_ );
    if ( outerParent_ )
      assert( outerParent_->Neighbor(outer_parent_face_id_) == innerParent_ );
    
    // creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
      
 } // end (constructor that infers face from higher-dimensional parent elements)







/**
    @note used by MeshManager::AddFace
*/
template<uint32_t dim>
Face<dim>::Face( const FiniteElementManager& fem_manager,
                 const FiniteVolumeStencilManager<dim>* fvm_manager,
                 Element<dim>* const inner_parent,
                 Element<dim>* const outer_parent,
                 uint32_t inner_parent_face_id,
                 uint32_t outer_parent_face_id,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,Face>( fem_manager.E( inner_parent->FE()->ElementTypeOfFace( inner_parent_face_id ) ) ),
    FiniteVolumePolicy<dim,Face>( (fvm_manager) ? fvm_manager->Stencil( inner_parent->FE()->ElementTypeOfFace( inner_parent_face_id ) ) : nullptr ),
    idx_(NULL_IDX),
    innerParent_(inner_parent),
    outerParent_(outer_parent),
    inner_parent_face_id_(inner_parent_face_id),
    outer_parent_face_id_(outer_parent_face_id)
 {
    assert( innerParent_ != nullptr );
    assert( outerParent_ != nullptr );
    if ( innerParent_ == outerParent_ ) {
         cerr <<"\n\nFace<"<< dim <<">(ctor: of face between parents): supplied pointers point to the same element: ";
         cerr << inner_parent->Idx();
         inner_parent->Out();
         innerParent_ = nullptr;
         outerParent_ = nullptr;
         assert( innerParent_ != outerParent_ );
         return;
      }
    assert( innerParent_->Neighbor(inner_parent_face_id_) == outerParent_ );
    assert( outerParent_->Neighbor(outer_parent_face_id_) == innerParent_ );

    // resizing the node and neighbor vectors
    node_connector_.resize(this->FE()->Nodes(),nullptr);
    face_connector_.resize(this->FE()->Faces(),nullptr);
    
    // assigning the nodes
    uint32_t i_node{0};
    for ( const auto& i : inner_parent->FE()->NodesOfFace( inner_parent_face_id_ ) )
      node_connector_[i_node++] = inner_parent->N(i);

    // creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
      
 } // end (constructor that infers face from higher-dimensional parent elements)









/**
    constructs Face object as an exact match of the target face of the supplied inner parent element
 
    @note boundary element from which face is constructed may have multiple 
    boundary faces, therefore the nth_boundary face variable is required

    @note used by MeshManager::AddBoundaryFace
*/
template<uint32_t dim>
Face<dim>::Face( Element<dim>& e,
                 csmp::FiniteElement* FE_type_for_face,
                 const FiniteVolumeStencilManager<dim>* fvm_manager,
                 uint32_t boundary_face,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::Face>(FE_type_for_face),
    FiniteVolumePolicy<dim,csmp::Face>( (fvm_manager) ? fvm_manager->Stencil(FE_type_for_face->ElementType()) : nullptr ),
    idx_(NULL_IDX),
    node_connector_(FE_type_for_face->Nodes(),nullptr),
    face_connector_(FE_type_for_face->Neighbors(),nullptr),
    innerParent_(&e),
    outerParent_(nullptr),
    inner_parent_face_id_(boundary_face)
 {
    assert( boundary_face < e.Faces() );
    if ( e.Neighbor(boundary_face) != nullptr ) outerParent_ = e.Neighbor(boundary_face);
    if constexpr ( dim == 2 ) assert( e.IsSurface() );
    if constexpr ( dim == 3 ) assert( e.IsVolume() );

    // 1. creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );

    // 2. assigning nodes to face in the same order as the face nodes
    //    of the inner parent element
    uint32_t n_count{0U};
    // nodes of the Element object from wich this Face is constructed
    for ( const auto& j : e.FE()->NodesOfFace( boundary_face ) ) {
         assert( e.N(j) != nullptr );
         node_connector_[n_count++] = e.N(j);
      }
   
 } // end constructor







/**
    Constructs model-edge line-element face connected with the two surface  elements of model boundary sharing its nodes
    
    The higher-dimensional Faces  have to share an edge with the line element.
    
    @param edge_nodes contains node pointers in the sequence in which they appear on the edge of
    the inner volumetric element that is adjacent to the Face
    
    @attention in this case, the higher-dimensional neighbors of the Face do not share faces,
    but edges with the lower-dimensional element. Therefore the face_node_id's are not assigned.

    @note used by MeshManager::AddEdgeFace
*/
template<uint32_t dim>
Face<dim>::Face( csmp::FiniteElement* FE_type_of_boundary_face,
                 const FiniteVolumeStencilManager<dim>* fvm_manager,
                 Element<dim>* const parent_of_face1,
                 Element<dim>* const parent_of_face2,
                 uint32_t parent_elmt1_segm_id,
                 uint32_t parent_elmt2_segm_id,
                 const std::vector<Node<dim>*>&  edge_nodes,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )
 : FiniteElementPolicy<dim,csmp::Face>(FE_type_of_boundary_face),
   FiniteVolumePolicy<dim,csmp::Face>( (fvm_manager) ? fvm_manager->Stencil(FE_type_of_boundary_face->ElementType()) : nullptr ),
   idx_(NULL_IDX),
   innerParent_(parent_of_face1),
   outerParent_(parent_of_face2),
   inner_parent_face_id_(parent_elmt1_segm_id),
   outer_parent_face_id_(parent_elmt2_segm_id),
   node_connector_(edge_nodes),
   face_connector_(2U,nullptr)
 {
    assert( innerParent_ != nullptr );
    // assert( outerParent_ != nullptr );
   
    // creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, ip );
    else
        this->ResizePropertyStorage( ep );
      
 } // end edge constructor






/**
    Re-constructor with storage but without connectivity.
    
    @note used to reconstruct model from VSet::VData / native binary file (MeshManager::Initialise)
*/
template<uint32_t dim>
Face<dim>::Face( size_t index,
                 csmp::FiniteElement* f,
                 const csmp::FiniteVolumeStencil<dim>* fvs,
                 const LocalVariables& ep,
                 const IntegrationPointVariables& ip )

  : FiniteElementPolicy<dim,csmp::Face>(f),
    FiniteVolumePolicy<dim,csmp::Face>(fvs),
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
template<uint32_t dim>
Face<dim>::Face( const Face<dim>& fc )
  : FiniteElementPolicy<dim,csmp::Face>(fc.FE()),
    FiniteVolumePolicy<dim,csmp::Face>(fc.FV()),
    idx_(fc.idx_),
    face_connector_(fc.face_connector_),
    node_connector_( fc.node_connector_),
    innerParent_(fc.innerParent_),
    outerParent_(fc.outerParent_),
    inner_parent_face_id_(fc.inner_parent_face_id_),
    outer_parent_face_id_(fc.outer_parent_face_id_)
  {
    assert( this->FE() != nullptr /* detect unitialized element*/ );
    assert( !face_connector_.empty() /* detect unitialized element*/ );
    // variable storage: call of initialization function
    this->LVS( fc.LVS() );
  }




/// move constructor
template<uint32_t dim>
Face<dim>::Face( Face<dim>&& fc )
  : FiniteElementPolicy<dim,csmp::Face>(fc.FE()),
    FiniteVolumePolicy<dim,csmp::Face>(fc.FV()),
    idx_(move(fc.idx_)),
    inner_parent_face_id_(move(fc.inner_parent_face_id_)),
    outer_parent_face_id_(move(fc.outer_parent_face_id_)),
    innerParent_(move(fc.innerParent_)),
    outerParent_(move(fc.outerParent_)),
    node_connector_(move(fc.node_connector_)),
    face_connector_(move(fc.face_connector_))
 {
    this->LVS( move(fc.LVS()) );
    
//    cerr <<"\nFace: called move constructor.";
 }




/**
     Before deleting the Face make sure that:
     
    // disconnecting the neighbor faces that are connected to this element
    if ( !face_connector_.empty() )
      for ( auto& it : face_connector_ )
        if ( it != nullptr  && !it->face_connector_.empty() )
          it->Unassign( this );
*/
template<uint32_t dim>
Face<dim>::~Face()
 {
 }




/**
    private as there should never be a need to use it
    and the results are undefined in terms of the pointers to parent and neighbor elements
    
    @note assignment operator is used in container::find(Face) operations; why?
*/
template<uint32_t dim>
Face<dim>&  Face<dim>::operator=( const Face<dim>& fc )
 {
    if ( &fc != this ) {
        FiniteElementPolicy<dim,csmp::Face>::Assign(fc.FE());
        FiniteVolumePolicy<dim,csmp::Face>::AssignFiniteVolume(fc.FV());
        idx_                  = fc.idx_;
        face_connector_       = fc.face_connector_;
        node_connector_       = fc.node_connector_;
        innerParent_          = fc.innerParent_; // problematic pointer assignment
        outerParent_          = fc.outerParent_;
        inner_parent_face_id_ = fc.inner_parent_face_id_;
        outer_parent_face_id_ = fc.outer_parent_face_id_;
        this->LVS( fc.LVS() );
      }
      
//    cerr <<"\nFace::operator=  called assignment operator.";
    
    return *this;
 }




/**
    @note a temporary value cannot be equivalent to lvalue.
*/
template<uint32_t dim>
Face<dim>&  Face<dim>::operator=( Face<dim>&& fc )
 {
    if ( &fc != this ) {
        FiniteElementPolicy<dim,csmp::Face>::Assign(fc.FE());
        FiniteVolumePolicy<dim,csmp::Face>::AssignFiniteVolume(fc.FV());
        idx_                  = move( fc.idx_ );
        inner_parent_face_id_ = move( fc.inner_parent_face_id_ );
        outer_parent_face_id_ = move( fc.outer_parent_face_id_ );
        innerParent_          = move( fc.innerParent_ );
        outerParent_          = move( fc.outerParent_ );
        face_connector_       = move(fc.face_connector_);
        node_connector_       = move(fc.node_connector_);
        this->LVS( move(fc.LVS()) );
     }
    
//    cerr <<"\nFace::operator=  called move-assignment operator.";

    return *this;
 }




// CONSTRUCTION PROCESS

/**
    Connect Face object with its nodes.
*/
template<uint32_t dim>
void Face<dim>::Assign( uint32_t i, csmp::Node<dim>* const nd_ptr )
 {
    assert( i < node_connector_.size() );
    node_connector_[i] = nd_ptr;
 }



template<uint32_t dim>
void Face<dim>::Assign( uint32_t i, csmp::Face<dim>* const f_ptr )
 {
    assert( i < face_connector_.size() );
    face_connector_[i] = f_ptr;
 }


/**
      Removes the first occurrence of the neighbor pointed to by f_ptr from the neighbor array.
      Any further occurrence is ignored.
      
      @note call twice if there are multple occurrences.
      
      @note if called with a nullptr argument method will return false.
*/
template<uint32_t dim>
bool Face<dim>::Unassign( const Face<dim>* const f_ptr )
{
  if ( f_ptr == nullptr ) return false;
	for ( auto i{0U}; i < face_connector_.size(); ++i )
		if ( f_ptr == face_connector_[i] )
      {
         face_connector_[i] = nullptr;
         return true;
      }
	return false;
  
} // end Unassign



template<uint32_t dim>
uint32_t  Face<dim>::Nodes() const
{
	return static_cast<uint32_t>(node_connector_.size());
}



template<uint32_t dim>
uint32_t  Face<dim>::Neighbors() const
{
	return static_cast<uint32_t>(face_connector_.size());
}


template<uint32_t dim>
uint32_t  Face<dim>::ConnectedNeighbors() const
{
	size_t nulls(0);
	for (auto& f : face_connector_)
		if ( f == nullptr ) nulls++;
	return static_cast<uint32_t>(face_connector_.size() - nulls);
}


template<uint32_t dim>
uint32_t  Face<dim>::Faces() const
{
	return static_cast<uint32_t>(face_connector_.size());
}


/**
    Connects face to the higher-dimensional parent elements which it may be sandwiched between.
    Node pointers are reassigned if they were already pointing to other nodes.
    
    @note method assumes that the Face object already has valid node pointers
    
    @attention the nodes must be assigned to the Face before this method is called
 
    @attention The node numbering of the Face determines the direction of its normal.
    Here the convention is assumed that the first parent element is that on the inside
    of the Face with regard to the outward pointing normal and the second element is on
    the outside. It follows that the node sequences are the same.
*/
template<uint32_t dim>
void Face<dim>::Assign( Element<dim>* const innerElement, Element<dim>* const outerElement )
  {
    // 1. argument checks and assignments
    // ----------------------------------
    assert( innerElement != nullptr ); // inner element must be defined
    if constexpr ( dim == 3U ) assert( innerElement->IsVolume() );
    if constexpr ( dim == 2U ) assert( innerElement->IsSurface() );
    innerParent_ = innerElement;
    
    if ( outerElement != nullptr ) { // outer element is defined if face is in model interior
         if constexpr ( dim == 3U ) assert( outerElement->IsVolume() );
         if constexpr ( dim == 2U ) assert( outerElement->IsSurface() );
         outerParent_ = outerElement;
      }
    
    // 2. finding the face of the inner element which corresponds to this Face and checking its node numbering
    // -------------------------------------------------------------------------------------------------------
    // creating a search key for Face
    set<Node<dim>*>  face_nds;
    const size_t face_nodes(Nodes());
    for ( auto i{0U}; i<face_nodes; ++i ) {
         // are the nodes there?
         assert( this->N(i) != nullptr );
         face_nds.insert( this->N(i) );
      }

    // 2.1 matching the lower-dimensional face to a face of inner higher-dimensional element
    // -------------------------------------------------------------------------------------
    if ( ((dim == 3U) && this->IsSurface() ) || // Face is either a triangle or a quadrilateral in 3D
         ((dim == 2U) && this->IsLine()) )      // Face is a line in 2D
      {
         // searching the matching Face of the inner parent element
         bool             matching_face_found(false);
         const size_t     faces(innerParent_->Faces());
         set<Node<dim>*>  parent_nds;
        
         // inner parent
         for ( auto i{0U}; i<faces; ++i ) {
              for ( const auto& j : innerParent_->FE()->NodesOfFace(i) )
                parent_nds.insert( innerParent_->N(j) );
              // checking
              if ( face_nds == parent_nds ) {
                    matching_face_found = true;
                    inner_parent_face_id_ = i;
                    break;
                }
              parent_nds.clear();
           }         

         // outer parent if any
         if ( outerElement != nullptr ) {
             for ( auto i{0U}; i<faces; ++i ) {
                  for ( const auto& j : outerParent_->FE()->NodesOfFace(i) )
                    parent_nds.insert( outerParent_->N(j) );
                  // checking
                  if ( face_nds == parent_nds ) {
                        matching_face_found = true;
                        outer_parent_face_id_ = i;
                        break;
                    }
                  parent_nds.clear();
               }
           }

         return;
      }
    
    // 3. matching the dimension-2 face with an edge of the higher-dimensional element
    // -------------------------------------------------------------------------------
    // (if the higher-dimensional element has 2 more dimensions that the Face, a comparison
    //  with its edges needs to be performed)
    // searching for the matching Segment (Edge) of the inner parent element
    bool             matching_segment_found(false);
    const size_t     segments(innerParent_->Segments());
    vector<uint32_t>   nodes_of_segm;
    set<Node<dim>*>  parent_nds;
    
    // inner parent
    for ( auto i{0U}; i<segments; ++i ) {
         innerParent_->FE()->NodesOfSegment( i, nodes_of_segm );
         for ( size_t j{0U}; j<nodes_of_segm.size(); ++j )
           parent_nds.insert( innerParent_->N( nodes_of_segm[j] ) );
         // checking wether line element matches edge dim+2 element
         if ( face_nds == parent_nds ) {
               inner_parent_face_id_ = i;
               matching_segment_found = true;
               break;
           }
         parent_nds.clear();
      }

    // outer parent
    if ( outerElement != nullptr ) {
        for ( auto i{0U}; i<segments; ++i ) {
             outerParent_->FE()->NodesOfSegment( i, nodes_of_segm );
             for ( size_t j{0U}; j<nodes_of_segm.size(); ++j )
               parent_nds.insert( outerParent_->N( nodes_of_segm[j] ) );
             // checking wether line element matches edge dim+2 element
             if ( face_nds == parent_nds ) {
                   outer_parent_face_id_ = i;
                   matching_segment_found = true;
                   break;
               }
             parent_nds.clear();
          }
      }
      
 } // end assign




template<uint32_t dim>
void Face<dim>::Unassign( const csmp::Node<dim>* const nd_ptr )
  {
    if ( nd_ptr == nullptr ) return;
    for ( auto i{0U}; i < node_connector_.size(); i++ )
        if ( nd_ptr == node_connector_[i] ) {
            node_connector_[i] = nullptr;
            break;
          }
  }




/// iterator to the element nodes
template<uint32_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Face<dim>::NodesBegin() const
{
  return node_connector_.begin();
}

template<uint32_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Face<dim>::NodesEnd() const
{
  return node_connector_.end();
}

template<uint32_t dim>
typename std::vector<Face<dim>*>::const_iterator  Face<dim>::NeighborsBegin() const
{
  return face_connector_.begin();
}

template<uint32_t dim>
typename std::vector<Face<dim>*>::const_iterator  Face<dim>::NeighborsEnd() const
{
  return face_connector_.end();
}




// VISITOR

/// SKM temporary patch to get Visitor objects to work on boundaries
template<uint32_t dim>
void Face<dim>::Accept( csmp::Visitor<dim>& vis )
 {
    if ( vis.ApplicationTarget() == FACE ) {
         vis.Visit(this);
         return;
      }
    if ( vis.ApplicationTarget() == NODE ) {
         uint32_t n_nodes{ this->Nodes() };
         for ( auto i{0U}; i< n_nodes; i++ )
             this->N(i)->Accept( vis );
         return;
      }
    throw logic_error("Face<dim>::Accept: target of visitation unresolved.");

 } // end Accept




  

// ACCESSORS



/// returns the current index of this face assuming that a meaningful value was assigned earlier
template<uint32_t dim>
void  Face<dim>::Idx( size_t idx_to_assign ) const
  {
    idx_ = idx_to_assign;
  }

template<uint32_t dim>
size_t   Face<dim>::Idx() const
  {
    return idx_;
  }


template<uint32_t dim>
typename std::vector<csmp::Face<dim>*>&  Face<dim>::NeighborElementVector()
  {
    return face_connector_;
  }


template<uint32_t dim>
csmp::Node<dim>* const Face<dim>::N( uint32_t n ) const
  {
     assert( n < Nodes() );
     return node_connector_[n];
  }


/**
    watch out if there is no neighbor this returns a NULL pointer
*/
template<uint32_t dim>
csmp::Face<dim>* const Face<dim>::Neighbor( uint32_t n ) const
 {
    assert( n < Neighbors() );
    return face_connector_[n];
 }


/**
    Returns either the INSIDE or OUTSIDE elements
    connected to this face if they exist.
*/
template<uint32_t dim>
Element<dim>* const Face<dim>::Parent( INTERFACE_SIDE side ) const
 {
    assert( side != MIDDLE );
    assert( side != 0 );
    assert( innerParent_ != nullptr );
    return ( side == INSIDE ) ? innerParent_ : outerParent_;
 }


/**
    inner parent is always initialised for a valid face
*/
template<uint32_t dim>
Element<dim>* const Face<dim>::InnerParent() const
 {
    assert( innerParent_ != nullptr );
    return innerParent_;
 }


/**
    @attention outer parent will not be initialised if the Face lies on the outside boundary of the model
*/
template<uint32_t dim>
Element<dim>* const Face<dim>::OuterParent() const
 {
    return outerParent_;
 }




/// return face ID of inner parent element
template<uint32_t dim>
uint32_t  Face<dim>::InnerParentFaceID() const
{
  if ( innerParent_ == nullptr ) return NULL_IDX;
  return inner_parent_face_id_;
}




/// return face ID of outer parent element
template<uint32_t dim>
uint32_t  Face<dim>::OuterParentFaceID() const
{
  if ( outerParent_ == nullptr ) return NULL_IDX;
  return outer_parent_face_id_;
}




/// return face ID of inner or outer parent element depending on index
template<uint32_t dim>
uint32_t  Face<dim>::ParentFaceID( INTERFACE_SIDE side ) const
{
  if ( side == INSIDE )
  {
    assert( innerParent_ != nullptr );
    return inner_parent_face_id_;
  }
  else if ( side == OUTSIDE )
  {
    assert( outerParent_ != nullptr );
    return outer_parent_face_id_;
  }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.Note( WARNING, "csmp::Face<dim>::ParentFaceID:", "Interface 'side' could not be determined." );

  return inner_parent_face_id_;
}




/// assign face ID of inner or outer parent element depending on their local face numbering
template<uint32_t dim>
void  Face<dim>::ParentFaceID( INTERFACE_SIDE side, uint32_t idx )
{
  if ( side == INSIDE ) {
      assert( innerParent_ != nullptr );
      inner_parent_face_id_ = idx;
      return;
    }
  
  if ( side == OUTSIDE ) {
      assert( outerParent_ != nullptr );
      outer_parent_face_id_ = idx;
      return;
    }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.Note( WARNING, "csmp::Face<dim>::ParentFaceID:", "Interface 'side' could not be determined." );

} // end ParentFaceID(assignment)










// GEOMETRY

/**
     Returns a vector (not a unit normal) connecting the barycenter of the Face
     with that of the higher-dimensional parent element on its inside.
*/
template<uint32_t dim>
void Face<dim>::VectorToInnerElementBaryCenter( VectorVariable<dim>& vectorToInner ) const
{
   assert( Parent(INSIDE) != nullptr );
   vectorToInner = Parent(INSIDE)->BaryCenter() - this->BaryCenter();
}


/// returns area (3D) or length (2D) of face
template<uint32_t dim>
double Face<dim>::Area() const
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

The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<uint32_t dim>
void  Face<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
  {
    const auto n_nodes{ Nodes() };
    XY.Resize( n_nodes, dim );
    for ( auto i{0U}; i<n_nodes; ++i )
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
template<uint32_t dim>
Point<dim>  Face<dim>::BaryCenter() const
  {
    Point<dim>  pt(N(0U)->Coordinate());
    const auto  n_nodes{ Nodes() };
    for ( auto i{1U}; i<n_nodes; ++i )
      pt += N(i)->Coordinate();

    return pt / static_cast<double>(Nodes());
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
template<uint32_t dim>
double  Face<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
  {
    double fMinTemp( static_cast<double>( DBL_MAX) );
    double fMaxTemp( static_cast<double>(-DBL_MAX) );

    // this normalisation is necessary because the vector variable
    // being any physical quantity may have any magnitude
    const double fMagnitudeOfDirection(vecDirection.Length());
    // avoid division by zero
    assert( fMagnitudeOfDirection >= numeric_limits<double>::epsilon() );

    const auto n_nodes{ Nodes() };
    for ( auto i{0U}; i<n_nodes; ++i ) {
        // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
        double fTemp(vecDirection.DotProduct( N(i)->Coordinate() ));
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
template<uint32_t dim>
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
    const auto  n_nodes{ Nodes() };
    V.resize(n_nodes);

    for ( auto i{0U}; i<n_nodes; i++ )
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
template<uint32_t dim>
void  Face<dim>::Out() const
 {
    cout <<"\n\n\nFace<"<< dim <<">::Out: number: "<< idx_;

    string str;
    cout <<" ("<< parseFiniteElementType(this->FE_Type()) <<")";
    if ( outerParent_ == nullptr ) cout <<", face is located at model boundary.";
    cout << endl;

    cout <<"\nInternal data: "<< endl;

    cout <<"\n\tconnected nodes with boundary flags:  ";
    for ( auto i{0U}; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i)->AtBoundary());
         cout << N(i)->Idx() <<":"<< str <<"  ";
      }
    cout << endl;

    cout <<"\n\tconnected neighbor Face types / boundary flags:\n";
    for ( auto i{0U}; i<this->Neighbors(); i++ )
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
         cout <<"\t"<<"inside higher-dim parent Element: "<< this->innerParent_->Idx();
         cout  <<" ("<< parseFiniteElementType(this->Parent(INSIDE)->FE_Type()) <<")"<< endl;
      }
    else cout <<"\tnone.\n";
    if ( this->outerParent_ != nullptr ) {
         cout <<"\t"<<"outside higher-dim parent Element: "<< this->outerParent_->Idx();
         cout <<" ("<< parseFiniteElementType(this->Parent(OUTSIDE)->FE_Type()) <<")"<< endl;
      }
    else cout <<"\tnone.\n";

 } // end Out


// ==========================================================================================================
// TO DEPRECATE
// ==========================================================================================================

/**
  SKM revised 8/02/22
*/
template<uint32_t dim>
bool  Face<dim>::operator==( const Face<dim>& fc ) const
 {
	 if (&fc != this) {
       // do they have the same higher-dimensional neighbors
       if ( innerParent_ != fc.innerParent_ || outerParent_ != fc.outerParent_ ) return false;
       // are these neighbors adjacent to the same faces
       if ( inner_parent_face_id_ != fc.inner_parent_face_id_ ||
            outer_parent_face_id_ != fc.outer_parent_face_id_ ) return false;
       // do they have the same nodes
       if ( node_connector_ != fc.node_connector_ ) return false;
       // do they have the same neighbors
       if ( face_connector_ != fc.face_connector_ ) return false;
     }
    return true;
 }




template<uint32_t dim>
void* Face<dim>::operator new( size_t size )
  {
//      std::cout<< "\nFace<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

template<uint32_t dim>
void Face<dim>::operator delete( void* p )
  {
//     std::cout<< "\nFace<"<< dim <<">: called overloaded delete operator.\n";
     free(p);
     p = nullptr;
  }


// DEPRECATED

/*  RECREATE IF NEEDED USING THE NODE TO PARENT ELEMENT CAPABILITY
template<uint32_t dim>
size_t  Face<dim>::ParentNodeNumber( size_t n ) const
 {
    assert( n < parent_element_node_ids_.size() );
    return parent_element_node_ids_[n];
 }
*/





/**
     Reports the local number (0..faces-1) of the face of the inner higher-dimensional parent element that
     matches the nodes of the Face object (also in terms of the sequence of these nodes).
     
     @return local face number or max() uint if no index could be found.

     @attention assumes that the Face has valid nodes and its inner parent element is connected
*/
/*
template<uint32_t dim>
size_t Face<dim>::ParentFaceNumber(INTERFACE_SIDE side) const
 {

    Element<dim>* parent;
    switch (side){
    case INSIDE: parent = innerParent_;
      break;
    case OUTSIDE: parent = outerParent_;
      break;
    case MIDDLE: throw csmp::Exception(ERROR, "Face::ParentFaceNumber" , "Face does not have a MIDDLE parent!");
      break;
    }

    assert(parent != nullptr);
    assert( !node_connector_.empty() );
   
     // 1. creating a unique key from the nodes of the Face
    set<uint32_t>  face_key;
    for ( auto nit=node_connector_.begin(); nit!=node_connector_.end(); ++nit ) {
         assert( (*nit) != nullptr );
         face_key.insert( (*nit)->Idx() );
      }
    
    // 2. creating face keys for the inner parent element and trying to match them
    //    with the one created for the current face
    vector<uint32_t> fnids;
    set<uint32_t>    face_key_n;
    const size_t faces(parent->Faces());
    for ( auto i{0U}; i<faces; ++i ) {
         parent->FE()->NodesOfFace( i, fnids );
         for ( size_t j{0U}; j<fnids.size(); ++j )
           face_key_n.insert( parent->N( fnids[j] )->Idx() );
         // has a matching face been found?
         if ( face_key == face_key_n ) return i;
         else face_key_n.clear();
      }

    return std::numeric_limits<uint32_t>::max();
 }
*/


// explicits
template class Face<1U>;
template class Face<2U>;
template class Face<3U>;


} // end namespace csmp
