#include "InterFace.h"
#include "Element.h"
#include "Face.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "Visitor.h"
#include "variableOperations.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {


/** New! SKM 29/7/2022: constructs  InterFace using the nodes and their numbering from the InterFace's higher-dimensional neighbors.
    A check is performed to ascertain that the nodes are indeed collocated.
    If not, the outside nodes are rotated until a match is obtained or a FATAL_ERROR is reported.
 
    @attention constructor expects that the Nodes have already been multiplicated and turned into manifolds elsewhere.
               Nodes of face of outside element are taken and rotated to be matching in the "reverse" order to Inside Node.
               Note: "reverse" is in quotation marks because it's not quite reversed if you consider MidPoint Nodes

   @attention outside nodes are rotated so that the last corner node of Inside matches with the first corner node on the outside
    this is needed to ensure MatchingN(i,INSIDE) == MatchingN(i,OUTSIDE)
    this breaks the idea that outside nodes match with the nodes of the face of outside element (they may be rotated )
    Therefore N(i,OUTSIDE) != OuterParent->N( OuterParent->NodesOfFace( outer_parent_face_id)[i] )
*/
template<uint32_t dim>
InterFace<dim>::InterFace( csmp::Element<dim>& elmt,
                           csmp::Element<dim>* inner_parent,
                           csmp::Element<dim>* outer_parent,
                           uint32_t adjacent_face_of_inner_element,
                           uint32_t adjacent_face_of_outer_element,
                           const LocalVariables&  interface_props,
                           const IntegrationPointVariables&  interface_integration_point_props )
  : FiniteElementPolicy<dim,csmp::InterFace>{ elmt.FE() },
    FiniteVolumePolicy<dim,csmp::InterFace>{ elmt.FV() },
    idx_{ numeric_limits<size_t>::max() },
    innerParent_{ inner_parent },
    outerParent_{ outer_parent },
    middleElement_{ nullptr },
    node_connector_( elmt.Nodes() * 2, nullptr ),
    interface_connector_( elmt.Neighbors(), nullptr ),
    inner_parent_face_id_{ adjacent_face_of_inner_element },
    outer_parent_face_id_{ adjacent_face_of_outer_element },
    current_side_{ INSIDE }
{
   assert( elmt.FE() != nullptr );
   assert( innerParent_ != nullptr );
   assert( outerParent_ != nullptr );
   assert( innerParent_->Neighbor(inner_parent_face_id_) == outerParent_ );
   assert( outerParent_->Neighbor(outer_parent_face_id_) == innerParent_ );

   // setting up storage for the local variables
   if ( this->UsesLocalCoordinates() ) {
        this->ResizePropertyStorage( interface_props, interface_integration_point_props );
     }
   else this->ResizePropertyStorage( interface_props );

   // 1. assigning the nodes to the new InterFace
   // -------------------------------------------
   // INSIDE nodes ---------
   uint32_t count{0U};
   for ( const auto& nit : innerParent_->FE()->NodesOfFace( inner_parent_face_id_ ) )
     Assign( count++, innerParent_->N(nit), INSIDE );

   // OUTSIDE nodes --------
   count = 0U;
   for ( const auto& nit : outerParent_->FE()->NodesOfFace( outer_parent_face_id_ ) )
     Assign( count++, outerParent_->N(nit), OUTSIDE );

   // Outside nodes are then rotated so that the last corner node of Inside matches with the first corner node on the outside
   // this is needed to ensure MatchingN(i,INSIDE) == MatchingN(i,OUTSIDE)
   // this breaks the idea that outside nodes match with the nodes of the face of outside element (they may be rotated )
   // Therefore N(i,OUTSIDE) != OuterParent->N( OuterParent->NodesOfFace( outer_parent_face_id)[i] )
   InitialiseNodeVector();

   // 2. checking that the nodes are collocated
   // ----------------------------------------------------------------------
#ifdef DEBUG
   if ( !AreNodesCollocated() )
     throw csmp::Exception( ERROR, "InterFace(costum contructor", "supplied interface nodes are not collocated");
#endif

   // 3. detaching the higher-dimensional element neighbors from one another
   // ----------------------------------------------------------------------
   innerParent_->Unassign( outerParent_ );
   outerParent_->Unassign( innerParent_ );

 } // end complete custom constructor (Element)




/**
   Contructs complete InterFace using the Face nodes as inside nodes; outside nodes ready to use in opposite order are supplied by node-pointer vector'
   
   @note it takes care of assigning the higher dimensional neighbor elements
   @note assigns outside nodes to InterFace also changing these nodes on the higher-dimensional outside element
   @note detaches the neighbor connection between the higher dimensional elements

 */
template<uint32_t dim>
InterFace<dim>::InterFace( csmp::Face<dim>* fptr,
               const LocalVariables&  interface_props,
               const IntegrationPointVariables&  interface_integration_point_props,
               vector<Node<dim>*> outside_nodes )
  : FiniteElementPolicy<dim,csmp::InterFace>{ fptr->FE() },
    FiniteVolumePolicy<dim,csmp::InterFace>{ fptr->FV() },
    innerParent_{ fptr->InnerParent() },
    outerParent_{ fptr->OuterParent() },
    middleElement_{ nullptr },
    node_connector_( fptr->Nodes() * 2, nullptr ),
    interface_connector_( fptr->Neighbors(), nullptr ),
    idx_{ numeric_limits<size_t>::max() },
    inner_parent_face_id_{ fptr->InnerParentFaceID() },
    outer_parent_face_id_{ fptr->OuterParentFaceID() },
    current_side_{ INSIDE }
{
   assert( fptr != nullptr );
   assert( outerParent_ != nullptr ); // interfaces must have neighbors on all sides
   assert( outside_nodes.size() == fptr->Nodes() );

   // setting up storage for the local variables
   if ( this->UsesLocalCoordinates() ) {
        this->ResizePropertyStorage( interface_props, interface_integration_point_props );
     }
   else this->ResizePropertyStorage( interface_props );

#ifdef DEBUG
   for ( const auto& nit : outside_nodes ) assert( nit != nullptr );
#endif

   // 1. assigning the nodes to the new InterFace
   // -------------------------------------------
   const auto n_nodes{ outside_nodes.size() };
   for ( uint32_t i{0U}; i< n_nodes; i++ ) {
        Assign( i, fptr->N(i), INSIDE );
        Assign( i, outside_nodes[i], OUTSIDE );
     }
     
   // 2. replacing the nodes on the outside element with the new outside nodes
   // ------------------------------------------------------------------------
   uint32_t n_count{0U};
   for ( const auto& i : outerParent_->FE()->NodesOfFace( outer_parent_face_id_ ) )
     outerParent_->Assign( i, outside_nodes[ n_count++ ] );
     
   // 3. detaching the higher-dimensional element neighbors from one another
   // ----------------------------------------------------------------------
   innerParent_->Unassign( outerParent_ );
   outerParent_->Unassign( innerParent_ );

 } // end complete custom constructor (Face)





/// custom constructor which also builds variable storage; used in most cases
template<uint32_t dim>
InterFace<dim>::InterFace( csmp::FiniteElement* f,
                           const csmp::FiniteVolumeStencil<dim>* fvs,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::InterFace>{ f },
    FiniteVolumePolicy<dim,csmp::InterFace>{ fvs },
    innerParent_{ nullptr },
    outerParent_{ nullptr },
    middleElement_{ nullptr },
    node_connector_( f->Nodes() * 2, nullptr ),
    interface_connector_( f->Neighbors(), nullptr ),
    idx_{ numeric_limits<size_t>::max() },
    current_side_{ INSIDE }
{
   assert( f != nullptr );

   if ( this->UsesLocalCoordinates() )
     this->ResizePropertyStorage( ep, ip );
   else
     this->ResizePropertyStorage( ep );
}


template<uint32_t dim>
InterFace<dim>::InterFace( size_t index,
                           csmp::FiniteElement* f,
                           const csmp::FiniteVolumeStencil<dim>* fvs,
                           const LocalVariables& ep,
                           const IntegrationPointVariables& ip )
  : FiniteElementPolicy<dim,csmp::InterFace>{ f },
    FiniteVolumePolicy<dim,csmp::InterFace>{ fvs },
    innerParent_{ nullptr },
    outerParent_{ nullptr },
    middleElement_{ nullptr },
    node_connector_( f->Nodes() * 2, nullptr ),
    interface_connector_( f->Neighbors(), nullptr ),
    idx_{ index },
    current_side_{ INSIDE }
{
   assert( f != nullptr );

   if ( this->UsesLocalCoordinates() )
     this->ResizePropertyStorage( ep, ip );
   else
     this->ResizePropertyStorage( ep );
}


/// copy constructor
template<uint32_t dim>
InterFace<dim>::InterFace( const InterFace<dim>& ifc )
  : FiniteElementPolicy<dim,csmp::InterFace>{ ifc },
    FiniteVolumePolicy<dim,csmp::InterFace>{ ifc.FV() },
    LocalVariableStorage<dim,csmp::InterFace>{ ifc },
    innerParent_{ ifc.innerParent_ },
    outerParent_{ ifc.outerParent_ },
    middleElement_{ ifc.middleElement_ },
    node_connector_{ ifc.node_connector_ },
    interface_connector_{ ifc.interface_connector_ },
    idx_{ ifc.idx_ },
    inner_parent_face_id_{ ifc.inner_parent_face_id_ },
    outer_parent_face_id_{ ifc.outer_parent_face_id_ },
    current_side_{ ifc.current_side_ }
{
  assert( !interface_connector_.empty() );
}



/// move constructor
template<uint32_t dim>
InterFace<dim>::InterFace( InterFace<dim>&& ifc )
  : FiniteElementPolicy<dim, csmp::InterFace>{ ifc },
    FiniteVolumePolicy<dim, csmp::InterFace>{ ifc.FV() },
    LocalVariableStorage<dim,csmp::InterFace>{ ifc },
    innerParent_{ ifc.innerParent_ },
    outerParent_{ ifc.outerParent_ },
    middleElement_{ ifc.middleElement_ },
    node_connector_{ ifc.node_connector_ },
    interface_connector_{ ifc.interface_connector_ },
    idx_{ ifc.idx_ },
    inner_parent_face_id_{ ifc.inner_parent_face_id_ },
    outer_parent_face_id_{ ifc.outer_parent_face_id_ },
    current_side_{ ifc.current_side_ }
 {
    assert( !interface_connector_.empty() ); // detected unitialized element

    ifc.AssignFiniteElementNullPtr();
    ifc.AssignFiniteVolumeNullPtr();
 }




/// assignment
template<uint32_t dim>
InterFace<dim>&  InterFace<dim>::operator=( const InterFace<dim>& ifc )
{
  if ( &ifc != this )
    {
      FiniteElementPolicy<dim, csmp::InterFace>::Assign( ifc.FE() );
      FiniteVolumePolicy<dim, csmp::InterFace>::AssignFiniteVolume( ifc.FV() );
      LocalVariableStorage<dim,csmp::InterFace>::LVS( ifc.LVS() );
      idx_ = ifc.idx_;
      node_connector_ = ifc.node_connector_;
      interface_connector_ = ifc.interface_connector_;
      middleElement_ = ifc.middleElement_;
      current_side_ = ifc.current_side_;
      innerParent_ = ifc.innerParent_;
      outerParent_ = ifc.outerParent_;
      inner_parent_face_id_ = ifc.inner_parent_face_id_;
      outer_parent_face_id_ = ifc.outer_parent_face_id_;
    }
  return *this;
}



template<uint32_t dim>
InterFace<dim>&  InterFace<dim>::operator=( InterFace<dim>&& ifc )
{
  assert( &ifc != this );

  FiniteElementPolicy<dim, csmp::InterFace>::Assign( ifc.FE() );
  FiniteVolumePolicy<dim, csmp::InterFace>::AssignFiniteVolume( ifc.FV() );
  LocalVariableStorage<dim,csmp::InterFace>::LVS( ifc.LVS() );
  idx_                  = ifc.idx_;
  inner_parent_face_id_ = ifc.inner_parent_face_id_;
  outer_parent_face_id_ = ifc.outer_parent_face_id_;
  innerParent_          = ifc.innerParent_;
  outerParent_          = ifc.outerParent_;
  middleElement_        = ifc.middleElement_;
  interface_connector_  = std::move( ifc.interface_connector_ );
  node_connector_       = std::move( ifc.node_connector_ );
  current_side_         = ifc.current_side_;

  ifc.AssignFiniteElementNullPtr();
  ifc.AssignFiniteVolumeNullPtr();

  return *this;
}





/**
    WARNING: this operator is used specifically in the process of creation of particular SplitBoundary.
    Therefore only important infromation for that process is taken into account in order to distinguish two InterFace's.
     That must be reference to inner and outer parent Elements, inner and outer parent face ID's
*/
template<uint32_t dim>
bool  InterFace<dim>::operator==( const InterFace<dim>& ifc ) const
{
  if ( &ifc != this )
    if ( innerParent_ != ifc.innerParent_ ||
         outerParent_ != ifc.outerParent_ ||
         inner_parent_face_id_ != ifc.inner_parent_face_id_ ||
         outer_parent_face_id_ != ifc.outer_parent_face_id_ )
      return false;

  return true;
}



template<uint32_t dim>
void* InterFace<dim>::operator new( size_t size )
  {
//      cout<< "\nInterFace<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

template<uint32_t dim>
void InterFace<dim>::operator delete( void* p )
  {
//     cout<< "\nInterFace<"<< dim <<">: called overloaded delete operator.\n";
     free(p);
     p = nullptr;
  }


/// iterator to the element nodes
template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator  InterFace<dim>::NodesBegin() const noexcept
{
   if ( current_side_ == INSIDE ) return node_connector_.begin();
   return next( node_connector_.begin(), this->FE()->Nodes() );
}

template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator  InterFace<dim>::NodesEnd() const noexcept
{
  if ( current_side_ == INSIDE ) return next( node_connector_.begin(), this->FE()->Nodes() );
  return node_connector_.end();
}



template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator InterFace<dim>::CornerNodesBegin() const noexcept
{
   if ( current_side_ == INSIDE ) return node_connector_.begin();
   return next( node_connector_.begin(), this->FE()->Nodes() );
}


template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator InterFace<dim>::CornerNodesEnd() const noexcept
{
  if ( current_side_ == INSIDE ) return next( node_connector_.begin(), this->FE()->Nodes() );
  return next(node_connector_.begin(),this->FE()->Nodes() + this->FE()->CornerNodes());
}



template<uint32_t dim>
typename vector<InterFace<dim>*>::const_iterator  InterFace<dim>::NeighborsBegin() const noexcept
{
  return interface_connector_.begin();
}

template<uint32_t dim>
typename vector<InterFace<dim>*>::const_iterator  InterFace<dim>::NeighborsEnd() const noexcept
{
  return interface_connector_.end();
}



// VISITOR
template<uint32_t dim>
void InterFace<dim>::Accept( csmp::Visitor<dim>& vis )
{
  if ( vis.ApplicationTarget() == INTER_FACE ) {
      vis.Visit( this );
      return;
    }
  throw csmp::Exception( ERROR, "InterFace<dim>::Accept", "Target of visitation unresolved (nodes not an option)." );
  
} // end Accept


/**
    Assigning neighbor InterFace objects on either side of the InterFace.
*/
template<uint32_t dim>
void InterFace<dim>::Assign( uint32_t i, InterFace<dim>* const ifc_ptr ) // neighbor interface
{
  assert( i < interface_connector_.size() );
  interface_connector_[i] = ifc_ptr;
  
} // end InterFace<dim>::Assign(neighbor)



template<uint32_t dim>
void InterFace<dim>::Assign( Element<dim>* const base_elmt )
{
  assert( base_elmt != nullptr );
  middleElement_ = base_elmt;
}


/**
Connects face to the higher-dimensional elements which it is sandwiched between.
These elements are stored in the local pair  of pointers called parents_.

@attention The node numbering of the Face determines the direction of its normal.
Here the convention is assumed that the first parent element is that on the inside
of the InterFace with regard to the outward pointing normal and the second element is on
the outside.

Same method as in Face, but - in addition - connects InterFace to the multiplicated nodes
on either side.
*/
template<uint32_t dim>
void InterFace<dim>::Assign( Element<dim>* const inner_elmt, uint32_t inner_local_face_id,
                             Element<dim>* const outer_elmt, uint32_t outer_local_face_id )
{
  assert( inner_elmt != nullptr );
  assert( outer_elmt != nullptr );
  innerParent_ = inner_elmt;
  outerParent_ = outer_elmt;
  assert( inner_local_face_id < innerParent_->Faces() );
  assert( outer_local_face_id < outerParent_->Faces() );
  inner_parent_face_id_ = inner_local_face_id;
  outer_parent_face_id_ = outer_local_face_id;
}




template<uint32_t dim>
void InterFace<dim>::Assign( Element<dim>* const parentElement, uint32_t faceId, INTERFACE_SIDE side )
{
  assert( parentElement != nullptr );

  if ( side == INSIDE ) {
      innerParent_ = parentElement;
      inner_parent_face_id_ = faceId;
      return;
    }
  
  if ( side == OUTSIDE ) {
      outerParent_ = parentElement;
      outer_parent_face_id_ = faceId;
      return;
    }
  
  if ( side == MIDDLE ) {
      middleElement_ = parentElement;
      throw csmp::Exception( ERROR, "csmp::InterFace<dim>::Assign( Element, face_id, side )", "Wrong method to assign intervening element!" );
    }

} // end Assign





/// connects interFace to its higher-dimensional neighbors and nodes establishing connections by itself
template<uint32_t dim>
void InterFace<dim>::AssignElementsAndNodes( Element<dim>* const inner_elmt,
                                             Element<dim>* const outer_elmt )
 {
     assert( inner_elmt != nullptr );
     assert( outer_elmt != nullptr );
     
     // 0. higher-dimensional elements adjacent to InterFace
     innerParent_ = inner_elmt;
     outerParent_ = outer_elmt;
     
     // 1. numbers of the shared faces
     pair<uint32_t, uint32_t> shared_faces = SharedElementFacesAndFaceIDs();
     inner_parent_face_id_ = shared_faces.first;
     outer_parent_face_id_ = shared_faces.second;
     
     // 2. connect nodes
     InitialiseNodeVector();
     
 } // end AssignElementsAndNodes




/**
    Assigns nodes to node_connector_ vector of the InterFace

    @note the nodes on the inner side of the interface are the first stored in the node connector
    vector, the outer ones follow in this single node vector.

    @author SKM 16/6/2016
*/
template<uint32_t dim>
void InterFace<dim>::Assign( uint32_t n_local, Node<dim>* nptr, INTERFACE_SIDE side )
{
   assert( nptr != nullptr );
   assert( this->FE() != nullptr );
   assert( !node_connector_.empty() );
   
   const auto finite_element_nodes( this->FE()->Nodes() );
   assert( finite_element_nodes * 2 == node_connector_.size() );
   assert( n_local < finite_element_nodes );

   if ( side == INSIDE ) {
        node_connector_[n_local] = nptr;
        return;
     }

   if ( side == OUTSIDE ) {
        assert( n_local + finite_element_nodes < node_connector_.size() );
        node_connector_[ n_local + finite_element_nodes ] = nptr;
        return;
     }

  if ( side == MIDDLE )
    throw csmp::Exception( ERROR, "csmp::InterFace<dim>::Assign( uint32_t, Node, INTERFACE_SIDE )",
                          "Called wrong method to assign intervening element!" );

} // end Assign node pointers




/**
    Unassigns the neighbor face, setting the pointer in the 'face_connector' vector to null
    @return whether removal was successful
*/
template<uint32_t dim>
bool InterFace<dim>::Unassign( const InterFace<dim>* const if_ptr )
  {
    if ( if_ptr == nullptr ) return false;
    const auto n_nbors{ static_cast<uint32_t>(interface_connector_.size()) };
    for ( uint32_t i{0U}; i < n_nbors; i++ )
      if ( if_ptr == interface_connector_[i] ) {
          interface_connector_[i] = nullptr;
          return true;
        }
    return false;
  }


template<uint32_t dim>
void InterFace<dim>::UnassignNeighbor( uint32_t nbor )
  {
     assert( nbor < interface_connector_.size() );
     interface_connector_[nbor] = nullptr;
  }




/**
    return the number of the Interface object neighbors which are connected with the InterFace and not null.
*/
template<uint32_t dim>
uint32_t InterFace<dim>::ConnectedNeighbors() const noexcept
{
  uint32_t connections = static_cast<uint32_t>(interface_connector_.size());
  for ( auto& f : interface_connector_ )
    if ( f == nullptr ) connections--;
  return connections;
}



  /// node_connector_.size() = total nodes on both sides of InterFace
template<uint32_t dim>
uint32_t  InterFace<dim>::Nodes() const noexcept
 { 
    assert( this->FE()!=nullptr ); 
    assert( (this->FE()->Nodes()*2) == node_connector_.size() );
    return static_cast<uint32_t>(node_connector_.size());
 }


  
/**
   Returns the total number of neighbors of the InterFace , which is equivalent
   to the number of interface finite element faces although either side of the Interface
   has different nodes.
*/
template<uint32_t dim>
uint32_t  InterFace<dim>::Neighbors() const noexcept
  {
     assert( this->FE()!=nullptr ); 
     assert( this->FE()->Neighbors() == interface_connector_.size() );
     return static_cast<uint32_t>(interface_connector_.size());
  }



/**
    Again this number is the same as that of the number of faces of the underlying finite element
    inspite of the fact that the InterFace has two sides.
*/
template<uint32_t dim>
uint32_t  InterFace<dim>::Faces() const noexcept
 { 
    assert( this->FE()!=nullptr ); 
    assert( this->FE()->Faces() == interface_connector_.size() );
    return static_cast<uint32_t>(interface_connector_.size());
 }



/**
    Finds the local numbers of the faces of the higher-dimensional element that will be connected by this interface;
    
    uses point coordinates that must be matched across the interface to find the nodes.
*/
template<uint32_t dim>
pair<uint32_t,uint32_t>  InterFace<dim>::SharedElementFacesAndFaceIDs()
{
  assert( innerParent_ != nullptr );
  assert( outerParent_ != nullptr );

  // 1. building search maps that we will use to find the shared interfaces
  //    key=pointset   face iD
  map<set<Point<dim> >, pair<INTERFACE_SIDE,uint32_t> >   inner_elmt_faces, outer_elmt_faces;
  // inner parent element
  const auto ifaces(innerParent_->Faces());
  for ( auto face{0U}; face < ifaces; ++face ) {
      set<Point<dim> >  face_key;
      for ( const auto& j : innerParent_->CornerNodesOfFace(face) )
        face_key.insert( j->Coordinate() );
      inner_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
    }
  // outer parent element
  const uint32_t ofaces(outerParent_->Faces());
  for ( uint32_t face = 0U; face < ofaces; ++face ) {
    set<Point<dim> >  face_key;
    for ( const auto& j : outerParent_->CornerNodesOfFace(face) )
      face_key.insert( j->Coordinate() );
    outer_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
  }

  // 2. finding the shared faces
  bool found( false );
  uint32_t inner_face_id, outer_face_id;
  for ( auto& inner_face : inner_elmt_faces ) {
    for ( auto& outer_face : outer_elmt_faces ) {
      // compares the sets of the point coordinates of potentially opposing faces
      if ( inner_face.first == outer_face.first ) {
        inner_face_id = inner_face.second.second;
        outer_face_id = outer_face.second.second;
        found = true;
        break;
      }
    }
    if ( found ) break;
  }
  return make_pair( inner_face_id, outer_face_id );

} // end SharedElementFacesAndFaceIDs




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
template<uint32_t dim>
void InterFace<dim>::InitialiseNodeVector()
{
  assert( innerParent_ );
  assert( outerParent_ );
  
  // resizing the nodevector
  if ( node_connector_.empty() ) {
      if constexpr ( dim == 1U ) node_connector_.resize( 2U );
      else node_connector_.resize( this->FE()->Nodes() * 2U );
    }
  node_connector_.shrink_to_fit();

   const auto n_face_nodes_inner{ innerParent_->FE()->NodesPerFace(inner_parent_face_id_) };
   const auto n_face_nodes_outer{ outerParent_->FE()->NodesPerFace(outer_parent_face_id_) };
   assert( n_face_nodes_inner == n_face_nodes_outer );
   
   //     4.1 Inside face: straightforward assignment from face indices
   uint32_t node_count{0U};
   for( const auto& i : innerParent_->FE()->NodesOfFace( inner_parent_face_id_) )
     Assign( node_count++, innerParent_->N(i), INSIDE );

   //    4.2 Outside face: find correct circular permutation of face indices
   //    which will preserve Node collocation
   auto nids = outerParent_->FE()->NodesOfFace( outer_parent_face_id_ );
   node_count = 0U;

   const auto n_corner_nodes_outer{ this->FE()->CornerNodes() };

// TODO: SKM comment: is there no better way than to rely on Point comparitor which can be unsafe
//       perhaps use lexicographical compare operator< of std::array
   while( outerParent_->N( nids[0] )->Coordinate() != N( n_corner_nodes_outer - 1 , INSIDE )->Coordinate() && node_count < n_face_nodes_outer ) {
         rotate( nids.begin(), nids.begin()+1, nids.end() );
         node_count++;
     }
   if ( node_count == nids.size() )
     throw Exception( ERROR, "InterFace::InitialiseNodeVector", "No circular permutation found from INSIDE face to OUTSIDE" );

   node_count = 0U;
   for ( auto i : nids )
     Assign( node_count++, outerParent_->N(i), OUTSIDE );

} // end InitialiseNodeVector






template<uint32_t dim>
void  InterFace<dim>::CurrentSide( INTERFACE_SIDE side_to_assign ) noexcept
{
  current_side_ = side_to_assign;
}


// ACCESSORS
template<uint32_t dim>
void  InterFace<dim>::Idx( size_t idx_to_assign ) const noexcept
{
  idx_ = idx_to_assign;
}

template<uint32_t dim>
size_t  InterFace<dim>::Idx() const noexcept
{
  return idx_;
}

template<uint32_t dim>
INTERFACE_SIDE  InterFace<dim>::CurrentSide() const noexcept
{
  return current_side_;
}


/**

@param offset offset of test_function operand in the equation with the primary variable indicated by top_key
Usage

@code
    std::vector<size_t> eq_idx_vec;

    cell.ActiveEquationIndices(
        DOF_indexes_,
        prop_key,
        offset,
        [&](size_t e) { eq_idx_vec.push_back(e); }
    );

@endcode
*/
template<uint32_t dim>
template<class Inserter>
void InterFace<dim>::ActiveEquationIndices( const csmp::Index& top_key,
                                            const vector<size_t>& DOF_indexes, size_t offset,
                                            Inserter& eq_idx_vec ) const
{
    // MAX_LOCAL_DOF bounds our stack allocation. 
    // 128 is sufficient for standard higher-order volumetric elements.
    constexpr size_t MAX_LOCAL_DOF = 128; 
    std::array<size_t, MAX_LOCAL_DOF> local_eqs;
    size_t local_count = 0;

    // Lambda to safely push to the local stack buffer
    auto push = [&](size_t global_dof) {
        size_t eq = DOF_indexes[global_dof];
        if (eq != NULL_IDX) {
            // Optional: assert(local_count < MAX_LOCAL_DOF) for debug builds
            local_eqs[local_count++] = eq;
        }
    };
    
    // 1. Gather all equations for this element
    for (const Node<dim>* node : node_connector_)
    {
        const size_t base = node->Idx();

        switch (top_key.type)
        {
            case SCALAR:
                push(base + offset);
                break;

            case VECTOR:
                for (uint32_t i = 0; i < dim; ++i)
                    push(base * dim + i + offset);
                break;

            case TENSOR:
                for (uint32_t i = 0; i < dim; ++i)
                    for (uint32_t j = 0; j < dim; ++j)
                        push(base * dim * dim + i * dim + j + offset);
                break;

            case ARRAY:
            case FLAGGEDARRAY:
                for (uint32_t i = 0; i < top_key.dataDepth; ++i)
                    push(base * top_key.dataDepth + i + offset);
                break;
        }
    }
    
    // 2. Sort and deduplicate equations specific to this element
    auto begin_it = local_eqs.begin();
    auto end_it = begin_it + local_count;
    
    std::sort(begin_it, end_it);
    auto unique_end = std::unique(begin_it, end_it);

    // 3. Push the strictly unique equations to the global inserter
    for (auto it = begin_it; it != unique_end; ++it) {
        eq_idx_vec.push_back(*it); 
    }
    
} // end ActiveEquationIndices

template void InterFace<1>::ActiveEquationIndices( const csmp::Index&, const vector<size_t>&, size_t, vector<size_t>& ) const;
template void InterFace<2>::ActiveEquationIndices( const csmp::Index&, const vector<size_t>&, size_t, vector<size_t>& ) const;
template void InterFace<3>::ActiveEquationIndices( const csmp::Index&, const vector<size_t>&, size_t, vector<size_t>& ) const;



template<uint32_t dim>
vector<csmp::InterFace<dim>*>&  InterFace<dim>::NeighborElementVector()
{
  return interface_connector_;
}



/**
      Performs test using the distanceTo operator on the node points
*/
template<uint32_t dim>
bool  InterFace<dim>::AreNodesCollocated(double tolerance) const
 {
    const uint32_t n_nodes{ this->FE()->Nodes() };
    for ( uint32_t i{0U}; i<n_nodes; i++ )
      if ( !approximatelyEqual( distance( this->MatchingN(i,INSIDE)->Coordinate(), this->MatchingN(i,OUTSIDE)->Coordinate() ), 0., tolerance ) )
        return false;
      
    return true;
 }



/**
    END_POINT is a  classifier that applies on the perimeter of SplitBoundary objects (perimeter InterFaces)
    terminating within models where INSIDE and OUTSIDE nodes are identical.
    
    @param n_local either a node on the inside or on the outside of the interface,
    which is the same as that on the opposite side if the node is on the perimeter (3D) or at a free-standing end point (2D) of a SplitBoundary.
*/
template<uint32_t dim>
bool  InterFace<dim>::IsEndPointNode( uint32_t n_local ) const
 {
    assert( n_local < node_connector_.size() );
    const auto n_nodes = this->FE()->Nodes();
    
    // if this is an inside node
    if ( n_local < n_nodes ) {
         // assuming that all the node pointers are valid
         assert( this->N(n_local) != nullptr );
         assert( this->N(n_nodes - n_local) != nullptr );
         if ( this->N(n_local) ==  this->N(n_nodes - n_local) )
           return true;
      }

    // if this is an outside node
    if ( n_local >= n_nodes ) {
         // assuming that all the node pointers are valid
         assert( this->N(n_local) != nullptr );
         assert( this->N(n_local - n_nodes) != nullptr );
         if ( this->N(n_local) ==  this->N(n_local - n_nodes) )
           return true;
      }

    return false;
    
 } // end




/**
    Returns pointers to the nodes on either side of the Interface.

    @section input Input Arguments

    An integer from 0...n-1, where n is the number of nodes per face of the Element.
    These are numbered counterclockwise looking from the outside into the face
    of parent element 1.

    @param side  side refers to the first or second parent element.

    @section implementation Implementation
    
    @attention since the nodes match the face of of the adjacent higher-dimensional elements,
    they are numbered like these within the node container. It follows that the inside nodes in the
    node connector are in normal order, but the ones for the outside are in reverse order starting
    with the last node. This is taking into account when they are returned by this method.

    A range check is performed.

    @return A pointer to the Target node.
*/
template<uint32_t dim>
csmp::Node<dim>* const InterFace<dim>::N( uint32_t n, INTERFACE_SIDE side ) const
{
  assert( n < this->FE()->Nodes() );
  if ( side == INSIDE ) return node_connector_[n];

  if ( side == OUTSIDE ) {
      // the number of nodes on a single side of the interface
      // SKM fix to pass InterFace_Test
      // const uint32_t outside_idx = static_cast<uint32_t>(node_connector_.size()) - 1U - n;
      const uint32_t outside_idx = n + this->FE()->Nodes();
      return node_connector_[outside_idx];
    }

  assert( side == MIDDLE );
  if ( middleElement_ != nullptr )
    return middleElement_->N( n );

  throw csmp::Exception( ERROR, "InterFace<dim>::N( local_id, side ) const", "Base Element does not exist!" );

  return nullptr;
}




/**
    Access to nodes on the CURRENT_SIDE of the interface. Taken from member current_side_
    @param n_local must be the local node index counting from 0 to the number of nodes on one side of the element
*/
template<uint32_t dim>
csmp::Node<dim>* const InterFace<dim>::N( uint32_t n_local ) const noexcept
{
  assert( n_local < this->FE()->Nodes() );
  return this->N(n_local,current_side_);
}



/**
    Returns pointers to the node which match with the INSIDE ordering of the interface

    @section input Input Arguments

    An integer from 0...n-1, where n is the number of nodes per face of the Element adjacent to the InterFace.
    The nodes on the outside match with (are collocated with) the nodes on the INSIDE, and therefore they no
    longer reflect the numbering given by the outer parent elemnts Face.

    @param side  side refers to the first or second parent element.

    @section implementation Implementation

    @attention since the nodes match the face nodes of of the adjacent higher-dimensional elements,
    they are numbered like these within the node container. It follows that the inside nodes in the
    node connector are in normal order, but the ones for the outside are in reverse order starting
    with the last node. Consequently, this method traverses the outside nodes in a reverse order, in order
    to output Nodes on the outside which are matched with nodes on the inside.

    @warning The ordering on the outside is INCONSISTENT with the ordering given from N(i,outside)!

    A range check is performed.

    @return A pointer to the Target node.
*/
template<uint32_t dim>
csmp::Node<dim>* const InterFace<dim>::MatchingN( uint32_t n, INTERFACE_SIDE side ) const
{
  assert( n < this->FE()->Nodes() );
  if ( side == INSIDE ) return node_connector_[n];

  uint32_t const cn_nodes = this->FE()->CornerNodes();
  uint32_t const fe_nodes = this->FE()->Nodes();
  if ( side == OUTSIDE ) {
    if (  n  < cn_nodes ){
      //Traverse nodes backwards from the last corner node
      const uint32_t outside_idx = fe_nodes + cn_nodes - 1 - n;
      assert(outside_idx >= fe_nodes );
      return node_connector_[outside_idx];
    }

    //Then we are on the midside nodes
    uint32_t one{1}, md_nodes = this->FE()->MidSideNodes();
    if (n < cn_nodes + md_nodes ){
     //Traverse the midside  nodes in reverse, but starting one node before the last node
      const uint32_t outside_idx = fe_nodes + cn_nodes + md_nodes - 1 - static_cast<uint32_t>(one % md_nodes) - (n-cn_nodes);
      assert(outside_idx >= fe_nodes );
      return node_connector_[outside_idx];
    } else {
      //this is a barycentric node
      const uint32_t outside_idx = n + fe_nodes;
      return node_connector_[outside_idx];
    }
  }

  assert( side == MIDDLE );
  if ( middleElement_ != nullptr )
    return middleElement_->N( n );

  throw csmp::Exception( ERROR, "InterFace<dim>::N( local_id, side ) const", "Intervening Element does not exist!" );

  return nullptr;
}







/**
    Returns the equal dimensional neighbor of the InterFace which also is an interface element.
*/
template<uint32_t dim>
csmp::InterFace<dim>* const InterFace<dim>::Neighbor( uint32_t n ) const noexcept
{
  assert( interface_connector_.size() == this->Neighbors() );
  assert( n < this->Neighbors() );
  return interface_connector_[n];
}





/**
    Returns which local node number the current node corresponds to in the given
    parent element.
 
   @note costly method, use judiciously.
*/
template<uint32_t dim>
uint32_t  InterFace<dim>::ParentNodeNumber( uint32_t n, INTERFACE_SIDE side ) const
{
  assert( n < node_connector_.size() );

  if ( side == INSIDE ) {
    // finding the parent element that corresponds to inside element
    for ( auto i{0U}; i<node_connector_[n]->Parents(); ++i )
      if ( node_connector_[n]->Parent( i ) == innerParent_ )
        return node_connector_[n]->ParentNodeNumber( i );
  }
  else if ( side == OUTSIDE ) {
    for ( auto i{0U}; i<node_connector_[n]->Parents(); ++i )
      if ( node_connector_[n]->Parent( i ) == outerParent_ )
        return node_connector_[n]->ParentNodeNumber( i );
  }
  else if ( side == MIDDLE )
    return n;

  throw csmp::Exception( ERROR, "InterFace<dim>::ParentNodeNumber:",
                         "Node does not seem to be connected to parent element." );
  return 0;
}




// watch out if there is no base element this returns a nullptr pointer
template<uint32_t dim>
Element<dim>* const InterFace<dim>::Parent( INTERFACE_SIDE side ) const noexcept
{
  if ( side == INSIDE )
    return innerParent_;
  else if ( side == OUTSIDE )
    return outerParent_;
  return middleElement_;
}

template<uint32_t dim>
Element<dim>* const InterFace<dim>::InnerParent() const noexcept
{
  return innerParent_;
}
template<uint32_t dim>
Element<dim>* const InterFace<dim>::OuterParent() const noexcept
{
  return outerParent_;
}

template<uint32_t dim>
Element<dim>*  InterFace<dim>::InterveningElement() noexcept
{
  return middleElement_;
}

template<uint32_t dim>
const Element<dim>*  InterFace<dim>::InterveningElement() const noexcept
{
  return middleElement_;
}


/// return face ID of inner parent element
template<uint32_t dim>
uint32_t  InterFace<dim>::InnerParentFaceID() const noexcept
{
  if ( innerParent_ == nullptr ) return NULL_IDX32U;
  return inner_parent_face_id_;
}


/// return face ID of outer parent element
template<uint32_t dim>
uint32_t  InterFace<dim>::OuterParentFaceID() const noexcept
{
  if ( outerParent_ == nullptr ) return NULL_IDX32U;
  return outer_parent_face_id_;
}


/// return face ID of inner or outer parent element depending on index
template<uint32_t dim>
uint32_t  InterFace<dim>::ParentFaceID( INTERFACE_SIDE side ) const
{
  if ( side == INSIDE ) {
      assert( innerParent_ != nullptr );
      return inner_parent_face_id_;
    }
  else if ( side == OUTSIDE ) {
      assert( outerParent_ != nullptr );
      return outer_parent_face_id_;
    }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.Note( WARNING, "csmp::InterFace<dim>::ParentFaceID:", "Interface 'side' could not be determined." );

  return inner_parent_face_id_;
}




/// assign face ID of inner or outer parent element depending on their local face numbering 
template<uint32_t dim>
void  InterFace<dim>::ParentFaceID( INTERFACE_SIDE side, uint32_t idx )
{
  if ( side == INSIDE ) {
      assert( innerParent_ != nullptr );
      inner_parent_face_id_ = idx;
      return;
    }
  else if ( side == OUTSIDE ) {
      assert( outerParent_ != nullptr );
      outer_parent_face_id_ = idx;
      return;
    }

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.Note( WARNING, "csmp::InterFace<dim>::ParentFaceID:", "Interface 'side' could not be determined." );

} // end ParentFaceID(assignment)






// GEOMETRY

inline double triangleArea( const Point<1U>&, const Point<1U>&, const Point<1U>& ) {
    cerr <<"\nInterFace<1>::Area() called triangleArea = function stub that is meaningless in 1D.\n";
    return 1.;
 }

inline Point<1U> normalOfTriangle( const Point<1U>&, const Point<1U>&, const Point<1U>& ) {
   cerr <<"\nInterFace<1>::UnitNormal() called normalOfTriangle = function stub that is meaningless in 1D.\n";
   return 1.;
}

/**
    Computes the interface area from scratch and not using the CoordinateMatrix / InterFace FEM machinery.
    The area is computed taking the side of the interface into account.
    Thus, it will give different results for INSIDE and OUTSIDE if the nodes on either side no longer match.
*/
template<uint32_t dim>
double InterFace<dim>::Area( INTERFACE_SIDE side ) const
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );

   if ( side == MIDDLE ) {
        // Volume() internally computes the coordinate matrix for the middle element
        if ( middleElement_ != nullptr ) return middleElement_->Volume();
        csmp_error.Note( ERROR, "InterFace<dim>::Area:", "InterFace FE type not recognized." );
        return numeric_limits<double>::signaling_NaN();
     }
   else {
        // this uses the Coordinate matrix initialised here
        NodeCoordinateMatrix( this->FE()->XY, side );
        return this->FE()->Volume();
     }

   // error
   return numeric_limits<double>::signaling_NaN();
}



/**
    Returns the outward-pointing unit normal to the InterFace pointing in the direction of the OUTSIDE.
      Per default, the normal is sourced from the MIDDLE element.
      If it is not there, the coordinate matrix is initalized with the averages of inside and outside node coordinates.
*/
template<uint32_t dim>
csmp::Point<dim>  InterFace<dim>::UnitNormal() const
{
   if ( middleElement_ ) return middleElement_->UnitNormal();
   // fall-back option
   BisectorCoordinateMatrix();
   return Point<dim>( this->FE()->UnitNormal() );
}




/**
   returns the normal pointing from the inside to the outside higher-dimensional Element of the InterFace, calculated for bisector plane.
*/
template<uint32_t dim>
csmp::Point<dim>  InterFace<dim>::UnitNormal( INTERFACE_SIDE side ) const
 {
    if ( side == MIDDLE && middleElement_ != nullptr )
      return middleElement_->UnitNormal();
    
    NodeCoordinateMatrix( this->FE()->XY, side );
    return Point<dim>( this->FE()->UnitNormal() );

 } // end UnitNormal




/**
Initialises VectorVariable with vector between node pair - returns false if overlap, true if distant
*/
template<uint32_t dim>
double InterFace<dim>::NodeSpacing( uint32_t n ) const
  {
    const uint32_t n_nodes{ this->FE()->Nodes() };
    assert( n < n_nodes );
    return distance( MatchingN(n,INSIDE)->Coordinate(), MatchingN(n,OUTSIDE)->Coordinate() );
  }





/**

Returns a matrix with 'nodes'-rows and 'coordinate-directions' columns.
This matrix defines the positions of the elements nodes for
the finite-element matrix assembly. Since the number of element nodes
may vary among different elements types, the number of rows in XY may
also vary from element to element.

@param XY is the DenseMatrix<DM_MIN> class object (value type fT) stored in FiniteElement.
This matrix is dynamically resized if necessary but must have been constructed with a finite size
before passing it to CoordinateMatrix().

The node coordinates are returned into the supplied matrix.

@section application Application

Method will be called by FiniteElementPolicy to initialise XY matrix inside the finite element

@attention when INTERFACE_SIDE == MIDDLE, the node locations on either side of the interface
are used to find mid-points.

*/
template<uint32_t dim>
void  InterFace<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
{
   if ( current_side_ == MIDDLE ){
     BisectorCoordinateMatrix(XY);
     return;
   }

   NodeCoordinateMatrix( XY, current_side_ );

} // end NodeCoordinateMatrix


template<uint32_t dim>
void  InterFace<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY, INTERFACE_SIDE side ) const
{
 
  if ( side == MIDDLE ){
      BisectorCoordinateMatrix(XY);
      return;
  }


  const auto n_nodes( this->FE()->Nodes() );
  XY.Resize( n_nodes, dim );

  for ( uint32_t i{0U}; i<n_nodes; ++i )
    XY.AssignRow( i, N( i, side )->Coordinate() );

} // end NodeCoordinateMatrix


/**
       Interface bisector plane.
*/
template<uint32_t dim>
void  InterFace<dim>::BisectorCoordinateMatrix() const
{
  const auto n_nodes( this->FE()->Nodes() );
  this->FE()->XY.Resize( n_nodes, dim );

  for ( uint32_t i{0U}; i<n_nodes; ++i ) {
       const Point<dim> mid_point = (this->MatchingN( i, INSIDE )->Coordinate() + this->MatchingN( i, OUTSIDE )->Coordinate()) / 2.;
       this->FE()->XY.AssignRow( i, mid_point );
    }

} // end CoordinateMatrix


/**
       Interface bisector with XY Dense matrix to be overwritten.
*/
template<uint32_t dim>
void  InterFace<dim>::BisectorCoordinateMatrix(DenseMatrix<DM_MIN>& XY) const
{
  const auto n_nodes( this->FE()->Nodes() );
  XY.Resize( n_nodes, dim );

  for ( uint32_t i{0U}; i<n_nodes; ++i ) {
       const Point<dim> mid_point = (this->MatchingN( i, INSIDE )->Coordinate() + this->MatchingN( i, OUTSIDE )->Coordinate()) / 2.;
       XY.AssignRow( i, mid_point );
    }

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
Point<dim>  InterFace<dim>::BaryCenter() const
{
  Point<dim>  pt; // initialised to zero

  // all the nodes on both sides
  for ( const auto& n : node_connector_ )
    pt += n->Coordinate();

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
double  InterFace<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
{
  double fMinTemp( static_cast<double>(DBL_MAX) );
  double fMaxTemp( static_cast<double>(-DBL_MAX) );

  // this normalisation is necessary because the vector variable
  // being any physical quantity may have any magnitude
  const double fMagnitudeOfDirection( vecDirection.Length() );
  // avoid division by zero
  assert( fMagnitudeOfDirection >= numeric_limits<double>::epsilon() );

  for ( const auto& n : node_connector_ ) {
    // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
    double fTemp( vecDirection.DotProduct( n->Coordinate() ) );
    fTemp /= fMagnitudeOfDirection;

    // update minimum value
    fMinTemp = min( fMinTemp, fTemp );
    // update maximum value
    fMaxTemp = max( fMaxTemp, fTemp );
  }

  //substract magnitudes
  return fMaxTemp - fMinTemp;
}




/**
    returns property values at the nodes.
    
         @attention the order of the nodes corresponds to the numbering of the face of the corresponding higher dimensional element 
         on the in- or outside. Unless the return_nodes_outside_that_match_inside == true
*/
template<uint32_t dim>
template< class Var>
void  InterFace<dim>::NodePropertyVector( const csmp::Index& idx, vector<Var>& V, INTERFACE_SIDE side ) const
{
  if ( idx.place != NODE ) {
    cerr << "\nInterFace<" << dim;
    cerr << ">::NodePropertyVector: Requested property ";
    cerr << "is not placed on the nodes; property Index: " << endl;
    idx.Out();
    return;
  }

  // resizing V if necessary
  const auto  n_nodes( this->FE()->Nodes() );
  V.resize( n_nodes );

  for ( uint32_t i{0U}; i<n_nodes; i++ )
    N( i, side )->Read( idx, V[i] ); // OUTSIDE nodes correspond to ordering of face of higher dim parent
}


// scalar
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
// vector
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, vector<VectorVariable<1U> >&, INTERFACE_SIDE) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, vector<VectorVariable<2U> >&, INTERFACE_SIDE) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, vector<VectorVariable<3U> >&, INTERFACE_SIDE) const;
// tensor
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, vector<TensorVariable<1U> >&, INTERFACE_SIDE) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, vector<TensorVariable<2U> >&, INTERFACE_SIDE) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, vector<TensorVariable<3U> >&, INTERFACE_SIDE) const;
// array
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
// flagged array
template void  InterFace<1U>::NodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::NodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::NodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;


/**
    returns property values at the nodes.

         @attention the order of the nodes corresponds to the numbering of the face of the corresponding higher dimensional element
         on the in- or outside. Unless the return_nodes_outside_that_match_inside == true
*/
template<uint32_t dim>
template< class Var>
void  InterFace<dim>::MatchingNodePropertyVector( const csmp::Index& idx, vector<Var>& V, INTERFACE_SIDE side ) const
{
  if ( idx.place != NODE ) {
    cerr << "\nInterFace<" << dim;
    cerr << ">::NodePropertyVector: Requested property ";
    cerr << "is not placed on the nodes; property Index: " << endl;
    idx.Out();
    return;
  }

  // resizing V if necessary
  const auto  n_nodes( this->FE()->Nodes() );
  V.resize( n_nodes );

  for ( uint32_t i{0U}; i<n_nodes; i++ )
    MatchingN( i, side )->Read( idx, V[i] );              //OUTSIDE nodes correspond to ordering of face of higher dim parent
}




// scalar
template void  InterFace<1U>::MatchingNodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::MatchingNodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::MatchingNodePropertyVector( const csmp::Index&, vector<ScalarVariable>&, INTERFACE_SIDE) const;
// vector
template void  InterFace<1U>::MatchingNodePropertyVector( const csmp::Index&, vector<VectorVariable<1U> >&, INTERFACE_SIDE) const;
template void  InterFace<2U>::MatchingNodePropertyVector( const csmp::Index&, vector<VectorVariable<2U> >&, INTERFACE_SIDE) const;
template void  InterFace<3U>::MatchingNodePropertyVector( const csmp::Index&, vector<VectorVariable<3U> >&, INTERFACE_SIDE) const;
// tensor
template void  InterFace<1U>::MatchingNodePropertyVector( const csmp::Index&, vector<TensorVariable<1U> >&, INTERFACE_SIDE) const;
template void  InterFace<2U>::MatchingNodePropertyVector( const csmp::Index&, vector<TensorVariable<2U> >&, INTERFACE_SIDE) const;
template void  InterFace<3U>::MatchingNodePropertyVector( const csmp::Index&, vector<TensorVariable<3U> >&, INTERFACE_SIDE) const;
// array
template void  InterFace<1U>::MatchingNodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::MatchingNodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::MatchingNodePropertyVector( const csmp::Index&, vector<ArrayVariable>&, INTERFACE_SIDE) const;
// flagged array
template void  InterFace<1U>::MatchingNodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<2U>::MatchingNodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;
template void  InterFace<3U>::MatchingNodePropertyVector( const csmp::Index&, vector<FlaggedArrayVariable>&, INTERFACE_SIDE) const;









// SCREEN OUTPUT
template<uint32_t dim>
void  InterFace<dim>::Out() const
{
  cout << "\nInterFace<" << dim << ">::Out: number: " << idx_;
  cout <<": "<< parseFiniteElementType( this->FE_Type() );
  cout << "\n\tconnected inside nodes with boundary flags:   ";
  string str;
  for ( uint32_t i{0U}; i<this->FE()->Nodes(); i++ ) {
    str = parseBoundary( N( i, INSIDE )->AtBoundary() );
    cout << N( i, INSIDE )->Idx() << ":" << str << "  ";
  }
  cout << "\n\tconnected outside nodes with boundary flags:  ";
  for ( uint32_t i{0U}; i<this->FE()->Nodes(); i++ ) {
    str = parseBoundary( N( i, OUTSIDE )->AtBoundary() );
    cout << N( i, OUTSIDE )->Idx() << ":" << str << "  ";
  }
  cout << endl;

  cout << "\n\tconnected neighbor InterFace types / boundary flags: ";
  for ( uint32_t i{0U}; i<this->Neighbors(); i++ )
    if ( Neighbor( i ) != nullptr ) {
        cout << "\t\t" << Idx() << ": ";
        cout << parseFiniteElementType( Neighbor( i )->FE_Type() ) << ": ";
        //str = parseBoundary(Neighbor(i)->AtBoundary());
        //cout << str;
        cout << endl;
      }
    else cout << "none ";
    cout << endl;

    Point<dim>  pt( this->BaryCenter() );

    if ( dim == 1U )
      cout << "\n\tBarycentre at (xyz): " << pt[0] << endl;
    else if ( dim == 2U )
      cout << "n\tBarycentre at (xyz): " << pt[0] << ", " << pt[1] << endl;
    else
      cout << "\n\tBarycentre at (xyz): " << pt[0] << ", " << pt[1] << ", " << pt[2] << endl;

    const auto ipoints( this->IntegrationPoints() );
    if ( ipoints > 0U ) {
      cout << "\n\tStorage sites for IntegrationPoint properties: " << ipoints << endl;
    }

    cout <<"\n\t"<<"Parent (higher-dimensional) Element objects:\n";
    if ( innerParent_ != nullptr ) {
         cout <<"\t\t"<<"inside "<< dim <<"-dimensional parent Element: "<< this->innerParent_->Idx();
         cout  <<": "<< parseFiniteElementType(this->Parent(INSIDE)->FE_Type()) << endl;
      }
    else cout <<"\t\tnone.\n";
    if ( this->outerParent_ != nullptr ) {
         cout <<"\t\t"<<"outside "<< dim <<"-dimensional parent Element: "<< this->outerParent_->Idx();
         cout <<": "<< parseFiniteElementType(this->Parent(OUTSIDE)->FE_Type()) << endl;
      }
    else cout << "\t\tnone.\n";
    if ( this->middleElement_ != nullptr ) {
         cout <<"\t\t"<<"intervening "<< dim-1 <<"-dimensional Element: "<< this->middleElement_->Idx();
         cout <<": "<< parseFiniteElementType(this->Parent(MIDDLE)->FE_Type()) << endl;
      }
    else cout << "\t\tnone.\n";

    cout << "\n\tUnit normal:  ";
    Point<dim> un = UnitNormal();
    for ( uint32_t i{0U}; i<dim; i++ ) cout << un[i] << ", ";
    cout << endl;

} // end Out


template class InterFace<1U>;
template class InterFace<2U>;
template class InterFace<3U>;

} // end namespace csmp
