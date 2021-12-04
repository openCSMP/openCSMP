#include "Node.h"
#include "Element.h"
#include "Visitor.h"
#include "NodeManifold.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
    default constructor
    
    Node coordinate is initialised to origin (0,0,0)
*/
template<size_t dim>
Node<dim>::Node()
    : idx_(ULONG_MAX),
      at_boundary_(NOT)
  {
  }



/**
    Initialises everything except for parent element related vectors.
*/
template<size_t dim>
Node<dim>::Node( size_t idx, const Point<dim>& pt, const LocalVariables& lvs, BOX_BOUNDARY boundary_flag )
 : LocalVariableStorage<dim,Node>(lvs),
   xyz_(pt),
   idx_(idx),
   at_boundary_(boundary_flag)
 {
 }



/**
    Copy constructor also copies the pointer assignments (!).
*/
template<size_t dim>
Node<dim>::Node( const Node<dim>& nd )
  : xyz_(nd.xyz_), idx_(nd.idx_),
    parent_node_indexes_(nd.parent_node_indexes_),
    parent_element_pointers_(nd.parent_element_pointers_),
    at_boundary_(nd.at_boundary_)
  {
    this->LVS( nd.LVS() );
  }



/**
    Move constructor also copies the pointer assignments (!).
*/
template<size_t dim>
Node<dim>::Node( Node<dim>&& nd )
  : xyz_{ move(nd.xyz_) },
    idx_{nd.idx_},
    parent_node_indexes_{ move(nd.parent_node_indexes_) },
    parent_element_pointers_{ move(nd.parent_element_pointers_) },
    at_boundary_{nd.at_boundary_}
  {
    this->LVS( move(nd.LVS()) );
  }




template<size_t dim>
Node<dim>::~Node()
 {
    if ( manifold_ != nullptr )
      manifold_->Remove( this );
      
//    cerr <<"\nNode "<< Idx() <<": called destructor.";
 }




template<size_t dim>
Node<dim>& Node<dim>::operator=( const Node<dim>& nd )
 {
    if ( &nd != this ) 
      {
         idx_                     = nd.idx_;
         at_boundary_             = nd.at_boundary_;
         parent_node_indexes_     = nd.parent_node_indexes_;
         parent_element_pointers_ = nd.parent_element_pointers_; 
         xyz_                     = nd.xyz_; 
         this->LVS( move( nd.LVS() ) );
      }
    return *this;
 }




/**
    Move assignment, relying on that similar operators exist for the nodes components.
    
    @note this assumes that the supplied node is a temporary.
*/
template<size_t dim>
Node<dim>& Node<dim>::operator=( Node<dim>&& nd )
 {
    assert( this != &nd );
    xyz_                     = nd.xyz_;
    idx_                     = nd.idx_;
    at_boundary_             = nd.at_boundary_;
    parent_node_indexes_     = move( nd.parent_node_indexes_ );
    parent_element_pointers_ = move( nd.parent_element_pointers_ );
    this->LVS( nd.LVS() );
 
    return *this;
 }




/**
    This comparison operator was designed specifically for the creation of particular Region, Boundary and SplitBoundary objects.
    Hence, only information important for this process is taken into account, distinguishing 2 Nodes from each other.
    That must be coordinate and parent elements

    @author Roman Manasipov
    @date   2014
*/
template<size_t dim>
bool  Node<dim>::operator==( const Node<dim>& nd )
 {
    // TODO: fix Node comparitor
    cerr <<"Node<"<< dim <<">: called operato== Node comparitor, extremely costly and of questionable value\n";
    if ( &nd != this )
    {
        if( Coordinate() == nd.Coordinate() )
        {
            if( parent_element_pointers_.size() == nd.parent_element_pointers_.size() )
            {
                if( parent_element_pointers_.empty() )
                {
                    if( idx_ == nd.idx_ )
                        return true;
                    return false;
                }

                std::set<Element<dim>*> elmts1(parent_element_pointers_.begin(), parent_element_pointers_.end());
                std::set<Element<dim>*> elmts2(nd.parent_element_pointers_.begin(), nd.parent_element_pointers_.end());
                std::vector<Element<dim>*> elmts_intersect;
                std::set_intersection( elmts1.begin(), elmts1.end(),
                                       elmts2.begin(), elmts2.end(),
                                       std::back_inserter(elmts_intersect) );
                if( elmts_intersect.size() == parent_element_pointers_.size() )
                    return true;

                return false;
            }
            return false;
        }
        return false;
    }
    return true;
 }
 
 
 // private operators

template<size_t dim>
void* Node<dim>::operator new( size_t size )
  {
//      std::cout<< "\nNode<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

template<size_t dim>
void Node<dim>::operator delete( void* p )
  {
//     std::cout<< "\nNode<"<< dim <<">: called overloaded delete operator.\n";
     free(p);
     p = nullptr;
  }


/**
 
Assigns parent element pointer and remembers which node (local idx)
of the parent element this node is.

@param pnode is the node number in the parent element of the assigned node.
@param element is pointer to the parent element that gets assigned.

@section arguments Input Arguments

The local number of the parent element which is going to be assigned.  
The local node number of this node in the parent element.  
A pointer to the parent Element which shall be added.  
*/
template<size_t dim>
void Node<dim>::Assign( size_t pnode, Element<dim>* element )
 {
    assert( parent_node_indexes_.size() == parent_element_pointers_.size() );
    assert( pnode <= FIFTY );
    
    for ( size_t parent{0}; parent<parent_node_indexes_.size(); parent++ )
      if ( parent_node_indexes_[parent] == NOT_INITIALIZED ) {
             parent_node_indexes_[parent]     = static_cast<ONE_BYTE_NUMBER>(pnode);
             parent_element_pointers_[parent] = element;
             return;
          }

 } // end Assign




/**
    Sets the pointer to the target element to zero, and the corresponding node number to NOT_INITIALIZED.
*/
template<size_t dim>
bool Node<dim>::Unassign( Element<dim>* element )
 {
    assert( parent_node_indexes_.size() ==  parent_element_pointers_.size() );
    for ( size_t parent(0); parent < Parents(); ++parent )
      if ( Parent(parent) == element )
        {
          parent_node_indexes_[parent]     = NOT_INITIALIZED;
          parent_element_pointers_[parent] = nullptr;
          return true;
        }
    return false;
    
 } // end Unassign




/**
    Removes parent elements pointers that were set to nullptr.
*/
template<size_t dim>
void Node<dim>::EraseNullPointerParents()
 {
    assert( parent_node_indexes_.size() ==  parent_element_pointers_.size() );
    const size_t  n_parents{Parents()};
    size_t        n_new_parents{0};

    for ( size_t parent{0}; parent < n_parents; ++parent ) {
         if ( parent_element_pointers_[parent] == nullptr )
           parent_node_indexes_[parent] = NOT_INITIALIZED;
         else n_new_parents++;
      }

    // if nullptr parents were detected, the parent storage needs to be rebuild
    if ( n_new_parents < n_parents ) {
        parent_element_pointers_.erase( remove( parent_element_pointers_.begin(),
                                                parent_element_pointers_.end(), nullptr ),
                                                parent_element_pointers_.end() );

        parent_node_indexes_.erase( remove( parent_node_indexes_.begin(),
                                            parent_node_indexes_.end(), NOT_INITIALIZED ),
                                            parent_node_indexes_.end() );
      }
    
 } // end EraseNullPointerParents




template<size_t dim>
void  Node<dim>::Accept( csmp::Visitor<dim>& v )
{
    v.Visit(this);
}


/**

   If the parent elements are of the same dimension as the node,
   then there are exactly as many elements as are nodes.
   
   @attention calculation is costly, compute only once per node loop!
   
   @test correct SKM (unit test exists)
*/
template<size_t dim>
size_t  Node<dim>::Neighbors() const
 {
    size_t  node_neighbors(parent_node_indexes_.size());
    
    // subtracting number of lower dimensional parent elements as these node sharing
    // but otherwise inconsequential elements would lead to a wrong node count
    if constexpr ( dim == 3U ) {
         for ( typename vector<Element<dim>*>::const_iterator
               it=parent_element_pointers_.begin(); it!=parent_element_pointers_.end(); ++it )
           if ( (*it)->IsSurfaceElement() || (*it)->IsLineElement() ) node_neighbors--;
         return node_neighbors;
      }
    if constexpr ( dim == 2U ) {
         for ( typename vector<Element<dim>*>::const_iterator
               it=parent_element_pointers_.begin(); it!=parent_element_pointers_.end(); it++ )
           if ( (*it)->IsLineElement() ) node_neighbors--;
         return node_neighbors;
      }
      
    // 1D version
    return node_neighbors;
 }


/**
    Implements node connectivity graph.
 
    Counter-clockwise node-numbering convention of the parent 
    elements is used.
    
    @test OK - unit test, see NodeNeighborConnectivity_Test
*/
template<size_t dim>
Node<dim>*  Node<dim>::Neighbor( size_t neighbor_node ) const
 {
    assert( neighbor_node < Neighbors() );
    assert( Parent( neighbor_node ) != nullptr );
   
    size_t next_node(ParentNodeNumber(neighbor_node) + 1U);
    if ( next_node == Parent( neighbor_node )->Nodes() ) next_node = 0U;
    
    return parent_element_pointers_[ neighbor_node ]->N( next_node );
 }






// ACCESSORS


template<size_t dim>
double  Node<dim>::operator[]( size_t i ) const { return xyz_[i]; }

template<size_t dim>
double&  Node<dim>::operator[]( size_t i ) { return xyz_[i]; }

template<size_t dim>
double&  Node<dim>::operator()( size_t i ) { return xyz_[i]; }

template<size_t dim>
Point<dim>  Node<dim>::Coordinate() const { return xyz_; }

template<size_t dim>
 void Node<dim>::Coordinate( const Point<dim>& p ) { xyz_=p; }



/**

The dynamic storage for the parent element data is resized preserving the
existing entries. The values of potential new elements are set to zero.

@param n The desired new size of the storage.
*/
template<size_t dim>
void  Node<dim>::ResizeParentStorage( size_t n )
  {
     parent_node_indexes_.resize( n, NOT_INITIALIZED );
     std::vector<ONE_BYTE_NUMBER>( parent_node_indexes_ ).swap( parent_node_indexes_ );
     parent_element_pointers_.resize( n, nullptr );
     std::vector<Element<dim>*>( parent_element_pointers_ ).swap( parent_element_pointers_ );
  }



template<size_t dim>
void  Node<dim>::EraseParents()
  {
     parent_element_pointers_.clear();
     parent_node_indexes_.clear();
  }



template<size_t dim>
size_t   Node<dim>::Parents() const
  { return parent_element_pointers_.size(); }



template<size_t dim>
size_t   Node<dim>::ParentNodeNumber( size_t parent_element_number ) const
  {
     assert( parent_element_number < parent_node_indexes_.size() );
     return parent_node_indexes_[ parent_element_number ];
  }


template<size_t dim>
Element<dim>*  Node<dim>::Parent( size_t parent_element_number ) const
  {
     assert( parent_element_number < parent_element_pointers_.size() );
     return parent_element_pointers_[ parent_element_number ];
  }


    /// checks whether Element is a parent of the node
template<size_t dim>
bool  Node<dim>::IsParent( const Element<dim>* const eptr ) const
  {
     return binary_search( parent_element_pointers_.begin(),
                           parent_element_pointers_.end(), eptr );
  }




/**
     sorts parent vector for searching and eliminates potential nullpointers
*/
template<size_t dim>
void  Node<dim>::UpdateParents() {
     assert( parent_element_pointers_.size() == parent_node_indexes_.size() );

     // sorting
     // creating indices
     vector<int> indices(parent_element_pointers_.size());
     iota( indices.begin(), indices.end(), 0 );
    
     // organising indices in the order that 'b' will have once it is sorted
     sort( indices.begin(), indices.end(),
           [&]( int i, int j ) -> bool {
                assert( parent_element_pointers_[i] != nullptr );
                return parent_element_pointers_[i] < parent_element_pointers_[j];
             }
        );
  
    // sorting the vectors
    sort( parent_element_pointers_.begin(), parent_element_pointers_.end() );
  
    // extra vector needed for tempory
    vector<ONE_BYTE_NUMBER> temp{ parent_node_indexes_.size() };
    int n{0};
    for ( auto i : indices ) temp[n++] = parent_node_indexes_[i];
    parent_node_indexes_ = temp;

  } // end UpdateParents





template<size_t dim>
void  Node<dim>::Idx( size_t idx_to_assign ) const
 { idx_ = idx_to_assign; }


template<size_t dim>
size_t   Node<dim>::Idx() const
 { return idx_; }



/**
   Set or return the global boundary flag of the Node.
*/
template<size_t dim>
void Node<dim>::AtBoundary( BOX_BOUNDARY b ) { at_boundary_ = b; }


template<size_t dim>
BOX_BOUNDARY  Node<dim>::AtBoundary() const { return at_boundary_; }

/// Set or return the coordinates of the current node.
template<size_t dim>
void            Node<dim>::x( double xc )  { xyz_[0u] = xc; }

template<size_t dim>
void            Node<dim>::y( double yc )  { xyz_[1u] = yc; }

template<size_t dim>
void            Node<dim>::z( double zc )  { xyz_[2u] = zc; }


template<size_t dim>
double          Node<dim>::x() const { return xyz_[0u]; }

template<size_t dim>
double          Node<dim>::y() const { return xyz_[1u]; }

template<size_t dim>
double          Node<dim>::z() const { return xyz_[2u]; }


// MANIFOLDS

/// access to manifold if any; returns nullptr if the node is not a manifold
template<size_t dim>
bool Node<dim>::IsManifold() const { return (manifold_ == nullptr); }


template<size_t dim>
NodeManifold<dim>* const Node<dim>::Manifold() const { return manifold_; }


template<size_t dim>
void Node<dim>::Assign( NodeManifold<dim>* const md )
 {
    assert( md != nullptr );
    manifold_ = md;
 }



// OUTPUT

/**
 
Outputs internal data of the node. 
*/
template<size_t dim>
void Node<dim>::Out() const
 {
    cout <<"\n\nNode<"<< dim <<">: "<< idx_;
    if ( at_boundary_ != NOT ) {
         string str(parseBoundary(at_boundary_));
         cout <<", Boundary flag: "<< str;
      }
    cout <<", Coordinates: "<< xyz_;
    cout << endl;
#ifndef NDEBUG
    if ( parent_node_indexes_.size() > 0u ) {
         cout <<"\nElement objects sharing the node / node position therein:\n"<< endl;
         for ( size_t i=0U; i<Parents(); i++ ) {
              if ( Parent(i) == NULL ) cout <<"NONE (null pointer) ";
              else Parent(i)->Out();
              cout <<"(node "<< ParentNodeNumber(i) <<"), ";
           }
         cout << endl;
      } 
#endif
    if ( manifold_ != nullptr ) {
         cout <<"\nconnected nodes: ";
         manifold_->Out();
      }
    
 } // end Out
 

template class Node<1U>;
template class Node<2U>;
template class Node<3U>;
 
 
 
 // NON-MEMBER FUNCTIONS
 
/**
   returns parent elements shared by nodes of face, inner side is reported first; outer next else application: give nodes of lower-dimensional face to find element on either side
   
   @todo method probably sill contains a large number of redundant operations.
*/
template<size_t dim>
pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( typename vector<Node<dim>*>::const_iterator first,
                                                               typename vector<Node<dim>*>::const_iterator last )
 {
    assert( first != last );
    // create sets of the parent elements of the face nodes checking which ones are shared
    const typename vector<Node<dim>*>::const_iterator nodesEnd{last};
    typename vector<Node<dim>*>::const_iterator nit{first};

    // creating a set of the parent elements of the first node
    assert( (*nit)->Parents() > 0 );
    const size_t n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;
    for ( size_t i{0}; i<n_parents; ++i ) {
         assert( (*nit)->Parent(i) != nullptr );
         if constexpr ( dim == 3 ) if ( !(*nit)->Parent(i)->IsVolumeElement() ) continue;
         if constexpr ( dim == 2 ) if ( !(*nit)->Parent(i)->IsSurfaceElement() ) continue;
         shared_parents.insert( (*nit)->Parent(i) );
      }
      
    // advancing the node iterator
    nit++;

    // searching for shared parent elements in the following nodes
    while ( nit != nodesEnd ) {
         set<Element<dim>*> temp;
         const size_t n_parents{(*nit)->Parents()};
         for ( size_t i{0}; i<n_parents; ++i ) {
              assert( (*nit)->Parent(i) != nullptr );
              if constexpr ( dim == 3 ) if ( !(*nit)->Parent(i)->IsVolumeElement() ) continue;
              if constexpr ( dim == 2 ) if ( !(*nit)->Parent(i)->IsSurfaceElement() ) continue;
              if ( shared_parents.find( (*nit)->Parent(i) ) != shared_parents.end() )
                temp.insert( (*nit)->Parent(i) );
           }
         shared_parents = temp;
         nit++;
      }
      
    // drawing the results together
    assert( !shared_parents.empty() );
    if ( shared_parents.size() == 1 ) return make_pair( (*shared_parents.begin()), nullptr );
    
    // if two parent elements were found, the one on the inside needs to be determined
    assert( shared_parents.size() == 2 );
    // initial guess
    pair<Element<dim>*,Element<dim>*> result( (*shared_parents.begin()), (*shared_parents.rbegin()) );
    
    // 3D case where face is either a triangle or a quadrilateral
    if constexpr ( dim == 3 ) {
         vector<Node<3>*> face_nodes( first, last );
         assert( face_nodes.size() >= 3 );
         // getting normal to face from the first 3 node coordinates
         Point<3> vec1(face_nodes[0]->Coordinate() - face_nodes[1]->Coordinate()); // cw
         Point<3> vec2(face_nodes[2]->Coordinate() - face_nodes[1]->Coordinate()); // ccw
         Point<3> nrml =  crossProduct( vec2, vec1 );
         // checking whether a vector from the second element barycentre to the first parent element center yields a negative or positive dot product
         Point<3> bvec = (*shared_parents.rbegin())->BaryCenter() - (*shared_parents.begin())->BaryCenter();
         // using dot-product to find inner element: if normal is pointing toward barycentre of first element, initial order needs to be reversed
         if ( dotProduct( nrml, bvec ) < 0. ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }
    
    // 2D case where the face is line and the non-existing normal points out of the plane
    if constexpr ( dim == 2 ) {
         vector<Node<2>*> face_nodes( first, last );
         assert( face_nodes.size() == 2 );
         Point<2> vec(face_nodes[1]->Coordinate() - face_nodes[0]->Coordinate()); // line element node numbering
         // rotating this line clockwise to get the normal
         Point<2> nrml( -vec[1] /* -y */, vec[0] /* x */ );
         // checking whether a vector from the faces barycentre to the parent element center yields a negative or positive dot product
         Point<2> bvec = (*shared_parents.rbegin())->BaryCenter() - (*shared_parents.begin())->BaryCenter();
         // using the dot-product to find inner element, if the normal is pointing toward first barycentre, initial order needs to be reversed
         if ( dotProduct( nrml, bvec ) < 0. ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }
 
     // in a 1D model faces coincide with nodes and have just a single node
     if constexpr ( dim == 1 ) {
         // the inside element is that for which the node is node 2
         const size_t parent_element{0};
         if ( (*first)->ParentNodeNumber( parent_element ) == 0 ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }

    return result;
    
 } // end parentElementsSharedByFace

template pair<Element<3>*,Element<3>*>  parentElementsSharedByFace( typename vector<Node<3>*>::const_iterator, typename std::vector<Node<3>*>::const_iterator );
template pair<Element<2>*,Element<2>*>  parentElementsSharedByFace( typename vector<Node<2>*>::const_iterator, typename std::vector<Node<2>*>::const_iterator );
template pair<Element<1>*,Element<1>*>  parentElementsSharedByFace( typename vector<Node<1>*>::const_iterator, typename std::vector<Node<1>*>::const_iterator );




/**
       Returns pointer to element with the supplied face nodes.
       For use in the case where only one parent is expected, for instance, when the Face is at a model boundary.
       
       @return pair of Element and the face or segment of the element that has the same nodes.
*/
template<size_t dim>
pair<Element<dim>*,size_t>  parentElement( typename vector<Node<dim>*>::const_iterator first,
                                           typename vector<Node<dim>*>::const_iterator last )
 {
    // ascertain that there are multiple nodes
    assert( first != last );
    
    // 1. create set of parent elements that share all face nodes
    // ----------------------------------------------------------
    const typename vector<Node<dim>*>::const_iterator nodesEnd{last};
    typename vector<Node<dim>*>::const_iterator nit{first};
    nit++;

    // creating set of parent elements shared by first and second node
    assert( nit != nodesEnd );
    assert( (*nit)->Parents() > 0 );
    const size_t n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;
    for ( size_t i{0}; i<n_parents; ++i ) {
         assert( (*nit)->Parent(i) != nullptr );
         if constexpr ( dim == 3 ) if ( !(*nit)->Parent(i)->IsVolumeElement() ) continue;
         if constexpr ( dim == 2 ) if ( !(*nit)->Parent(i)->IsSurfaceElement() ) continue;
         if ( (*first)->IsParent( (*nit)->Parent(i) ) )
           shared_parents.insert( (*nit)->Parent(i) );
      }
      
    // advancing the node pointer
    nit++;

    // deleting elements that aren't parents of the remaining nodes
    while ( nit != nodesEnd ) {
         for ( auto it=shared_parents.begin(); it!=shared_parents.end(); ) {
              if ( !(*nit)->IsParent(*it) )
                it = shared_parents.erase( it );
              else it++;
           }
         nit++;
      }
      
    // verifying that the results are as expected
    assert( !shared_parents.empty() );
    if (  shared_parents.size() > 1 ) {
         ErrorHandler::Instance().notice( ERROR, "parentElement", "more than one element was found",
                                         "this may be the case for a lower-dimensional element inside the model; use other function");

         return make_pair( (*shared_parents.begin()), UINT_MAX );
      }
    
    
    // 2. find the element's face that matches the nodes
    // -------------------------------------------------
    Element<dim>*  eptr{ (*shared_parents.begin()) };
    vector <Node<dim>*> face_nodes( first, last );
    sort( face_nodes.begin(), face_nodes.end() );
    
    // are the corner nodes of the faces contained in the input node pointer range?
    const size_t n_faces{ eptr->Faces() };
    for ( size_t face{0}; face < n_faces; ++face ) {
         vector<size_t> fnids = eptr->FE()->CornerNodesOfFace(face);
         const size_t n_cnr_nodes{ fnids.size() };
         bool all_nodes_are_contained{true};
         for ( size_t j{0}; j<n_cnr_nodes; ++j )
           if ( !binary_search( face_nodes.begin(), face_nodes.end(), eptr->N( fnids[j] ) ) ) {
                all_nodes_are_contained = false;
                break;
             }
         // when the face has been found the result is returned
         if ( all_nodes_are_contained )
           return make_pair( eptr, face );
             
      }
    
     // 3. if none of the faces contains all nodes, perhaps a segment will
    // -------------------------------------------------------------------
    const size_t n_segments{ eptr->Segments() };
    for ( size_t segm{0}; segm < n_segments; ++segm ) {
         vector<size_t> snids;
         eptr->FE()->NodesOfSegment( segm, snids );
         const size_t n_segm_nodes{ snids.size() };
         bool all_nodes_are_contained{true};
         for ( size_t j{0}; j<n_segm_nodes; ++j )
           if ( !binary_search( face_nodes.begin(), face_nodes.end(), eptr->N( snids[j] ) ) ) {
                all_nodes_are_contained = false;
                break;
             }
         // when the face has been found the result is returned
         if ( all_nodes_are_contained )
           return make_pair( eptr, segm );
      }
   
    return make_pair( nullptr, UINT_MAX );
    
 } // end parentElement

template pair<Element<3>*,size_t> parentElement( typename vector<Node<3>*>::const_iterator,
                                                 typename vector<Node<3>*>::const_iterator );
template pair<Element<2>*,size_t> parentElement( typename vector<Node<2>*>::const_iterator,
                                                 typename vector<Node<2>*>::const_iterator );
template pair<Element<1>*,size_t> parentElement( typename vector<Node<1>*>::const_iterator,
                                                 typename vector<Node<1>*>::const_iterator );



template<size_t dim>
void printParents( const Node<dim>* const nptr )
 {
    assert( nptr != nullptr );
    assert( nptr->Parents() > 0 );
    const size_t n_parents{nptr->Parents()};
    set<Element<dim>*> parents;
    
    cout <<"\nNode "<< nptr->Idx();
    for ( size_t i{0}; i<n_parents; ++i ) {
         if ( nptr->Parent(i) == nullptr ) cout <<" null";
         else {
              cout <<" "<< parseAbbreviated_FE_Type( nptr->Parent(i)->FE_Type() );
              cout <<":"<< nptr->Parent(i)->Idx() <<"(n"<< nptr->ParentNodeNumber(i) <<")";
           }
         parents.insert( nptr->Parent(i) );
      }
      
    const size_t n_duplicates = n_parents - parents.size();
    if ( n_duplicates > 0 )
      cout <<" parent vector contains "<< n_duplicates << " duplicates.";
      
    cout << endl;
      
 } // end printParents

template void printParents( const Node<3>* const );
template void printParents( const Node<2>* const );
template void printParents( const Node<1>* const );

} // end namespace csmp

