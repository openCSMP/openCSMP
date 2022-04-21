#include "Node.h"
#include "Element.h"
#include "Visitor.h"
#include "NodeManifold.h"
#include "ErrorHandler.h"
#include "MeshManagementUtilities.h"

using namespace std;

namespace csmp {

/**
    default constructor
    
    Node coordinate is initialised to origin (0,0,0)
*/
template<uint32_t dim>
Node<dim>::Node()
    : idx_(numeric_limits<size_t>::max()),
      manifold_(nullptr),
      at_boundary_(NOT)
  {
  }



/**
    Initialises everything except for parent element related vectors.
*/
template<uint32_t dim>
Node<dim>::Node( size_t idx, const Point<dim>& pt, const LocalVariables& lvs, BOX_BOUNDARY boundary_flag )
 : LocalVariableStorage<dim,Node>(lvs),
   xyz_(pt),
   idx_(idx),
   manifold_(nullptr),
   at_boundary_(boundary_flag)
 {
 }



/**
    Copy constructor also copies the pointer assignments (!).
    
    @attention when copy constructing manifold nodes, make sure to add this new Node to it
*/
template<uint32_t dim>
Node<dim>::Node( const Node<dim>& nd )
  : xyz_(nd.xyz_), idx_(nd.idx_),
    parent_element_pointers_(nd.parent_element_pointers_),
    neighbor_node_pointers_(nd.neighbor_node_pointers_),
    manifold_(nd.manifold_),
    parent_node_indexes_(nd.parent_node_indexes_),
    at_boundary_(nd.at_boundary_)
  {
    this->LVS( nd.LVS() );
    // NB: if this is a manifold, the new Node must be added to it,
    //     but this can only be done once the node has been constructed
  }



/**
    Move constructor also copies the pointer assignments (!).
*/
template<uint32_t dim>
Node<dim>::Node( Node<dim>&& nd )
  : xyz_{ move(nd.xyz_) },
    idx_{ move(nd.idx_) },
    parent_element_pointers_{ move(nd.parent_element_pointers_) },
    neighbor_node_pointers_{ move(nd.neighbor_node_pointers_) },
    manifold_{ move(nd.manifold_) },
    parent_node_indexes_{ move(nd.parent_node_indexes_) },
    at_boundary_{ move(nd.at_boundary_) }
  {
    this->LVS( move(nd.LVS()) );
  }




template<uint32_t dim>
Node<dim>::~Node()
 {
//    if ( manifold_ != nullptr ) manifold_->Remove( this );
//    cerr <<"\nNode "<< Idx() <<": called destructor.";
 }




template<uint32_t dim>
Node<dim>& Node<dim>::operator=( const Node<dim>& nd )
 {
    if ( &nd != this ) 
      {
         xyz_                     = nd.xyz_;
         idx_                     = nd.idx_;
         at_boundary_             = nd.at_boundary_;
         parent_element_pointers_ = nd.parent_element_pointers_;
         neighbor_node_pointers_  = nd.neighbor_node_pointers_;
         manifold_                = nd.manifold_;
         parent_node_indexes_     = nd.parent_node_indexes_;
         this->LVS( move( nd.LVS() ) );
      }
    return *this;
 }




/**
    Move assignment, relying on that similar operators exist for the nodes components.
    
    @note this assumes that the supplied node is a temporary.
*/
template<uint32_t dim>
Node<dim>& Node<dim>::operator=( Node<dim>&& nd )
 {
    assert( this != &nd );
    xyz_                     = nd.xyz_;
    idx_                     = nd.idx_;
    at_boundary_             = nd.at_boundary_;
    parent_element_pointers_ = move( nd.parent_element_pointers_ );
    neighbor_node_pointers_  = move( nd.neighbor_node_pointers_ );
    manifold_                = move( nd.manifold_ );
    parent_node_indexes_     = move( nd.parent_node_indexes_ );
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
template<uint32_t dim>
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

template<uint32_t dim>
void* Node<dim>::operator new( size_t size )
  {
//      std::cout<< "\nNode<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

template<uint32_t dim>
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
template<uint32_t dim>
void Node<dim>::Assign( uint32_t pnode, Element<dim>* element )
 {
    assert( parent_node_indexes_.size() == parent_element_pointers_.size() );
    assert( pnode <= FIFTY );
    
    for ( uint32_t parent{0}; parent<parent_node_indexes_.size(); parent++ )
      if ( parent_node_indexes_[parent] == NOT_INITIALIZED ) {
           parent_node_indexes_[parent]     = static_cast<ONE_BYTE_NUMBER>(pnode);
           parent_element_pointers_[parent] = element;
           return;
        }

 } // end Assign




/**
    Sets the pointer to the argument element to zero, and the corresponding node number to NOT_INITIALIZED.
    
    @note the argument element is not deleted. To do this use EraseNullPointerParents().
*/
template<uint32_t dim>
bool Node<dim>::Unassign( Element<dim>* element )
 {
    assert( parent_node_indexes_.size() ==  parent_element_pointers_.size() );
    for ( uint32_t parent(0); parent < Parents(); ++parent )
      if ( Parent(parent) == element )
        {
          parent_node_indexes_[parent]     = NOT_INITIALIZED;
          parent_element_pointers_[parent] = nullptr;
          return true;
        }
    return false;
    
 } // end Unassign




/**
    Removes parent elements pointers, but only if these were set to nullptr before.
*/
template<uint32_t dim>
void Node<dim>::EraseNullPointerParents()
 {
    assert( parent_node_indexes_.size() ==  parent_element_pointers_.size() );
    const auto  n_parents{Parents()};
    uint32_t    n_new_parents{0};

    for ( auto parent{0}; parent < n_parents; ++parent ) {
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




template<uint32_t dim>
void  Node<dim>::Accept( csmp::Visitor<dim>& v )
{
// TODO: implement dispatching of visitor to neighbor nodes
    v.Visit(this);
}


// node neighbor functionality

/// initialises the corner-node to neighbor corner node pointer vector
template<uint32_t dim>
void  Node<dim>::Assign( std::set<Node<dim>*>& neighbor_nodes )
 {
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
 }
 
 
 
template<uint32_t dim>
void  Node<dim>::Assign( std::vector<Node<dim>*>& neighbor_nodes, bool sort_neighbors )
 {
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
    if ( sort_neighbors )
      sort( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end() );
 }



template<uint32_t dim>
void  Node<dim>::AssignPropertyValuesFrom( const Node<dim>& nd )
  {
     // copies the property values
     this->LVS( nd.LVS() );
  }



/**
    Re-establishes the node neighbors of the node using the parent element connectivity.
    Since lower-dimensional elements share the nodes with the higher dimensional ones, they are ignored.
    Only corner nodes are considered in the first pass.
    Subsequently midside nodes which exist in higher-order elements are detected and their neighbors are assigned.
*/
template<uint32_t dim>
uint32_t  Node<dim>::AssignNodeNeighbors()
 {
    set<Node<dim>*>  current_nbors;
    const auto       n_parents{ Parents() };
    
    for ( auto i{0U}; i < n_parents; ++i )
      // only considering parent elements with the same dimensions as the model
      if ( Parent(i) && Parent(i)->IsEquidimensional() ) {
           for ( auto& j :  Parent(i)->CornerNodesConnectedTo( ParentNodeNumber(i) ) )
             current_nbors.insert( j );
        }
    
    neighbor_node_pointers_.assign( current_nbors.begin(), current_nbors.end() );
    neighbor_node_pointers_.shrink_to_fit();
    
    return static_cast<uint32_t>(neighbor_node_pointers_.size());
 }




/**
   Sorts vector and removes duplicates and nullptrs.
   The vector is trimmed so that size matches capacity.
*/
template<uint32_t dim>
uint32_t  Node<dim>::UpdateNeighbors()
 {
    sort( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end() );
    
    neighbor_node_pointers_.erase( unique( neighbor_node_pointers_.begin(),
                                           neighbor_node_pointers_.end() ),
                                   neighbor_node_pointers_.end() );
                                   
    neighbor_node_pointers_.erase( remove( neighbor_node_pointers_.begin(),
                                           neighbor_node_pointers_.end(), nullptr ),
                                   neighbor_node_pointers_.end() );
                                   
    neighbor_node_pointers_.shrink_to_fit();
    
    return static_cast<uint32_t>(neighbor_node_pointers_.size());
 }


 
template<uint32_t dim>
bool  Node<dim>::IsNeighbor( const Node<dim>* const nptr ) const
 {
    return binary_search( neighbor_node_pointers_.begin(),
                          neighbor_node_pointers_.end(), nptr );
 }


template<uint32_t dim>
void  Node<dim>::AddNeighbor( Node<dim>* neighbor_node )
 {
    // grow the vector in small increments only
    if ( neighbor_node_pointers_.size() == neighbor_node_pointers_.capacity() )
      neighbor_node_pointers_.reserve( neighbor_node_pointers_.size() + 2 );
      
    neighbor_node_pointers_.push_back( neighbor_node );
    sort( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end() );
 }


/// just moves the unwanted element to the end of the vector, use UpdateNeighbors to shrink vector to new size
template<uint32_t dim>
void  Node<dim>::RemoveNeighbor( const Node<dim>* const neighbor_node )
 {
    remove( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end(), neighbor_node );
 }
 
 




/**

   If the parent elements are of the same dimension as the node,
   then there are exactly as many elements as are nodes.
   
   @attention calculation is costly, compute only once per node loop!
   
   @test correct SKM (unit test exists)
*/
template<uint32_t dim>
uint32_t  Node<dim>::Neighbors() const
 {
   return static_cast<uint32_t>(neighbor_node_pointers_.size());
 }


/**
    Implements node connectivity graph, giving access to all the Node objects that the Node is connected to via Element, Face or InterFace edges.
    Node neighbors are enlisted in an order that is determined by sorting the points to them (this makes the neighbor vector searchable).
    
    @param neighbor_node the nth Node in the sortes vector of pointers to nodes stored in the node.
    
    @attention If the node is a manifold (topologically co-located with other nodes that can be accessed looping over the branches of the manifold,
    then only those Nodes are neighbors who are on the same contiguous mesh patch as the Node, i.e., not shared with the other nodes in the manifold.
    
    @author SKM
    @date 8/10/2021
*/
template<uint32_t dim>
Node<dim>*  Node<dim>::Neighbor( uint32_t neighbor_node ) const
 {
    assert( neighbor_node < neighbor_node_pointers_.size() );
    return neighbor_node_pointers_[ neighbor_node ];
 }






// ACCESSORS


template<uint32_t dim>
double  Node<dim>::operator[]( uint32_t i ) const { return xyz_[i]; }

template<uint32_t dim>
double&  Node<dim>::operator[]( uint32_t i ) { return xyz_[i]; }

template<uint32_t dim>
double&  Node<dim>::operator()( uint32_t i ) { return xyz_[i]; }

template<uint32_t dim>
Point<dim>  Node<dim>::Coordinate() const { return xyz_; }

template<uint32_t dim>
 void Node<dim>::Coordinate( const Point<dim>& p ) { xyz_=p; }



/**

The dynamic storage for the parent element data is resized preserving the
existing entries. The values of potential new elements are set to zero.

@param n The desired new size of the storage.
*/
template<uint32_t dim>
void  Node<dim>::ResizeParentStorage( uint32_t n )
  {
     // while this does not release the memory, it is essential to zap previous content of parent storage
     parent_node_indexes_.clear();
     parent_element_pointers_.clear();
     // resize, initialising contained pointers to null and indices to NOT_INITIALIZED
     parent_node_indexes_.resize( n, NOT_INITIALIZED );
     std::vector<ONE_BYTE_NUMBER>( parent_node_indexes_ ).swap( parent_node_indexes_ );
     parent_element_pointers_.resize( n, nullptr );
     std::vector<Element<dim>*>( parent_element_pointers_ ).swap( parent_element_pointers_ );
  }



template<uint32_t dim>
void  Node<dim>::EraseParents()
  {
     parent_element_pointers_.clear();
     parent_node_indexes_.clear();
  }



template<uint32_t dim>
uint32_t   Node<dim>::Parents() const
  { return static_cast<uint32_t>(parent_element_pointers_.size()); }



template<uint32_t dim>
uint32_t  Node<dim>::ParentNodeNumber( uint32_t parent_element_number ) const
  {
     assert( parent_element_number < parent_node_indexes_.size() );
     return parent_node_indexes_[ parent_element_number ];
  }


template<uint32_t dim>
Element<dim>*  Node<dim>::Parent( uint32_t parent_element_number ) const
  {
     assert( parent_element_number < parent_element_pointers_.size() );
     return parent_element_pointers_[ parent_element_number ];
  }


    /// checks whether Element is a parent of the node
template<uint32_t dim>
bool  Node<dim>::IsParent( const Element<dim>* const eptr ) const
  {
     return binary_search( parent_element_pointers_.begin(),
                           parent_element_pointers_.end(), eptr );
  }




/**
     sorts parent vector for searching.
     
     @attention parent vector must not contain any nullptrs.
*/
template<uint32_t dim>
void  Node<dim>::SortParents() {
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
    for ( auto& i : indices ) temp[n++] = parent_node_indexes_[i];
    parent_node_indexes_ = temp;

  } // end UpdateParents





template<uint32_t dim>
void  Node<dim>::Idx( size_t idx_to_assign ) const
 {
    idx_ = idx_to_assign;
 }


template<uint32_t dim>
size_t   Node<dim>::Idx() const
 {
    return idx_;
 }



/**
   Set or return the global boundary flag of the Node.
*/
template<uint32_t dim>
void Node<dim>::AtBoundary( BOX_BOUNDARY b ) { at_boundary_ = b; }


template<uint32_t dim>
BOX_BOUNDARY  Node<dim>::AtBoundary() const { return at_boundary_; }

/// Set or return the coordinates of the current node.
template<uint32_t dim>
void            Node<dim>::x( double xc )  { xyz_[0u] = xc; }

template<uint32_t dim>
void            Node<dim>::y( double yc )  { xyz_[1u] = yc; }

template<uint32_t dim>
void            Node<dim>::z( double zc )  { xyz_[2u] = zc; }


template<uint32_t dim>
double          Node<dim>::x() const { return xyz_[0u]; }

template<uint32_t dim>
double          Node<dim>::y() const { return xyz_[1u]; }

template<uint32_t dim>
double          Node<dim>::z() const { return xyz_[2u]; }


// MANIFOLDS

/// access to manifold if any; returns nullptr if the node is not a manifold
template<uint32_t dim>
bool Node<dim>::IsManifold() const { return (manifold_ != nullptr); }


template<uint32_t dim>
NodeManifold<dim>* const Node<dim>::Manifold() const { return manifold_; }


template<uint32_t dim>
void Node<dim>::Assign( NodeManifold<dim>& nmf )
 {
    manifold_ = &nmf;
 }



// OUTPUT

/**
 
Outputs internal data of the node. 
*/
template<uint32_t dim>
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
         for ( auto i{0U}; i<Parents(); i++ ) {
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
template<uint32_t dim>
pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( typename vector<Node<dim>*>::const_iterator first,
                                                               typename vector<Node<dim>*>::const_iterator last )
 {
    assert( first != last );
    // create sets of the parent elements of the face nodes checking which ones are shared
    const typename vector<Node<dim>*>::const_iterator nodesEnd{last};
    typename vector<Node<dim>*>::const_iterator nit{first};

    // creating a set of the parent elements of the first node
    assert( (*nit)->Parents() > 0 );
    const auto n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;
    for ( auto i{0U}; i<n_parents; ++i ) {
         assert( (*nit)->Parent(i) != nullptr );
         if constexpr ( dim == 3 ) if ( !(*nit)->Parent(i)->IsVolume() ) continue;
         if constexpr ( dim == 2 ) if ( !(*nit)->Parent(i)->IsSurface() ) continue;
         shared_parents.insert( (*nit)->Parent(i) );
      }
      
    // advancing the node iterator
    nit++;

    // searching for shared parent elements in the following nodes
    while ( nit != nodesEnd ) {
         set<Element<dim>*> temp;
         const auto parents{(*nit)->Parents()};
         for ( auto i{0U}; i<parents; ++i ) {
              assert( (*nit)->Parent(i) != nullptr );
              if constexpr ( dim == 3 ) if ( !(*nit)->Parent(i)->IsVolume() ) continue;
              if constexpr ( dim == 2 ) if ( !(*nit)->Parent(i)->IsSurface() ) continue;
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
         const uint32_t parent_element{0};
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
       
       @param first and last are iterators to the nodes of the Face.
       
       @return pair of Element and the face or segment of the element that has the same nodes.
       
       @attention method works only if all the parent element pointers that the node stores are valid.
*/
template<uint32_t dim>
pair<Element<dim>*,size_t>  parentElement( typename vector<Node<dim>*>::const_iterator first,
                                           typename vector<Node<dim>*>::const_iterator last )
 {
    // ascertain that there are multiple nodes
    assert( first != last );
    
    // 1. create set of higher-dimensional parent elements that share all face nodes
    // -----------------------------------------------------------------------------
    const typename vector<Node<dim>*>::const_iterator nodesEnd{last};
    typename vector<Node<dim>*>::const_iterator       nit{first};

    // for all equidimensional parents elements of first node, select the ones that also are parents of the other nodes
    assert( nit != nodesEnd );
    assert( (*nit)->Parents() > 0 );
    const auto         n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;

    // creating set of parent elements shared by first and second node
    for ( auto i{0U}; i<n_parents; ++i )
      if ( (*nit)->Parent(i) != nullptr ) {
           if constexpr ( dim == 3U ) if ( !(*nit)->Parent(i)->IsVolume() ) continue;
           if constexpr ( dim == 2U ) if ( !(*nit)->Parent(i)->IsSurface() ) continue;
           // is this parent also one of the first node
           bool parent_to_all{true};
           for ( auto nit2=first; nit2!=nodesEnd; ++nit2 )
             if ( !(*nit2)->IsParent( (*nit)->Parent(i) ) ) {
                  parent_to_all = false;
                  break;
               }
           if ( parent_to_all )
             shared_parents.insert( (*nit)->Parent(i) );
        }
      
    // verifying that the results are as expected
    if (  shared_parents.empty() ) {
          for ( ; first!=last; ++first )
            printParents( (*first) );
          ErrorHandler::Instance().notice( ERROR, "parentElement", "no suitable parent element was found" );

          return make_pair( (*shared_parents.begin()), numeric_limits<size_t>::max() );
       }
    if (  shared_parents.size() > 1 ) {
    
// DEBUGGING - visualising the discovered higher dimensional elements
Element<dim>* elmt1 = (*shared_parents.begin());
elmt1->Idx( 1 );
elmt1->CoordinateMatrix();
DenseMatrix<DM_MIN>  DATA1( 1, elmt1->Nodes() );
for ( int i{0}; i<elmt1->Nodes(); ++i ) DATA1(0,i) = static_cast<double>(elmt1->N(i)->AtBoundary());
elmt1->FE()->OutputNodeDataToVTK( "parent_elmt", "node_flag", DATA1 );
Element<dim>* elmt2 = (*shared_parents.rbegin());
elmt2->Idx( 2 );
elmt2->CoordinateMatrix();
DenseMatrix<DM_MIN>  DATA2( 1, elmt2->Nodes() );
for ( int i{0}; i<elmt2->Nodes(); ++i ) DATA2(0,i) = static_cast<double>(elmt2->N(i)->AtBoundary());
elmt2->FE()->OutputNodeDataToVTK( "parent_elmt", "node_flag", DATA2 );
// if there are two elements, are they overlapping?
if ( interPenetrating<dim>( elmt1, elmt2 ) )
  ErrorHandler::Instance().notice( ERROR, "parentElement", "more than one element was found",
                                         "and they are interpenetrating (=partially or fully overlapping)");

         ErrorHandler::Instance().notice( ERROR, "parentElement", "more than one element was found",
                                         "this may be the case for a lower-dimensional element inside the model; use other function");

         return make_pair( (*shared_parents.begin()), numeric_limits<size_t>::max() );
      }
    
    
    // 2. find the element's face that matches the nodes
    // -------------------------------------------------
    Element<dim>*  eptr{ (*shared_parents.begin()) };
    vector <Node<dim>*> face_nodes( first, last );
    sort( face_nodes.begin(), face_nodes.end() );
    
    // are the corner nodes of the faces contained in the input node pointer range?
    const auto n_faces{ eptr->Faces() };
    for ( auto face{0}; face < n_faces; ++face ) {
         vector<uint32_t> fnids = eptr->FE()->CornerNodesOfFace(face);
         const auto n_cnr_nodes{ fnids.size() };
         bool all_nodes_are_contained{true};
         for ( uint32_t j{0}; j<n_cnr_nodes; ++j )
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
    const auto n_segments{ eptr->Segments() };
    for ( uint32_t segm{0}; segm < n_segments; ++segm ) {
         vector<uint32_t> snids;
         eptr->FE()->NodesOfSegment( segm, snids );
         const auto n_segm_nodes{ snids.size() };
         bool all_nodes_are_contained{true};
         for ( auto j{0U}; j<n_segm_nodes; ++j )
           if ( !binary_search( face_nodes.begin(), face_nodes.end(), eptr->N( snids[j] ) ) ) {
                all_nodes_are_contained = false;
                break;
             }
         // when the face has been found the result is returned
         if ( all_nodes_are_contained )
           return make_pair( eptr, segm );
      }
   
    return make_pair( nullptr, numeric_limits<uint32_t>::max() );
    
 } // end parentElement

template pair<Element<3>*,size_t> parentElement( typename vector<Node<3>*>::const_iterator,
                                                 typename vector<Node<3>*>::const_iterator );
template pair<Element<2>*,size_t> parentElement( typename vector<Node<2>*>::const_iterator,
                                                 typename vector<Node<2>*>::const_iterator );
template pair<Element<1>*,size_t> parentElement( typename vector<Node<1>*>::const_iterator,
                                                 typename vector<Node<1>*>::const_iterator );



template<uint32_t dim>
void printNeighbors( const Node<dim>* const nptr )
 {
    assert( nptr != nullptr );
    const auto n_nbors{nptr->Neighbors()};
    assert( nptr->Parents() > 0 );
    
    cout <<"\nnode "<< nptr->Idx() <<":";
    for ( auto i{0U}; i<n_nbors; ++i ) {
         const Node<dim>* const nd_nbor = nptr->Neighbor(i);
         if ( nd_nbor == nullptr ) cout <<" NULL";
         else cout <<" "<< nd_nbor->Idx(); // <<":"<< parseBoundary( nd_nbor->AtBoundary() );
      }
      
    cout <<" ("<< parseBoundary( nptr->AtBoundary() ) <<")";
      
 } // end printParents

template void printNeighbors( const Node<3>* const );
template void printNeighbors( const Node<2>* const );
template void printNeighbors( const Node<1>* const );





template<uint32_t dim>
void printParents( const Node<dim>* const nptr )
 {
    assert( nptr != nullptr );
    assert( nptr->Parents() > 0 );
    const auto n_parents{nptr->Parents()};
    set<Element<dim>*> parents;
    
    cout <<"\nNode "<< nptr->Idx();
    for ( auto i{0U}; i<n_parents; ++i ) {
         if ( nptr->Parent(i) == nullptr ) cout <<" NULL";
         else {
              cout <<" "<< parseAbbreviated_FE_Type( nptr->Parent(i)->FE_Type() );
              cout <<":"<< nptr->Parent(i)->Idx() <<"(n"<< nptr->ParentNodeNumber(i) <<")";
           }
         parents.insert( nptr->Parent(i) );
      }
      
    const auto n_duplicates = n_parents - parents.size();
    if ( n_duplicates > 0 )
      cout <<" parent vector contains "<< n_duplicates << " duplicates.";
      
    cout << endl;
      
 } // end printParents
 
 
template void printParents( const Node<3>* const );
template void printParents( const Node<2>* const );
template void printParents( const Node<1>* const );






// calculates the size of the Node excluding the stored variables
template<uint32_t dim>
size_t sizeOf( const Node<dim>* const nptr )
  {
    size_t total_size = sizeof( *nptr ); // padded static store of object
    // dynamic allocation
    total_size += nptr->Parents() * sizeof( Element<dim>* );
    total_size += nptr->Parents() * sizeof( ONE_BYTE_NUMBER );
    // + local variable storage
    
    return total_size;

  } // end sizeOf

template size_t sizeOf( const Node<3>* const );
template size_t sizeOf( const Node<2>* const );
template size_t sizeOf( const Node<1>* const );



} // end namespace csmp

