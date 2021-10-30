#include "Node.h"
#include "Element.h"
#include "Visitor.h"
#include "NodeManifold.h"

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
   idx_(idx),
   xyz_(pt),
   at_boundary_(boundary_flag)
 {
 }



/**
    Copy constructor also copies the pointer assignments (!).
*/
template<size_t dim>
Node<dim>::Node( const Node<dim>& nd )
  : idx_(nd.idx_),
    at_boundary_(nd.at_boundary_),
    xyz_(nd.xyz_),
    parent_node_indexes_(nd.parent_node_indexes_),
    parent_element_pointers_(nd.parent_element_pointers_)
  {
    this->LVS( nd.LVS() );
  }



/**
    Move constructor also copies the pointer assignments (!).
*/
template<size_t dim>
Node<dim>::Node( Node<dim>&& nd )
  : idx_{nd.idx_},
    at_boundary_{nd.at_boundary_},
    xyz_{nd.xyz_},
    parent_node_indexes_{nd.parent_node_indexes_},
    parent_element_pointers_{nd.parent_element_pointers_}
  {
    this->LVS( move(nd.LVS()) );
  }




template<size_t dim>
Node<dim>::~Node()
 {
    if ( manifold_ != nullptr )
      manifold_->Remove( this );
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
         this->LVS( nd.LVS() );
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
    idx_                     = nd.idx_;
    at_boundary_             = nd.at_boundary_;
    parent_node_indexes_     = nd.parent_node_indexes_;
    parent_element_pointers_ = nd.parent_element_pointers_;
    xyz_                     = nd.xyz_;
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

@section arguments Input Arguments

The local number of the parent element which is going to be assigned.  
The local node number of this node in the parent element.  
A pointer to the parent Element which shall be added.  
*/
template<size_t dim>
void Node<dim>::Assign( size_t pnode, Element<dim>* element )
 {
    assert( element != nullptr );
    assert( parent_node_indexes_.size() == parent_element_pointers_.size() );
    assert( pnode <= FIFTY );
    
    for ( size_t parent=0U; parent<parent_node_indexes_.size(); parent++ )
      if ( parent_node_indexes_[parent] == NOT_INITIALIZED ) {
           parent_node_indexes_[parent]     = static_cast<ONE_BYTE_NUMBER>(pnode);
           parent_element_pointers_[parent] = element;
           return;
        }

 } // end Assign




/// Removes the provided element as parent and returns true, false if not found
template<size_t dim>
bool Node<dim>::Unassign( Element<dim>* element )
 {
    assert( parent_node_indexes_.size() ==  parent_element_pointers_.size() );
    for( size_t parent(0); parent < Parents(); ++parent )
      if( element == Parent(parent) )
        {
          parent_node_indexes_.erase( parent_node_indexes_.begin()+parent );
          parent_element_pointers_.erase( parent_element_pointers_.begin()+parent );
          // free memory
          parent_node_indexes_.swap( parent_node_indexes_ );
          parent_element_pointers_.swap( parent_element_pointers_ );
          return true;
        }
    return false;
 } // end Unassign


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
    if ( dim == 3U ) {
         for ( typename vector<Element<dim>*>::const_iterator
               it=parent_element_pointers_.begin(); it!=parent_element_pointers_.end(); it++ )
           if ( (*it)->FE()->IsSurfaceElement() or (*it)->FE()->IsLineElement() ) node_neighbors--; 
         return node_neighbors;
      }
    if ( dim == 2U ) {
         for ( typename vector<Element<dim>*>::const_iterator
               it=parent_element_pointers_.begin(); it!=parent_element_pointers_.end(); it++ )
           if ( (*it)->FE()->IsLineElement() ) node_neighbors--; 
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
   
    size_t next_node(ParentNodeNumber(neighbor_node) + 1U);
    if ( next_node == Parent( neighbor_node )->Nodes() ) next_node = 0U;
    
    return parent_element_pointers_[ neighbor_node ]->N( next_node );
 }






// ACCESSORS


template<size_t dim>
double64  Node<dim>::operator[]( size_t i ) const { return xyz_[i]; }

template<size_t dim>
double64&  Node<dim>::operator[]( size_t i ) { return xyz_[i]; }

template<size_t dim>
double64&  Node<dim>::operator()( size_t i ) { return xyz_[i]; }

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
void            Node<dim>::x( double64 xc )  { xyz_[0u] = xc; }

template<size_t dim>
void            Node<dim>::y( double64 yc )  { xyz_[1u] = yc; }

template<size_t dim>
void            Node<dim>::z( double64 zc )  { xyz_[2u] = zc; }


template<size_t dim>
double64          Node<dim>::x() const { return xyz_[0u]; }

template<size_t dim>
double64          Node<dim>::y() const { return xyz_[1u]; }

template<size_t dim>
double64          Node<dim>::z() const { return xyz_[2u]; }


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
   returns parent elements shared by face, inner side is reported first; outer next else application: give nodes of lower-dimensional face to find element on either side
*/
template<size_t dim>
pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( const vector<Node<dim>*>& face_nodes )
 {
    assert( !face_nodes.empty() );
    // create sets of the parent elements of the face nodes checking which ones are shared
    const typename vector<Node<dim>*>::const_iterator nodesEnd{face_nodes.end()};
    typename vector<Node<dim>*>::const_iterator nit{face_nodes.begin()};

    // creating a set of the parent elements of the first node
    assert( (*nit)->Parents() > 0 );
    const size_t n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;
    for ( size_t i{0}; i<n_parents; ++i ) {
         assert( (*nit)->Parent(i) != nullptr );
         if ( dim == parseFiniteElementDimension( (*nit)->Parent(i)->FE()->ElementType() ) )
           shared_parents.insert( (*nit)->Parent(i) );
      }
    nit++;

    // searching for shared parent elements in the following nodes
    while ( nit != nodesEnd ) {
         set<Element<dim>*> temp;
         const size_t n_parents{(*nit)->Parents()};
         for ( size_t i{0}; i<n_parents; ++i ) {
              assert( (*nit)->Parent(i) != nullptr );
              if ( dim == parseFiniteElementDimension( (*nit)->Parent(i)->FE()->ElementType() ) )
                if ( shared_parents.find( (*nit)->Parent(i) ) != shared_parents.end() )
                  temp.insert( (*nit)->Parent(i) );
           }
         // shared parents is narrowed down further
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
    if constexpr ( dim == 3 ) {
         // getting normal to face from the first 3 node coordinates
         assert( face_nodes.size() >= 3 );
         Point<3> vec1(face_nodes[0]->Coordinate() - face_nodes[1]->Coordinate()); // cw
         Point<3> vec2(face_nodes[2]->Coordinate() - face_nodes[1]->Coordinate()); // ccw
         Point<3> nrml =  crossProduct( vec2, vec1 );
         // getting the barycentre of the face (just considering 3 nodes, assuming it is reasonably planar)
         Point<3> fbarycentre = face_nodes[0]->Coordinate() + face_nodes[1]->Coordinate() + face_nodes[2]->Coordinate();
         fbarycentre /= 3.;
         // checking whether a vector from the faces barycentre to the parent element center yields a negative or positive dot product
         Point<3> ebarycentre =(*shared_parents.begin())->BaryCenter();
         Point<3> vec3(ebarycentre - fbarycentre);
         // using dot-product to find inner element: if normal is pointing toward barycentre of first element, initial order needs to be reversed
         if ( dotProduct( nrml, vec3 ) > 0. ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }
    if constexpr ( dim == 2 ) {
         // 2D case where the face is line and the non-existing normal points out of the plane
         assert( face_nodes.size() == 2 );
         Point<2> vec(face_nodes[1]->Coordinate() - face_nodes[0]->Coordinate()); // line element node numbering
         // rotating this line clockwise to get the normal
         Point<2> nrml( -vec[1] /* -y */, vec[0] /* x */ );
         // getting the barycentre of the face (just considering 3 nodes, assuming it is reasonably planar)
         Point<2> fbarycentre = face_nodes[0]->Coordinate() + face_nodes[1]->Coordinate();
         fbarycentre /= 2.;
         // checking whether a vector from the faces barycentre to the parent element center yields a negative or positive dot product
         Point<2> ebarycentre = (*shared_parents.begin())->BaryCenter();
         Point<2> vec1(ebarycentre - fbarycentre);
         // using the dot-product to find inner element, if the normal is pointing toward barycentre, initial order needs to be reversed
         if ( dotProduct( nrml, vec1 ) > 0. ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }
 
     // in a 1D model faces coincide with nodes and have just a single node
     if constexpr ( dim == 1 ) {
         assert( face_nodes.size() == 1 );
         // the inside element is that for which the node is node 2
         const size_t parent_element{0};
         if ( face_nodes[0]->ParentNodeNumber( parent_element ) == 0 ) {
              auto swap     = result.second;
              result.second = result.first;
              result.first  = swap;
           }
      }

    return result;
    
 } // end parentElementsSharedByFace

template pair<Element<3>*,Element<3>*>  parentElementsSharedByFace( const vector<Node<3>*>& );
template pair<Element<2>*,Element<2>*>  parentElementsSharedByFace( const vector<Node<2>*>& );
template pair<Element<1>*,Element<1>*>  parentElementsSharedByFace( const vector<Node<1>*>& );


} // end namespace csmp

