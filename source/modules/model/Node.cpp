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
      at_boundary_(NOT),
      BREP_entity_(MESH_VERTEX)
  {
  }



/**
    Initialises everything except for parent element related vectors.
*/
template<uint32_t dim>
Node<dim>::Node( size_t idx,
                 const Point<dim>& pt,
                 const LocalVariables& lvs,
                 BOX_BOUNDARY boundary_flag, TOPOTYPE topotype )
 : LocalVariableStorage<dim,Node>(lvs),
   xyz_(pt),
   idx_(idx),
   manifold_(nullptr),
   at_boundary_(boundary_flag),
   BREP_entity_(topotype)
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
    at_boundary_(nd.at_boundary_),
    BREP_entity_(nd.BREP_entity_)
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
  : xyz_{ std::move(nd.xyz_) },
    idx_{ std::move(nd.idx_) },
    parent_element_pointers_{ std::move(nd.parent_element_pointers_) },
    neighbor_node_pointers_{ std::move(nd.neighbor_node_pointers_) },
    manifold_{ std::move(nd.manifold_) },
    parent_node_indexes_{ std::move(nd.parent_node_indexes_) },
    at_boundary_{ std::move(nd.at_boundary_) },
    BREP_entity_{ std::move(nd.BREP_entity_)}
  {
    this->LVS( std::move(nd.LVS()) );
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
         BREP_entity_             = nd.BREP_entity_;
         parent_element_pointers_ = nd.parent_element_pointers_;
         neighbor_node_pointers_  = nd.neighbor_node_pointers_;
         manifold_                = nd.manifold_;
         parent_node_indexes_     = nd.parent_node_indexes_;
         this->LVS( std::move( nd.LVS() ) );
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
    BREP_entity_             = std::move( nd.BREP_entity_ );
    parent_element_pointers_ = std::move( nd.parent_element_pointers_ );
    neighbor_node_pointers_  = std::move( nd.neighbor_node_pointers_ );
    manifold_                = std::move( nd.manifold_ );
    parent_node_indexes_     = std::move( nd.parent_node_indexes_ );
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
    cerr <<"Node<"<< dim <<">: called operator==(Node) comparitor, extremely costly and of questionable value\n";
    if ( &nd != this )
      {
          if ( BREP_entity_ != nd.BREP_entity_ ) return false;
          if ( Coordinate() == nd.Coordinate() )
            {
                if ( parent_element_pointers_.size() == nd.parent_element_pointers_.size() )
                {
                    if ( parent_element_pointers_.empty() )
                      {
                          if( idx_ == nd.idx_ )
                              return true;
                          return false;
                      }

                    set<Element<dim>*> elmts1(parent_element_pointers_.begin(), parent_element_pointers_.end());
                    set<Element<dim>*> elmts2(nd.parent_element_pointers_.begin(), nd.parent_element_pointers_.end());
                    vector<Element<dim>*> elmts_intersect;
                    set_intersection( elmts1.begin(), elmts1.end(),
                                      elmts2.begin(), elmts2.end(),
                                      back_inserter(elmts_intersect) );
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
void  Node<dim>::Assign( set<Node<dim>*>& neighbor_nodes )
 {
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
 }
 
 
 
template<uint32_t dim>
void  Node<dim>::Assign( vector<Node<dim>*>& neighbor_nodes, bool sort_neighbors )
 {
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
    if ( sort_neighbors )
      sort( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end() );
 }



    /// copies property values from the argument node to the current node
template<uint32_t dim>
void  Node<dim>::CopyPropertyValuesFrom( Node<dim>& nd )
 {
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
    
    for ( uint32_t i{0U}; i < n_parents; ++i )
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


/**
    Moves  unwanted element to the end of the vector before erasing it.
    @attention Use UpdateNeighbors to shrink vector to new size.
*/
template<uint32_t dim>
void Node<dim>::RemoveNeighbor( const Node<dim>* const neighbor_node )
 {
    auto it = remove( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end(), neighbor_node );
//    neighbor_node_pointers_.erase( remove( neighbor_node_pointers_.begin(),
//                                           neighbor_node_pointers_.end(), neighbor_node ),
//                                   neighbor_node_pointers_.end() );
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


/**
   If node is part of line elements (2D) or surface elements (3D), method returns a unit normal that represents the average of the normals of the connected elements.
   
   @return false if 1) the node is not located on the inside of a patch of lower dimensional elements, 2) not properly initialised or 3) adjacent element normals are pointing in opposite directions
*/
template<uint32_t dim>
bool Node<dim>::UnitNormal( Point<dim>& avg_nrml ) const
 {
    // 1. verifying that the node indeed lies on an internal surface
    // 1.1 finding the surface elements connected to the node and their normals
    uint32_t    dim_1_elmt_count{0U};
    Point<dim>  surf_nrml; // initialised to zero
    avg_nrml = 0.;

    if constexpr ( dim == 3U ) {
        for ( auto i{0U}; i<Parents(); i++ )
          if ( Parent(i) && Parent(i)->IsSurface() ) {
              if ( dim_1_elmt_count >= 1 && dotProduct( surf_nrml, Parent(i)->UnitNormal() ) < 0. ) {
                   Out();
                   ErrorHandler::Instance().Note( ERROR, "Node<3U>::UnitNormal",
                                                         "surface element normals point into opposite directions");
                   return false;
                }
              surf_nrml = Parent(i)->UnitNormal();
              avg_nrml += surf_nrml;
              dim_1_elmt_count++;
           }
       }

    if constexpr ( dim == 2U ) {
        for ( auto i{0U}; i<Parents(); i++ )
          if ( Parent(i) && Parent(i)->IsLine() ) {
              if ( dim_1_elmt_count >= 1  && dotProduct( surf_nrml, Parent(i)->UnitNormal() ) < 0. ) {
                   Out();
                   ErrorHandler::Instance().Note( ERROR, "Node<2U>::UnitNormal",
                                                         "line element normals point into opposite directions");
                   return false;
                }
              surf_nrml = Parent(i)->UnitNormal();
              avg_nrml += surf_nrml;
              dim_1_elmt_count++;
           }
       }

    if ( dim_1_elmt_count <= 2 ) {
         ErrorHandler::Instance().Note( ERROR, "Node<dim>::UnitNormal",
                                               "supplied node does not lie in the interior of a surface");
         return false;
      }
    // 1.2 obtaining average normal orientation by normalisation
    avg_nrml /= static_cast<double>(dim_1_elmt_count);
    
    return true;

 } // end UnitNormal



// OUTPUT

/**
 
Outputs internal data of the node. 
*/
template<uint32_t dim>
void Node<dim>::Out() const
 {
    cout <<"\n\nNode<"<< dim <<">: "<< idx_;
    if ( BREP_entity_ != MESH_VERTEX ) {
         string str(parseTopology(BREP_entity_));
         cout <<", topologic role: "<< str;
      }
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
   Returns  the element(s) that  share the face nodes which are supplied to the function via node iterators.
   The inside element is reported first; the outer one next. A nullptr is returned second if only the inside element is found.
   Since the function uses node numbers, it does not depend on a valid element neighbor connectivity.
   
   Application: give nodes of lower-dimensional face to find element on either side
   
   @todo method probably sill contains a large number of redundant operations.
   
   @attention in the case of quadratic FEs, this function must be given only the corner nodes of the element face
*/
template<uint32_t dim>
pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( typename vector<Node<dim>*>::const_iterator first,
                                                               typename vector<Node<dim>*>::const_iterator last )
 {
    assert( first != last );
    // 1. creating sets of the parent elements of the face nodes, testing which ones are shared
    // ----------------------------------------------------------------------------------------
    //const auto nodes_of_face = distance(first,last); // 1 for line, 2-for surface, 3 or 4 for volume
    const auto nodesEnd{last};
    auto       nit{first};

    // creating a set of the parent elements of the first node
    assert( (*nit)->Parents() > 0 );
    const auto n_parents{(*nit)->Parents()};
    set<Element<dim>*> shared_parents;
    // first node
    for ( uint32_t i{0U}; i<n_parents; ++i ) {
         if ( (*nit)->Parent(i) == nullptr ) continue;
         // only considering equi-dimensional elements
         if constexpr ( dim == 3U ) if ( !(*nit)->Parent(i)->IsVolume() ) continue;
         if constexpr ( dim == 2U ) if ( !(*nit)->Parent(i)->IsSurface() ) continue;
         shared_parents.insert( (*nit)->Parent(i) );
      }
      
    // advancing the node iterator
    nit++;

    // searching for shared parent elements in subsequent nodes
    while ( nit != nodesEnd ) {
         set<Element<dim>*> temp;
         const auto parents{(*nit)->Parents()};
         for ( uint32_t i{0U}; i<parents; ++i ) {
              if ( (*nit)->Parent(i) == nullptr ) continue;
              if constexpr ( dim == 3U ) if ( !(*nit)->Parent(i)->IsVolume() ) continue;
              if constexpr ( dim == 2U ) if ( !(*nit)->Parent(i)->IsSurface() ) continue;
              if ( shared_parents.find( (*nit)->Parent(i) ) != shared_parents.end() )
                temp.insert( (*nit)->Parent(i) );
           }
         shared_parents = temp;
         if ( shared_parents.size() == 2U ) break;
         nit++;
      }
      
    // 2. analysing the results and determining the inside element if there are two
    // ----------------------------------------------------------------------------
    if ( shared_parents.empty() ) {
         cerr << endl << endl <<"input nodes: ";
         while ( first != last ) { cerr << (*first)->Idx() <<" "; first++; }
         throw csmp::Exception( ERROR, "parentElementsSharedByFace", "no shared parent elements found" );
      }
       
    // if there is only one element that shares the nodes it must be located on the inside of a boundary
    if ( shared_parents.size() == 1U ) return make_pair( (*shared_parents.begin()), nullptr );
    
    // once two parent elements were found, the one on the inside needs to be determined
    // initial guess
    pair<Element<dim>*,Element<dim>*> result( (*shared_parents.begin()), (*shared_parents.rbegin()) );
    
    // 3D case where face is either a triangle or a quadrilateral
    if constexpr ( dim == 3U ) {
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
    if constexpr ( dim == 2U ) {
         vector<Node<2>*> face_nodes( first, last );
         assert( face_nodes.size() == 2U );
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
     if constexpr ( dim == 1U ) {
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
          ErrorHandler::Instance().Note( ERROR, "parentElement", "no suitable parent element was found" );
          return make_pair( nullptr, numeric_limits<size_t>::max() );
       }
    if (  shared_parents.size() > 1 ) {
    
// DEBUGGING - visualising the discovered higher dimensional elements
#ifdef NODE_DEBUG
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
  ErrorHandler::Instance().Note( WARNING, "parentElement", "more than one element was found",
                                         "and they are interpenetrating (=partially or fully overlapping)");
#endif
         ErrorHandler::Instance().Note( WARNING, "parentElement", "more than one element was found",
                                         "this may be the case for a lower-dimensional element inside the model; use other function");
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



/**
    Determines role of Node (simple mesh node vs. geometric constraint), using BOX_BOUNDARY flagging and parent element connectivity.
    Returns the inferred topo type classifier.
    
    @attention full diagnostics are possible only if the mesh contains line elements for each (edge) curve and triangles or quadrilaterals for each surface of the original BREP.
    
    @return geometric classifier of the node point.
*/
template<uint32_t dim>
TOPOTYPE checkModelPartThatNodeBelongsTo( const Node<dim>* const node )
 {
   // establishing the geometric attributes of the nodes
   const bool at_external_boundary = ( node->AtBoundary() != NOT && node->AtBoundary() != INTERNAL ) ? true : false;
   const bool BREP_defining_node   = ( at_external_boundary == true || node->AtBoundary() == INTERNAL ) ? true : false;
   
   // using parent elements to determine whether the node is situated on a line or surface
   int part_of_line{0}, part_of_surface{0};
   const auto n_parent_elmts{ node->Parents() };
   if ( BREP_defining_node ) {
       // ! works only if the topology defining cells are elements because node has no connection to Face or InterFace
       if ( n_parent_elmts > 0U ) {
           for ( uint32_t i{0U}; i<n_parent_elmts; i++ )
             if ( node->Parent(i) ) {
                if ( node->Parent(i)->IsLine() )
                  part_of_line++;
                else if ( node->Parent(i)->IsSurface() )
                  part_of_surface++;
             }
         }
     }
 
   switch( node->Attribute() ) {
        // if the point serves the sole purpose of discretisation, away from any BREP entities
        case MESH_VERTEX: if ( !BREP_defining_node ) return MESH_VERTEX;
        // essential point of the input geometry such as a line intersection, two surfaces touching etc.
        case INTERSECTION_POINT: if ( BREP_defining_node && part_of_line > 2U ) return INTERSECTION_POINT;
        // includes end points of lines
        case PERIMETER_POINT:
               if constexpr ( dim == 2U ) if ( BREP_defining_node && part_of_line == 1U )
                 return PERIMETER_POINT;
               if constexpr ( dim == 3U ) if ( BREP_defining_node && part_of_line == 2U )
                 return PERIMETER_POINT;
        // point where a line touches the outside boundary of a model
        case EXTERIOR_POINT: if ( at_external_boundary && part_of_line == 0U ) return EXTERIOR_POINT;
        case INTERIOR_LINE:
               if constexpr ( dim == 2U ) if ( BREP_defining_node && !at_external_boundary && part_of_line >= 1U )
                 return INTERIOR_LINE;
               if constexpr ( dim == 3U ) if ( BREP_defining_node && !at_external_boundary && part_of_line >= 2U )
                 return INTERIOR_LINE;
        case PERIMETER_LINE:
               if constexpr ( dim == 3U ) if ( BREP_defining_node && !at_external_boundary && part_of_line >= 2U )
                 return INTERIOR_LINE;
        case EXTERIOR_LINE:
               if ( BREP_defining_node && at_external_boundary && part_of_line >= 1U ) return EXTERIOR_LINE;
        case INTERSECTION_LINE:
               if constexpr ( dim == 3U ) if ( BREP_defining_node && !at_external_boundary &&
                                               part_of_line >= 1U && part_of_surface >= 2U )
                 return INTERSECTION_LINE;
        case INTERIOR_SURFACE:
               if constexpr ( dim == 3U ) if ( BREP_defining_node && !at_external_boundary &&
                                               part_of_line == 0U && part_of_surface == 1U )
                 return INTERIOR_SURFACE;
        case PERIMETER_SURFACE:
               if constexpr ( dim == 3U ) if ( BREP_defining_node && at_external_boundary &&
                                               part_of_line == 0U && part_of_surface == 1U )
                 return PERIMETER_SURFACE;
        case EXTERIOR_SURFACE:
               if constexpr ( dim == 3U ) if ( BREP_defining_node && at_external_boundary &&
                                               part_of_line == 0U && part_of_surface == 1U )
                 return EXTERIOR_SURFACE;
        default:
          cerr <<"\nconsistencyCheck(NodeManifold): TOPOTYPE of Node could not be resolved."<< endl;
      
      } // end switch
      
    return MESH_VERTEX;
          
 } // end checkModelPartThatNodeBelongsTo

template TOPOTYPE checkModelPartThatNodeBelongsTo( const Node<3U>* const );
template TOPOTYPE checkModelPartThatNodeBelongsTo( const Node<2U>* const );
template TOPOTYPE checkModelPartThatNodeBelongsTo( const Node<1U>* const );







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



/**
    For a node that lies on an internal surface, method finds it volumetric (3D) or surface (2D) parent elements on the inside or outside of this lower dimensional feature.
    functions throws if assumptions are not met, i.e., the node does not lie in the interior of a lower dimensional feature.
    
    @attention the inside outside classification is based on a geometric average of the normals of the surface elements that the node forms part of.
    If any of these normals diverges by more that 90o from the others the result is inconclusive
    
    @param node on the interior of a surface inside of the model
    
    @return vectors of pointers to the equidimensional elements on the inside and the outside of the surface that the node lies in the interior of.
    
*/
template<uint32_t dim>
pair<vector<Element<dim>*>,vector<Element<dim>*>>  parentElementsAdjacentTo( const Node<dim>* const node )
 {
    // 0. checking prerequities
    assert( node != nullptr );
    assert( node->AtBoundary() == INTERNAL || node->AtBoundary() == NOT );
    assert( node->Parents() >= 2U ); // must be initialised
 
    // 1. verifying that the node lies on an internal surface and determining unit normal to it by averaging
    Point<dim>  avg_nrml;
    if ( !node->UnitNormal( avg_nrml ) )
      throw csmp::Exception( ERROR, "parentElementsAdjacentTo",
                            "supplied node does not lie in the interior of a surface");
    
    vector<Element<dim>*> inside_elmt_ptrs, outside_elmt_ptrs;
    uint32_t              e_count{ 0U };
    
    // 2. for the volumetric parent elements of the node, determine which side of the surface they lie on
    //   (assumption: if dot product between node and barycentre of these elements is negative they lie on the inside)
    for ( auto i{0U}; i<node->Parents(); i++ )
      if ( node->Parent(i) && node->Parent(i)->IsEquidimensional() ) {
           // find distance between node and barycentre of volumetric/surface element
           Point<dim> bctr_vec = node->Parent(i)->BaryCenter() - node->Coordinate();
           // inside elements are found
           if ( dotProduct( avg_nrml, bctr_vec ) < 0. )
             inside_elmt_ptrs.push_back( node->Parent(i) );
           else
             outside_elmt_ptrs.push_back( node->Parent(i) );
           e_count++;
       }
       
    if ( e_count >= 1 )
    return make_pair( inside_elmt_ptrs, outside_elmt_ptrs );
    
    return pair<vector<Element<dim>*>,std::vector<Element<dim>*>>{};
     
 } // end parentElementsInsideAndOutsideOfSurface

template pair<vector<Element<3>*>,vector<Element<3>*>>  parentElementsAdjacentTo( const Node<3>* const );
template pair<vector<Element<2>*>,vector<Element<2>*>>  parentElementsAdjacentTo( const Node<2>* const );




/**
    For the rectangular of cube-shaped bounding box defined by its min/max corner coordinates, determine whether the Node is contained
    @code
      if ( !(p.x < box.left || p.x > box.right || p.y > box.bottom || p.y < box.top || ... )
    @endcode
*/
template<>
bool isWithinBoundingBox( const Point<3U>& pmin, const Point<3U>& pmax, const Node<3U>* const nptr )
 {
    assert( nptr != nullptr );
    Point<3U> pt = nptr->Coordinate();
    if ( !(pt[0] < pmin[0] || pt[0] > pmax[0] || pt[1] < pmin[1] || pt[1] > pmax[1] || pt[2] < pmin[2] || pt[2] > pmax[2]) ) return true;
    return false;
 }
template<>
bool isWithinBoundingBox( const Point<2U>& pmin, const Point<2U>& pmax, const Node<2U>* const nptr )
 {
    assert( nptr != nullptr );
    Point<2U> pt = nptr->Coordinate();
    if ( !( pt[0] < pmin[0] || pt[0] > pmax[0] || pt[1] < pmin[1] || pt[1] > pmax[1] ) ) return true;
    return false;
 }
template<>
bool isWithinBoundingBox( const Point<1U>& pmin, const Point<1U>& pmax, const Node<1U>* const nptr )
 {
    assert( nptr != nullptr );
    Point<1U> pt = nptr->Coordinate();
    if ( !( pt[0] < pmin[0] || pt[0] > pmax[0] ) ) return true;
    return false;
 }



} // end namespace csmp

