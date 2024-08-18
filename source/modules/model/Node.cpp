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
    assert( element != nullptr );
    
    // looking for a free slot in existing parent vectors
    const auto parents{ parent_node_indexes_.size() };
    for ( uint32_t parent{0u}; parent < parents; parent++ )
      if ( parent_node_indexes_[parent] == NOT_INITIALIZED ) {
           parent_node_indexes_[parent]     = static_cast<ONE_BYTE_NUMBER>(pnode);
           parent_element_pointers_[parent] = element;
           return;
        }
        
    // if the parent element storage must be extended
    parent_node_indexes_.push_back( static_cast<ONE_BYTE_NUMBER>(pnode) );
    parent_element_pointers_.push_back( element );

 } // end Assign




/**
    Sets the pointer to the argument element to zero, and the corresponding node number to NOT_INITIALIZED.
    
    @note the argument element is not deleted. To do this use EraseNullPointerParents().
*/
template<uint32_t dim>
bool Node<dim>::Unassign( Element<dim>* element )
 {
    if ( element == nullptr ) return false;
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
     parent_node_indexes_.shrink_to_fit();
     parent_element_pointers_.resize( n, nullptr );
     parent_element_pointers_.shrink_to_fit();
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



// TODO: refactor to a more efficient design, perhaps vector<pair<size_t,eptr>
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
      Find all parent elements that contain the supplied range of nodes.
      
  @return subset of parent elements that share all the supplied nodes
*/
template<uint32_t dim>
vector<Element<dim>*> parentElementsContaining( typename vector<Node<dim>*>::const_iterator first,
                                                typename vector<Node<dim>*>::const_iterator last )
 {
    // the supplied range of nodes must contain at least two nodes
    assert( next(first,1) != last );
    
    // sorting parent element vectors for the intersection algorithm
    (*first)->SortParents();
    (*next(first,1))->SortParents();

    vector<Element<dim>*>  elmts_with_all_nodes, isect;
    // store the elements shared between the first and the second node in 'isect'
    set_intersection( (*first)->ParentElementsBegin(), (*first)->ParentElementsEnd(),
                      (*next(first,1))->ParentElementsBegin(), (*next(first,1))->ParentElementsEnd(),
                      back_inserter(elmts_with_all_nodes) );
    first++;
    first++;
    
    // find the shared parent elements for all the supplied nodes
    while ( first != last ) {
         (*first)->SortParents();
         set_intersection( elmts_with_all_nodes.begin(), elmts_with_all_nodes.end(),
                           (*first)->ParentElementsBegin(), (*first)->ParentElementsEnd(),
                           back_inserter(isect) );
         elmts_with_all_nodes = isect;
         isect.clear();
         first++;
      }
      
    return elmts_with_all_nodes;
      
  } // end parentElementsContaining

template vector<Element<3U>*> parentElementsContaining( typename vector<Node<3U>*>::const_iterator,
                                                        typename vector<Node<3U>*>::const_iterator );
template vector<Element<2U>*> parentElementsContaining( typename vector<Node<2U>*>::const_iterator,
                                                        typename vector<Node<2U>*>::const_iterator );
template vector<Element<1U>*> parentElementsContaining( typename vector<Node<1U>*>::const_iterator,
                                                        typename vector<Node<1U>*>::const_iterator );




 
/**
   Returns  the dim-dimensional element(s) that  share the face corner nodes supplied to the function via node iterators.
   The inside element is reported first; the outer one next. A nullptr is returned second if only the inside element is found.
   Since the function uses node numbers, it does not depend on a valid element neighbor connectivity.
   
   Application: give nodes of lower-dimensional face to find element on either side
   
   @attention THE SUPPLIED NODES MUST BE IN CORRECT ORDER because they are used to identify the inside element
   
   @attention in the case of quadratic FEs, this function must be given only the corner nodes of the element face
*/
template<uint32_t dim>
pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t>>  parentElements( typename vector<Node<dim>*>::const_iterator first,
                                                                                 typename vector<Node<dim>*>::const_iterator last )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    assert( first != last );
    
    // 1. Find all dim-dimensional parent elements that contain the supplied range of nodes
    // ------------------------------------------------------------------------------------
    vector<Element<dim>*>  shared_parents = parentElementsContaining<dim>( first, last );
    const auto remaining_cells = eraseLowerDimensionalOrInvalidCells( shared_parents );
      
    // 2. analysing the results
    // ----------------------------------------------------------------------------
    if ( remaining_cells == 0U ) {
         cerr << endl << endl <<"input nodes: ";
         auto first1{ first };
         while ( first1 != last ) { cerr << (*first1)->Idx() <<" "; first1++; }
         csmp_error.Note( ERROR, "parentElements", "no shared parent elements found; returning null pointers" );
         return make_pair( make_pair( nullptr, numeric_limits<uint32_t>::max()), make_pair( nullptr, numeric_limits<uint32_t>::max() ) );
      }
    // if there is only one element that shares the nodes it must be located on the inside of a boundary
    else if ( remaining_cells == 1U ) {
         csmp_error.Note( WARNING, "parentElements", "only one shared parent element found which is returned element" );
         auto face_id = faceWithCornerNodes( (*shared_parents.begin()), first, last );
         return make_pair( make_pair( (*shared_parents.begin()), face_id ), make_pair( nullptr, numeric_limits<uint32_t>::max() ) );
      }
    else if ( remaining_cells > 2U ) {
         csmp_error.Note( ERROR, "parentElements", "found more than 2 parent elements ?!, returning first two" );
         return make_pair( make_pair( (*shared_parents.begin()), numeric_limits<uint32_t>::max() ),
                           make_pair( (*next(shared_parents.begin(),1)), numeric_limits<uint32_t>::max() ) );
      }
    
    // 3. finding the shared faces of the elements
    // -------------------------------------------------------------------------------
    // (assuming that the neighbor connectivity between the elements works)
    pair<uint32_t,uint32_t> shared_faces = findAdjacentFacesFromNeighbors( shared_parents[0], shared_parents[1] );
    if ( shared_faces.first == numeric_limits<uint32_t>::max() ||
         shared_faces.second == numeric_limits<uint32_t>::max() ) {
         csmp_error.Note( ERROR, "parentElements", "could not determine face indices of share faces" );
         shared_faces = findAdjacentFacesFromNodes( shared_parents[0], shared_parents[1] );
         assert( shared_faces.first != numeric_limits<uint32_t>::max() );
         assert( shared_faces.second != numeric_limits<uint32_t>::max() );
      }

    // 4. determining the inside element which must have the same node ordering as the one supplied
    // --------------------------------------------------------------------------------------------
    vector<Node<dim>*> cnr_node_vec( first, last ), eface_node_vec1;
    for ( const auto& nd : shared_parents[0]->FE()->CornerNodesOfFace( shared_faces.first ) )
      eface_node_vec1.push_back( shared_parents[0]->N(nd) );
    
    // if there is a match the elements are returned in their current order
    if ( cnr_node_vec == eface_node_vec1 )
      return make_pair( make_pair( shared_parents[0], shared_faces.first ), make_pair( shared_parents[1], shared_faces.second ) );
      
    // else the second lot of face nodes are compared
    vector<Node<dim>*> eface_node_vec2;
    for ( const auto& nd : shared_parents[1]->FE()->CornerNodesOfFace( shared_faces.second ) )
      eface_node_vec2.push_back( shared_parents[1]->N(nd) );
    if ( cnr_node_vec == eface_node_vec2 )
      return make_pair( make_pair( shared_parents[1], shared_faces.second ), make_pair( shared_parents[0], shared_faces.first ) );
      
    // 5. hopefully we never get here: diagnostics
    // -------------------------------------------
#ifndef NDEBUG
    cout <<"\n"<<"parentElements: failed to match nodes:";
    for ( const auto& nit : cnr_node_vec ) cout <<" "<< nit->Idx();
    cout <<", with face nodes of shared elements:";
    cout <<"\n\t"<< shared_parents[0]->Idx() <<": " << parseFiniteElementType( shared_parents[0]->FE_Type() );
    cout <<": nodes of shared face: "<< shared_faces.first <<": ";
    for ( const auto& nit : eface_node_vec1 ) cout << nit->Idx() <<" ";
    cout <<"\n\t"<< shared_parents[1]->Idx() <<": " << parseFiniteElementType( shared_parents[1]->FE_Type() );
    cout <<": nodes of shared face: "<< shared_faces.second <<": ";
    for ( const auto& nit : eface_node_vec2 ) cout << nit->Idx() <<" ";
    cout <<"\n"<<"trying node rotation now..."<< endl;
#endif

    // 6. attempting the same after rotating the corner nodes
    // ------------------------------------------------------
    // first element first
    for ( uint32_t i{0u}; i<cnr_node_vec.size(); ++i ) {
         rotate( eface_node_vec1.begin(), eface_node_vec1.begin()+1, eface_node_vec1.end() );
         if ( cnr_node_vec == eface_node_vec1 )
           return make_pair( make_pair( shared_parents[0], shared_faces.first ), make_pair( shared_parents[1], shared_faces.second ) );
      }
    // second element
    for ( uint32_t i{0u}; i<cnr_node_vec.size(); ++i ) {
         rotate( eface_node_vec2.begin(), eface_node_vec2.begin()+1, eface_node_vec2.end() );
         if ( cnr_node_vec == eface_node_vec2 )
           return make_pair( make_pair( shared_parents[1], shared_faces.second ), make_pair( shared_parents[0], shared_faces.first ) );
      }

#ifndef NDEBUG
    cout <<"\n"<<"parentElements: failed to match nodes even after rotation:";
    for ( const auto& nit : cnr_node_vec ) cout <<"\n\t\t"<< nit->Idx() <<": "<< nit->Coordinate();
    cout <<"\n\t"<<"with face nodes of shared elements:";
    cout <<"\n\t\t"<< shared_parents[0]->Idx() <<": " << parseFiniteElementType( shared_parents[0]->FE_Type() );
    cout <<"\n\t\t"<<"nodes of shared face: "<< shared_faces.first <<": ";
    for ( const auto& nit : eface_node_vec1 ) cout <<"\n\t\t\t"<< nit->Idx() <<": "<< nit->Coordinate();
    cout <<"\n\t\t"<< shared_parents[1]->Idx() <<": " << parseFiniteElementType( shared_parents[1]->FE_Type() );
    cout <<"\n\t\t"<<"nodes of shared face: "<< shared_faces.second <<": ";
    for ( const auto& nit : eface_node_vec2 ) cout <<"\n\t\t\t"<< nit->Idx() <<": "<< nit->Coordinate();
    cout << endl;
    csmp_error.Note( WARNING, "parentElements", "could not match face nodes on inside and outside" );
    cout <<"\n"<<"trying matching of node sets now..."<< endl;
#endif

    // 7. returning those faces which have the same collection of nodes
    // ----------------------------------------------------------------
    set<Node<dim>*> key_cnr_nodes( cnr_node_vec.begin(), cnr_node_vec.end() ),
                    key_e1_nodes( eface_node_vec1.begin(), eface_node_vec1.end() ),
                    key_e2_nodes( eface_node_vec2.begin(), eface_node_vec2.end() );
                    
    if ( key_cnr_nodes == key_e1_nodes && key_cnr_nodes == key_e2_nodes )
      return make_pair( make_pair( shared_parents[0], shared_faces.first ), make_pair( shared_parents[1], shared_faces.second ) );

    // 8. Oh dear! - was the wrong face identified?
    // ----------------------------------------------------------------
#ifndef NDEBUG
    const auto e_matched      = ( key_cnr_nodes != key_e1_nodes ) ? 0 : 1;
    const auto e_not_found    = ( key_cnr_nodes != key_e1_nodes ) ? 1 : 0;
    const auto e_problem_face = ( key_cnr_nodes != key_e1_nodes ) ? shared_faces.second : shared_faces.first;
    cout <<"\n"<<"parentElements: only Element "<< shared_parents[e_matched]->Idx() <<" contains the nodes:\n\t\t        ";
    for ( const auto& nit : cnr_node_vec ) cout <<" "<< nit->Idx();
    cout <<"\n\t"<<"Element "<< shared_parents[e_not_found]->Idx() <<": ";
    cout << parseFiniteElementType( shared_parents[e_not_found]->FE_Type() ) <<": face: "<< e_problem_face <<" cannot be matched.";
    cout <<" The nodes of its faces are:";
    for ( uint32_t face{0u}; face < shared_parents[e_not_found]->Faces(); ++face ) {
         cout <<"\n\t\t"<<"face: "<< face <<": ";
         for ( auto nd : shared_parents[e_not_found]->FE()->NodesOfFace(face) )
           cout << shared_parents[e_not_found]->N(nd)->Idx() <<" ";
      }
    cout << endl;
    csmp_error.Note( ERROR, "parentElements", "can you see any of the faces matching nodes that have been missed?" );
#endif

    throw csmp::Exception( ERROR, "parentElements", "found 1 face-matched element but for element 2 the nodes do not match");
    
    // whatever has been found
    return make_pair( make_pair( shared_parents[0], shared_faces.first ), make_pair( shared_parents[1], shared_faces.second ) );
    
 } // end parentElements

template pair<pair<Element<3>*,uint32_t>,pair<Element<3>*,uint32_t>>  parentElements( typename vector<Node<3>*>::const_iterator, typename std::vector<Node<3>*>::const_iterator );
template pair<pair<Element<2>*,uint32_t>,pair<Element<2>*,uint32_t>>  parentElements( typename vector<Node<2>*>::const_iterator, typename std::vector<Node<2>*>::const_iterator );
template pair<pair<Element<1>*,uint32_t>,pair<Element<1>*,uint32_t>>  parentElements( typename vector<Node<1>*>::const_iterator, typename std::vector<Node<1>*>::const_iterator );


/* CUTOUT CODE DUPLICATING OTHER FUNCTIONALITY

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
*/






/**
    Eliminates surface and line elements from cell range in 3D, and line elements from cell range in 2D.
    
    @param elmt_pointers range to elements of potentially different types
    @return number of remaining cells
    
    @attention access to entire vector is necessary because erase() is a member  function of vector
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t eraseLowerDimensionalOrInvalidCells( vector<CELL<dim>*>& elmt_pointers )
 {
    if ( elmt_pointers.empty() ) return 0ul;
    
    // eliminating potential lower-dimensional elements or nullprts from result vector
    auto new_end = remove_if( elmt_pointers.begin(), elmt_pointers.end(),
                              []( const Element<dim>* eit )
                               {
                                  if constexpr ( dim == 3U )
                                    return ( eit==nullptr || !eit->IsVolume() );
                                  else if constexpr ( dim == 2U )
                                    return ( eit==nullptr || !eit->IsSurface() );
                                  else
                                    return ( eit==nullptr );
                               } );
                               
    // deleting cells beyond the vectors new end
    elmt_pointers.erase( new_end, elmt_pointers.end() );

    return elmt_pointers.size();
    
 } // end eraseLowerDimensionalOrInvalidCells

template size_t eraseLowerDimensionalOrInvalidCells( vector<Element<3>*>& );
template size_t eraseLowerDimensionalOrInvalidCells( vector<Element<2>*>& );
template size_t eraseLowerDimensionalOrInvalidCells( vector<Element<1>*>& );



/**
       Erases cells from the vector that have a different dimension than the current cell.
       @return the number of cells remaining in the vector
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t eraseDifferentDimensionalOrInvalidCells( const CELL<dim>* const cptr, vector<CELL<dim>*>& elmt_pointers )
 {
    if ( elmt_pointers.empty() ) return 0ul;
    if constexpr ( dim == 1U ) return 0ul;
    
    // TODO: there should be a virtual function in FiniteElement that returns cellshape!
    const CELL_SHAPE shape = parseFiniteElementDimension( cptr->FE_Type() );
    
    // eliminating potential lower-dimensional elements or nullprts from result vector
    auto new_end = remove_if( elmt_pointers.begin(), elmt_pointers.end(),
                              [shape]( const Element<dim>* eit )
                               {
                                  if ( eit==nullptr ) return true;
                                  if constexpr ( dim == 3U ) {
                                       if ( shape == VOLUME  && !eit->IsVolume() ) return true;
                                       if ( shape == SURFACE && !eit->IsSurface() ) return true;
                                       if ( shape == LINE    && !eit->IsLine() ) return true;
                                    }
                                  else if constexpr ( dim == 2U ) {
                                       if ( shape == SURFACE && !eit->IsSurface() ) return true;
                                       if ( shape == LINE    && !eit->IsLine() ) return true;
                                    }
                                  return false;
                               } );
                               
    // deleting cells beyond the vectors new end
    elmt_pointers.erase( new_end, elmt_pointers.end() );

    return elmt_pointers.size();
    
 } // end eraseDifferentDimensionalOrInvalidCells

template size_t eraseDifferentDimensionalOrInvalidCells( const Element<3>* const, vector<Element<3>*>& );
template size_t eraseDifferentDimensionalOrInvalidCells( const Element<2>* const, vector<Element<2>*>& );
template size_t eraseDifferentDimensionalOrInvalidCells( const Element<1>* const, vector<Element<1>*>& );




/**
     Finds the number of the cell face which consists of the supplied range of nodes.
     
     @return the local face number 0..faces-1, or uint32_t::max
*/
template<uint32_t dim, template<uint32_t> class CELL>
uint32_t faceWithCornerNodes( const CELL<dim>* const cell_ptr,
                              typename vector<Node<dim>*>::const_iterator first,
                              typename vector<Node<dim>*>::const_iterator last )
  {
     assert( cell_ptr != nullptr );
     
     set<Node<dim>*> node_set( first, last );
     for ( uint32_t face{0u}; face<cell_ptr->Faces(); ++face )
       if ( cell_ptr->CornerNodesOfFace(face) == node_set ) return face;
        
     return numeric_limits<uint32_t>::max();
     
  } // end faceWithCornerNodes
      
template uint32_t faceWithCornerNodes( const Element<3U>* const,
                                       typename vector<Node<3U>*>::const_iterator,
                                       typename vector<Node<3U>*>::const_iterator );
template uint32_t faceWithCornerNodes( const Element<2U>* const,
                                       typename vector<Node<2U>*>::const_iterator,
                                       typename vector<Node<2U>*>::const_iterator );
template uint32_t faceWithCornerNodes( const Element<1U>* const,
                                       typename vector<Node<1U>*>::const_iterator,
                                       typename vector<Node<1U>*>::const_iterator );
     
      

/**
       Returns pointer to the highest dimensional element that has a face with with the supplied range of nodes.
       When more than a single element is found, this may mean that the Element found is not on a model boundary.
       If so, the method selects the element with the least number of neighbors as the return value.
       
       @param first and last are iterators to the nodes of the Face.
       
       @return pair of Element and  the number of its face which shares the supplied nodes.
       
       @attention method works only if all the parent element pointers that the node stores are valid.
       @attention method only returns single element even if there are multiple ones that need covering
       but there can always only be a single element that matches the input lower-dim element that has the supplied nodes.
*/
template<uint32_t dim>
pair<Element<dim>*,uint32_t> parentElement( typename vector<Node<dim>*>::const_iterator first,
                                            typename vector<Node<dim>*>::const_iterator last )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    // 0. Find all parent elements that contain the supplied range of nodes
    // --------------------------------------------------------------------
    vector<Element<dim>*>  elmts_with_all_nodes = parentElementsContaining<dim>( first, last );

    // 1. If more than a single parent element shares all nodes
    // --------------------------------------------------------
    if ( elmts_with_all_nodes.size() > 1 )
      {
#ifndef NDEBUG
         csmp_error.Note( WARNING, "parentElement()", "more than one element containing all nodes found.");
#endif
         // 1.1 trying erasing potential lower-dimensional elements or nullprts from result vector
         if ( eraseLowerDimensionalOrInvalidCells( elmts_with_all_nodes ) == 0 ) {
             if  constexpr( dim == 3 ) {
                   csmp_error.Note( ERROR, "parentElement",
                                           "no volume element containing all of the nodes in the supplied range found");
                   return make_pair( nullptr, numeric_limits<uint32_t>::max() );
                }
             else if  constexpr( dim == 2 ) {
                   csmp_error.Note( ERROR, "parentElement",
                                           "no surface element containing all of the nodes in the supplied range found");
                   return make_pair( nullptr, numeric_limits<uint32_t>::max() );
                }
           }
         // 1.2 extra diagnostics if more than a single higher-dimensional element was found this function is unsuitable,
         //     but diagnostics will be offered before returning
         if ( elmts_with_all_nodes.size() == 2 )
           {
                 string node_numbers;
                 while ( first != last ) { node_numbers += to_string( (*first)->Idx() ); node_numbers +=","; first++; }
                 // checking whether the nodes belong to a face inside of the model
                 auto faces = findAdjacentFacesFromNeighbors( elmts_with_all_nodes[0], // FASTER
                                                              elmts_with_all_nodes[1] );
                                                              
                cout <<"\n"<<"parentElement: supplied nodes "<< node_numbers <<" lie on the model inside between "<< endl;
                cout <<"\t"<< parseAbbreviated_FE_Type(elmts_with_all_nodes[0]->FE_Type());
                cout <<":"<< elmts_with_all_nodes[0]->Idx() <<" face:"<< faces.first;
                cout <<" and "<< parseAbbreviated_FE_Type(elmts_with_all_nodes[1]->FE_Type());
                cout <<":"<< elmts_with_all_nodes[1]->Idx() <<" face:"<< faces.second << endl;
                cout.flush();
                   
                csmp_error.Note( WARNING, "parentElement: two valid parents found:",
                                "use parentElements() to handle model-interior Face with 2 equidimensional parent elements");

                return make_pair( nullptr, numeric_limits<uint32_t>::max() );
            }
            
      } // end more than one parent case

    assert( elmts_with_all_nodes.size() == 1U );
    
    
    // 2. Find the parent element's face that consists of the nodes
    // ------------------------------------------------------------
    uint32_t face_id = faceWithCornerNodes( elmts_with_all_nodes[0], first, last );
                              
     // when the face has been found the result is returned
     if ( face_id != numeric_limits<uint32_t>::max() )
       return make_pair( elmts_with_all_nodes[0], face_id );
 
   csmp_error.Note( ERROR, "parentElement", "could not find face shared with higher dimensional element");

   // nothing useful
   return make_pair( nullptr, numeric_limits<uint32_t>::max() );
    
 } // end parentElement1

template pair<Element<3>*,uint32_t> parentElement( typename vector<Node<3>*>::const_iterator,
                                                   typename vector<Node<3>*>::const_iterator );
template pair<Element<2>*,uint32_t> parentElement( typename vector<Node<2>*>::const_iterator,
                                                   typename vector<Node<2>*>::const_iterator );
template pair<Element<1>*,uint32_t> parentElement( typename vector<Node<1>*>::const_iterator,
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
    // neighbor nodes
    total_size += nptr->Neighbors() * sizeof( Node<dim>* );
    // manifold
    total_size += sizeof( NodeManifold<dim>* );
    // coordinates
    total_size += sizeof( Point<dim> );
    // idx
    total_size += sizeof(size_t);
    // flags
    total_size += sizeof(BOX_BOUNDARY);
    total_size += sizeof(TOPOTYPE);
    // + local variable storage
    total_size += sizeof(nptr->LVS());
    
    return total_size;

  } // end sizeOf

template size_t sizeOf( const Node<3>* const );
template size_t sizeOf( const Node<2>* const );
template size_t sizeOf( const Node<1>* const );





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

