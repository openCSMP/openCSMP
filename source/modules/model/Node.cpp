// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Node.h"
#include "Element.h"
#include "NodeParentElementVector.h"
#include "Visitor.h"
#include "NodeManifold.h"
#include "ErrorHandler.h"
#include "meshManagementUtilities.h"
#include "ConvexPolygon.h"

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
    parents_(nd.parents_),
    neighbor_node_pointers_(nd.neighbor_node_pointers_),
    manifold_(nd.manifold_),
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
Node<dim>::Node( Node<dim>&& nd ) noexcept
  : xyz_{ std::move(nd.xyz_) },
    idx_{ nd.idx_ },
    parents_{ std::move(nd.parents_) },
    neighbor_node_pointers_{ std::move(nd.neighbor_node_pointers_) },
    manifold_{ nd.manifold_ },
    at_boundary_{ nd.at_boundary_ },
    BREP_entity_{ nd.BREP_entity_ }
  {
    this->LVS( nd.LVS() );
    
    nd.manifold_ = nullptr;
  }




template<uint32_t dim>
Node<dim>& Node<dim>::operator=( const Node<dim>& nd )
 {
    if ( &nd != this ) 
      {
         xyz_                     = nd.xyz_;
         idx_                     = nd.idx_;
         parents_                 = nd.parents_;
         neighbor_node_pointers_  = nd.neighbor_node_pointers_;
         manifold_                = nd.manifold_;
         at_boundary_             = nd.at_boundary_;
         BREP_entity_             = nd.BREP_entity_;
         
         this->LVS( nd.LVS() );
      }
    return *this;
 }



/**
    Move assignment, relying on that similar operators exist for the nodes components.
    
    @note this assumes that the supplied node is a temporary.
*/
template<uint32_t dim>
Node<dim>& Node<dim>::operator=( Node<dim>&& nd ) noexcept
 {
    assert( this != &nd );
    xyz_                     = nd.xyz_;
    idx_                     = nd.idx_;
    at_boundary_             = nd.at_boundary_;
    BREP_entity_             = nd.BREP_entity_ ;
    parents_                 = std::move(nd.parents_);
    neighbor_node_pointers_  = std::move(nd.neighbor_node_pointers_);
    manifold_                = nd.manifold_;

    this->LVS( nd.LVS() );
    
    nd.manifold_ = nullptr;
 
    return *this;
 }






/**
    Tries to distinguish 2 nodes from one another.
    
    This comparison operator was designed specifically for the creation of particular Region, Boundary and SplitBoundary objects.
    Hence, only information important for this process is taken into account, distinguishing 2 Nodes from each other.
    That must be coordinate and parent elements
*/
template<uint32_t dim>
bool  Node<dim>::operator==( const Node<dim>& nd ) const noexcept
 {
    if ( &nd == this ) return true;
    if ( idx_ != nd.idx_ ) return false;
    if ( BREP_entity_ != nd.BREP_entity_ ) return false;
    // collocation is not a criterion:  if ( distance( Coordinate(), nd.Coordinate() ) > numeric_limits<double>::epsilon() ) return false;
    if ( AtBoundary() != nd.AtBoundary() ) return false;
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
    assert( pnode <= FIFTY );
    assert( element != nullptr );
    parents_.Assign( element, pnode );

 } // end Assign




/**
    Sets the pointer to the argument element to zero, and the corresponding node number to NOT_INITIALIZED.
    
    @note the argument element is not deleted. To do this use EraseNullPointerParents().
*/
template<uint32_t dim>
void Node<dim>::Unassign( const Element<dim>* const element )
 {
    parents_.Unassign( element );
    
 } // end Unassign




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
    assert( neighbor_nodes.count(this) == 0 );
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
 }
 
 
 
template<uint32_t dim>
void  Node<dim>::Assign( vector<Node<dim>*>& neighbor_nodes, bool sort_neighbors )
 {
    assert(find(neighbor_nodes.begin(), neighbor_nodes.end(), this) == neighbor_nodes.end());
    neighbor_node_pointers_.assign( neighbor_nodes.begin(), neighbor_nodes.end() );
    if ( sort_neighbors )
      sort( neighbor_node_pointers_.begin(), neighbor_node_pointers_.end() );
 }



    /// copies property values from the argument node to the current node
template<uint32_t dim>
void  Node<dim>::CopyPropertyValuesFrom( Node<dim>& nd ) noexcept
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
    // needed by unique() see below
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
bool  Node<dim>::IsNeighbor( const Node<dim>* const nptr ) const noexcept
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
    Moves  unwanted Node  to the end of the vector before erasing it.
    @attention Use UpdateNeighbors to shrink vector to new size.
*/
template<uint32_t dim>
void Node<dim>::RemoveNeighbor( const Node<dim>* const neighbor_node )
 {
    neighbor_node_pointers_.erase( remove( neighbor_node_pointers_.begin(),
                                           neighbor_node_pointers_.end(), neighbor_node ),
                                   neighbor_node_pointers_.end() );
 }
 
 
template<uint32_t dim>
void Node<dim>::EraseNeighbors()
 {
    vector<Node<dim>*>().swap(neighbor_node_pointers_);
 }




/**

   If the parent elements are of the same dimension as the node,
   then there are exactly as many elements as are nodes.
   
   @attention calculation is costly, compute only once per node loop!
   
   @test correct SKM (unit test exists)
*/
template<uint32_t dim>
uint32_t  Node<dim>::Neighbors() const noexcept
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
Node<dim>*  Node<dim>::Neighbor( uint32_t neighbor_node ) const noexcept
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
Point<dim>  Node<dim>::Coordinate() const  noexcept { return xyz_; }

template<uint32_t dim>
 void Node<dim>::Coordinate( const Point<dim>& p )  noexcept { xyz_=p; }





template<uint32_t dim>
void  Node<dim>::EraseParents()
  {
     parents_.Erase();
  }


template<uint32_t dim>
uint32_t   Node<dim>::Parents() const noexcept
  { return static_cast<uint32_t>(parents_.Size()); }


template<uint32_t dim>
uint32_t  Node<dim>::ParentNodeNumber( uint32_t parent_element_number ) const noexcept
  {
     assert( parent_element_number < parents_.Size() );
     return parents_.LocalNodeNumber( parent_element_number );
  }


template<uint32_t dim>
Element<dim>* Node<dim>::Parent( uint32_t parent_element_number ) noexcept
  {
     assert( parent_element_number < parents_.Size() );
     return parents_.ParentElement( parent_element_number );
  }

template<uint32_t dim>
const Element<dim>* const  Node<dim>::Parent( uint32_t parent_element_number ) const noexcept
  {
     assert( parent_element_number < parents_.Size() );
     return parents_.ParentElement( parent_element_number );
  }


    /// checks whether Element is a parent of the node
template<uint32_t dim>
bool  Node<dim>::IsParent( const Element<dim>* const eptr ) const noexcept
  {
     return parents_.IsParent(eptr);
  }







template<uint32_t dim>
void  Node<dim>::Idx( size_t idx_to_assign ) const noexcept
 {
    idx_ = idx_to_assign;
 }


template<uint32_t dim>
size_t   Node<dim>::Idx() const noexcept
 {
    return idx_;
 }



/**
   Set or return the global boundary flag of the Node.
*/
template<uint32_t dim>
void Node<dim>::AtBoundary( BOX_BOUNDARY b ) noexcept { at_boundary_ = b; }


template<uint32_t dim>
BOX_BOUNDARY  Node<dim>::AtBoundary() const noexcept { return at_boundary_; }

/// Set or return the coordinates of the current node.
template<uint32_t dim>
void            Node<dim>::x( double xc ) noexcept  { xyz_[0u] = xc; }

template<uint32_t dim>
void            Node<dim>::y( double yc ) noexcept  { xyz_[1u] = yc; }

template<uint32_t dim>
void            Node<dim>::z( double zc ) noexcept  { xyz_[2u] = zc; }


template<uint32_t dim>
double          Node<dim>::x() const noexcept { return xyz_[0u]; }

template<uint32_t dim>
double          Node<dim>::y() const noexcept { return xyz_[1u]; }

template<uint32_t dim>
double          Node<dim>::z() const noexcept { return xyz_[2u]; }


// MANIFOLDS

/// access to manifold if any; returns nullptr if the node is not a manifold
template<uint32_t dim>
bool Node<dim>::IsManifold() const noexcept { return (manifold_ != nullptr); }


template<uint32_t dim>
NodeManifold<dim>* const Node<dim>::Manifold() const noexcept { return manifold_; }

// node cannot be member of multiple manifolds at the same time
template<uint32_t dim>
void Node<dim>::Assign( NodeManifold<dim>& nmf ) noexcept
 {
    manifold_ = &nmf;
 }

// node cannot be member of multiple manifolds at the same time
template<uint32_t dim>
void Node<dim>::Disconnect() noexcept
 {
    manifold_ = nullptr;
 }


/**
   If node is part of line elements (2D) or surface elements (3D), method returns a unit normal that represents the average of the normals of the connected elements.
   
   @return false if 1) the node is not located on the inside of a patch of lower dimensional elements, 2) not properly initialised or 3) adjacent element normals are pointing in opposite directions
*/
template<uint32_t dim>
bool Node<dim>::UnitNormal( Point<dim>& avg_nrml ) const noexcept
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






/**
    Computes normal to Node from its neighbor nodes, but ignoring neighbors that are not contained with the supplied Node vector  (as in computer graphics)
    Supplying this vector allows to narrow the calculation to nodes that lie on a line or surface.
    
    @note The method was created to help check the consistency of the node numbering of lower-dimensional elements representing surfaces prior to converting them into boundaries.
    
    @param first iterator to first node in the restricted node range that shall be considered (so that a geometric element can be captured)
    @param last iterator to last node in the restricted node range that shall be considered (so that a geometric element can be captured)
    @return either the vertex normal or a Point object initialised with signalling NaN
    
    @attention std::find() is used so that the supplied Node range does not have to be sorted (or may consist of 2 sorted ranges, like in a model subdomain)
    
    Explanation of algorithm:  for each vertex, calculate the normal of the plane formed by the two edges entering and leaving that vertex.
    More formally, given vertices 𝐯1,𝐯2,…𝐯𝑛 with counterclockwise winding, define the normal at the 𝑖th vertex as:
    
         𝐧𝑖=(𝐯𝑖−𝐯𝑖−1)×(𝐯𝑖+1−𝐯𝑖)
         
    (where the indices wrap around). The algorithm is implemented using the MJL library in the  discretisation/ directory
    
    @author SKM
    @date 12/10/24
*/
template<uint32_t dim>
Point<dim> Node<dim>::VertexNormal( typename vector<Node<dim>*>::const_iterator first,
                                    typename vector<Node<dim>*>::const_iterator last ) const noexcept
 {
//     throw csmp::Exception( ERROR, "Node<dim>::VertexNormal", "method fails if node ordering on surface is not correct" );
 
     // 0. creating a subset of the node neighbor coordinates which lie on the feature of interest
     vector<Point<dim>> nbors_on_feature;
     nbors_on_feature.reserve( Neighbors() );

     for ( auto nit=NeighborsBegin(); nit!=NeighborsEnd(); ++nit )
       if ( find( first, last, (*nit) ) != last )
         nbors_on_feature.push_back( (*nit)->Coordinate() );
 
     // at least two neighbor nodes are needed to perform the normal construction (but may not be enough)
     if ( nbors_on_feature.size() <= 2 )
       return Point<dim>( numeric_limits<double>::signaling_NaN() );
       
     // 1. Order the vertices found so that they form an anti-clockwise convex polygon around the node/vertex of interest
     ConvexPolygon<dim> nbor_polygon( nbors_on_feature.begin(), nbors_on_feature.end() );

     // 2. compute the normal as an average of the spokes
     //    (where Point a is the coordinate of the node)
     return nbor_polygon.UnitNormal( Coordinate() );
     
  } // end VertexNormal




    /// which geometric part of the discretisation of the initial boundary representation (BREP) of the model geometry the node belongs to
 template<uint32_t dim>
 TOPOTYPE Node<dim>::Attribute() const noexcept {
     return BREP_entity_;
  }
    
template<uint32_t dim>
void Node<dim>::Attribute( TOPOTYPE geom_feature ) noexcept {
    BREP_entity_ = geom_feature;
 }




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
    if ( parents_.Size() > 0u ) {
         cout <<"\nElement objects sharing the node / node position therein:\n"<< endl;
         parents_.Out();
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
    
    // create sorted parent element vectors for the intersection algorithm
    // vector node 1
    vector<Element<dim>*>  elmts_with_all_nodes1( (*first)->Parents() );
    for ( uint32_t i{0u}; i<(*first)->Parents(); ++i ) elmts_with_all_nodes1[i] = (*first)->Parent(i);
    sort( elmts_with_all_nodes1.begin(), elmts_with_all_nodes1.end() );
    // vector node 2
    vector<Element<dim>*>  elmts_with_all_nodes2( (*next(first,1))->Parents() ), isect;
    for ( uint32_t i{0u}; i<(*next(first,1))->Parents(); ++i ) elmts_with_all_nodes2[i] = (*next(first,1))->Parent(i);
    sort( elmts_with_all_nodes2.begin(), elmts_with_all_nodes2.end() );
    
    // store the elements shared between the first and the second node in 'isect'
    set_intersection( elmts_with_all_nodes1.begin(), elmts_with_all_nodes1.end(),
                      elmts_with_all_nodes2.begin(), elmts_with_all_nodes2.end(),
                      back_inserter(isect) );
    first++;
    first++;
    
    // find the shared parent elements for all the supplied nodes
    while ( first != last ) {
         // vector 1
         elmts_with_all_nodes1 = isect;
         isect.clear();
         // vector 2
         elmts_with_all_nodes2.resize( (*first)->Parents() );
         for ( uint32_t i{0u}; i<(*first)->Parents(); ++i ) elmts_with_all_nodes2[i] = (*first)->Parent(i);
         sort( elmts_with_all_nodes2.begin(), elmts_with_all_nodes2.end() );
         set_intersection( elmts_with_all_nodes1.begin(), elmts_with_all_nodes1.end(),
                           elmts_with_all_nodes2.begin(), elmts_with_all_nodes2.end(),
                           back_inserter(isect) ); // custom comparitor used for the sets
         first++;
      }
      
    return isect;
      
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
    eface_node_vec1.reserve( cnr_node_vec.size() );
    for ( const auto& nd : shared_parents[0]->FE()->CornerNodesOfFace( shared_faces.first ) )
      eface_node_vec1.push_back( shared_parents[0]->N(nd) );
    
    // if there is a match the elements are returned in their current order
    if ( cnr_node_vec == eface_node_vec1 )
      return make_pair( make_pair( shared_parents[0], shared_faces.first ), make_pair( shared_parents[1], shared_faces.second ) );
      
    // else the second lot of face nodes are compared
    vector<Node<dim>*> eface_node_vec2;
    eface_node_vec2.reserve( cnr_node_vec.size() );
    for ( const auto& nd : shared_parents[1]->FE()->CornerNodesOfFace( shared_faces.second ) )
      eface_node_vec2.push_back( shared_parents[1]->N(nd) );
    if ( cnr_node_vec == eface_node_vec2 )
      return make_pair( make_pair( shared_parents[1], shared_faces.second ), make_pair( shared_parents[0], shared_faces.first ) );
      
    // interim diagnostics before cell matching by rotation is attempted
    // -----------------------------------------------------------------
#ifdef DEBUG_NODE_FUNCTIONS
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
    const size_t e_matched      = ( key_cnr_nodes != key_e1_nodes ) ? 0 : 1;
    const size_t e_not_found    = ( key_cnr_nodes != key_e1_nodes ) ? 1 : 0;
    const size_t e_problem_face = ( key_cnr_nodes != key_e1_nodes ) ? shared_faces.second : shared_faces.first;
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
    set<const Element<dim>*> parents;
    
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

