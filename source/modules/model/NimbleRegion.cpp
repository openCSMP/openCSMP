#include "NimbleRegion.h"
#include "ModelSubDomain.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "PDE_Integrator.h"

using namespace std;

//#define NIMBLE_REGION_DEBUG

namespace csmp {

/**
   construction from nodes only, relying on existing node-parent-element connectivity to identify elements
 
   the assumption is made that the node pointers are unique.
 
   the element vector is also sorted for searching.
 
   the elements are identified by looping over the parents of the nodes,
   then additional nodes are added where necessary to capture a halo for the pressure equation.
 
*/
template<uint32_t dim>
NimbleRegion<dim>::NimbleRegion( const PropertyDatabase<dim>& p,
                                 typename vector<Node<dim>*>::const_iterator first,
                                 typename vector<Node<dim>*>::const_iterator last )
 : ModelSubDomain<dim,Element>( "NimbleRegion", p ),
   n_interior_nodes_(0U)
 {
    //cout <<"\nNimbleRegion<"<< dim <<">(constructor): building region from ";
    //cout << distance(first,last) <<" nodes.\n";
   
    // universal function that is most efficient when the region needs to be build from scratch
    Update3( first, last );

 } // end (constructur)








/**
    As above, but perimeter nodes are identified via neighborhood relations of elements.
*/
template<uint32_t dim>
void NimbleRegion<dim>::Update( typename vector<Node<dim>*>::const_iterator first,
                                typename vector<Node<dim>*>::const_iterator last )
 {
    const size_t n_target_nodes( distance(first,last));
    if (  n_target_nodes == 0U )
      throw csmp::Exception( ERROR, "NimbleRegion<dim>::Update", "Supplied node range is empty; nothing was done." );
 
    // creating unique node and element sets for the nimble regions
    // the node set will also be used to find the perimeter nodes
    set<Element<dim>*>   elmts;
    const auto nodes_end(last);
    for ( auto nit(first); nit!=nodes_end; ++nit ) {
         const auto parent_elements((*nit)->Parents());
         for ( auto i{0U}; i<parent_elements; ++i ) {
              assert( (*nit)->Parent(i) != nullptr );
              elmts.insert( (*nit)->Parent(i) );
           }
      }
    // replacing element vector by the elements pointers contained in the set

    elements_.assign( elmts.begin(), elmts.end() );
    // ? not needed because already sorted:
    sort( elements_.begin(), elements_.end() );
   
    // check that the elements only have corner nodes, i.e. linear interpolation
    if ( elements_[0U]->FE()->Interpolation() != 1 )
      throw csmp::Exception( FATAL_ERROR, "NimbleRegion(constructor):",
                             "this implementation works only for straight-sided elements with linear interpolation functions.");


    // finding the boundary nodes, as those nodes that are on element faces that have no or only an outside-region neighbor
    set<Node<dim>*> perimeter_nodes;
   
    for ( auto& it : elements_ ) {
         const auto neighbors(it->Neighbors());
         for ( auto i{0U}; i<neighbors; ++i )
            if ( it->Neighbor(i) == nullptr ||
                !binary_search( elements_.begin(), elements_.end(), it->Neighbor(i) ) ) {
                 // collecting the perimeter nodes into a set
                 for ( const auto& j : it->FE()->NodesOfFace(i) )
                   perimeter_nodes.insert( it->N(j) );
              }
            //else cerr <<"i ";
      }
   
    // putting the nodes into sorted vector with interior and perimeter ranges
    // destructing existing node and element members, but keeping potential storage
    if ( !nodes_.empty() ) nodes_.clear();
    nodes_.reserve(n_target_nodes);

#ifdef NIMBLE_REGION_DEBUG
    size_t  counter(0U);
#endif
   
    auto nit(first);
    while ( nit != nodes_end ) {
#ifdef NIMBLE_REGION_DEBUG
         // debugging: labeling the input nodes continuously
         (*nit)->Idx( counter++ );
#endif
         // storing interior nodes
         if ( (*nit)->AtBoundary() == NOT ) nodes_.push_back( (*nit) );
         nit++;
      }
   
    // making the interior node vector searchable by sorting
    n_interior_nodes_ = nodes_.size();
    sort( nodes_.begin(), nodes_.end() );


// SHOULD NOT BE NECESSARY (but else, the interior nodes are overwritten by superfluous perimeter nodes???
for ( auto nit2(first); nit2!=nodes_end; ++nit2 )
  perimeter_nodes.erase( (*nit2) );

#ifdef NIMBLE_REGION_DEBUG
    // debugging: labeling the perimeter nodes continuously
    for ( typename set<Node<dim>*>::iterator nit=perimeter_nodes.begin(); nit!=perimeter_nodes.end(); ++nit )
      (*nit)->Idx( counter++ );
#endif

    // inserting the sorted perimeter nodes at the end of the node vector so that it now consists of 2 sorted ranges
    nodes_.reserve( n_interior_nodes_ + perimeter_nodes.size() );
    const typename set<Node<dim>*>::iterator perimeter_nodes_end(perimeter_nodes.end());
    for ( typename set<Node<dim>*>::iterator nit2=perimeter_nodes.begin(); nit2!=perimeter_nodes_end; ++nit2 )
      // if the node is located on the model boundary, it is a perimeter node
      nodes_.push_back( (*nit2) );

 } // end Update




/**
    Removes all storage and recreates region from scratch.

    construction from nodes only, relying on existing node-parent-element connectivity to identify elements
 
    the assumption is made that the node pointers are unique.
 
    the element vector is also sorted for searching.
 
    the elements are identified by looping over the parents of the nodes,
    then additional nodes are added where necessary to capture a halo for the pressure equation.
 
    @todo PERIMETER IDENTIFICATION DOES NOT WORK FOR DISCONTIGOUS REGIONS YET
*/
template<uint32_t dim>
void NimbleRegion<dim>::Update2( typename vector<Node<dim>*>::const_iterator first,
                                      typename vector<Node<dim>*>::const_iterator last )
 {
    const size_t n_target_nodes( distance(first,last));
    if (  n_target_nodes == 0U )
      throw csmp::Exception( ERROR, "NimbleRegion<dim>::Update2", "Supplied node range is empty; nothing was done." );
 
    // destructing existing node and element members, but keeping potential storage
    if ( !nodes_.empty() ) nodes_.clear();
    if ( !elements_.empty() ) elements_.clear();
    nodes_.reserve(n_target_nodes);

#ifdef NIMBLE_REGION_DEBUG
    size_t  counter(0U);
#endif
   
    // recording the 'interior' nodes, excluding nodes at model boundary
    while ( first != last ) {
#ifdef NIMBLE_REGION_DEBUG
         // debugging: labeling the input nodes continuously
         (*first)->Idx( counter++ );
#endif
         // storing interior nodes
         if ( (*first)->AtBoundary() == NOT ) nodes_.push_back( (*first) );
         first++;
      }
   
    // making the interior node vector searchable by sorting
    n_interior_nodes_ = nodes_.size();
    sort( nodes_.begin(), nodes_.end() );
   
    // creating unique node and element sets for the nimble regions
    // the node set will also be used to find the perimeter nodes
    set<Element<dim>*>  elmts;
    set<Node<dim>*>  perimeter_nodes;
    const auto nodes_end(nodes_.end());
    for ( auto& nit : nodes_ )
       for ( auto i{0U}; i<nit->Parents(); ++i ) {
            assert( nit->Parent(i) != nullptr );
            auto eit = elmts.insert( nit->Parent(i) );
            // if this is element was not encountered before, we keep track of its nodes
            if ( eit.second ) {
                const auto nodes(nit->Parent(i)->Nodes());
                for ( auto j{0U}; j<nodes; ++j )
                    // but only those nodes that are not contained in the interior
                    if ( !binary_search( nodes_.begin(), nodes_.end(), nit->Parent(i)->N(j) ) ) {
                         // storing the unique perimeter nodes
                         perimeter_nodes.insert( nit->Parent(i)->N(j) );
#ifdef NIMBLE_REGION_DEBUG
                          // debugging
                         if ( perimeter_nodes.find(nit->Parent(i)->N(j)) != perimeter_nodes.end() )
                           nit->Parent(i)->N(j)->Idx( counter++ );
#endif
                      }
                }
         }
    // transfer element set to a sorted vector
    elements_.assign( elmts.begin(), elmts.end() );
    // sort( elements_.begin(), elements_.end() ); // needed?
   
    // check that the elements only have corner nodes, i.e. linear interpolation
    if ( elements_[0U]->FE()->Interpolation() != 1 )
      throw csmp::Exception( FATAL_ERROR, "NimbleRegion(constructor):",
                             "this implementation works only for straight-sided elements with linear interpolation functions.");

    // inserting the sorted perimeter nodes at the end of the node vector so that it now consists of 2 sorted ranges
    nodes_.reserve( n_interior_nodes_ + perimeter_nodes.size() );
    const typename set<Node<dim>*>::const_iterator perimeter_nodes_end(perimeter_nodes.end());
    for ( auto nit=perimeter_nodes.begin(); nit!=perimeter_nodes_end; ++nit )
      // if the node is located on the model boundary, it is a perimeter node
      nodes_.push_back( (*nit) );

 } // end Update2




/**
    combine update and update2.
*/
template<uint32_t dim>
void NimbleRegion<dim>::Update3( typename vector<Node<dim>*>::const_iterator first,
                                      typename vector<Node<dim>*>::const_iterator last )
 {
    const size_t n_target_nodes( distance(first,last));
    if (  n_target_nodes == 0U )
      throw csmp::Exception( ERROR, "NimbleRegion<dim>::Update", "Supplied node range is empty; nothing was done." );

    // destructing existing node members, but keeping potential storage
    if ( !nodes_.empty() ) nodes_.clear();
    nodes_.reserve(n_target_nodes);

#ifdef NIMBLE_REGION_DEBUG
    size_t  counter(0U);
#endif
   
    auto       nit(first);
    const auto nodes_end(last);
    
    while ( nit != nodes_end ) {
#ifdef NIMBLE_REGION_DEBUG
         // debugging: labeling the input nodes continuously
         (*nit)->Idx( counter++ );
#endif
         // storing interior nodes
         /*if ( (*nit)->AtBoundary() == NOT )*/ nodes_.push_back( (*nit) );
         nit++;
    }
   
    // making the interior node vector searchable by sorting
    n_interior_nodes_ = nodes_.size();
    
    sort( nodes_.begin(), nodes_.end() );
    
    // creating unique node and element sets for the nimble regions
    // the node set will also be used to find the perimeter nodes
    set<Element<dim>*>  elmts;
    set<Node<dim>*>  perimeter_nodes;
    //for ( typename vector<Node<dim>*>::iterator nit=nodes_.begin(); nit!=nodes_end; ++nit )
    for ( auto& nit2 : nodes_ )
       for ( auto i{0U}; i<nit2->Parents(); ++i ) {
            assert( nit2->Parent(i) != nullptr );
            auto eit = elmts.insert( nit2->Parent(i) );
            // if this is element was not encountered before, we keep track of its nodes
            if ( eit.second ) {
                const auto nodes(nit2->Parent(i)->Nodes());
                for ( auto j{0U}; j<nodes; ++j )
                    // but only those nodes that are not contained in the interior
                    if ( !binary_search( nodes_.begin(), nodes_.end(), nit2->Parent(i)->N(j) ) ) {
                         // storing the unique perimeter nodes
                         perimeter_nodes.insert( nit2->Parent(i)->N(j) );
                      }
                }
         }
    // transfer element set to a sorted vector
    elements_.assign( elmts.begin(), elmts.end() );
    sort( elements_.begin(), elements_.end() ); // needed?
    
    // check that the elements only have corner nodes, i.e. linear interpolation
    if ( elements_[0U]->FE()->Interpolation() != 1 )
      throw csmp::Exception( FATAL_ERROR, "NimbleRegion(constructor):",
                             "this implementation works only for straight-sided elements with linear interpolation functions.");
    

    // finding addtional boundary nodes, as those nodes that are on element faces that have no or only an outside-region neighbor
    for ( auto& it : elements_ ) {
         const auto neighbors(it->Neighbors());
         //cout<<"    neighbors = "<<neighbors<<endl;
         for ( auto i{0U}; i<neighbors; ++i )
            if ( it->Neighbor(i) == nullptr ||
                !binary_search( elements_.begin(), elements_.end(), it->Neighbor(i) ) ) {
                 for ( const auto& j : it->FE()->NodesOfFace(i) )
                   perimeter_nodes.insert( it->N(j) );
              }
            //else cerr <<"i ";
      }

    // SHOULD NOT BE NECESSARY (but else, the interior nodes are overwritten by superfluous perimeter nodes???
    for ( auto nit2(first); nit2!=nodes_end; ++nit2 )
        perimeter_nodes.erase( (*nit2) );

#ifdef NIMBLE_REGION_DEBUG
    // debugging: labeling the perimeter nodes continuously
    for ( typename set<Node<dim>*>::iterator nit=perimeter_nodes.begin(); nit!=perimeter_nodes.end(); ++nit )
      (*nit)->Idx( counter++ );
#endif

    // inserting the sorted perimeter nodes at the end of the node vector so that it now consists of 2 sorted ranges
    nodes_.reserve( n_interior_nodes_ + perimeter_nodes.size() );
    const typename set<Node<dim>*>::iterator perimeter_nodes_end(perimeter_nodes.end());
    for ( auto nit2=perimeter_nodes.begin(); nit2!=perimeter_nodes_end; ++nit2 )
      // if the node is located on the model boundary, it is a perimeter node
      nodes_.push_back( (*nit2) );
 } // end Update





// USER METHODS

// retrieving information

template<uint32_t dim>
size_t NimbleRegion<dim>::Nodes() const
 {
     return nodes_.size();
 }


template<uint32_t dim>
size_t NimbleRegion<dim>::InteriorNodes() const
 {
     return n_interior_nodes_;
 }


template<uint32_t dim>
size_t NimbleRegion<dim>::PerimeterNodes() const
 {
     return nodes_.size() - n_interior_nodes_;
 }


template<uint32_t dim>
size_t NimbleRegion<dim>::Cells() const
 {
     return elements_.size();
 }


/**
    Consecutive numbering, interior nodes first, perimeter nodes second.
*/
template<uint32_t dim>
size_t NimbleRegion<dim>::RenumberNodes() const
 {
     size_t node_number(0U);
     for ( auto& nit : nodes_ )
       nit->Idx( node_number++ );
   
     return node_number;
 }


// accessors
template<uint32_t dim>
Element<dim>* const NimbleRegion<dim>::E( size_t n )
 {
    assert( n<elements_.size() );
    return elements_[n];
 }


template<uint32_t dim>
const Element<dim>* const NimbleRegion<dim>::E( size_t n ) const
 {
    assert( n<elements_.size() );
    return elements_[n];
 }

 
template<uint32_t dim>
Node<dim>* const NimbleRegion<dim>::N( size_t n )
 {
    assert( n < nodes_.size() );
    return nodes_[n];
 }


template<uint32_t dim>
const Node<dim>* const NimbleRegion<dim>::N( size_t n ) const
 {
    assert( n < nodes_.size() );
    return nodes_[n];
 }
 



/**
    When you do not want to call the destructor yet.
*/
template<uint32_t dim>
void NimbleRegion<dim>::Clear() {
	nodes_.clear();
	elements_.clear();
	n_interior_nodes_ = 0U;
}




/**
    Deleting everything including the storage.
*/
template<uint32_t dim>
void NimbleRegion<dim>::Erase()
{
    nodes_.erase( nodes_.begin(), nodes_.end() );
    elements_.erase( elements_.begin(), elements_.end() );
    n_interior_nodes_ = 0U;
}





/**
    Print the NimbleRegion to screen.
*/
template<uint32_t dim>
void NimbleRegion<dim>::Out() const
 {
    cout <<"\nNimbleRegion<dim>::Out: "<< nodes_.size() <<" nodes, "<< elements_.size() <<" elements.\n";
    cout <<"\t"<< n_interior_nodes_ <<" interior nodes.\n";
    cout <<"\nindices of member nodes:\ninterior: ";
    for ( auto nit=nodes_.begin(); nit!=next(nodes_.begin(),n_interior_nodes_); nit++ )
      cout << (*nit)->Idx() <<" ";
    cout <<"\nperimeter: ";
    for ( auto nit=next(nodes_.begin(),n_interior_nodes_); nit!=nodes_.end(); nit++ )
      cout << (*nit)->Idx() <<" ";
    cout <<"\n";
   
    cout <<"\nindices of member elements (not renumbered):\n";
    for ( auto it : elements_ )
      cout << it->Idx() <<" ";
    cout <<"\n";
   
    cout << endl;

 } // end Out




// EXPERIMENTAL METHODS

/**
    Assuming that there was no remeshing, the NimbleRegion may only have grown at its perimeter
    and former interior nodes will remain interior nodes.
 
    Where new nodes are added, these and former perimeter nodes may become interior nodes
    because elements that they belong to will also be added, introducing additional nodes.
    It follows that the interior - perimeter ordering has to be revisited.
 
    New elements will be added in the process.
*/
/*
template<uint32_t dim>
void NimbleRegion<dim>::Grow( typename vector<Node<dim>*>::iterator first,
                              typename vector<Node<dim>*>::iterator last )
{
   // inserting the new (exterior) nodes at the end of the nodes vector
   copy( first, last, inserter(nodes_, nodes_.end()));
  
   // checking whether former perimeter nodes have become interior nodes,
   // looking into the region from the new perimeter nodes
   // TODO:
  
  
} // end Grow
*/


template class NimbleRegion<1>;
template class NimbleRegion<2>;
template class NimbleRegion<3>;

} // end namespace csmp
