#include "NimbleRegion.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "PDE_Integrator.h"

using namespace std;

namespace csmp {


/**
   construction from nodes only, relying on existing node-parent-element connectivity to identify elements
 
   the assumption is made that the node pointers are unique.
 
   the element vector is also sorted for searching.
 
   the elements are identified by looping over the parents of the nodes,
   then additional nodes are added where necessary to capture a halo for the pressure equation.
 
*/
template<size_t dim>
NimbleRegion<dim>::NimbleRegion( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last )
 : n_interior_nodes_(0U)
 {
    cout <<"\nNimbleRegion<"<< dim <<">(constructor): building region from ";
    cout << distance(first,last) <<" nodes.\n";
    Rebuild( first, last );

 } // end (constructur)




/**
    Removes all storage and recreates region from scratch.

   construction from nodes only, relying on existing node-parent-element connectivity to identify elements
 
   the assumption is made that the node pointers are unique.
 
   the element vector is also sorted for searching.
 
   the elements are identified by looping over the parents of the nodes,
   then additional nodes are added where necessary to capture a halo for the pressure equation.
*/
template<size_t dim>
void NimbleRegion<dim>::Rebuild( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last )
 {
    if ( distance(first,last) == 0U )
      throw csmp::Exception( ERROR, "NimbleRegion<dim>::Rebuild", "Supplied node range is empty" );
 
    // allocating a reasonable amount of memory
    if ( !nodes_.empty() ) nodes_.clear();
    if ( !elements_.empty() ) elements_.clear();
    nodes_.reserve(distance(first,last));
   
    // assigning the 'interior' nodes, i.e. the ones without a halo, and excluding nodes at the model boundary
    while ( first != last ) {
         if ( (*first)->AtBoundary() == NOT ) nodes_.push_back( (*first) );
         first++;
      }
   
    // making the interior node vector searchable by sorting
    n_interior_nodes_ = nodes_.size();
    sort( nodes_.begin(), nodes_.end() );
   
    // making a set for finding the perimeter nodes, identifying a unique set of member elements as parents of interior nodes
    set<CellType*>  elmts;
    set<Node<dim>*>  marginal_nodes;
    for ( typename vector<Node<dim>*>::const_iterator nit=nodes_.begin(); nit!=nodes_.end(); ++nit )
       for ( size_t i=0; i<(*nit)->Parents(); ++i ) {
            pair<typename set<CellType*>::iterator,bool> eit = elmts.insert( (*nit)->Parent(i) );
            // if this is an element that was not added before, we keep track of its nodes
            if ( eit.second )
              for ( size_t j=0U; j<(*nit)->Parent(i)->Nodes(); ++j )
                marginal_nodes.insert( (*nit)->Parent(i)->N(j) );
         }
    // transfer creating sorted vector
    elements_.assign( elmts.begin(), elmts.end() );
    sort( elements_.begin(), elements_.end() ); // needed?
   
    // check that the elements only have corner nodes, i.e. linear interpolation
    if ( elements_[0U]->FE()->Interpolation() != 1 )
      throw csmp::Exception( FATAL_ERROR, "NimbleRegion(constructor):",
                             "this implementation works only for straight-sided elements with linear interpolation functions.");
   
    // separating the perimeter nodes into the corresponding vector
    const typename vector<Node<dim>*>::iterator interior_nodes_end(next(nodes_.begin(),n_interior_nodes_));
    for ( typename set<Node<dim>*>::const_iterator nit=marginal_nodes.begin(); nit!=marginal_nodes.end(); ++nit ) {
          // if the node is located on the model boundary, it is a perimeter node
          if ( (*nit)->AtBoundary() != NOT ) nodes_.push_back( (*nit) );
          else {
               // if the node is not contained in the interior, it is a perimeter / halo node
               if ( !binary_search( nodes_.begin(), interior_nodes_end, (*nit) ) ) {
                     nodes_.push_back( (*nit) );
                 }
            }
      }
    // sorting the perimeter nodes
    sort( next(nodes_.begin(),n_interior_nodes_), nodes_.end() );

 } // end Rebuild


/*

    /// adds multiple nodes and potential extra elements to region, does not remove any nodes or elements
template<size_t dim>

void NimbleRegion<dim>::AddNodes( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last )
 {
 }
 
    /// removes nodes and elements that might have become disconnected from the region
template<size_t dim>

void NimbleRegion<dim>::RemoveNodes(std::vector<Node<dim>*>& )
 {
 }
*/



// USER METHODS

// retrieving information

template<size_t dim>

size_t NimbleRegion<dim>::Nodes() const
 {
     return nodes_.size();
 }


template<size_t dim>

size_t NimbleRegion<dim>::InteriorNodes() const
 {
     return n_interior_nodes_;
 }


template<size_t dim>

size_t NimbleRegion<dim>::PerimeterNodes() const
 {
     return nodes_.size() - n_interior_nodes_;
 }


template<size_t dim>

size_t NimbleRegion<dim>::Elements() const
 {
     return elements_.size();
 }


/**
    Consecutive numbering, interior nodes first, perimeter nodes second.
*/
template<size_t dim>

size_t NimbleRegion<dim>::RenumberNodes() const
 {
     size_t node_number(0U);
     const typename vector<Node<dim>*>::const_iterator nodes_end(nodes_.end());
     for ( typename vector<Node<dim>*>::const_iterator nit=nodes_.begin(); nit!=nodes_end; ++nit )
       (*nit)->Idx( node_number++ );
   
     return node_number;
 }


// accessors
template<size_t dim>
typename NimbleRegion<dim>::CellType* const NimbleRegion<dim>::E( size_t n )
 {
    assert( n<elements_.size() );
    return elements_[n];
 }


template<size_t dim>
const typename NimbleRegion<dim>::CellType* const NimbleRegion<dim>::E( size_t n ) const
 {
    assert( n<elements_.size() );
    return elements_[n];
 }

 
template<size_t dim>

Node<dim>* const NimbleRegion<dim>::N( size_t n )
 {
    assert( n < nodes_.size() );
    return nodes_[n];
 }


template<size_t dim>

const Node<dim>* const NimbleRegion<dim>::N( size_t n ) const
 {
    assert( n < nodes_.size() );
    return nodes_[n];
 }
 

// iterators
template<size_t dim>

typename std::vector<csmp::Node<dim>*>::iterator  NimbleRegion<dim>::NodesBegin() { return nodes_.begin(); }

template<size_t dim>

typename std::vector<csmp::Node<dim>*>::iterator  NimbleRegion<dim>::NodesEnd() { return nodes_.end(); }

template<size_t dim>

typename std::vector<csmp::Node<dim>*>::const_iterator  NimbleRegion<dim>::NodesBegin() const  { return nodes_.begin(); }

template<size_t dim>

typename std::vector<csmp::Node<dim>*>::const_iterator  NimbleRegion<dim>::NodesEnd() const  { return nodes_.end(); }

template<size_t dim>

typename std::vector<Element
<dim>*>::iterator  NimbleRegion<dim>::ElementsBegin() { return elements_.begin(); }

template<size_t dim>

typename std::vector<Element
<dim>*>::iterator  NimbleRegion<dim>::ElementsEnd() { return elements_.end(); }

template<size_t dim>

typename std::vector<Element
<dim>*>::const_iterator NimbleRegion<dim>::ElementsBegin() const { return elements_.begin(); }

template<size_t dim>

typename std::vector<Element
<dim>*>::const_iterator NimbleRegion<dim>::ElementsEnd() const { return elements_.end(); }


/**
    Assuming that there was no remeshing, the region can only grow at its perimeter.
    Thus all former interior nodes will remain interior nodes.
 
    Where the new nodes are added, former perimeter nodes may become interior nodes.
 
    New elements will be added in the process.
*/

template<size_t dim>

void NimbleRegion<dim>::AddNodes( typename std::vector<Node<dim>*>::iterator first,
                                       typename std::vector<Node<dim>*>::iterator last )
{
   // inserting the new nodes at the end of the vector
	 std::copy( first, last, std::inserter(nodes_, nodes_.end()));
	 set<Node<dim>*> set(nodes_.begin(), nodes_.end());
	 nodes_.assign(set.begin(), set.end());

}



// LUAT's METHODS - COMMENTED OUT - NOT SURE WHAT THESE METHODS DO BECAUSE THERE IS NO DOCUMENTATION

/*
template<size_t dim>

void NimbleRegion<dim>::CollectRegionElements() {
	for (auto& nIter = nodes_.cbegin(); nIter != nodes_.cend(); ++nIter) {
		for (size_t i(0); i < (*nIter)->Parents(); ++i) {
			elements_.push_back((*nIter)->Parent(i));
		}
	}
	set<CellType*> set(elements_.begin(), elements_.end());
	elements_.assign(set.begin(), set.end());
}

// use it if one wants to add large number of node
template<size_t dim>

void NimbleRegion<dim>::CreatePerimeters() {
	set<Node<dim>*> internalNode(nodes_.begin(), nodes_.end());
	set<Node<dim>*> p_set;
	set<Node<dim>*>::iterator nIter;
	for (auto& eIter = elements_.cbegin(); eIter != elements_.cend(); ++eIter) {
		for (size_t i(0); i < (*eIter)->Nodes(); ++i) {
			nIter = internalNode.find((*eIter)->N(i));
			if (nIter == internalNode.end()) {
				p_set.insert((*eIter)->N(i));
			}
		}
	}
	perimeters_nodes_.assign(p_set.begin(), p_set.end());
	std::copy(begin(perimeters_nodes_), end(perimeters_nodes_), std::back_inserter(nodes_));
}
*/

template<size_t dim>

void NimbleRegion<dim>::Erase() {
	nodes_.clear();
	elements_.clear();
	n_interior_nodes_ = 0U;
}



template<size_t dim>

typename std::vector<csmp::Node<dim>*>::const_iterator NimbleRegion<dim>::PerimeterNodesBegin() const
  {
     return next( nodes_.begin(), n_interior_nodes_ );
  }


template<size_t dim>

typename std::vector<csmp::Node<dim>*>::iterator NimbleRegion<dim>::PerimeterNodesBegin()
  {
     return next( nodes_.begin(), n_interior_nodes_ );
  }



template class NimbleRegion<1>;
template class NimbleRegion<2>;
template class NimbleRegion<3>;

} // end namespace csmp
