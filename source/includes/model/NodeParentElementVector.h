//
//  NodeParentElementVector.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 22/8/2024.
//

#ifndef CSMP_NODE_PARENT_ELEMENT_VECTOR_H
#define CSMP_NODE_PARENT_ELEMENT_VECTOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Element;

/**
  @pre values must be unique, indices not
  @pre one cannot search by index
*/
template<uint32_t dim>
class NodeParentElementVector {
  public:
    /// to pre-reserve some elements to control over memory growth
    explicit NodeParentElementVector( size_t capacity );
    NodeParentElementVector() = default;
    ~NodeParentElementVector() = default;
    NodeParentElementVector( const NodeParentElementVector<dim>& ) = default;
    NodeParentElementVector( NodeParentElementVector<dim>&& ) = default;
    NodeParentElementVector<dim>& operator=( const NodeParentElementVector<dim>& ) = default;
    
    /// to avoid unneccasary reallocations of the storage in case the desired number of entries is known before multiple Assign's are called
    void Reserve( size_t expected_size );
    
    /// assign parent-element pointer - parent-element node-number pairs to the storage
    void Assign( Element<dim>* const, uint32_t index );
    
    /// erase parent element from vector
    void Unassign( const Element<dim>* const );
    
    /// to prune excess capacity
    void ShrinkToFit();
    
    /// get rid of current storage contents
    void Erase();
    
    // SortParents() - not needed because linear search is faster than binary_search for <50 elmts
    
    /// find integer associated with given value
    uint32_t LocalNodeNumber( uint32_t vector_entry ) const;
    uint32_t LocalNodeNumber( const Element<dim>* const ) const;
    
    /// retrieving parent element of the Node by vector index (implemented with find_if)
    Element<dim>* ParentElement( uint32_t vector_entry );
    const Element<dim>* const ParentElement( uint32_t vector_entry ) const;
   
    /// is contained in the vector?
    bool IsParent( const Element<dim>* const ) const;

    /// parents
    uint32_t Parents() const;
    size_t   Size() const;
    size_t   Capacity() const { return values_.capacity(); }
 
     /// iterators
    typename std::vector<std::pair<Element<dim>*,short> >::const_iterator ParentsBegin() const;
    typename std::vector<std::pair<Element<dim>*,short> >::const_iterator ParentsEnd() const;
    
    /// storage requirements
    size_t SizeOf() const;
    
    void Out() const;

  private:
    //                   value pointer, local integer
    std::vector<std::pair<Element<dim>*,short> > values_;  ///< parent element pointers
};



} // end csmp


#endif /* CSMP_NODE_PARENT_ELEMENT_VECTOR_H */
