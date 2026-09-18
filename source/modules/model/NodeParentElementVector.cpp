// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  NodeParentElementVector.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 23/8/2024.
//

#include "NodeParentElementVector.h"
#include "Element.h"

using namespace std;

namespace csmp {

// INLINE FUNCTIONS

/// to have some control over memory growth
template<uint32_t dim>
NodeParentElementVector<dim>::NodeParentElementVector( size_t size ) { values_.reserve(size); }
    

/// erase parent element from vector
template<uint32_t dim>
void NodeParentElementVector<dim>::Reserve( size_t size )
 {
    values_.reserve(size);
 }




template<uint32_t dim>
void NodeParentElementVector<dim>::Assign( Element<dim>* const elmt, uint32_t index )
 {
    // rejecting nullpointer
    if ( elmt == nullptr ) return;
    // consistency check
    assert( index < elmt->Nodes() );
    // linear search to check whether value already exists
    if ( !values_.empty() ) {
         auto it = find_if( values_.begin(), values_.end(),
                            [elmt]( const pair<Element<dim>*,uint32_t>& val ) -> bool
                                         { return ( elmt == val.first ); } );
         if ( it == values_.end() )
           values_.push_back( make_pair( elmt, static_cast<short>(index) ) );
      }
    else
      values_.push_back( make_pair( elmt, static_cast<short>(index) ) );
 }



/// erase parent element from vector
template<uint32_t dim>
void NodeParentElementVector<dim>::Unassign( const Element<dim>* const elmt )
 {
    // linear search to check whether value already exists
    values_.erase( remove_if( values_.begin(), values_.end(),
                                   [elmt]( const pair<Element<dim>*,short>& val ) -> bool
                                     { return ( elmt == val.first ); }
                                 ), values_.end() );
  }



/// erase parent element from vector
template<uint32_t dim>
void NodeParentElementVector<dim>::Erase()
 {
    //  vector<T>().swap(x); - swaps memory with an empty vector
    vector<pair<Element<dim>*,short> >().swap(values_);
    // clear is not guaranteed to release the memory
    // values_.clear();
 }



/// find integer associated with given value
template<uint32_t dim>
uint32_t NodeParentElementVector<dim>::LocalNodeNumber( uint32_t vector_entry ) const
  { assert(vector_entry < values_.size()); return static_cast<uint32_t>(values_[vector_entry].second); }


/**
       Handcoded linear search.
*/
template<uint32_t dim>
uint32_t NodeParentElementVector<dim>::LocalNodeNumber( const Element<dim>* const eptr ) const
 {
    for ( const auto& it : values_ )
      if ( eptr == it.first )
        return static_cast<uint32_t>(it.second);

    // if not found
    return numeric_limits<uint32_t>::max();
 }


/* profiled as slower than LocalNodeNumber()
template<uint32_t dim>
inline uint32_t NodeParentElementVector<dim>::LocalNodeNumber1( const Element<dim>* const eptr ) const
 {
    auto it = find_if( values_.begin(), values_.end(),
                            [eptr](const pair<Element<dim>* const,short>& elmt ){ return ( eptr == elmt.first ); } );
    // if not found
    if ( it == values_.end() ) return numeric_limits<uint32_t>::max();
    return static_cast<uint32_t>((*it).second);
 }
*/
 


/// retrieving parent element of the Node by vector index
template<uint32_t dim>
Element<dim>* NodeParentElementVector<dim>::ParentElement( uint32_t vector_entry )
  { assert(vector_entry<values_.size()); return values_[vector_entry].first; }

template<uint32_t dim>
const Element<dim>* const NodeParentElementVector<dim>::ParentElement( uint32_t vector_entry ) const
  { assert(vector_entry<values_.size()); return values_[vector_entry].first; }


template<uint32_t dim>
bool NodeParentElementVector<dim>::IsParent( const Element<dim>* const eptr ) const
 {
    for ( const auto& it : values_ )
      if ( eptr == it.first )
        return true;

    return false;
 }

/* profiled slower version
template<uint32_t dim>
inline bool NodeParentElementVector<dim>::IsParent( Element<dim>* const elmt ) const
 {
    // linear search to check whether Element is contained in the vector
    auto it = find_if( values_.begin(), values_.end(),
                            [elmt]( const pair<Element<dim>*,uint32_t>& val ) -> bool
                              { return ( elmt == val.first ); }
                          );
    if ( it != values_.end() ) return true;
    return false;
 }
*/


template<uint32_t dim>
typename vector<pair<Element<dim>*,short> >::const_iterator NodeParentElementVector<dim>::ParentsBegin() const
  { return values_.begin(); }

template<uint32_t dim>
typename vector<pair<Element<dim>*,short> >::const_iterator NodeParentElementVector<dim>::ParentsEnd() const
  { return values_.end(); }



template<uint32_t dim>
uint32_t NodeParentElementVector<dim>::Parents() const { return static_cast<uint32_t>(values_.size()); }

template<uint32_t dim>
size_t NodeParentElementVector<dim>::Size() const { return values_.size(); }
    

template<uint32_t dim>
size_t NodeParentElementVector<dim>::SizeOf() const
 { return sizeof(*this) + values_.size() * sizeof(pair<Element<dim>*,short>); }
    

/// to prune excess capacity from vector
template<uint32_t dim>
void NodeParentElementVector<dim>::ShrinkToFit() { values_.shrink_to_fit(); }
    

template<uint32_t dim>
void NodeParentElementVector<dim>::Out() const {
   cout <<"\n"<<"NodeParentElementVector: size: "<< values_.size() <<": of Element - local-node-number pairs:"<< endl;
   for ( const auto& it : values_ )
        cout <<"\t\t"<< it.first->Idx() <<": "<< parseFiniteElementType( it.first->FE_Type() ) <<": "<< it.second << endl;
}


template class NodeParentElementVector<3U>;
template class NodeParentElementVector<2U>;
template class NodeParentElementVector<1U>;

} // end csmp
