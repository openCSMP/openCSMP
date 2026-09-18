// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "MeshRenumberVisitor.h"
#include "Node.h"
#include "Element.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
MeshRenumberVisitor<dim>::MeshRenumberVisitor()
    : e_counter_(0U), n_counter_(0U)
  { 
  }



// MeshRenumberVisitor Methods ====================

template<uint32_t dim>
void MeshRenumberVisitor<dim>::Visit( Element<dim>* eptr ) 
 { 
    eptr->Idx( e_counter_++ );

    for ( typename std::vector<csmp::Node<dim>*>::const_iterator
          nit=eptr->NodesBegin(); nit!=eptr->NodesEnd(); nit++ )
      {
         pair<typename set<csmp::Node<dim>*>::iterator,bool>  np(node_ptrs_.insert(*nit));
         if ( np.second ) (*nit)->Idx( n_counter_++ );
      }
    
 }


template<uint32_t dim>
size_t MeshRenumberVisitor<dim>::VisitedElements() const
 { return e_counter_; }
 

template<uint32_t dim>
size_t MeshRenumberVisitor<dim>::VisitedNodes() const
 { return n_counter_; }


template<uint32_t dim>
void MeshRenumberVisitor<dim>::Reset()
 { e_counter_=0U; n_counter_=0U; node_ptrs_.clear(); }


template class MeshRenumberVisitor<1U>;
template class MeshRenumberVisitor<2U>;
template class MeshRenumberVisitor<3U>;

} // end namespace csmp
