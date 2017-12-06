//
//  SourceTermRHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "SourceTermRHS.h"
#include "Model.h"
#include "Node.h"
#include "Element.h"

using namespace std;

namespace csmp {

    /// note that this source variable may be an element or a nodal one
template<size_t dim>
SourceTermRHS<dim>::SourceTermRHS( const Model<dim>& model, const char* source_variable, ACCUMULATION_MODE accumulation_mode )
 : VectorOperator<dim>(accumulation_mode),
   src_key_(model.Database().StorageKey(source_variable))
 {
 }


template<size_t dim>
void SourceTermRHS<dim>::AccumulateFV( const Node<dim>* nptr, std::vector<double64>& rhs ) const
 {
    assert( src_key_.place == NODE );
    assert( src_key_.type  == SCALAR );
    rhs[ nptr->Idx() ] += nptr->Read( src_key_ );
 }
  

    /// accumulates distributed values of the source term on the finite volume
template<size_t dim>
void SourceTermRHS<dim>::AccumulateStencil( const Element<dim>* eptr, std::vector<double64>& rhs ) const
 {
 }


} // end csmp
