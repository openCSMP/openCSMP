//
//  SourceTermRHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "SourceTermRHS.h"
#include "Node.h"
#include "Element.h"

using namespace std;

namespace csmp {

    /// note that this source variable may be an element or a nodal one
template<size_t dim>
SourceTermRHS<dim>::SourceTermRHS( const csmp::INDEX<SCALAR,ELEMENT>& source_variable_key,
                                   const csmp::INDEX<SCALAR,NODE>& nodal_source_variable_key,
                                   const csmp::INDEX<SCALAR,NODE>& transported_variable_key  )
 : esrc_key_(source_variable_key),
   nsrc_key_(nodal_source_variable_key),
   adv_key_(transported_variable_key)
 {
 }


template<size_t dim>
void SourceTermRHS<dim>::AccumulateFiniteVolume( const Node<dim>& n, std::vector<double64>& rhs ) const
 {
    rhs[ n.Idx() ] += n.Read( nsrc_key_ ) * n.Read( adv_key_ );
 }
  

    /// accumulates distributed values of the source term on the finite volume
template<size_t dim>
void SourceTermRHS<dim>::AccumulateStencil( const Element<dim>& e, std::vector<double64>& rhs ) const
 {
    double64 esource = e.Read( esrc_key_ );
    if ( fabs(esource) < numeric_limits<double64>::epsilon() ) return;
    
    const size_t nodes = e.Nodes();
    esource /= static_cast<double64>(nodes);
    for ( size_t i=0U; i<nodes; ++i )
      rhs[ e.N(i)->Idx() ] += esource * e.N(i)->Read( adv_key_ );
 }

template class SourceTermRHS<1U>;
template class SourceTermRHS<2U>;
template class SourceTermRHS<3U>;

} // end csmp
