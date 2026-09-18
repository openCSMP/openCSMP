// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NT_op.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NT_op<dim,CELL>::NT_op( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test)
 {
    MathOperatorRHS<dim,CELL>::Name("NT_op", oper, test );
 
     if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR || 
          MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != NODE ) 
     {
        throw csmp::Exception( FATAL_ERROR, "MathOperatorRHS->NT_op<dim>::(constructor)", 
                               oper, "Operand must be a scalar variable placed on the nodes." );
     }
 }





/** Reads scalar node data for further processing
by the ComputeContribution() method.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NT_op<dim,CELL>::GetOperands( const CELL<dim>& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), M_ );
   }



/**
 
Accumulates nodal vector forces from nodes which are flagged Neumann with
regard to the target variable. The Operand is therefore a node
variable.  

@section arguments Input Arguments 

The element from which accumulation into the global matrix takes place
and the time-increment over which the force will be applied.  

The result is returned into the MathOperator vector.


@section application Application

To compute nodal forces acting on the boundary of a model.  
 
@todo (1) Not tested yet, maybe, the loads must be on the midside nodes 

*/
template<uint32_t dim, template<uint32_t> class CELL>
void NT_op<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
   MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() * dim );
   fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );
  
   for ( auto i{0}; i<e.Nodes(); i++ )
     if ( M_[i].Flag() == NEUMANN )
       // contributions must be divided by number of elements which share the node
       // to avoid multiple accumulation
       MathOperatorRHS<dim,CELL>::RHS[i] = M_[i]() / static_cast<double>(e.N(i)->Parents());
     
} // end ComputeContribution


template class NT_op<1U>;
template class NT_op<2U>;
template class NT_op<3U>;

template class NT_op<1U,Face>;
template class NT_op<2U,Face>;
template class NT_op<3U,Face>;

} // csmp
