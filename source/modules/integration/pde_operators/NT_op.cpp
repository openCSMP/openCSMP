#include "NT_op.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {



template<size_t dim,class SIMPLEX>
NT_op<dim,SIMPLEX>::NT_op( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test)
 {
    MathOperatorRHS<dim>::Name("NT_op", oper, test );
 
     if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR || 
          MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE ) 
     {
        throw csmp::Exception( FATAL_ERROR, "MathOperatorRHS->NT_op<dim>::(constructor)", 
                        oper, "Operand must be a scalar variable placed on the nodes." );
     }
 }





/** Reads scalar node data for further processing
by the ComputeContribution() method.  
*/
template<size_t dim,class SIMPLEX>
void NT_op<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), M_ );
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
template<size_t dim,class SIMPLEX>
void NT_op<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
   MathOperatorRHS<dim>::RHS.resize( e.Nodes() * dim );
   fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
  
   for ( size_t i=0U; i<e.Nodes(); i++ )
     if ( M_[i].Flag() == NEUMANN )
       // contributions must be divided by number of elements which share the node
       // to avoid multiple accumulation
       MathOperatorRHS<dim>::RHS[i] = M_[i].Value() / static_cast<double64>(e.N(i)->Parents());
     
} // end ComputeContribution


template class NT_op<1U,Element<1U> >;
template class NT_op<2U,Element<2U> >;
template class NT_op<3U,Element<3U> >;

template class NT_op<1U,Face<1U> >;
template class NT_op<2U,Face<2U> >;
template class NT_op<3U,Face<3U> >;

} // csmp
