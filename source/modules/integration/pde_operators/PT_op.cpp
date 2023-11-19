#include "PT_op.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
PT_op<dim,CELL>::PT_op( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test)
 {
    MathOperatorRHS<dim,CELL>::Name("PT_op", oper, test );
 
     if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != VECTOR || 
          MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != NODE ) 
       throw csmp::Exception( FATAL_ERROR, "MathOperatorRHS->PT_op<dim>::(constructor):",
                              oper, "Operand must be a vector variable placed on the nodes." );
 }





/**
 
Reads node data (=vector properties) for further processing
by the ComputeContribution() method.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void PT_op<dim,CELL>::GetOperands( const CELL<dim>& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), NODAL_FORCE );
   }



/**
 
Accumulates nodal vector forces from nodes which are flagged Dirichlet with
regard to the target variable. The Operand is therefore a node
variable.  

@section arguments Input Arguments 

The element from which accumulation into the global matrix takes place
and the time-increment over which the force will be applied.  

The result gets stored into the MathOperator right-hand vector.

@section application Application

To compute nodal forces acting on the boundary of a model.
 
*/
template<uint32_t dim, template<uint32_t> class CELL>
void PT_op<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
   MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() * dim );
   fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );
  
   const size_t nodes(e.Nodes());
  
   for ( uint32_t i{0u}; i<nodes; i++ )
     for ( uint32_t j{0U}; j<dim; j++ )
         // forces must be divided by number of elements they will be accumulated from
         // to avoid multiple accumulation
         MathOperatorRHS<dim,CELL>::RHS[ i * dim + j ] = NODAL_FORCE[i][j] / static_cast<double>(e.N(i)->Parents());

} // end ComputeContribution
     

template class PT_op<1U>;
template class PT_op<2U>;
template class PT_op<3U>;

template class PT_op<1U,Face>;
template class PT_op<2U,Face>;
template class PT_op<3U,Face>;

} // csmp
