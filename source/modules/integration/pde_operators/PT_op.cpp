#include "PT_op.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
PT_op<dim,SIMPLEX>::PT_op( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test)
 {
    MathOperatorRHS<dim>::Name("PT_op", oper, test );
 
     if ( MathOperatorRHS<dim>::MaterialOperandType() != VECTOR || 
          MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE ) 
       throw csmp::Exception( CSMP_FATAL_ERROR, "MathOperatorRHS->PT_op<dim>::(constructor):",
                              oper, "Operand must be a vector variable placed on the nodes." );
 }





/**
 
Reads node data (=vector properties) for further processing
by the ComputeContribution() method.  
*/
template<size_t dim,class SIMPLEX>
void PT_op<dim,SIMPLEX>::GetOperands(  SIMPLEX& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), NODAL_FORCE );
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
template<size_t dim,class SIMPLEX>
void PT_op<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
   MathOperatorRHS<dim>::RHS.resize( e.Nodes() * dim );
   fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
  
   const size_t nodes(e.Nodes());
  
   for ( size_t i=0U; i<nodes; i++ )
     for ( size_t j=0U; j<dim; j++ )
         // forces must be divided by number of elements they will be accumulated from
         // to avoid multiple accumulation
         MathOperatorRHS<dim>::RHS[ i * dim + j ] = NODAL_FORCE[i][j] / static_cast<double64>(e.N(i)->Parents());

} // end ComputeContribution
     

template class PT_op<1U,Element<1U> >;
template class PT_op<2U,Element<2U> >;
template class PT_op<3U,Element<3U> >;

template class PT_op<1U,Face<1U> >;
template class PT_op<2U,Face<2U> >;
template class PT_op<3U,Face<3U> >;

} // csmp
