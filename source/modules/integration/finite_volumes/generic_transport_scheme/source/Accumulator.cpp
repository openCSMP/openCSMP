#include "Accumulator.h"
#include "IntegralEquation.h"
#include "ImplicitTransport.h"
#include "Region.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {


/**
  Should use information from equation to define computational problem (solution variables etc.)
*/
template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::SetUp() 
{
   // sizing the linear algebraic system
   User()->ComputationDomain().UpdateMemberIndexes();
   User()->LinearSystem().Resize( User()->ComputationDomain().Nodes() );
   
   // filling the solution vector with the current value of the transport variable
   size_t node_counter(0U);
   for ( typename vector<Node<dim>*>::const_iterator
         nit=User()->ComputationDomain().NodesBegin(); nit!=User()->ComputationDomain().NodesEnd(); ++nit )
     User()->LinearSystem().X[node_counter++] = (*nit)->Read( User()->Notation.key_C0 );
   
}


/**
   Default sequence
     1. pore volume in diagonal of LHS
     2. upwinded flows coupling diagonal with off-diagonal and outgoing flows into diagonal
     3.  pore volume times concentration into rhs
     4. sources and sinks into rhs
     
    Inside of the domain and for the halo stencils, accumulation by stencil is used, while truncated boundary FVs need to be accumulated finite volume by finite volume (only certain terms)
*/
template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::Accumulate( double time_increment )
 {
    ErrorHandler& error_handler( ErrorHandler::Instance() );

    // multipass process (processing all operators with execution level 1, which should be at the beginning of ordered map)
    for ( typename GoverningEquation<dim>::MatrixOperatorConstIterator
             et=User()->LHS_OperatorsBegin(); et!= User()->LHS_OperatorsEnd(); ++et ) {
            // ascertaining that operations are performed in meaningful sequence
            int execution_level = (*et).first.ExecutionLevel();
            if ( execution_level == 1 && ( (*et).first.OperationType() == MULTIPLY || (*et).first.OperationType() == DIVIDE ) ) {
                 (*et).first.Out();
                 cerr << typeid( (*et).second ).name();
                 error_handler.Note( ERROR, "Accumulator<dim,SUBDOMAIN,USER>::Accumulate",
                                      "mutiplication/division on empty matrix has no effect; nothing was done.");
                 continue;  
              }
            // imposing operation onto MatrixOperator (time-increment was set before)
            double factor(1.);
            switch( (*et).first.OperationType() ) {
                 case SUBTRACT: factor *= -1.;
                   break;
                 case MULTIPLY: error_handler.Note( ERROR, "Accumulator<dim,SUBDOMAIN,USER>::Accumulate", "MULTIPLY not implemented yet" );
                      // TODO: mutiply corresponding matrix entries with the ones from the matrix?
                   break;
                 case DIVIDE: error_handler.Note( ERROR, "Accumulator<dim,SUBDOMAIN,USER>::Accumulate", "DIVIDE not implemented yet" );
                   break;
                 default: factor = 1.; // = ADD
              }
            if ( (*et).first.MultiplyWithTimeIncrement() ) factor *= time_increment;
            (*et).second->Factor( factor );
            
            // processing interior of solution domain
            AccumulateByStencil( User()->ComputationDomain().CellsBegin(), 
                                 User()->ComputationDomain().CellsEnd(), 
                                 (*et).second, User()->LinearSystem().LHS );  
                                   
            // processing potential halo of solution domain
            if ( User()->HasHaloStencils() )
              AccumulateByStencil( User()->HaloCellsBegin(), User()->HaloCellsEnd(), 
                                   (*et).second, 
                                   User()->LinearSystem().LHS );    
       }

    // boundary of solution domain
 }





    /// compensate (+) balance at inflow boundaries, (-) balance at outflow boundaries, and any potential divergence of flow at no-flow boundaries              
template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::BalanceFlowsThroughTruncatedBoundaryFiniteVolumes( SparseMatrix& lhs, 
                                                                               vector<double>& rhs,
                                                                               double time_increment )
 {
    // 1. compensate outflows by adding the (-) flux balances to the diagonal of the solution matrix
    
    // 2. compensate inflows by subtracting the (+) flux balance at truncated boundary FVs to the RHS
    
    // 3. for no-flow boundaries make sure that the flux-balance where there are no sources and sinks is equal to zero 
    
 } // end BalanceFlowsThroughTruncatedBoundaryFiniteVolumes






/**
    Accumulates lefthand-side = matrix A of the equation A x = b
*/
template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::AccumulateByStencil( typename vector<Element<dim>*>::const_iterator first, 
                                                           typename vector<Element<dim>*>::const_iterator last,
                                                           const MatrixOperator<dim>* const mat_op,
                                                           SparseMatrix& lhs ) const
 {
   const typename vector<Element<dim>*>::const_iterator cells_end(last);
   // accumulation
   while( first != cells_end ) {
        mat_op->AccumulateStencil( *(*first), lhs );
        first++;
     }

 } // end AccumulateByStencil



template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::AccumulateByStencil( typename vector<Element<dim>*>::const_iterator first, 
                                                           typename vector<Element<dim>*>::const_iterator last,
                                                           const VectorOperator<dim>* const vec_op,
                                                           vector<double>& rhs ) const
 {
   const typename vector<Element<dim>*>::const_iterator cells_end(last);
   // accumulation
   while( first != cells_end ) {
        vec_op->AccumulateStencil( *(*first), rhs );
        first++;
     }

 } // end AccumulateByStencil




template<uint32_t dim, template<uint32_t> class USER>
void Accumulator<dim,USER>::AccumulateByFiniteVolume( typename std::vector<Node<dim>*>::const_iterator nit,
                                                                typename std::vector<Node<dim>*>::const_iterator end,
                                                                const MatrixOperator<dim>* const mat_op, ///< modified to contain operation information
                                                                SparseMatrix& lhs ) const
{
   const typename std::vector<Node<dim>*>::const_iterator  nodes_end(end);
   while( nit != nodes_end ) {
        mat_op->AccumulateFiniteVolume( *(*nit), lhs );
        nit++;
     } 
  
} // end AccumulateByFiniteVolume


template<uint32_t dim,template<uint32_t> class USER>
void Accumulator<dim,USER>::AccumulateByFiniteVolume( typename std::vector<Node<dim>*>::const_iterator nit,
                                                                typename std::vector<Node<dim>*>::const_iterator end,
                                                                const VectorOperator<dim>* const vec_op, ///< modified to contain operation information
                                                                std::vector<double>& rhs ) const
{
   const typename std::vector<Node<dim>*>::const_iterator nodes_end(end);
   while( nit != nodes_end ) {
        vec_op->AccumulateFiniteVolume( *(*nit), rhs );
        nit++;
     } 
  
} // end AccumulateByFiniteVolume

  
template class Accumulator<3U,ImplicitTransport>;

} // end csmp

