#include "Integral_NT_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
Integral_NT_op_N_dV<dim,SIMPLEX>::Integral_NT_op_N_dV( const PropertyDatabase<dim>& pref,
                                                       const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    INN(3,3)
 {
    MathOperatorRHS<dim>::Name("Integral_NT_op_N_dV", oper, test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_NT_op_N_dV<dim>::(constructor)", 
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_op_N_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }





/** Reads the Operand values from the elements.
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_op_N_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
     // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

  // reading Young's modulus (must be an element variables)
   e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), sc );
}




/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_op_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());

    // consistent formulation
    if ( !MathOperatorRHS<dim>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( INN );
         fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0 );
         // the matrix is contracted into a vector
         for ( size_t i=0; i<e.Nodes(); i++ ) 
           for ( size_t j=0; j<e.Nodes(); j++ ) 
             MathOperatorRHS<dim>::RHS[i] += INN(i,j) * sc();
      }
    // lumped formulation  
    else
      {
         double64 res = (e.Volume() * sc()) / static_cast<double64>(e.Nodes());
         fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), res );
      }

} // end ComputeContribution



template class Integral_NT_op_N_dV<1U,Element<1U> >;
template class Integral_NT_op_N_dV<2U,Element<2U> >;
template class Integral_NT_op_N_dV<3U,Element<3U> >;

template class Integral_NT_op_N_dV<1U,Face<1U> >;
template class Integral_NT_op_N_dV<2U,Face<2U> >;
template class Integral_NT_op_N_dV<3U,Face<3U> >;

} // csmp
