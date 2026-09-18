// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_NT_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_NT_op_N_dV<dim,CELL>::Integral_NT_op_N_dV( const PropertyDatabase<dim>& pref,
                                                    const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    INN(3,3)
 {
    MathOperatorRHS<dim,CELL>::Name("Integral_NT_op_N_dV", oper, test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_NT_op_N_dV<dim>::(constructor)", 
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_op_N_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }





/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
     // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

  // reading Young's modulus (must be an element variables)
   e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), sc );
}




/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());

    // consistent formulation
    if ( !MathOperatorRHS<dim,CELL>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( INN );
         fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
         // the matrix is contracted into a vector
         for ( auto i{0U}; i<e.Nodes(); i++ ) 
           for ( auto j{0U}; j<e.Nodes(); j++ ) 
             MathOperatorRHS<dim,CELL>::RHS[i] += INN(i,j) * sc();
      }
    // lumped formulation  
    else
      {
         double res = (e.Volume() * sc()) / static_cast<double>(e.Nodes());
         fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), res );
      }

} // end ComputeContribution



template class Integral_NT_op_N_dV<1U,Element>;
template class Integral_NT_op_N_dV<2U,Element>;
template class Integral_NT_op_N_dV<3U,Element>;

template class Integral_NT_op_N_dV<1U,Face>;
template class Integral_NT_op_N_dV<2U,Face>;
template class Integral_NT_op_N_dV<3U,Face>;

} // csmp
